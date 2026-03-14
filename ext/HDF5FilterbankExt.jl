module HDF5FilterbankExt

import Blio.Filterbank: Header, fil2h5, mmap, h52fil
import Base.write
import HDF5: write_attribute

if isdefined(Base, :get_extension)
    using HDF5: File, Dataset, h5open, attributes, create_dataset, datatype,
        write_dataset, Filters, get_chunk
else
    using ..HDF5: File, Dataset, h5open, attributes, create_dataset, datatype,
        write_dataset, Filters, get_chunk
end

"""
    fil2h5(fbname[, h5name]; kwargs...)
    fil2h5(fbname[, h5name]; copy=true, chunk, filters, kwargs...)

Convert a SIGPROC filterbank file into an HDF5 file.

The first form without the `copy` keyword argument (or with `copy=false`) will
create an HDF5 file with an external `data` dataset (i.e. it will point to the
data that lives inside the filterbank file).  The second form with `copy=true`
will copy the data from the filterbank file into the HDF5 file.  In this case,
the data will be chunked and compressed according the `chunk` and `filters`
keyword arguments.

In both forms, `h5name` will default to `fbname` with any `.fil` extension
removed and a `.h5` extension added.  An exception will be thrown if `h5name`
already exists.  Additional `kwargs` may be passed to provide additional (and/or
corrected) attributes for the `data` dataset.

The `chunk` keyword argument should be passed as a tuple of 3 integers.  The
default value is `(512,1,16)`.

The `filters` argument should be passed as an `HDF5.Filter` object or a Vector
of `HDF5.Filter` objects.  The default for `filters` is `[Shuffle(),
Deflate(3)]`, which uses filters that are built-in the HDF5 library.  Another
option (of many) is to pass `filters=BitshuffleFilter(compressor=:lz4)` to use
the bitshuffle filter from the `H5Zbitshuffle` package, which you would need to
add as a dependency of your project because Blio does not depend on it.
"""
function fil2h5(fbname, h5name=replace(fbname, r"\.fil$"=>"") * ".h5";
    copy=false, chunk=(512,1,16), filters=[Filters.Shuffle(), Filters.Deflate(3)],
    kwargs...
)
    ispath(h5name) && error("$h5name already exists")

    fbh = open(io->read(io, Header), fbname)

    h5open(h5name, "w") do h5
        # For blimpy compatability
        write_attribute(h5, "CLASS", "FILTERBANK")
        write_attribute(h5, "VERSION", "1.0")

        # Create `data` dataset
        size = (fbh[:nchans], fbh[:nifs], fbh[:nsamps])
        kwchunk = isempty(chunk) ? () : (; chunk)
        data = if copy
            fbd = open(io->mmap(io, fbh), fbname)
            ds = create_dataset(h5, "data", Float32, size; kwchunk..., filters)
            write_dataset(ds, datatype(fbd), fbd)
            ds
        else
            external = (fbname, fbh[:header_size], fbh[:data_size])
            create_dataset(h5, "data", Float32, size; external)
        end

        # Write fbh as attributes of `data` datset
        write_attribute(data, fbh)
        # Allow additional (or corrected) attributes from `kwargs`
        for (k,v) in kwargs
            write_attribute(data, string(k), v)
        end
    end

    h5name
end

"""
    h52fil(h5name[, fbname]; kwargs...)

Convert FBH5 HDF5 file `h5name` into SIGPROC filterbank file `fbname`.

`fbname` will default to `h5name` with any `.(h5|fbh5|hdf5)` extension
removed and a `.fil` extension added.  An exception will be thrown if `h5name`
already exists.  Additional `kwargs` may be passed to provide additional (and/or
corrected) SIGPROC Filterbank attributes.

Filterbank data that is `Int8` or `Float32` will be copied directly.  Other data
types will be converted to `Float32`.  The `nbits` header item will be set
appropriately.

Any filters required to read the `data` dataset in `h5name` must be loaded
beforehand (e.g. `using H5Zbitshuffle` if `h5name` uses bitshuffle compression).
"""
function h52fil(h5name, fbname=replace(h5name, r"\.(h5|fbh5|hdf5)$"=>"") * ".fil";
    kwargs...
)
    ispath(fbname) && error("$fbname already exists")

    # Open file, make Filterbank.Header, merge in kwargs, get `data` dataset
    h5open(h5name) do h5
        fbh = Header(h5)
        merge!(fbh, kwargs)
        data::Dataset = h5["data"]
        intype = eltype(data)
        outtype = intype ∈ (UInt8, Int8, Float32) ? intype : Float32

        # Ensure nbits matches outtype
        fbh[:nbits] = 8*sizeof(outtype)

        # Copy data in batches of size (nchans, nifs, chunksize[3])
        nchans = fbh[:nchans]
        nifs = fbh[:nifs]
        ntpb = get_chunk(data)[3]
        ibuf = Array{intype}(undef, nchans, nifs, ntpb)
        obuf = (intype === outtype) ? ibuf : Array{outtype}(undef, nchans, nifs, ntpb)

        nbatches, nleftover = divrem(size(data,3), ntpb)

        open(fbname, "w") do io
            # Write header
            write(io, fbh; quiet=true)

            # Copy batches of data
            i=1
            for b in 1:nbatches
                copyto!(ibuf, data, :, :, i:i+ntpb-1)
                if ibuf !== obuf
                    copyto!(obuf, ibuf)
                end
                write(io, obuf)
                i += ntpb
            end

            # Copy leftovers
            ivw = view(ibuf, :, :, 1:nleftover)
            ovw = view(obuf, :, :, 1:nleftover)

            copyto!(ivw, data, :, :, i:i+nleftover-1)
            if ibuf !== obuf
                copyto!(ovw, ivw)
            end
            write(io, ovw)
        end
    end

    fbname
end

"""
    Header(h5ds::HDF5.Dataset)::Filterbank.Header
    Header(h5::HDF5.File)::Filterbank.Header

Create `Filterbank.Header` object from dataset `h5ds` or `h5["data"]`.

The names/types of the attributes are not checked for Filterbank validity.  The
`data_size` and `nsamps` attributes will be added to the `Header` object if they
do not exist as attributes.
"""
function Header(h5ds::Dataset)::Header
    fbh = Header()
    attrs = attributes(h5ds)
    for k in keys(attrs)
        # Skip DIMENSION_LABELS, but keep all others (even non-SIGPROC ones)
        k == "DIMENSION_LABELS" && continue
        fbh[Symbol(k)] = attrs[k][]
    end
    # Add :data_size and :nsamps if not present
    if !haskey(fbh, :data_size)
        fbh[:data_size] = sizeof(eltype(h5ds)) * prod(size(h5ds))
    end
    if !haskey(fbh, :nsamps)
        fbh[:nsamps] = size(h5ds, ndims(h5ds))
    end
    fbh
end

function Header(h5::File)::Header
    Header(h5["data"])
end

"""
    write_attribute(h5ds::HDF5.Dataset, fbh::Filterbank.Header)
    write_attribute(h5::HDF5.File, fbh::Filterbank.Header)

Write `fbh` contents as attributes of dataset `h5ds` or `h5["data"]`.

The `data_size` and `nsamps` fields will not be written since they can be
inferred from the dimensions of the dataset.
"""
function write_attribute(h5ds::Dataset, fbh::Header)
    for (k,v) in fbh
        # Don't copy these synthetic attributes
        k == :data_size && continue
        k == :nsamps && continue
        write_attribute(h5ds, string(k), v)
    end
end

function write_attribute(h5::File, fbh::Header)
    write_attribute(h5["data"], fbh)
end

@deprecate write(h5ds::Dataset, fbh::Header) write_attribute(h5ds::Dataset, fbh::Header)
@deprecate write(h5ds::File, fbh::Header) write_attribute(h5ds::File, fbh::Header)

end # module HDF5FilterbankExt
