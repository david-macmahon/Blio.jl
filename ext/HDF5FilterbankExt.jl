module HDF5FilterbankExt

import Blio.Filterbank: Header, fil2h5, mmap
import Base.write
import HDF5: write_attribute

if isdefined(Base, :get_extension)
    using HDF5: File, Dataset, h5open, attributes, create_dataset, datatype,
        write_dataset, Filters
else
    using ..HDF5: File, Dataset, h5open, attributes, create_dataset, datatype,
        write_dataset, Filters
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

    h5 = h5open(h5name, "w")
    # For blimpy compatability
    write_attribute(h5, "CLASS", "FILTERBANK")
    write_attribute(h5, "VERSION", "1.0")

    # Create `data` dataset
    size = (fbh[:nchans], fbh[:nifs], fbh[:nsamps])
    data = if copy
        fbd = open(io->mmap(io, fbh), fbname)
        ds = create_dataset(h5, "data", Float32, size; chunk, filters)
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

    close(h5)

    h5name
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
