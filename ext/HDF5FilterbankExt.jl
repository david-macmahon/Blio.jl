module HDF5FilterbankExt

import Blio.Filterbank: Header, fil2h5
import Base.write

if isdefined(Base, :get_extension)
    using HDF5: File, Dataset, h5open, attributes, create_dataset, write_attribute
else
    using ..HDF5: File, Dataset, h5open, attributes, create_dataset, write_attribute
end

function fil2h5(fbname, h5name="$fbname.h5"; kwargs...)
    ispath(h5name) && error("$h5name already exists")

    fbh = open(io->read(io, Header), fbname)
    h5 = h5open(h5name, "w")

    # For blimpy compatability
    fattrs = attributes(h5)
    fattrs["CLASS"] = "FILTERBANK"
    fattrs["VERSION"] = "1.0"

    size = (fbh[:nchans], fbh[:nifs], fbh[:nsamps])
    external = (fbname, fbh[:header_size], fbh[:data_size])
    data = create_dataset(h5, "data", Float32, size; external)

    attrs = attributes(data)
    for k in keys(fbh)
        attrs[string(k)] = fbh[k]
    end
    # Allow additional (or corrected) attributes from `kwargs`
    for k in keys(kwargs)
        attrs[string(k)] = kwargs[k]
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
    write(h5ds::HDF5.Dataset, fbh::Filterbank.Header)
    write(h5::HDF5.File, fbh::Filterbank.Header)

Write `fbh` contents as attributes of dataset `h5ds` or `h5["data"]`.

The `data_size` and `nsamps` fields will not be written since they can be
inferred from the dimensions of the dataset.
"""
function write(h5ds::Dataset, fbh::Header)
    for (k,v) in fbh
        # Don't copy these synthetic attributes
        k == :data_size && continue
        k == :nsamps && continue
        write_attribute(h5ds, string(k), v)
    end
end

function write(h5::File, fbh::Header)
    write(h5["data"], fbh)
end

end # module HDF5FilterbankExt
