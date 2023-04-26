module HDF5FilterbankExt

import Blio.Filterbank: Header, fil2h5

if isdefined(Base, :get_extension)
    import HDF5: h5open, attributes, create_dataset
else
    import ..HDF5: h5open, attributes, create_dataset
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

end # module HDF5FilterbankExt
