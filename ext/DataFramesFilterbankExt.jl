module DataFramesFilterbankExt

import Blio.Filterbank.Header

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, Tables, push!
else
    import ..DataFrames: DataFrame, Tables, push!
end

function DataFrame(v::Vector{Header})
    DataFrame(Tables.dictcolumntable(v))
end

function push!(df::DataFrame, fbh::Header)
    push!(df, getfield(fbh, :dict), cols=:union)
end

end # module DataFramesFilterbankExt
