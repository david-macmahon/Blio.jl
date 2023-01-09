module DataFramesGuppiRawExt

import Blio.GuppiRaw.Header
import Blio.GuppiRaw.load

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, Tables, push!, select!
else
    import ..DataFrames: DataFrame, Tables, push!, select!
end

function DataFrame(v::Vector{Header})
    df = DataFrame(Tables.dictcolumntable(v))
    select!(df, sort(names(df)))
end

function push!(df::DataFrame, grh::Header)
    push!(df, getfield(grh, :dict), cols=:union)
end

"""
  load(io::IO ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])
  load(fn::AbstractString ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])

Equivalent to `load(io; headers=DataFrame(), datablocks)` except that the
returned `headers` DataFrame will have column names in sorted order.  Same for
the method taking `fn::AbstractString`.
"""
function load(io::IO, ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])
    headers = DataFrame()
    load(io; headers, datablocks)
    select!(headers, sort(names(headers))), datablocks
end

function load(rawname::AbstractString, ::Type{DataFrame};
              datablocks=Array{<:Complex{<:Integer}}[])
  open(rawname) do io
    load(io, DataFrame; datablocks)
  end
end

end # module DataFramesGuppiRawExt
