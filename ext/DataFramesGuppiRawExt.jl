module DataFramesGuppiRawExt

import Blio.GuppiRaw.Header
import Blio.GuppiRaw.load

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, push!
    using DataFrames: DataFrameRow, Tables, select!
else
    import ..DataFrames: DataFrame, push!
    using ..DataFrames: DataFrameRow, Tables, select!
end

# Outer-constructor to create Header from DataFrameRow
function Header(h::DataFrameRow)
    grh = Header()
    # If backend field exists, put it first
    haskey(h, :backend)  && (grh[:backend] = h[:backend])
    for (k,v) in zip(keys(h), values(h))
        k === :backend && continue
        v === nothing && continue
        ismissing(v) && continue
        grh[k] = v
    end
    grh
end

function DataFrame(v::Vector{Header})
    df = DataFrame(Tables.dictcolumntable(v))
    select!(df, sort(names(df)))
end

function push!(df::DataFrame, grh::Header)
    push!(df, getfield(grh, :dict), cols=:union)
end

"""
    load([predicate,] io::IO, ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])
    load([predicate,] rawname, ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])
    load([predicate,] rawnames, ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])

Equivalent to `load([predicate,] io; headers=DataFrame(), datablocks)` except
that the returned `headers` DataFrame will have column names in sorted order.
Same for the methods taking `rawname::AbstractString` and
`rawnames::AbstractVector{<:AbstractString}`.
"""
function load(predicate::Function, io::IO, ::Type{DataFrame};
              headers::DataFrame=DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
    load(predicate, io; headers, datablocks)
    select!(headers, sort(names(headers))), datablocks
end

function load(predicate::Function, rawname::AbstractString, ::Type{DataFrame};
              headers::DataFrame = DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
    open(rawname) do io
        load(predicate, io, DataFrame; headers, datablocks)
    end
end

function load(predicate::Function, rawnames::AbstractVector{<:AbstractString},
              ::Type{DataFrame};
              headers::DataFrame = DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
    for rawname in rawnames
        load(predicate, rawname, DataFrame; headers, datablocks)
    end
    (headers, datablocks)
end

function load(src, ::Type{DataFrame};
              headers::DataFrame=DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
    load(grh->true, src, DataFrame; headers, datablocks)
end

end # module DataFramesGuppiRawExt
