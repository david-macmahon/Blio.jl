"""
Module for interacting with GuppiRaw files.

See also:
[`GuppiRaw.Header`](@ref),
[`GuppiRaw.HeaderDataUnit`](@ref),
[`read!(io::IO, grh::GuppiRaw.Header)`](@ref),
[`read(io::IO, ::Type{GuppiRaw.Header})`](@ref),
[`write(io::IO, grh::GuppiRaw.Header)`](@ref),
[`Array(grh::GuppiRaw.Header, nchan::Int=0)`](@ref)
[`chanfreq(grh::GuppiRaw.Header, chan::Real)`](@ref)
[`chanfreqs(grh::GuppiRaw.Header, chans::AbstractRange)`](@ref)
[`ntime(GuppiRaw.Header)`](@ref)
[`resize_hdu(hdu::GuppiRaw.HeaderDataUnit)`](@ref)
"""
module GuppiRaw

export Header

using Printf
using OrderedCollections
import Blio

# Header record size
const HEADER_REC_SIZE = 80

# Maximum number of records in a header
const HEADER_MAX_RECS = 256

# Number of columns for numeric header values
const HEADER_NUMERIC_COLS = 23

# Used to detect END record
const END = map(Int, collect("END "))

"""
Type used to hold a GuppiRaw header
"""
struct Header <: AbstractDict{Symbol, Any}
  # OrderedDict that holds the header key=value pairs
  dict::OrderedDict{Symbol, Any}

  # Buffer used to read in header bytes (and then some)
  buf::Array{UInt8,2}

  # Inner constructor to ensure that dict isn't set directly
  function Header(maxrecs::Integer=HEADER_MAX_RECS)::Header
    # Create an OrderedDict
    dict = OrderedDict()
    # Create an Array to buffer header data
    buf = Array{UInt8}(undef, HEADER_REC_SIZE, maxrecs);
    # Initialize fields and return
    new(dict, buf)
  end
end

function Base.setindex!(h::Header, val::Any, key::Symbol)
  setindex!(getfield(h, :dict), val, key)
end

function Base.setindex!(h::Header, val::Any, key::AbstractString)
  setindex!(h, val, Symbol(lowercase(key)))
end

function Base.getindex(h::Header, key::Symbol)
  getindex(getfield(h, :dict), key)
end

function Base.getindex(h::Header, key::AbstractString)
  getindex(h, Symbol(lowercase(key)))
end

function Base.get(h::Header, key::Symbol, default=nothing)
  get(getfield(h, :dict), key, default)
end

function Base.get(h::Header, key::AbstractString, default=nothing)
  get(h, Symbol(lowercase(key)), default)
end

function Base.getproperty(h::Header, sym::Symbol)
  get(h, sym, nothing)
end

function Base.length(h::Header)
  length(getfield(h, :dict))
end

function Base.propertynames(h::Header)
  Tuple(keys(getfield(h, :dict)))
end

function Base.iterate(h::Header, state...)
  iterate(getfield(h, :dict), state...)
end

function Base.empty!(h::Header)
  empty!(getfield(h, :dict))
end

function Base.delete!(h::Header, key::Symbol)
  delete!(getfield(h, :dict), key)
end

# For Julia < 1.9.0
if !isdefined(Base, :get_extension)
  # Add DataFrame constructor for Vector{GuppiRaw.Header} if/when DataFrames is
  # imported.
  using Requires
end
@static if !isdefined(Base, :get_extension)
  function __init__()
    @require DataFrames="a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
      include("../ext/DataFramesGuppiRawExt.jl")
    end
  end
end

"""
Returns the number of channels per antenna (i.e. `obsnchan ÷ nants`).
Missing `nants` implies `nants == 1`.
"""
function antnchan(grh::GuppiRaw.Header)::Int
  @assert haskey(grh, :obsnchan) "header has no obsnchan field"
  obsnants = get(grh, :nants, 1)
  @assert typeof(obsnants) <: Int

  obsnchan = grh[:obsnchan]
  @assert typeof(obsnchan) <: Int
  @assert obsnchan % obsnants == 0 "nants must divide obsnchan"

  # Return number of channels per antenna
  obsnchan ÷ obsnants
end

"""
    read!(io::IO, grh::GuppiRaw.Header; skip_padding::Bool=true)::Bool

Read a GUPPI header from `io` and populate `grh`.

Returns `false` if not enough data remain in `io` for a header, otherwise
`true` is returned.  This can be used to loop through a GuppiRaw file.

If `skip_padding` is true and `grh[:directio]` is non-zero, padding bytes will
be skipped after the header such that the file position will be a multiple of
512 bytes.  This will leave the file positioned at the start of the data block
following the header.  If reading a file of concatenated headers (i.e. without
a data block after each header), use `skip_padding=false`.
"""
function Base.read!(io::IO, grh::GuppiRaw.Header;
                    skip_padding::Bool=true
                   )::Bool
  # Get buf from Header
  buf = getfield(grh, :buf)
  @assert size(buf, 1) == HEADER_REC_SIZE

  # If not enough bytes remaining (EOF), return false
  if filesize(io) - position(io) < sizeof(buf)
    return false
  end

  # Read bytes into buf
  read!(io, buf)

  # Find END record
  endidx = findfirst(c->buf[1:4,c] == END, 1:size(buf,2))

  if isnothing(endidx)
    error("GUPPI RAW header not found in $(length(buf)) bytes")
  end

  # Make grh empty
  empty!(grh)

  # Parse first endidx-1 records
  for i in 1:endidx-1
    rec = String(buf[:,i])
    k, v = split(rec, '=', limit=2)..., missing
    # Skip malformed records
    ismissing(v) && continue
    k = Symbol(lowercase(strip(k)))
    v = strip(v)
    if v[1] == '\''
      v = strip(v, [' ', '\''])
    elseif !isnothing(match(r"^[+-]?[0-9]+$", v))
      v = parse(Int, v)
    elseif !isnothing(tryparse(Float64, v))
      v = parse(Float64, v)
    end
    grh[k] = v
  end

  # Seek io to just after END rec
  skip(io, HEADER_REC_SIZE*endidx - sizeof(buf))

  # If directio exists and is non-zero, seek past padding
  if skip_padding && get(grh, :directio, 0) != 0
    skip(io, mod(-position(io), 512))
  end

  return true
end

"""
    read(io::IO, ::Type{GuppiRaw.Header};
         skip_padding::Bool=true)::GuppiRaw.Header

Create a `GuppiRaw::Header` object, call `read!()` to populate it, then return
it.
"""
function Base.read(io::IO, ::Type{Header}; skip_padding::Bool=true)::Header
  grh = Header()
  read!(io, grh, skip_padding=skip_padding)
  return grh
end

# write_header_item is an internal function whose methods are used to write
# different types of values as a GUPPI raw record.
function write_header_item(io::IO, kw::Symbol, val::Integer)
  s = uppercase(rpad(kw,8)[1:8]) * "= " * lpad(string(val), HEADER_NUMERIC_COLS)
  write(io, rpad(s, 80)[1:80])
end

# Not sure the best/appropriate formatting to use here,
# but this seems like a reasonable approach to start with.
function write_header_item(io::IO, kw::Symbol, val::AbstractFloat)
  valstr = @sprintf("%.17G", val)
  # If no decimal point was output, redo to force one
  if isnothing(findfirst('.', valstr))
    if isnothing(findfirst('E', valstr))
      valstr = @sprintf("%.1f", val)
    else
      valstr = @sprintf("%.1E", val)
    end
  end
  s = uppercase(rpad(kw,8)[1:8]) *  "= " * lpad(valstr, HEADER_NUMERIC_COLS)
  write(io, rpad(s, 80)[1:80])
end

function write_header_item(io::IO, kw::Symbol, val::AbstractIrrational)
  write_header_item(io, kw, Float64(val))
end

function write_header_item(io::IO, kw::Symbol, val::Rational)
  write_header_item(io, kw, Float64(val))
end

function write_header_item(io::IO, kw::Symbol, val::AbstractString)
  s = uppercase(rpad(kw,8)[1:8]) * "= '" * rpad(val, 8) * "'"
  write(io, rpad(s, 80)[1:80])
end

"""
    write(io::IO, grh::GuppiRaw.Header; skip_padding::Bool=true)

Write `grh` as a GUPPI header to `io`.  If `skip_padding` is true and
`grh[:directio]` is non-zero, padding bytes will be written after the header
such that the file position will be a multiple of 512 bytes.
"""
function Base.write(io::IO, grh::GuppiRaw.Header;
                    skip_padding::Bool=true
                   )
  bytes_written = 0

  for (k,v) in grh
    bytes_written += write_header_item(io, k, v)
  end

  bytes_written += write(io, rpad("END", 80));

  # Handle DIRECTIO padding
  if skip_padding && get(grh, :directio, 0) != 0
    padding = mod(-bytes_written, 512)
    if padding != 0
      bytes_written += write(io, zeros(UInt8, padding))
    end
  end

  return bytes_written
end

"""
   size(grh::GuppiRaw.Header[, dim])

Return a tuple containing the dimensions of the data block described by `grh`.
Optionally you can specify a dimension to just get the length of that
dimension.

If `NANTS` is 1 or unspecified, the returned tuple have three elements
corresponding to `(npol, ntime, obsnchan)`.

If `NANTS` is greater than 1, the returned tuple will have four elements
corresponding to `(npol, ntime, obsnchan÷nants, nants)`.
"""
function Base.size(grh::Header)
  @assert haskey(grh, :obsnchan) "header has no obsnchan field"

  obsnants = get(grh, :nants, 1)
  obsnchan = grh[:obsnchan]
  obsntime = Blio.ntime(grh)
  npol = get(grh, :npol, 1) < 2 ? 1 : 2

  if obsnants > 1
    @assert obsnchan % obsnants == 0 "obsnchn $obsnchan not divisible by nants $obsnants"
    dims = (npol, obsntime, obsnchan÷obsnants, obsnants)
  else
    dims = (npol, obsntime, obsnchan)
  end

  dims
end

Base.size(grh::Header, dim) = size(grh)[dim]

"""
Type alias for possible GuppiRaw data Arrays
"""
const RawArray = Union{Array{Complex{Int8}},Array{Complex{Int16}}}

"""
A HeaderDataUnit is a struct consisting of a Header and a RawArray
"""
struct HeaderDataUnit
  hdr::Header
  data::RawArray
end

"""
    HeaderDataUnit(grh::GuppiRaw.Header=GuppiRaw.Header())

Construct a HeaderDataUnit from a `grh`.  If `grh` is empty, the `data` field
will be `Complex{Int8}[]`, otherwise it will be `Array(grh)`.  Note that `grh`
is not copied, so when making multiple HeaderDataUnits be sure to use a
different Header instance for each one.
"""
function HeaderDataUnit(grh::Header=Header())::HeaderDataUnit
  data = isempty(grh) ? Complex{Int8}[] : Array(grh)
  HeaderDataUnit(grh, data)
end

# Show HDUs compactly
Base.show(io::IO, hdu::HeaderDataUnit) = print(io,
  "GuppiRaw.HeaderDataUnit(hdr: ", length(hdu.hdr), " items, ",
                         "data: ", typeof(hdu.data), size(hdu.data), ")"
)

"""
    Array(grh::GuppiRaw.Header, nchan::Int=0)::Array{Complex{Integer}}

Return an uninitialized 3 or 4 dimensional Array sized for `nchan` channels of
RAW data as specified by metadata in `header`, specifically the `BLOCSIZE`,
`NANTS`, `OBSNCHAN`, `NPOL`, and `NBITS` fields.  `nchan <= 0` implies all
channels.

The data type of the Array elements will be `Complex{Int8}` when `NBITS == 8`
or `Complex{Int16}` when `NBITS == 16`.

If `NANTS` is 1 or unspecified, the Array will be dimensioned as [npol, ntime,
nchan] to match the RAW data block layout.

If `NANTS` is greater than 1 and `nchan` is a multiple of `OBSNCHAN/NANTS`,
then the returned array will be dimensioned as [npol, ntime, obsnchan/nants,
nchan*nants/obschan] (i.e. it will have an extra antenna dimension).
"""
function Base.Array(grh::Header, nchan::Int=0)::RawArray
  @assert haskey(grh, :obsnchan) "header has no obsnchan field"

  obsnants = get(grh, :nants, 1)
  obsnchan = grh[:obsnchan]
  # antnchan is number of channels per antenna
  antnchan = GuppiRaw.antnchan(grh)

  dims = size(grh)

  # If a specific number of channels are requested
  if nchan > 0
    # If nants > 1 (i.e. length of dims is 4)
    if length(dims) == 4
      antnchan = dims[3]
      obsnants = dims[4]
      if nchan > obsnants * antnchan
        nchan = obsnants * antnchan
        @warn "nchan limited to $nchan"
      end
      if nchan % antnchan != 0
        @warn "nchan $nchan not divisible by antnchan $antnchan"
        dims = (dims[1:2]..., nchan)
      else
        obsnants = nchan ÷ antnchan
        dims = (dims[1:3]..., obsnants)
      end
    else
      # nants==1 case
      obsnchan = dims[3]
      if nchan > obsnchan
        nchan = obsnchan
        @warn "nchan limited to $nchan"
      end
      dims = (dims[1:2]..., nchan)
    end
  end

  nbits = get(grh, :nbits, 8)
  @assert nbits == 8 || nbits == 16 "unsupported nbits ($nbits)"
  eltype = nbits == 8 ? Int8 : Int16

  Array{Complex{eltype}}(undef, dims)
end

end # module GuppiRaw

export chanfreq
export chanfreqs
export ntime
export resize_hdu

"""
    chanfreq(grh::GuppiRaw.Header, chan::Real)::Float64

Returns the center frequency of the channel given by `chan` based on the
`obsfreq`, `chan_bw`, `obsnchan`, and `nants` fields of `grh`.  The first
channel in the file is considered to be channel 1 (i.e. `chan` is one-based).
"""
function chanfreq(grh::GuppiRaw.Header, chan::Real)::Float64
  @assert haskey(grh, :obsfreq) "header has no obsfreq field"
  @assert haskey(grh, :chan_bw) "header has no chan_bw field"

  # antnchan is number of channels per antenna
  antnchan = GuppiRaw.antnchan(grh)

  # Subtract 1 and modulo antnchan
  chan = (chan-1) % antnchan

  grh[:obsfreq] - antnchan*grh[:chan_bw]/2 + grh[:chan_bw]*(chan + 0.5)
end

"""
    chanfreqs(grh::GuppiRaw.Header,
              chans::AbstractRange=1:grh[:obsnchan]/grh[:nants])::AbstractRange

Returns the center frequencies of the channels given by `chans` based on the
`obsfreq`, `chan_bw`, `obsnchan`, and `nants` fields of `grh`.  The first
channel in the file is considered to be channel 1 (i.e. `chans` are
one-based).  Frequencies for channels beyond `obsnchan/nants` will be
returned if requested, but they are not valid.
"""
function chanfreqs(grh::GuppiRaw.Header)::AbstractRange
  # antnchan is number of channels per antenna
  antnchan = GuppiRaw.antnchan(grh)

  range(chanfreq(grh, 1), step=grh[:chan_bw], length=antnchan)
end

function chanfreqs(grh::GuppiRaw.Header,
                        chans::AbstractRange)::AbstractRange
  range(chanfreq(grh, first(chans)),
        step=grh[:chan_bw]*step(chans),
        length=length(chans)
       )
end

"""
    ntime(grh::GuppiRaw.Header)::Integer

Return the number of time samples per block
"""
function ntime(grh::GuppiRaw.Header)::Integer
  @assert haskey(grh, :blocsize) "header has no blocsize field"
  @assert haskey(grh, :obsnchan) "header has no obsnchan field"

  npol = get(grh, :npol, 1) < 2 ? 1 : 2
  nbits = get(grh, :nbits, 8)
  @assert nbits == 8 || nbits == 16 "unsupported nbits ($nbits)"

  nt, rem = divrem(8 * grh[:blocsize], 2 * grh[:obsnchan] * npol * nbits)
  @assert rem == 0

  return nt
end

"""
    resize_hdu(hdu::GuppiRaw.HeaderDataUnit)::GuppiRaw.HeaderDataUnit

Return a `HeaderDataUnit` whose `data` Array is appropriately typed and sized for
the data block described by `hdu.hdr`.  If `hdu.data` is already appropriately
typed and sized then `hdu` is returned, otherwise a new `HeaderDataUnit` is
returned with the same `Header` object as `hdu.hdr` but a new `data` Array that
is appropriately typed and sized according to `hdu.hdr`.
"""
function resize_hdu(hdu::GuppiRaw.HeaderDataUnit)::GuppiRaw.HeaderDataUnit
  nbits = get(hdu.hdr, :nbits, 8)
  @assert nbits == 8 || nbits == 16 "unsupported nbits ($nbits)"
  reimtype = nbits == 8 ? Int8 : Int16

  if isempty(hdu.hdr)
    if isempty(hdu.data) && eltype(hdu.data) == Complex{Int8}
      return hdu
    else
      return GuppiRaw.HeaderDataUnit(hdu.hdr, Complex{Int8}[])
    end
  elseif (size(hdu.hdr) == size(hdu.data)
      &&  eltype(hdu.data) == Complex{reimtype})
    return hdu
  end

  GuppiRaw.HeaderDataUnit(hdu.hdr, Array(hdu.hdr))
end
