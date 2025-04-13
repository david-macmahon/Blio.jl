module DimensionalDataFilterbankExt

import DimensionalData: Dim, dims
using Blio.Filterbank: Header, chanfreqs

"""
    dims(fbh::Filterbank.Header, pols=:auto) -> (Dim{:Frequency}, Dim{:Pol}, Dim{:Time})

Return a `Tuple{Dim{:Frequency}, Dim{:Pol}, Dim{:Time}}` suitable for wrapping
the data array described by `fbh` in a `DimArray`.  The frequency axis is fully
determined by the `:fch1`, `:foff`, and `:nchans` fields of `fbh`.  The time
axis contains the time of each sample relative to the first sample.  It is fully
determined by the `:tsamp` and `:nsamps` fields of `fbh`.  For standard
Filterbank headers the frequency and time axes will be in MHz and seconds, resp.

The polarization axis, `Dim{:Pol}`, is a categorical axis.  Filterbank headers
do not specify which polarization products are present; only the number of
polarization products is given in the `nifs` field of `fbh`.  The names of the
polarization products may be passed as the `pols` parameter, which can be either
a `Vector{Symbol}` to give the names explicitly or one of several Symbols that
specify pre-defined lists of common polarization product names (see table).  If
`:auto`, the default, is passed it will be treated as `:stokes` when
`fbh[:nifs]` is 1, otherwise it will be treated as `:xy`.

| pols::Symbol | Equivalent pols::Vector  |
|--------------|--------------------------|
| :auto        | (see text)               |
| :stokes      | [:I, :Q, :U, :V]         |
| :xy          | [:XX, :YY, :ReXY, :ImXY] |
| :yx          | [:YY, :XX, :ReYX, :ImYX] |
| :lr          | [:LL, :RR, :ReLR, :ImLR] |
| :rl          | [:RR, :LL, :ReRL, :ImRL] |
| :hv          | [:HH, :VV, :ReHV, :ImHV] |
| :vh          | [:VV, :VH, :ReVH, :ImVH] |
"""
function dims(fbh::Header, pols::AbstractVector{Symbol})
    nifs = fbh[:nifs]
    nifs > length(pols) && error("nifs $nifs > $(length(pols))")

    (
        Dim{:Frequency}(chanfreqs(fbh)),
        Dim{:Pol}(first(pols, nifs)),
        Dim{:Time}(range(0, step=fbh[:tsamp], length=fbh[:nsamps]))
    )
end

function dims(fbh::Header, pols::Symbol=:auto)
    if pols == :auto
        pols = fbh[:nifs] == 1 ? :stokes : :xy
    end

    polsvec = pols == :stokes   ? [:I,  :Q,  :U,    :V   ] :
              pols == :xy       ? [:XX, :YY, :ReXY, :ImXY] :
              pols == :lr       ? [:LL, :RR, :ReLR, :ImLR] :
              pols == :hv       ? [:HH, :VV, :ReHV, :ImHV] :
              pols == :yx       ? [:YY, :XX, :ReYX, :ImYX] :
              pols == :rl       ? [:RR, :LL, :ReRL, :ImRL] :
              pols == :vh       ? [:VV, :HH, :ReVH, :ImVH] :
              error("invalid pols $pols")

    dims(fbh, polsvec)
end

end # module DimensionalDataFilterbankExt
