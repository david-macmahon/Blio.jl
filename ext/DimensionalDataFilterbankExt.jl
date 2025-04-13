module DimensionalDataFilterbankExt

import DimensionalData: Dim, dims
using Blio.Filterbank: Header, chanfreqs

"""
    dims(fbh::Filterbank.Header, pols=:auto) -> (Dim{:Frequency}, Dim{:Pol}, Dim{:Time})

Return a `Tuple{Dim}` suitable for wrapping the data array described by `fbh` in
a `DimArray`.  `pols` can be `:stokes`, `:crosspol`, or `:auto`.  When `pol` is
`:stokes`, the polarization axis will be dimensioned using the first
`fbh[:nifs]` labels of `[:I, :Q, :U, :V]`.  When `pol` is `:crosspol` the
polarization axis will be dimensioned using the first `fbh[:nifs]` labels of
`[:XX, :YY, :ReXY, :ImXY]`.  When `pols` is `:auto`, the default, it will be
treated as `:stokes` if `fbh[:nifs]` is 1, otherwise it will be treated as
`:crosspol`.

For full control, `pols` can be passed an `AbstractVector{Symbol}` with a length
of at least `fbh[:nifs]`.
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
        pols = fbh[:nifs] == 1 ? :stokes : :crosspol
    end

    polsvec = pols == :stokes   ? [:I,  :Q,  :U,    :V   ] :
              pols == :crosspol ? [:XX, :YY, :ReXY, :ImXY] :
              error("invalid pols $pols")

    dims(fbh, polsvec)
end

end # module DimensionalDataFilterbankExt
