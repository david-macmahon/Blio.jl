# Blio.jl

[![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://david-macmahon.github.io/Blio.jl/stable/)
[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://david-macmahon.github.io/Blio.jl/dev/)

**Breakthrough Listen I/O for Julia**

Blio.jl provides I/O support for radio astronomy data file formats used by the
[Breakthrough Listen](https://breakthroughinitiatives.org/initiative/3) program.

## Supported Formats

- **SIGPROC Filterbank** (`.fil`) — reduced-rate filterbank spectra
- **GUPPI Raw** (`.raw`) — raw data from several radio telescope backends

## Features

- Read and write SIGPROC Filterbank headers and spectra
- Read and write GUPPI Raw headers and data blocks
- Memory map data arrays for efficient access to large files
- Convert between GuppiRaw headers and Filterbank headers
- Convert between `.fil` and `.h5` formats (requires [HDF5.jl](https://github.com/JuliaIO/HDF5.jl))
- Optionally load headers into `DataFrame`s (requires [DataFrames.jl](https://github.com/JuliaDataFrames/DataFrames.jl))
- Optionally wrap data arrays in `DimArray`s (requires [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl))

## Installation

Blio.jl is not yet registered in the General registry. You can install it
directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/david-macmahon/Blio.jl.git")
```

## License

Blio.jl is licensed under the BSD 2-Clause License. See [LICENSE](LICENSE) for details.
