"""
    IgBLAST

A Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

Prepare BLAST databases once with [`prepare`](@ref), then run IgBLAST inside
the `do` block. Temporary databases are deleted when the block returns.

# Exports
- `install_igblast`, `is_igblast_installed`
- `prepare`, `command`, `IgBLASTSession`
- `AbstractIgBLAST`, `IgBLASTn`, `IgBLASTp`
- `AbstractAuxiliary`, `NoAuxiliary`, `noauxiliary`, `AuxiliaryFile`
- `AbstractGermlines`, `VGermlines`, `VDJGermlines`, `germlines_for`

# Examples

```julia
using IgBLAST

dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
params = Dict("organism" => "human", "domain_system" => "imgt")

prepare(IgBLASTn, dbs; additional_params=params) do ig
    ig("query.fasta", "out.tsv")
end

# Time only igblastn (plain files; databases already prepared)
t = prepare(IgBLASTn, dbs; additional_params=params) do ig
    @elapsed run(command(ig, "query.fasta", "out.tsv"))
end
```
"""
module IgBLAST

using Artifacts
using CodecZlib
import Pkg: ensure_artifact_installed
using BioSequences
using FASTX

export install_igblast, is_igblast_installed
export prepare, command, IgBLASTSession
export AbstractIgBLAST, IgBLASTn, IgBLASTp
export AbstractAuxiliary, NoAuxiliary, noauxiliary, AuxiliaryFile
export AbstractGermlines, VGermlines, VDJGermlines, germlines_for

const IGBLAST_VERSION = "1.22.0"

include("types.jl")
include("paths.jl")
include("io.jl")
include("database.jl")
include("command.jl")
include("session.jl")
include("install.jl")

function __init__()
    artifact_toml = artifact_toml_path()
    if !isfile(artifact_toml) || !is_igblast_installed()
        @info "IgBLAST not found or not properly installed. Installing now..."
        install_igblast()
    else
        ensure_artifact_installed("IgBLAST", artifact_toml)
    end
end

end # module
