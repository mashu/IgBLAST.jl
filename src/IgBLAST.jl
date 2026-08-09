"""
    IgBLAST

A Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

The API is built around **multiple dispatch**: omit the auxiliary file, pass `nothing` /
`NoAuxiliary()`, or supply an `AuxiliaryFile` / path. Prefer `IgBLASTRunner` for repeated runs.

# Exports
- `install_igblast`, `is_igblast_installed`
- `run_igblast`, `IgBLASTRunner`
- `AbstractIgBLAST`, `IgBLASTn`, `IgBLASTp`
- `AbstractAuxiliary`, `NoAuxiliary`, `noauxiliary`, `AuxiliaryFile`
- `GermlineDatabases`

# Examples

```julia
using IgBLAST

install_igblast()

# No auxiliary file (CDR3 annotation not required)
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv";
            additional_params = Dict("organism" => "human", "domain_system" => "imgt"))

# Explicit nothing / NoAuxiliary
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", nothing, "out.tsv")
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", noauxiliary, "out.tsv")

# With auxiliary file
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", "human_gl.aux", "out.tsv")

# Callable runner
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner("query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv")
```
"""
module IgBLAST

using Artifacts
using CodecZlib
using TranscodingStreams
using Pkg.BinaryPlatforms
using ProgressMeter
import Pkg: ensure_artifact_installed
import Pkg.BinaryPlatforms: platform_key_abi
using BioSequences
using FASTX

export install_igblast, run_igblast, is_igblast_installed, IgBLASTRunner
export AbstractIgBLAST, IgBLASTn, IgBLASTp
export AbstractAuxiliary, NoAuxiliary, noauxiliary, AuxiliaryFile
export GermlineDatabases

const IGBLAST_VERSION = "1.22.0"

include("types.jl")
include("paths.jl")
include("fasta.jl")
include("database.jl")
include("command.jl")
include("progress.jl")
include("install.jl")
include("run.jl")

function __init__()
    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    if !isfile(artifact_toml) || !is_igblast_installed()
        @info "IgBLAST not found or not properly installed. Installing now..."
        install_igblast()
    else
        ensure_artifact_installed("IgBLAST", artifact_toml)
    end
end

end # module
