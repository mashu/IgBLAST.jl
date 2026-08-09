"""
    IgBLAST

A Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

Auxiliary data is **optional**: omit it, or pass a custom file only when needed.
Prefer [`IgBLASTRunner`](@ref) for repeated configured runs.

# Exports
- `install_igblast`, `is_igblast_installed`
- `run_igblast`, `IgBLASTRunner`
- `AbstractIgBLAST`, `IgBLASTn`, `IgBLASTp`
- `AbstractAuxiliary`, `NoAuxiliary`, `noauxiliary`, `AuxiliaryFile`
- `AbstractGermlines`, `VGermlines`, `VDJGermlines`, `GermlineDatabases`

# Examples

```julia
using IgBLAST

install_igblast()

# No auxiliary file
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv";
            additional_params = Dict("organism" => "human", "domain_system" => "imgt"))

# Custom auxiliary file when required
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", "human_gl.aux", "out.tsv")

# Typed germlines + runner
dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner("query.fasta", dbs, "out.tsv")
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
export AbstractGermlines, VGermlines, VDJGermlines, GermlineDatabases, germlines_for

const IGBLAST_VERSION = "1.22.0"

include("types.jl")
include("paths.jl")
include("fasta.jl")
include("database.jl")
include("command.jl")
include("progress.jl")
include("process.jl")
include("install.jl")
include("run.jl")

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
