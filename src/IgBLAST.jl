"""
    IgBLAST

A Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

This package provides a convenient interface to install and run IgBLAST, supporting both IgBLASTn and IgBLASTp variants.
It handles the installation of IgBLAST binaries, prepares input files, runs analyses, and monitors progress.

# Exports
- `install_igblast`: Function to install IgBLAST binaries
- `run_igblast`: Main function to run IgBLAST analyses
- `is_igblast_installed`: Function to check if IgBLAST is installed
- `IgBLASTn`: Type representing the nucleotide version of IgBLAST
- `IgBLASTp`: Type representing the protein version of IgBLAST

# Examples
```julia
using IgBLAST

# Install IgBLAST if not already installed
install_igblast()

# Run an IgBLASTn analysis
run_igblast(
    IgBLASTn,
    "query.fasta",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "auxiliary.txt",
    "output.txt",
    additional_params = Dict("organism" => "human", "domain_system" => "imgt")
)
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

export install_igblast, run_igblast, is_igblast_installed
export AbstractIgBLAST, IgBLASTn, IgBLASTp

const IGBLAST_VERSION = "1.22.0"

"""
    AbstractIgBLAST

An abstract type representing IgBLAST variants.
"""
abstract type AbstractIgBLAST end

"""
    IgBLASTn

A concrete type representing the nucleotide version of IgBLAST.
"""
struct IgBLASTn <: AbstractIgBLAST end

"""
    IgBLASTp

A concrete type representing the protein version of IgBLAST.
"""
struct IgBLASTp <: AbstractIgBLAST end

"""
    executable(::Type{T}) where T <: AbstractIgBLAST

Returns the executable name for a given IgBLAST variant.
"""
executable(::Type{IgBLASTn}) = "igblastn"
executable(::Type{IgBLASTp}) = "igblastp"

include("utils.jl")
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