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

## Nucleotide sequence assignment (IgBLASTn)

For nucleotide sequences, both query and database files should contain nucleotide sequences:

```julia
using IgBLAST

# Install IgBLAST if not already installed
install_igblast()

# Run an IgBLASTn analysis with nucleotide sequences
run_igblast(
    IgBLASTn,
    "query_nucleotide.fasta",  # Nucleotide query sequences
    "V_nucleotide.fasta",       # Nucleotide V gene database
    "D_nucleotide.fasta",       # Nucleotide D gene database
    "J_nucleotide.fasta",      # Nucleotide J gene database
    "auxiliary.txt",            # Auxiliary file (can be empty string "" if not needed)
    "output.txt",
    additional_params = Dict("organism" => "human", "domain_system" => "imgt")
)

# Run without auxiliary file (optional for assignments without CDR3 analysis)
run_igblast(
    IgBLASTn,
    "query_nucleotide.fasta",
    "V_nucleotide.fasta",
    "D_nucleotide.fasta",
    "J_nucleotide.fasta",
    "",  # Empty aux_file
    "output.txt",
    additional_params = Dict("organism" => "human", "domain_system" => "imgt")
)
```

## Protein sequence assignment (IgBLASTp)

For protein sequences, the query file must contain protein sequences (amino acids), while the database files should contain nucleotide sequences (which will be automatically translated to protein during database preparation):

```julia
# Run an IgBLASTp analysis with protein query sequences
run_igblast(
    IgBLASTp,
    "query_protein.fasta",      # Protein query sequences (must be amino acids, not nucleotides)
    "V_nucleotide.fasta",       # Nucleotide V gene database (will be translated to protein)
    "D_nucleotide.fasta",       # Nucleotide D gene database (will be translated to protein, but not used)
    "J_nucleotide.fasta",       # Nucleotide J gene database (will be translated to protein, but not used)
    "",                         # Auxiliary file not used by IgBLASTp
    "output.txt",
    additional_params = Dict("organism" => "human")
)
```

**Important notes:**
- For `IgBLASTn`: Both query and database files should contain nucleotide sequences. Returns V, D, and J assignments.
- For `IgBLASTp`: Query file must contain protein sequences (amino acids), database files should contain nucleotide sequences. **Only returns V assignments** (D and J databases are prepared but not used by IgBLASTp).
- The `aux_file` parameter can be an empty string `""` if not needed (useful for assignments without CDR3 analysis)
- `IgBLASTp` does not support the auxiliary file parameter
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