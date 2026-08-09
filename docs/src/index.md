```@meta
CurrentModule = IgBLAST
```

# IgBLAST

A Julia package for running [IgBLAST](https://github.com/mashu/IgBLAST.jl) (v1.22.0) on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

The API uses **multiple dispatch**: omit the auxiliary file, pass `nothing` / [`NoAuxiliary`](@ref) / [`noauxiliary`](@ref), or supply a path / [`AuxiliaryFile`](@ref). Prefer [`IgBLASTRunner`](@ref) for repeated runs.

## Installation

```julia
using Pkg
Pkg.add("IgBLAST")
```

## Quick start

```julia
using IgBLAST

install_igblast()

# No auxiliary file
run_igblast(
    IgBLASTn,
    "query.fasta",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "output.tsv";
    additional_params = Dict("organism" => "human", "domain_system" => "imgt"),
)

# Explicit no-aux alternatives
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", nothing, "output.tsv")
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", noauxiliary, "output.tsv")

# With auxiliary file
run_igblast(
    IgBLASTn,
    "query.fasta",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "human_gl.aux",
    "output.tsv";
    additional_params = Dict("organism" => "human", "domain_system" => "imgt"),
)

# Callable runner
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism" => "human"))
runner("query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv")
```

### IgBLASTp

Query sequences must be amino acids; germline FASTA files should be nucleotide (translated when building BLAST DBs). Only V assignments are returned; auxiliary data is ignored.

```julia
run_igblast(
    IgBLASTp,
    "query_protein.fasta",
    "V_nucleotide.fasta",
    "D_nucleotide.fasta",
    "J_nucleotide.fasta",
    "output.tsv";
    additional_params = Dict("organism" => "human"),
)
```

Empty string `""` for the auxiliary argument remains accepted for backward compatibility and means no auxiliary file.

## API

```@index
```

```@autodocs
Modules = [IgBLAST]
```
