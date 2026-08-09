```@meta
CurrentModule = IgBLAST
```

# IgBLAST

A Julia package for running [IgBLAST](https://github.com/mashu/IgBLAST.jl) (v1.22.0) on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

Auxiliary data is **optional**: omit it by default, or pass a custom [`AuxiliaryFile`](@ref) / path when needed.
Query FASTA and result files ending in `.gz` are handled automatically (decompress query / compress output).
Prefer [`IgBLASTRunner`](@ref) for repeated runs.

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

# Gzip query and/or output (autodetected)
run_igblast(
    IgBLASTn,
    "query.fasta.gz",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "output.tsv.gz",
)

# Custom auxiliary file
run_igblast(
    IgBLASTn,
    "query.fasta",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "human_gl.aux",
    "output.tsv",
)

# Typed germlines + runner
dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism" => "human"))
runner("query.fasta", dbs, "out.tsv")
```

### IgBLASTp

Only the V germline is used (D/J are not prepared):

```julia
run_igblast(
    IgBLASTp,
    "query_protein.fasta",
    "V_nucleotide.fasta",
    "output.tsv";
    additional_params = Dict("organism" => "human"),
)
```

## API

```@index
```

```@autodocs
Modules = [IgBLAST]
```
