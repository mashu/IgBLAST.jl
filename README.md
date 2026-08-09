# IgBLAST.jl
![igblast-logo-svg](https://github.com/user-attachments/assets/b5ceac6b-49cc-40a0-aa0a-f7ce0a494b62)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBLAST.jl/dev/)
[![Build Status](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBLAST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBLAST.jl)

A Julia package for running IgBLAST (v1.22.0) analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

## Features

- Automatic installation and management of IgBLAST binaries
- Support for both IgBLASTn and IgBLASTp
- Optional auxiliary file via multiple dispatch (`omit` / `nothing` / `NoAuxiliary` / path)
- Callable `IgBLASTRunner` for repeated configured runs
- Progress monitoring for long-running analyses

## Installation

```julia
using Pkg
Pkg.add("IgBLAST")
```

## Quick Start

```julia
using IgBLAST

install_igblast()
```

### Nucleotide assignment (IgBLASTn)

Auxiliary data is **optional**. Prefer omitting it when CDR3 annotation is not needed:

```julia
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

# Explicit alternatives
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
```

### Callable runner

```julia
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner("query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv")

runner_aux = IgBLASTRunner(IgBLASTn; aux="human_gl.aux")
runner_aux("query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv")
```

### Protein assignment (IgBLASTp)

Query sequences must be amino acids; germline FASTA files are nucleotide and are translated when building BLAST DBs. Only V assignments are returned:

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

**Notes:**
- `IgBLASTn`: nucleotide query and databases; returns V, D, and J.
- `IgBLASTp`: protein query; germline FASTA should be nucleotide; auxiliary data is ignored.
- Empty string `""` for `aux` remains accepted for backward compatibility and means no auxiliary file.

For more detail, see the [documentation](https://mashu.github.io/IgBLAST.jl/dev/).
