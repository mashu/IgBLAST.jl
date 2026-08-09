# IgBLAST.jl
![igblast-logo-svg](https://github.com/user-attachments/assets/b5ceac6b-49cc-40a0-aa0a-f7ce0a494b62)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBLAST.jl/dev/)
[![Build Status](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBLAST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBLAST.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

A Julia package for running IgBLAST (v1.22.0) analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

## Features

- Automatic installation and management of IgBLAST binaries
- Support for both IgBLASTn and IgBLASTp via multiple dispatch
- Optional auxiliary file — omit by default; supply a custom one only when needed
- Gzip autodetection for query FASTA (`.fasta.gz`) and result TSV (`.tsv.gz`)
- Typed germlines: `VDJGermlines` / `VGermlines`
- Callable `IgBLASTRunner` for repeated configured runs
- Progress monitoring for long-running analyses

## Supported platforms

| OS | Architecture | Status |
|----|--------------|--------|
| Linux | x86_64 | Supported |
| macOS | x86_64 | Supported |
| macOS | ARM (Apple Silicon) | Not supported (no NCBI ARM binary; Rosetta 2 / x86_64 Julia may work) |
| Windows | x86_64 | Not supported yet |

IgBLAST binaries are installed automatically from NCBI via Julia artifacts (v1.22.0).

## Installation

```julia
using Pkg
Pkg.add("IgBLAST")
```

## Quick Start

Binaries install automatically on first `using IgBLAST`. Then run analyses:

### Nucleotide assignment (IgBLASTn)

Auxiliary data is optional. Prefer omitting it unless you need a custom aux file.
Gzip compression for query/output is autodetected from a `.gz` suffix:

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

# Compressed query and/or output
run_igblast(
    IgBLASTn,
    "query.fasta.gz",
    "V.fasta",
    "D.fasta",
    "J.fasta",
    "output.tsv.gz";
    additional_params = Dict("organism" => "human", "domain_system" => "imgt"),
)

# Custom auxiliary file
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

# Typed germlines
dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
run_igblast(IgBLASTn, "query.fasta", dbs, "output.tsv")
```

### Callable runner

```julia
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner("query.fasta", "V.fasta", "D.fasta", "J.fasta", "out.tsv")

runner_aux = IgBLASTRunner(IgBLASTn; aux="human_gl.aux")
runner_aux("query.fasta", VDJGermlines("V.fasta", "D.fasta", "J.fasta"), "out.tsv")
```

### Protein assignment (IgBLASTp)

Query sequences must be amino acids; the V germline FASTA should be nucleotide (translated when building the BLAST DB). Only V is used:

```julia
run_igblast(
    IgBLASTp,
    "query_protein.fasta",
    "V_nucleotide.fasta",
    "output.tsv";
    additional_params = Dict("organism" => "human"),
)
```

**Notes:**
- `IgBLASTn`: nucleotide query and V/D/J databases.
- `IgBLASTp`: protein query; only V germline is prepared/used.
- Paths ending in `.gz` are treated as gzip for both query input and result output.
- Empty string `""` for aux remains accepted for backward compatibility.

## License

The **MIT License** in [`LICENSE`](./LICENSE) covers only this Julia package — the wrapper, packaging, and Julia source in this repository. It does **not** apply to the NCBI IgBLAST binaries that the package downloads and runs.

NCBI IgBLAST itself is distributed separately by NCBI. With the exception of certain third-party files summarized by NCBI, that software is a “United States Government Work” under the terms of the United States Copyright Act. It was written as part of the authors’ official duties as United States Government employees and thus cannot be copyrighted. This software is freely available to the public for use. The National Library of Medicine and the U.S. Government have not placed any restriction on its use or reproduction.

See the [NCBI IgBLAST documentation](https://ncbi.github.io/igblast/) for upstream details and any third-party notices bundled with the tool.
