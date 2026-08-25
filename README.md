# IgBLAST.jl
![igblast-logo-svg](https://github.com/user-attachments/assets/b5ceac6b-49cc-40a0-aa0a-f7ce0a494b62)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBLAST.jl/dev/)
[![Build Status](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBLAST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBLAST.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

A Julia package for running IgBLAST (v1.22.0) analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

## Features

- Automatic installation and management of IgBLAST binaries
- `prepare(...) do ig` prepares BLAST databases once; they are deleted when the block returns
- `command(ig, query, output)` is the raw NCBI command (plain files, no Julia gzip)
- Support for both IgBLASTn and IgBLASTp via multiple dispatch
- Optional auxiliary file
- Gzip autodetection for query FASTA (`.fasta.gz`) and result TSV (`.tsv.gz`) on the session call

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

Binaries install automatically on first `using IgBLAST`. Prepare databases once, then run:

### Nucleotide assignment (IgBLASTn)

```julia
using IgBLAST

dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
params = Dict("organism" => "human", "domain_system" => "imgt")

prepare(IgBLASTn, dbs; additional_params=params) do ig
    ig("query.fasta", "output.tsv")
end

# Gzip query and/or output (autodetected from `.gz`)
prepare(IgBLASTn, dbs; additional_params=params) do ig
    ig("query.fasta.gz", "output.tsv.gz")
end

# Custom auxiliary file
prepare(IgBLASTn, dbs; aux="human_gl.aux", additional_params=params) do ig
    ig("query.fasta", "output.tsv")
end
```

Path arguments are also accepted: `prepare(IgBLASTn, v, d, j; ...) do ig ... end`.

### Timing only the IgBLAST command

`makeblastdb` runs **before** the `do` block. Inside it, time NCBI IgBLAST writing a plain file — no Julia gzip, no extra copies, no chunking:

```julia
using IgBLAST

dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
params = Dict("organism" => "human", "domain_system" => "imgt")

seconds = prepare(IgBLASTn, dbs; additional_params=params, num_threads=8) do ig
    @elapsed run(command(ig, "query.fasta", "output.tsv"))
end

println("igblastn wall time: ", round(seconds; digits=3), " s")
```

Use **uncompressed** `.fasta` / `.tsv` for this. `command` is the raw `igblastn`/`igblastp` `Cmd`. For everyday use (including `.gz`), `ig("query.fasta", "output.tsv")` is enough.

Fair comparisons: one process, one query file, one output file; parallelize with `-num_threads`, not by splitting FASTA in Julia.

### Protein assignment (IgBLASTp)

Query sequences must be amino acids; the V germline FASTA should be nucleotide (translated when building the BLAST DB). Only V is used:

```julia
prepare(IgBLASTp, VGermlines("V_nucleotide.fasta");
        additional_params=Dict("organism" => "human")) do ig
    ig("query_protein.fasta", "output.tsv")
end
```

**Notes:**
- `IgBLASTn`: nucleotide query and V/D/J databases.
- `IgBLASTp`: protein query; only V germline is prepared/used.
- Paths ending in `.gz` are treated as gzip on `ig(query, output)`. `command` always expects uncompressed paths.
- NCBI IgBLAST does not read or write gzip; compression is only the Julia convenience path.

## License

The **MIT License** in [`LICENSE`](./LICENSE) covers only this Julia package — the wrapper, packaging, and Julia source in this repository. It does **not** apply to the NCBI IgBLAST binaries that the package downloads and runs.

NCBI IgBLAST itself is distributed separately by NCBI. With the exception of certain third-party files summarized by NCBI, that software is a “United States Government Work” under the terms of the United States Copyright Act. It was written as part of the authors’ official duties as United States Government employees and thus cannot be copyrighted. This software is freely available to the public for use. The National Library of Medicine and the U.S. Government have not placed any restriction on its use or reproduction.

See the [NCBI IgBLAST documentation](https://ncbi.github.io/igblast/) for upstream details and any third-party notices bundled with the tool.
