```@meta
CurrentModule = IgBLAST
```

# IgBLAST

A Julia package for running [IgBLAST](https://github.com/mashu/IgBLAST.jl) (v1.22.0) on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

[`prepare`](@ref) builds BLAST databases, then runs your `do` block with an [`IgBLASTSession`](@ref). The temporary databases are deleted when the block returns — the same resource pattern as `mktempdir` / `open`.

## Installation

```julia
using Pkg
Pkg.add("IgBLAST")
```

## Quick start

Binaries install automatically on first `using IgBLAST`.

```julia
using IgBLAST

dbs = VDJGermlines("V.fasta", "D.fasta", "J.fasta")
params = Dict("organism" => "human", "domain_system" => "imgt")

prepare(IgBLASTn, dbs; additional_params=params) do ig
    ig("query.fasta", "output.tsv")
end

# Gzip query and/or output (autodetected)
prepare(IgBLASTn, dbs; additional_params=params) do ig
    ig("query.fasta.gz", "output.tsv.gz")
end

# Custom auxiliary file
prepare(IgBLASTn, dbs; aux="human_gl.aux", additional_params=params) do ig
    ig("query.fasta", "output.tsv")
end
```

### Timing the IgBLAST command

`makeblastdb` runs before the `do` block. For uncompressed files, [`command`](@ref) is NCBI IgBLAST writing a plain file:

```julia
seconds = prepare(IgBLASTn, dbs; additional_params=params, num_threads=8) do ig
    @elapsed run(command(ig, "query.fasta", "output.tsv"))
end
```

Do not use `.gz` paths here — gzip is Julia-side and would be included in the clock.

### IgBLASTp

Only the V germline is used:

```julia
prepare(IgBLASTp, VGermlines("V_nucleotide.fasta");
        additional_params=Dict("organism" => "human")) do ig
    ig("query_protein.fasta", "output.tsv")
end
```

## API

```@index
```

```@autodocs
Modules = [IgBLAST]
```
