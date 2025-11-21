# IgBLAST.jl
![igblast-logo-svg](https://github.com/user-attachments/assets/b5ceac6b-49cc-40a0-aa0a-f7ce0a494b62)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBLAST.jl/dev/)
[![Build Status](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBLAST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBLAST.jl)

A Julia package for running IgBLAST (v1.22.0) analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

## Features

- Automatic installation and management of IgBLAST binaries
- Support for both IgBLASTn and IgBLASTp
- Easy-to-use interface with customizable parameters
- Progress monitoring for long-running analyses

## Installation

```julia
using Pkg
Pkg.add("IgBLAST")
```

## Quick Start

```julia
using IgBLAST

# Install IgBLAST (if not already installed)
install_igblast()
```

### Nucleotide Sequence Assignment (IgBLASTn)

For nucleotide sequences, both query and database files should contain nucleotide sequences:

```julia
# Run an IgBLASTn analysis with nucleotide sequences
run_igblast(
    IgBLASTn,
    "data/ERR4238106.fasta.gz",      # Nucleotide query sequences
    "data/Macaca_mulatta_V.fasta",   # Nucleotide V gene database
    "data/Macaca_mulatta_D.fasta",   # Nucleotide D gene database
    "data/Macaca_mulatta_J.fasta",   # Nucleotide J gene database
    "data/rhesus_monkey_gl.aux",     # Auxiliary file (can be "" if not needed)
    "ERR4238106.tsv",
    additional_params = Dict("organism" => "rhesus_monkey", "ig_seqtype" => "Ig")
)

# Run without auxiliary file (optional for assignments without CDR3 analysis)
run_igblast(
    IgBLASTn,
    "query_nucleotide.fasta",
    "V_nucleotide.fasta",
    "D_nucleotide.fasta",
    "J_nucleotide.fasta",
    "",  # Empty aux_file
    "output.tsv",
    additional_params = Dict("organism" => "human", "domain_system" => "imgt")
)
```

### Protein Sequence Assignment (IgBLASTp)

For protein sequences, the query file must contain protein sequences (amino acids), while the database files should contain nucleotide sequences (which will be automatically translated to protein during database preparation):

```julia
# Run an IgBLASTp analysis with protein query sequences
run_igblast(
    IgBLASTp,
    "query_protein.fasta",      # Protein query sequences (must be amino acids, not nucleotides)
    "V_nucleotide.fasta",       # Nucleotide V gene database (will be translated to protein)
    "D_nucleotide.fasta",       # Nucleotide D gene database (will be translated to protein)
    "J_nucleotide.fasta",       # Nucleotide J gene database (will be translated to protein)
    "",                         # Auxiliary file not used by IgBLASTp
    "output.tsv",
    additional_params = Dict("organism" => "human")
)
```

**Important notes:**
- For `IgBLASTn`: Both query and database files should contain nucleotide sequences. Returns V, D, and J assignments.
- For `IgBLASTp`: Query file must contain protein sequences (amino acids), database files should contain nucleotide sequences. **Only returns V assignments** (D and J databases are prepared but not used by IgBLASTp).
- The `aux_file` parameter can be an empty string `""` if not needed (useful for assignments without CDR3 analysis)
- `IgBLASTp` does not support the auxiliary file parameter

For more detailed information, please refer to the [documentation](https://mashu.github.io/IgBLAST.jl/dev/).
