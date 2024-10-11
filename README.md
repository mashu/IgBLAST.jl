# IgBLAST.jl
![igblast-logo-svg](https://github.com/user-attachments/assets/b5ceac6b-49cc-40a0-aa0a-f7ce0a494b62)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgBLAST.jl/dev/)
[![Build Status](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgBLAST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgBLAST.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgBLAST.jl)

A Julia package for running IgBLAST analyses on immunoglobulin (Ig) and T cell receptor (TCR) sequences.

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

For more detailed information, please refer to the [documentation](https://mashu.github.io/IgBLAST.jl/).
