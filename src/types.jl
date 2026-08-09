"""
    AbstractIgBLAST

Abstract supertype for IgBLAST program variants.
"""
abstract type AbstractIgBLAST end

"""
    IgBLASTn

Nucleotide IgBLAST (`igblastn`).
"""
struct IgBLASTn <: AbstractIgBLAST end

"""
    IgBLASTp

Protein IgBLAST (`igblastp`).
"""
struct IgBLASTp <: AbstractIgBLAST end

"""
    executable(::Type{<:AbstractIgBLAST})

Return the IgBLAST executable basename for the given variant.
"""
executable(::Type{IgBLASTn}) = "igblastn"
executable(::Type{IgBLASTp}) = "igblastp"

"""
    default_outfmt(::Type{<:AbstractIgBLAST})

Default `-outfmt` for the variant.
"""
default_outfmt(::Type{IgBLASTn}) = 19
default_outfmt(::Type{IgBLASTp}) = 7

# --- Molecule kind ---

"""
    AbstractMolecule

Molecule kind used when preparing BLAST databases.
"""
abstract type AbstractMolecule end

struct DNAMolecule <: AbstractMolecule end
struct ProteinMolecule <: AbstractMolecule end

"""
    molecule(::Type{<:AbstractIgBLAST})

Molecule kind required by germline databases for this variant.
"""
molecule(::Type{IgBLASTn}) = DNAMolecule()
molecule(::Type{IgBLASTp}) = ProteinMolecule()

# --- Auxiliary data (optional) ---

"""
    AbstractAuxiliary

Optional IgBLAST auxiliary (J-gene) annotation data.
Omit aux, or pass [`NoAuxiliary`](@ref) / [`noauxiliary`](@ref), unless a custom
[`AuxiliaryFile`](@ref) is required.
"""
abstract type AbstractAuxiliary end

"""
    NoAuxiliary

Sentinel meaning no auxiliary file is supplied.
"""
struct NoAuxiliary <: AbstractAuxiliary end

"""
    noauxiliary

Singleton [`NoAuxiliary`](@ref) instance.
"""
const noauxiliary = NoAuxiliary()

"""
    AuxiliaryFile{P<:AbstractString}

Path to a custom IgBLAST auxiliary data file.
"""
struct AuxiliaryFile{P<:AbstractString} <: AbstractAuxiliary
    path::P
end

"""
    supports_auxiliary(::Type{<:AbstractIgBLAST}) -> Bool

Whether this variant accepts `-auxiliary_data`.
"""
supports_auxiliary(::Type{<:AbstractIgBLAST}) = false
supports_auxiliary(::Type{IgBLASTn}) = true

"""
    normalize_auxiliary(aux) -> AbstractAuxiliary

Normalize `nothing`, empty string, path, or auxiliary values to [`AbstractAuxiliary`](@ref).
"""
normalize_auxiliary(::Nothing) = NoAuxiliary()
normalize_auxiliary(::NoAuxiliary) = NoAuxiliary()
normalize_auxiliary(aux::AuxiliaryFile) = aux

function normalize_auxiliary(path::AbstractString)
    isempty(path) && return NoAuxiliary()
    return AuxiliaryFile(path)
end

validate_inputs(::NoAuxiliary) = nothing

function validate_inputs(aux::AuxiliaryFile)
    isfile(aux.path) || throw(ArgumentError("Auxiliary file does not exist: $(aux.path)"))
    return nothing
end

# --- Germline databases ---

"""
    AbstractGermlines

Abstract germline FASTA collection used by an IgBLAST variant.
"""
abstract type AbstractGermlines end

"""
    VGermlines{P<:AbstractString}

V-only germline FASTA (used by [`IgBLASTp`](@ref)).
"""
struct VGermlines{P<:AbstractString} <: AbstractGermlines
    v::P
end

"""
    VDJGermlines{P<:AbstractString}

V, D, and J germline FASTA paths (used by [`IgBLASTn`](@ref)).
"""
struct VDJGermlines{P<:AbstractString} <: AbstractGermlines
    v::P
    d::P
    j::P
end

"""
    GermlineDatabases

Deprecated alias for [`VDJGermlines`](@ref).
"""
const GermlineDatabases = VDJGermlines

function validate_inputs(db::VGermlines)
    isfile(db.v) || throw(ArgumentError("V database file does not exist: $(db.v)"))
    return nothing
end

function validate_inputs(db::VDJGermlines)
    isfile(db.v) || throw(ArgumentError("V database file does not exist: $(db.v)"))
    isfile(db.d) || throw(ArgumentError("D database file does not exist: $(db.d)"))
    isfile(db.j) || throw(ArgumentError("J database file does not exist: $(db.j)"))
    return nothing
end

"""
    germlines_for(::Type{<:AbstractIgBLAST}, v, d, j)

Build the germline collection appropriate for the IgBLAST variant.
`IgBLASTp` keeps only V.
"""
germlines_for(::Type{IgBLASTn}, v::AbstractString, d::AbstractString, j::AbstractString) =
    VDJGermlines(v, d, j)

germlines_for(::Type{IgBLASTp}, v::AbstractString, ::AbstractString, ::AbstractString) =
    VGermlines(v)

germlines_for(::Type{IgBLASTp}, v::AbstractString) = VGermlines(v)

# --- Prepared BLAST databases (after makeblastdb) ---

"""
    AbstractPreparedDB

Prepared BLAST database paths produced for a specific IgBLAST variant.
"""
abstract type AbstractPreparedDB end

struct PreparedVDJ{P<:AbstractString} <: AbstractPreparedDB
    v::P
    d::P
    j::P
end

struct PreparedV{P<:AbstractString} <: AbstractPreparedDB
    v::P
end

# --- File encoding (query input / result output) ---

"""
    AbstractFileEncoding

Compression encoding inferred from a path (e.g. `.gz`).
"""
abstract type AbstractFileEncoding end

"""
    PlainEncoding

Uncompressed file encoding (default when the path does not end in `.gz`).
"""
struct PlainEncoding <: AbstractFileEncoding end

"""
    GzipEncoding

Gzip-compressed file encoding, used when a path ends in `.gz`.
"""
struct GzipEncoding <: AbstractFileEncoding end

"""
    file_encoding(path) -> AbstractFileEncoding

Autodetect encoding from the filename (`*.gz` → [`GzipEncoding`](@ref)).
"""
file_encoding(path::AbstractString) =
    endswith(lowercase(path), ".gz") ? GzipEncoding() : PlainEncoding()

# --- Output format (progress parsing) ---

abstract type AbstractOutputFormat end
struct AIRRFormat <: AbstractOutputFormat end
struct CommentedTabularFormat <: AbstractOutputFormat end

output_format(::Type{IgBLASTn}, outfmt::Integer) =
    Int(outfmt) == 19 ? AIRRFormat() : CommentedTabularFormat()

output_format(::Type{IgBLASTp}, ::Integer) = CommentedTabularFormat()
