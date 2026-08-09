"""
    prepare_db(makeblastdb, db_file, db_type, temp_dir, molecule)

Prepare a BLAST database from a FASTA file, dispatching on molecule kind.
"""
function prepare_db(
    makeblastdb::AbstractString,
    db_file::AbstractString,
    db_type::AbstractString,
    temp_dir::AbstractString,
    ::DNAMolecule,
)
    temp_db = joinpath(temp_dir, "$(db_type)_db.fasta")
    cp(db_file, temp_db; force=true)
    run(`$makeblastdb -in $temp_db -dbtype nucl -parse_seqids`)
    return temp_db
end

function prepare_db(
    makeblastdb::AbstractString,
    db_file::AbstractString,
    db_type::AbstractString,
    temp_dir::AbstractString,
    ::ProteinMolecule,
)
    temp_db = joinpath(temp_dir, "$(db_type)_db.fasta")
    open(FASTA.Writer, temp_db) do writer
        open(FASTA.Reader, db_file) do reader
            for record in reader
                dna_seq = LongDNA{4}(sequence(record))
                trim_length = length(dna_seq) - (length(dna_seq) % 3)
                dna_seq = dna_seq[1:trim_length]
                protein_seq = translate(dna_seq)
                new_record = FASTA.Record(FASTA.identifier(record), string(protein_seq))
                write(writer, new_record)
            end
        end
    end
    run(`$makeblastdb -in $temp_db -dbtype prot -parse_seqids`)
    return temp_db
end

"""
    prepare_databases(::Type{IgBLASTn}, makeblastdb, germlines, temp_dir)

Prepare V, D, and J nucleotide BLAST databases.
"""
function prepare_databases(
    ::Type{IgBLASTn},
    makeblastdb::AbstractString,
    germlines::VDJGermlines,
    temp_dir::AbstractString,
)
    mol = molecule(IgBLASTn)
    return PreparedVDJ(
        prepare_db(makeblastdb, germlines.v, "V", temp_dir, mol),
        prepare_db(makeblastdb, germlines.d, "D", temp_dir, mol),
        prepare_db(makeblastdb, germlines.j, "J", temp_dir, mol),
    )
end

"""
    prepare_databases(::Type{IgBLASTp}, makeblastdb, germlines, temp_dir)

Prepare only the V protein BLAST database.
"""
function prepare_databases(
    ::Type{IgBLASTp},
    makeblastdb::AbstractString,
    germlines::VGermlines,
    temp_dir::AbstractString,
)
    mol = molecule(IgBLASTp)
    return PreparedV(prepare_db(makeblastdb, germlines.v, "V", temp_dir, mol))
end

prepare_databases(
    ::Type{IgBLASTp},
    makeblastdb::AbstractString,
    germlines::VDJGermlines,
    temp_dir::AbstractString,
) = prepare_databases(IgBLASTp, makeblastdb, VGermlines(germlines.v), temp_dir)

"""
    stage_auxiliary(aux, temp_dir)

Copy an auxiliary file into `temp_dir`, or return [`NoAuxiliary`](@ref).
"""
stage_auxiliary(::NoAuxiliary, ::AbstractString) = NoAuxiliary()

function stage_auxiliary(aux::AuxiliaryFile, temp_dir::AbstractString)
    temp_aux = joinpath(temp_dir, "aux_file")
    cp(aux.path, temp_aux)
    return AuxiliaryFile(temp_aux)
end
