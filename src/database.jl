"""
    source_fasta(db_file, temp_dir, db_type, molecule)

FASTA path to feed `makeblastdb`, dispatching on molecule kind.
Nucleotide germlines are used in place; protein germlines are translated.
"""
source_fasta(db_file::AbstractString, ::AbstractString, ::AbstractString, ::DNAMolecule) =
    db_file

function source_fasta(
    db_file::AbstractString,
    temp_dir::AbstractString,
    db_type::AbstractString,
    ::ProteinMolecule,
)
    translated = joinpath(temp_dir, "$(db_type)_prot.fasta")
    open(FASTA.Writer, translated) do writer
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
    return translated
end

"""
    prepare_db(makeblastdb, db_file, db_type, temp_dir, molecule)

Build a BLAST database in `temp_dir`. Returns the database prefix.
"""
function prepare_db(
    makeblastdb::AbstractString,
    db_file::AbstractString,
    db_type::AbstractString,
    temp_dir::AbstractString,
    mol::AbstractMolecule,
)
    src = source_fasta(db_file, temp_dir, db_type, mol)
    prefix = joinpath(temp_dir, db_type)
    dt = blastdb_type(mol)
    run(`$makeblastdb -in $src -dbtype $dt -parse_seqids -out $prefix`)
    return prefix
end

"""
    prepare_databases(makeblastdb, germlines, temp_dir, molecule)

Prepare BLAST databases for a germline collection.
"""
function prepare_databases(
    makeblastdb::AbstractString,
    germlines::VDJGermlines,
    temp_dir::AbstractString,
    mol::DNAMolecule,
)
    return PreparedVDJ(
        prepare_db(makeblastdb, germlines.v, "V", temp_dir, mol),
        prepare_db(makeblastdb, germlines.d, "D", temp_dir, mol),
        prepare_db(makeblastdb, germlines.j, "J", temp_dir, mol),
    )
end

function prepare_databases(
    makeblastdb::AbstractString,
    germlines::VGermlines,
    temp_dir::AbstractString,
    mol::ProteinMolecule,
)
    return PreparedV(prepare_db(makeblastdb, germlines.v, "V", temp_dir, mol))
end
