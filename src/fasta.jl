"""
    count_fasta_sequences(filename::AbstractString)

Count FASTA records in `filename` by counting header lines.
"""
function count_fasta_sequences(filename::AbstractString)
    count = 0
    open(filename, "r") do file
        for line in eachline(file)
            startswith(line, '>') && (count += 1)
        end
    end
    return count
end

"""
    stage_query(input_file, output_file, encoding)

Copy or decompress a query FASTA into `output_file`, dispatching on encoding.
"""
function stage_query(
    input_file::AbstractString,
    output_file::AbstractString,
    ::PlainFASTA,
)
    cp(input_file, output_file)
    return output_file
end

function stage_query(
    input_file::AbstractString,
    output_file::AbstractString,
    ::GzipFASTA,
)
    open(GzipDecompressorStream, input_file) do io
        write(output_file, read(io))
    end
    return output_file
end

"""
    stage_query(input_file, output_file)

Infer encoding from the filename and stage the query.
"""
stage_query(input_file::AbstractString, output_file::AbstractString) =
    stage_query(input_file, output_file, query_encoding(input_file))

# Backward-compatible name
const copy_and_decompress = stage_query
