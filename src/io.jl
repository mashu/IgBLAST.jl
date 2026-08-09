"""
    count_fasta_sequences(filename::AbstractString)

Count FASTA records in `filename` by counting header lines.
Supports plain and gzip-compressed files (autodetected from the name).
"""
count_fasta_sequences(filename::AbstractString) =
    count_fasta_sequences(filename, file_encoding(filename))

function count_fasta_sequences(filename::AbstractString, ::PlainEncoding)
    count = 0
    open(filename, "r") do file
        for line in eachline(file)
            startswith(line, '>') && (count += 1)
        end
    end
    return count
end

function count_fasta_sequences(filename::AbstractString, ::GzipEncoding)
    count = 0
    open(GzipDecompressorStream, filename) do file
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
    ::PlainEncoding,
)
    cp(input_file, output_file)
    return output_file
end

function stage_query(
    input_file::AbstractString,
    output_file::AbstractString,
    ::GzipEncoding,
)
    open(GzipDecompressorStream, input_file) do input
        open(output_file, "w") do output
            write(output, input)
        end
    end
    return output_file
end

"""
    stage_query(input_file, output_file)

Infer encoding from the filename and stage the query for IgBLAST.
"""
stage_query(input_file::AbstractString, output_file::AbstractString) =
    stage_query(input_file, output_file, file_encoding(input_file))

"""
    igblast_output_path(user_output, temp_dir, encoding)

Path where IgBLAST should write. For gzip outputs this is a temporary plain file.
"""
igblast_output_path(user_output::AbstractString, ::AbstractString, ::PlainEncoding) =
    user_output

igblast_output_path(::AbstractString, temp_dir::AbstractString, ::GzipEncoding) =
    joinpath(temp_dir, "igblast_out")

"""
    finalize_output(igblast_path, user_output, encoding)

Copy or compress IgBLAST's plain output into the user-requested path.
"""
finalize_output(::AbstractString, ::AbstractString, ::PlainEncoding) = nothing

function finalize_output(
    igblast_path::AbstractString,
    user_output::AbstractString,
    ::GzipEncoding,
)
    open(igblast_path, "r") do input
        open(GzipCompressorStream, user_output, "w") do output
            write(output, input)
        end
    end
    return nothing
end
