"""
    ensure_parent_directory(path)

Create the parent directory of `path` when it is non-empty and missing.
"""
function ensure_parent_directory(path::AbstractString)
    dir = dirname(path)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    return nothing
end

"""
    stage_query(input_file, output_file, ::GzipEncoding)

Decompress a gzip FASTA into `output_file` for IgBLAST.
"""
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
    finalize_output(igblast_path, user_output, ::GzipEncoding)

Compress IgBLAST's plain output into the user-requested gzip path.
"""
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
