"""
    is_result_line(fmt, line, seen_header) -> (count_as_result, seen_header)

Dispatch progress-line classification on output format.
"""
function is_result_line(::CommentedTabularFormat, line::AbstractString, seen_header::Bool)
    startswith(line, '#') && return false, seen_header
    return true, seen_header
end

function is_result_line(::AIRRFormat, line::AbstractString, seen_header::Bool)
    startswith(line, '#') && return false, seen_header
    if !seen_header
        return false, true  # skip AIRR header row
    end
    return true, true
end

"""
    monitor_progress(output_file, total_sequences, done_channel, fmt)

Asynchronously update a progress bar while IgBLAST writes `output_file`.
"""
function monitor_progress(
    output_file::AbstractString,
    total_sequences::Integer,
    done_channel::Channel{Bool},
    fmt::AbstractOutputFormat,
)
    p = Progress(total_sequences; desc="Running IgBLAST... ", color=:green)
    return @async begin
        processed_sequences = 0
        last_position = 0
        seen_header = false
        while processed_sequences < total_sequences
            if isfile(output_file)
                open(output_file, "r") do f
                    seekend(f)
                    if position(f) > last_position
                        seek(f, last_position)
                        for line in eachline(f)
                            count_line, seen_header = is_result_line(fmt, line, seen_header)
                            if count_line
                                processed_sequences += 1
                                update!(p, processed_sequences)
                            end
                        end
                        last_position = position(f)
                    end
                end
            end
            sleep(0.1)
        end
        finish!(p)
        put!(done_channel, true)
    end
end
