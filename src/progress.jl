"""
    monitor_progress(output_file, total_sequences, done_channel)

Asynchronously update a progress bar while IgBLAST writes `output_file`.
Returns the monitoring `Task`.
"""
function monitor_progress(
    output_file::AbstractString,
    total_sequences::Integer,
    done_channel::Channel{Bool},
)
    p = Progress(total_sequences; desc="Running IgBLAST... ", color=:green)
    return @async begin
        processed_sequences = 0
        last_position = 0
        while processed_sequences < total_sequences
            if isfile(output_file)
                open(output_file, "r") do f
                    seekend(f)
                    if position(f) > last_position
                        seek(f, last_position)
                        for line in eachline(f)
                            if !startswith(line, '#')
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

"""
    run_process(cmd, output_file, total_sequences)

Run `cmd`, monitor progress against `output_file`, and throw on failure.
"""
function run_process(cmd::Cmd, output_file::AbstractString, total_sequences::Integer)
    done_channel = Channel{Bool}(1)
    monitor_task = monitor_progress(output_file, total_sequences, done_channel)

    out = IOBuffer()
    err = IOBuffer()
    process = run(pipeline(cmd, stdout=out, stderr=err); wait=false)

    while process_running(process)
        isready(done_channel) && break
        sleep(0.1)
    end

    wait(process)

    if !istaskdone(monitor_task)
        schedule(monitor_task, InterruptException(); error=true)
    end

    if process_running(process)
        kill(process)
    end

    success(process) || error("$(String(take!(err)))\nIgBLAST execution failed. Check the error message above.")
    return nothing
end
