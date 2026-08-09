"""
    run_process(cmd, output_file, total_sequences, fmt)

Run `cmd`, monitor progress against `output_file`, and throw on failure.
"""
function run_process(
    cmd::Cmd,
    output_file::AbstractString,
    total_sequences::Integer,
    fmt::AbstractOutputFormat,
)
    done_channel = Channel{Bool}(1)
    monitor_task = monitor_progress(output_file, total_sequences, done_channel, fmt)

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

    process_running(process) && kill(process)

    success(process) || error("$(String(take!(err)))\nIgBLAST execution failed. Check the error message above.")
    return nothing
end
