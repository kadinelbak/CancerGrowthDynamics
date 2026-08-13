using Pkg
Pkg.activate(@__DIR__)

using JSON3

struct NotebookValidationDisplay <: AbstractDisplay end
Base.display(::NotebookValidationDisplay, ::MIME"image/png", ::Vector{UInt8}) = nothing
Base.display(::NotebookValidationDisplay, ::MIME"text/html", ::Any) = nothing
pushdisplay(NotebookValidationDisplay())

notebook_path = joinpath(@__DIR__, "notebooks", "a2780_staged_model_comparison.ipynb")
output_path = joinpath(@__DIR__, "outputs", "reports", "a2780_staged_model_comparison.executed.ipynb")
ENV["A2780_NOTEBOOK_REUSE_OUTPUTS"] = get(ENV, "A2780_NOTEBOOK_REUSE_OUTPUTS", "true")
ENV["A2780_MAX_TIME_PER_FIT"] = get(ENV, "A2780_MAX_TIME_PER_FIT", "1.0")

notebook = JSON3.read(lstrip(read(notebook_path, String), ['\ufeff']), Dict{String,Any})
execution_count = Ref(0)

function result_summary(result)
    result === nothing && return "Cell executed successfully."
    rendered = try
        sprint(show, MIME("text/plain"), result; context = :limit => true)
    catch
        string(typeof(result))
    end
    return first(rendered, min(length(rendered), 4000))
end

for cell in notebook["cells"]
    get(cell, "cell_type", "") == "code" || continue
    execution_count[] += 1
    cell["execution_count"] = execution_count[]
    code = join(get(cell, "source", String[]))
    try
        result = Base.include_string(Main, code, notebook_path)
        cell["outputs"] = Any[
            Dict(
                "name" => "stdout",
                "output_type" => "stream",
                "text" => result_summary(result) * "\n",
            ),
        ]
    catch e
        bt = catch_backtrace()
        cell["outputs"] = Any[
            Dict(
                "ename" => string(typeof(e)),
                "evalue" => sprint(showerror, e),
                "output_type" => "error",
                "traceback" => split(sprint(showerror, e, bt), '\n'),
            ),
        ]
        mkpath(dirname(output_path))
        open(output_path, "w") do io
            JSON3.pretty(io, notebook)
        end
        rethrow()
    end
end

mkpath(dirname(output_path))
open(output_path, "w") do io
    JSON3.pretty(io, notebook)
end
println("executed_notebook=", output_path)
