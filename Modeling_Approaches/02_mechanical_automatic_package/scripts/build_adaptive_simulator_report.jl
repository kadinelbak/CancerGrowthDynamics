using MechanicalAutomaticModeling
using JSON3

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(PACKAGE_ROOT, "..", ".."))
const REPORT_NAME = "a2780_adaptive_treatment_simulations.html"
const CONFIG_JSON = "a2780_adaptive_simulator_config.json"
const CONFIG_JS = "a2780_adaptive_simulator_config.js"

function write_config(directory, config)
    mkpath(directory)
    open(joinpath(directory, CONFIG_JSON), "w") do io
        JSON3.pretty(io, config)
        println(io)
    end
    open(joinpath(directory, CONFIG_JS), "w") do io
        print(io, "window.A2780_ADAPTIVE_CONFIG = ")
        JSON3.write(io, config)
        println(io, ";")
    end
end

function build_report()
    report_directory = joinpath(PACKAGE_ROOT, "outputs", "reports")
    docs_directory = joinpath(REPOSITORY_ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports")
    template = joinpath(report_directory, REPORT_NAME)
    isfile(template) || error("report template is missing: $template")

    config = A2780AdaptiveAdapter.load_a2780_adaptive_config(PACKAGE_ROOT)
    write_config(report_directory, config)
    write_config(docs_directory, config)
    cp(template, joinpath(docs_directory, REPORT_NAME); force = true)

    println("Default Stage 3 model: ", config.winners.stage3.model)
    println("Default Stage 4 model: ", config.winners.stage4.model)
    println("Wrote report and Julia-produced configuration to:")
    println("  ", report_directory)
    println("  ", docs_directory)
    return config
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && build_report()
