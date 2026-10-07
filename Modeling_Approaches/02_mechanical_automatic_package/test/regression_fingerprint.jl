# SHA-256 of the files that define a run's configuration and data (test/make_reference.jl stores it,
# test/test_regression.jl requires it to match): a changed setting needs a new, deliberately accepted reference.
using SHA
const FINGERPRINT_FILES = ["stages.toml", "validation.toml", "control_config.toml", "data/sources.toml", "data/qc_rules.toml",
                           "data/schedules.toml", "data/measurements.csv"]
config_fingerprint() = Dict(f => bytes2hex(open(sha256, joinpath(normpath(joinpath(@__DIR__, "..")), f))) for f in FINGERPRINT_FILES)
