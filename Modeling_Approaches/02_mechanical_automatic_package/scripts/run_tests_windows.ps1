$ErrorActionPreference = "Stop"

$packageRoot = Split-Path -Parent $PSScriptRoot
Push-Location $packageRoot

try {
    & julia +1.10.4 --compiled-modules=no --startup-file=no --project=. test/runtests.jl
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}
finally {
    Pop-Location
}
