# Run a long bash job detached while asking Windows not to idle-sleep (a per-process request, released
# when this script exits; no power settings are changed).
#   powershell -File sensitivity/run_awake.ps1 -Script sensitivity/run_variants.sh [-Only "V4_ns V1_floor0"]
param([Parameter(Mandatory = $true)][string]$Script, [string]$Only = "")
Add-Type -Namespace Win32 -Name Power -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
$ES_CONTINUOUS = [uint32]2147483648; $ES_SYSTEM_REQUIRED = [uint32]1
[Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED) | Out-Null
$bash = 'C:\Program Files\Git\bin\bash.exe'
try {
    Set-Location (Split-Path $PSScriptRoot -Parent)
    if ($Only) { & $bash $Script ($Only -split ' ') } else { & $bash $Script }
} finally {
    [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS) | Out-Null
}
