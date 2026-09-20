param(
    [switch]$DryRun,
    [switch]$Verbose
)

$NodeCommand = Get-Command node -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
if (-not $NodeCommand) {
    [Console]::Error.WriteLine("Node.js was not found on PATH. Install Node.js to run update-funding.")
    exit 1
}

$NodeArgs = @((Join-Path $PSScriptRoot "update-funding.js"))
if ($DryRun) {
    $NodeArgs += "--dry-run"
}
if ($Verbose) {
    $NodeArgs += "--verbose"
}

& $NodeCommand.Source @NodeArgs
exit $LASTEXITCODE
