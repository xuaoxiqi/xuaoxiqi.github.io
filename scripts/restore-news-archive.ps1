[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9_-]+$')]
    [string]$Slug
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

Add-Type -AssemblyName System.IO.Compression.FileSystem

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$backupRoot = Join-Path $repoRoot 'news-backups'
$package = Join-Path $backupRoot ($Slug + '.zip')
$archivePath = Join-Path $repoRoot (Join-Path 'news-archive' ($Slug + '.html'))
if (-not (Test-Path -LiteralPath $package)) { throw "Cold archive not found: $package" }
if (Test-Path -LiteralPath $archivePath) {
    Write-Host "Archive is already restored: $archivePath"
    exit 0
}

$restoreRoot = Join-Path $backupRoot ('.restore-' + $Slug + '-' + [Guid]::NewGuid().ToString('N'))
$fullBackupRoot = [IO.Path]::GetFullPath($backupRoot).TrimEnd('\', '/') + [IO.Path]::DirectorySeparatorChar
$fullRestoreRoot = [IO.Path]::GetFullPath($restoreRoot)
if (-not $fullRestoreRoot.StartsWith($fullBackupRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Unsafe restore staging path: $fullRestoreRoot"
}

try {
    New-Item -ItemType Directory -Path $restoreRoot -Force | Out-Null
    $zip = [IO.Compression.ZipFile]::OpenRead($package)
    try {
        foreach ($entry in $zip.Entries) {
            $name = $entry.FullName.Replace('\', '/')
            if ($name.StartsWith('/') -or $name.Split('/') -contains '..') {
                throw "Unsafe path in archive package: $name"
            }
        }
    }
    finally { $zip.Dispose() }
    [IO.Compression.ZipFile]::ExtractToDirectory($package, $restoreRoot)
    $restoredArchive = Join-Path $restoreRoot (Join-Path 'news-archive' ($Slug + '.html'))
    if (-not (Test-Path -LiteralPath $restoredArchive)) {
        throw 'Restored package does not contain the expected archive page.'
    }
    Get-ChildItem -LiteralPath $restoreRoot -Force | ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $repoRoot -Recurse -Force
    }
    if (-not (Test-Path -LiteralPath $archivePath)) { throw 'Archive restore verification failed.' }
    Write-Host "Restored $Slug to $archivePath"
    Write-Host 'Run node scripts/archive-news.js to return restored files to cold storage.'
}
finally {
    if (Test-Path -LiteralPath $restoreRoot) { Remove-Item -LiteralPath $restoreRoot -Recurse -Force }
}
