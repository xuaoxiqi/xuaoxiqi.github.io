[CmdletBinding()]
param(
    [switch]$KeepExpanded,
    [switch]$SkipCoverOptimization
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$newsPath = Join-Path $repoRoot 'news.html'
$archiveRoot = Join-Path $repoRoot 'news-archive'
$imageRoot = Join-Path $repoRoot 'images\news'
$backupRoot = Join-Path $repoRoot 'news-backups'
$stagingRoot = Join-Path $backupRoot '.staging'
$manifestPath = Join-Path $backupRoot 'manifest.json'

function Assert-WithinRoot {
    param([string]$Path, [string]$Parent)
    $fullPath = [IO.Path]::GetFullPath($Path)
    $fullParent = [IO.Path]::GetFullPath($Parent).TrimEnd('\', '/') + [IO.Path]::DirectorySeparatorChar
    if (-not $fullPath.StartsWith($fullParent, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Unsafe path outside expected root: $fullPath"
    }
}

function Get-PropertyValue {
    param([string]$Block, [string]$Name)
    $pattern = '(?m)^    ' + [regex]::Escape($Name) + ': "([^"]*)"$'
    $match = [regex]::Match($Block, $pattern)
    if ($match.Success) { return $match.Groups[1].Value }
    return ''
}

function Get-Sha256 {
    param([string]$Path)
    $stream = [IO.File]::OpenRead($Path)
    $algorithm = [Security.Cryptography.SHA256]::Create()
    try {
        return ([BitConverter]::ToString($algorithm.ComputeHash($stream))).Replace('-', '').ToLowerInvariant()
    }
    finally {
        $algorithm.Dispose()
        $stream.Dispose()
    }
}

function Get-NewsItems {
    $source = [IO.File]::ReadAllText($newsPath)
    $frontMatter = [regex]::Match($source, '(?ms)^---\s*\r?\n(?<content>.*?)\r?\n---\s*(?:\r?\n|$)')
    if (-not $frontMatter.Success) { throw 'news.html front matter was not found.' }
    $blocks = [regex]::Matches($frontMatter.Groups['content'].Value, '(?ms)^  - date: "(?<date>[^"]*)"\r?\n(?<rest>.*?)(?=^  - date:|\z)')
    foreach ($match in $blocks) {
        $block = '  - date: "' + $match.Groups['date'].Value + '"' + [Environment]::NewLine + $match.Groups['rest'].Value
        $url = Get-PropertyValue $block 'url'
        if (-not $url) { continue }
        $slug = Get-PropertyValue $block 'archive'
        if (-not $slug) {
            $uri = [Uri]$url
            $slug = $uri.AbsolutePath.TrimEnd('/').Split('/')[-1]
        }
        [pscustomobject]@{
            Slug = $slug
            Title = Get-PropertyValue $block 'title'
            Source = Get-PropertyValue $block 'source'
            PublishedAt = Get-PropertyValue $block 'datetime'
            OriginalUrl = $url
            Cover = Get-PropertyValue $block 'image'
        }
    }
}

function Optimize-Cover {
    param([string]$Path)
    $extension = [IO.Path]::GetExtension($Path).ToLowerInvariant()
    if ($extension -notin @('.jpg', '.jpeg', '.png')) { return $false }

    Add-Type -AssemblyName System.Drawing
    $sourceImage = $null
    $bitmap = $null
    $graphics = $null
    $encoderParameters = $null
    $tempPath = Join-Path ([IO.Path]::GetDirectoryName($Path)) ('.' + [IO.Path]::GetFileNameWithoutExtension($Path) + '.optimized' + $extension)
    try {
        $sourceImage = [Drawing.Image]::FromFile($Path)
        $scale = [Math]::Min(1.0, [Math]::Min(960.0 / $sourceImage.Width, 960.0 / $sourceImage.Height))
        if ($scale -ge 1.0) { return $false }
        $width = [Math]::Max(1, [int][Math]::Round($sourceImage.Width * $scale))
        $height = [Math]::Max(1, [int][Math]::Round($sourceImage.Height * $scale))
        $bitmap = New-Object Drawing.Bitmap($width, $height)
        $graphics = [Drawing.Graphics]::FromImage($bitmap)
        $graphics.CompositingQuality = [Drawing.Drawing2D.CompositingQuality]::HighQuality
        $graphics.InterpolationMode = [Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
        $graphics.SmoothingMode = [Drawing.Drawing2D.SmoothingMode]::HighQuality
        if ($extension -in @('.jpg', '.jpeg')) { $graphics.Clear([Drawing.Color]::White) }
        else { $graphics.Clear([Drawing.Color]::Transparent) }
        $graphics.DrawImage($sourceImage, 0, 0, $width, $height)
        $graphics.Dispose()
        $graphics = $null
        $sourceImage.Dispose()
        $sourceImage = $null

        if ($extension -in @('.jpg', '.jpeg')) {
            $codec = [Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object MimeType -eq 'image/jpeg' | Select-Object -First 1
            $encoderParameters = New-Object Drawing.Imaging.EncoderParameters(1)
            $encoderParameters.Param[0] = [Drawing.Imaging.EncoderParameter]::new([Drawing.Imaging.Encoder]::Quality, [long]82)
            $bitmap.Save($tempPath, $codec, $encoderParameters)
        }
        else {
            $bitmap.Save($tempPath, [Drawing.Imaging.ImageFormat]::Png)
        }
        $bitmap.Dispose()
        $bitmap = $null
        if (-not (Test-Path -LiteralPath $tempPath) -or (Get-Item -LiteralPath $tempPath).Length -eq 0) {
            throw "Optimized cover is empty: $tempPath"
        }
        $check = [Drawing.Image]::FromFile($tempPath)
        $check.Dispose()
        [IO.File]::Delete($Path)
        [IO.File]::Move($tempPath, $Path)
        Write-Host "[cover] $([IO.Path]::GetFileName((Split-Path $Path -Parent))): ${width}x${height}"
        return $true
    }
    finally {
        if ($graphics) { $graphics.Dispose() }
        if ($bitmap) { $bitmap.Dispose() }
        if ($sourceImage) { $sourceImage.Dispose() }
        if ($encoderParameters) { $encoderParameters.Dispose() }
        if (Test-Path -LiteralPath $tempPath) { Remove-Item -LiteralPath $tempPath -Force }
    }
}

New-Item -ItemType Directory -Path $backupRoot -Force | Out-Null
New-Item -ItemType Directory -Path $stagingRoot -Force | Out-Null

$items = @(Get-NewsItems)
$itemBySlug = @{}
foreach ($item in $items) { $itemBySlug[$item.Slug] = $item }

$knownDates = @{}
if (Test-Path -LiteralPath $manifestPath) {
    $oldManifest = [IO.File]::ReadAllText($manifestPath) | ConvertFrom-Json
    foreach ($entry in @($oldManifest.archives)) { $knownDates[$entry.slug] = $entry.archivedAt }
}

$packagedDates = @{}
$packagedCount = 0
$optimizedCount = 0
$expanded = if (Test-Path -LiteralPath $archiveRoot) { @(Get-ChildItem -LiteralPath $archiveRoot -Filter '*.html' -File) } else { @() }
foreach ($archiveFile in $expanded) {
    $slug = $archiveFile.BaseName
    if (-not $itemBySlug.ContainsKey($slug)) {
        Write-Warning "Skipping unrecognized expanded archive: $($archiveFile.Name)"
        continue
    }
    $imageDirectory = Join-Path $imageRoot $slug
    $bodyFiles = @()
    if (Test-Path -LiteralPath $imageDirectory) {
        $bodyFiles = @(Get-ChildItem -LiteralPath $imageDirectory -File | Where-Object Name -NotMatch '^cover\.(?:avif|gif|jpe?g|png|svg|webp)$')
    }
    $archiveText = [IO.File]::ReadAllText($archiveFile.FullName)
    $dateMatch = [regex]::Match($archiveText, 'Archived (\d{4}-\d{2}-\d{2})\.')
    $packagedDates[$slug] = if ($dateMatch.Success) { $dateMatch.Groups[1].Value } else { [DateTime]::UtcNow.ToString('yyyy-MM-dd') }

    $stage = Join-Path $stagingRoot $slug
    Assert-WithinRoot $stage $stagingRoot
    if (Test-Path -LiteralPath $stage) { Remove-Item -LiteralPath $stage -Recurse -Force }
    $stageArchive = Join-Path $stage 'news-archive'
    $stageImages = Join-Path $stage (Join-Path 'images\news' $slug)
    New-Item -ItemType Directory -Path $stageArchive -Force | Out-Null
    Copy-Item -LiteralPath $archiveFile.FullName -Destination (Join-Path $stageArchive $archiveFile.Name)
    if ($bodyFiles.Count -gt 0) {
        New-Item -ItemType Directory -Path $stageImages -Force | Out-Null
        foreach ($file in $bodyFiles) { Copy-Item -LiteralPath $file.FullName -Destination (Join-Path $stageImages $file.Name) }
    }

    $destination = Join-Path $backupRoot ($slug + '.zip')
    $temporary = Join-Path $backupRoot ('.' + $slug + '.new.zip')
    if (Test-Path -LiteralPath $temporary) { Remove-Item -LiteralPath $temporary -Force }
    [IO.Compression.ZipFile]::CreateFromDirectory($stage, $temporary, [IO.Compression.CompressionLevel]::Optimal, $false)
    $zip = [IO.Compression.ZipFile]::OpenRead($temporary)
    try {
        $entries = @($zip.Entries | Where-Object { $_.Name } | ForEach-Object { $_.FullName.Replace('\', '/') })
        $expected = @('news-archive/' + $slug + '.html') + @($bodyFiles | ForEach-Object { 'images/news/' + $slug + '/' + $_.Name })
        foreach ($name in $expected) {
            if ($entries -notcontains $name) { throw "Package verification failed; missing $name" }
        }
    }
    finally { $zip.Dispose() }

    $previous = $destination + '.previous'
    if (Test-Path -LiteralPath $previous) { Remove-Item -LiteralPath $previous -Force }
    if (Test-Path -LiteralPath $destination) { Move-Item -LiteralPath $destination -Destination $previous }
    try {
        Move-Item -LiteralPath $temporary -Destination $destination
        if (Test-Path -LiteralPath $previous) { Remove-Item -LiteralPath $previous -Force }
    }
    catch {
        if ((Test-Path -LiteralPath $previous) -and -not (Test-Path -LiteralPath $destination)) {
            Move-Item -LiteralPath $previous -Destination $destination
        }
        throw
    }

    if (-not $KeepExpanded) {
        Remove-Item -LiteralPath $archiveFile.FullName -Force
        foreach ($file in $bodyFiles) { Remove-Item -LiteralPath $file.FullName -Force }
    }
    Remove-Item -LiteralPath $stage -Recurse -Force
    $packagedCount += 1
    Write-Host "[packed] $slug ($($expected.Count) files)"
}

if (-not $SkipCoverOptimization) {
    foreach ($item in $items) {
        $coverPath = $item.Cover -replace '^\./', ''
        if ($coverPath -notmatch '^images/news/') { continue }
        $absoluteCover = Join-Path $repoRoot ($coverPath -replace '/', '\')
        if ((Test-Path -LiteralPath $absoluteCover) -and (Optimize-Cover $absoluteCover)) { $optimizedCount += 1 }
    }
}

$manifestEntries = @()
foreach ($item in $items) {
    $package = Join-Path $backupRoot ($item.Slug + '.zip')
    if (-not (Test-Path -LiteralPath $package)) { continue }
    $zip = [IO.Compression.ZipFile]::OpenRead($package)
    try { $fileCount = @($zip.Entries | Where-Object { $_.Name }).Count }
    finally { $zip.Dispose() }
    $archivedAt = if ($packagedDates.ContainsKey($item.Slug)) { $packagedDates[$item.Slug] }
        elseif ($knownDates.ContainsKey($item.Slug)) { $knownDates[$item.Slug] }
        else { (Get-Item -LiteralPath $package).LastWriteTimeUtc.ToString('yyyy-MM-dd') }
    $manifestEntries += [ordered]@{
        slug = $item.Slug
        title = $item.Title
        source = $item.Source
        publishedAt = $item.PublishedAt
        originalUrl = $item.OriginalUrl
        archivedAt = $archivedAt
        package = './news-backups/' + $item.Slug + '.zip'
        cover = $item.Cover
        bytes = (Get-Item -LiteralPath $package).Length
        files = $fileCount
        sha256 = Get-Sha256 $package
    }
}

$manifest = [ordered]@{
    version = 1
    updatedAt = [DateTime]::UtcNow.ToString('yyyy-MM-ddTHH:mm:ssZ')
    archives = $manifestEntries
}
$json = $manifest | ConvertTo-Json -Depth 5
[IO.File]::WriteAllText($manifestPath, $json + [Environment]::NewLine, [Text.UTF8Encoding]::new($false))

if ((Test-Path -LiteralPath $stagingRoot) -and @(Get-ChildItem -LiteralPath $stagingRoot -Force).Count -eq 0) {
    Remove-Item -LiteralPath $stagingRoot -Force
}

Write-Host "Cold archive ready: $packagedCount packages created/refreshed, $($manifestEntries.Count) verified packages indexed, $optimizedCount covers resized."
