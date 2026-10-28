[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$OutputPath
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Drawing

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$stories = [IO.File]::ReadAllText((Join-Path $repoRoot '_data\news.json')) | ConvertFrom-Json

$items = @()
foreach ($story in $stories) {
    $primaryReport = $story.reports[0]
    $firstReport = $story.reports | Sort-Object datetime | Select-Object -First 1
    $items += [pscustomobject]@{
        Date = $firstReport.date
        Source = $primaryReport.source
        Title = $story.title
        Image = $primaryReport.image
    }
}

$columns = 3
$cellWidth = 470
$cellHeight = 320
$margin = 18
$rows = [Math]::Ceiling($items.Count / $columns)
$canvas = New-Object Drawing.Bitmap(($columns * $cellWidth), ([int]$rows * $cellHeight))
$graphics = [Drawing.Graphics]::FromImage($canvas)
$graphics.Clear([Drawing.Color]::FromArgb(242, 247, 247))
$graphics.InterpolationMode = [Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
$graphics.SmoothingMode = [Drawing.Drawing2D.SmoothingMode]::HighQuality
$labelFont = New-Object Drawing.Font('Microsoft YaHei', 10, [Drawing.FontStyle]::Regular)
$indexFont = New-Object Drawing.Font('Segoe UI', 10, [Drawing.FontStyle]::Bold)
$borderPen = New-Object Drawing.Pen([Drawing.Color]::FromArgb(195, 215, 218), 1)
$textBrush = New-Object Drawing.SolidBrush([Drawing.Color]::FromArgb(38, 70, 82))
$mutedBrush = New-Object Drawing.SolidBrush([Drawing.Color]::FromArgb(91, 112, 122))
$cardBrush = New-Object Drawing.SolidBrush([Drawing.Color]::White)

try {
    for ($index = 0; $index -lt $items.Count; $index += 1) {
        $item = $items[$index]
        $column = $index % $columns
        $row = [Math]::Floor($index / $columns)
        $x = $column * $cellWidth + $margin
        $y = $row * $cellHeight + $margin
        $card = New-Object Drawing.Rectangle($x, $y, ($cellWidth - 2 * $margin), ($cellHeight - 2 * $margin))
        $graphics.FillRectangle($cardBrush, $card)
        $graphics.DrawRectangle($borderPen, $card)

        $imagePath = Join-Path $repoRoot (($item.Image -replace '^\./', '') -replace '/', '\')
        $imageRect = New-Object Drawing.Rectangle(($x + 10), ($y + 10), ($card.Width - 20), 205)
        if ([IO.Path]::GetExtension($imagePath) -eq '.svg') {
            $graphics.DrawString('SVG text cover - preview in browser', $labelFont, $mutedBrush, $imageRect.X, ($imageRect.Y + 90))
        }
        elseif (Test-Path -LiteralPath $imagePath) {
            $image = [Drawing.Image]::FromFile($imagePath)
            try {
                $scale = [Math]::Min($imageRect.Width / $image.Width, $imageRect.Height / $image.Height)
                $drawWidth = [int][Math]::Round($image.Width * $scale)
                $drawHeight = [int][Math]::Round($image.Height * $scale)
                $drawX = $imageRect.X + [int](($imageRect.Width - $drawWidth) / 2)
                $drawY = $imageRect.Y + [int](($imageRect.Height - $drawHeight) / 2)
                $graphics.DrawImage($image, $drawX, $drawY, $drawWidth, $drawHeight)
            }
            finally { $image.Dispose() }
        }

        $graphics.DrawString(('{0:D2}' -f ($index + 1)), $indexFont, $mutedBrush, ($x + 12), ($y + 224))
        $label = $item.Source + ' · ' + $item.Date + [Environment]::NewLine + $item.Title
        $labelRect = New-Object Drawing.RectangleF(($x + 48), ($y + 222), ($card.Width - 60), 58)
        $graphics.DrawString($label, $labelFont, $textBrush, $labelRect)
    }
    $directory = [IO.Path]::GetDirectoryName([IO.Path]::GetFullPath($OutputPath))
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
    $canvas.Save($OutputPath, [Drawing.Imaging.ImageFormat]::Png)
}
finally {
    $cardBrush.Dispose()
    $mutedBrush.Dispose()
    $textBrush.Dispose()
    $borderPen.Dispose()
    $indexFont.Dispose()
    $labelFont.Dispose()
    $graphics.Dispose()
    $canvas.Dispose()
}

Write-Host "Cover audit sheet created: $OutputPath"
