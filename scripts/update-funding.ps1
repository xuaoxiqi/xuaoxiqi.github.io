param(
    [switch]$DryRun,
    [switch]$Verbose
)

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
$ConfigPath = Join-Path $Root "data\funding-projects.json"
$PublicationsPath = Join-Path $Root "publications.html"
$ResearchPath = Join-Path $Root "research.html"

$Config = Get-Content -Raw -Encoding UTF8 $ConfigPath | ConvertFrom-Json
$RankTokens = @($Config.rankTokens)
$Projects = @($Config.projects)

function Get-DefaultLabelFromStem {
    param([string]$Stem)
    return (($Stem -replace "^\d+_", "") -replace "_arXiv$", "")
}

function Normalize-Href {
    param([string]$Href)
    return (($Href -replace "^\./", "").Replace("/", [IO.Path]::DirectorySeparatorChar))
}

function Invoke-PdfToText {
    param([string]$PdfPath)
    $Command = if ($Config.pdftotextCommand) { $Config.pdftotextCommand } else { "pdftotext" }
    $Text = & $Command -layout $PdfPath - 2>$null
    if ($LASTEXITCODE -ne 0) {
        throw "pdftotext failed for $PdfPath"
    }
    return ($Text -join "`n")
}

function Get-FundingText {
    param([string]$Text)
    $MarkerRegex = [regex]::new("\b(?:Funding\s*[\.:]|Acknowledgements?\.?|Acknowledgments?\.?|ACKNOWLEDGMENTS)\b", [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    $Markers = $MarkerRegex.Matches($Text)
    if ($Markers.Count -eq 0) {
        return ""
    }

    $Start = $Markers[$Markers.Count - 1].Index
    $Length = [Math]::Min(3500, $Text.Length - $Start)
    $Chunk = $Text.Substring($Start, $Length)
    $StopRegex = [regex]::new("\b(?:Declaration of interests|Conflict of interest|Author ORCIDs|Data availability|AUTHOR DECLARATIONS|DATA AVAILABILITY|APPENDIX|Appendix|REFERENCES|References|Supplementary materials)\b", [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    $Tail = if ($Chunk.Length -gt 1) { $Chunk.Substring(1) } else { "" }
    $Stop = $StopRegex.Match($Tail)
    if ($Stop.Success) {
        $Chunk = $Chunk.Substring(0, $Stop.Index + 1)
    }
    return $Chunk
}

function Find-TokenIndex {
    param(
        [string]$Text,
        $Token
    )
    $Best = $null
    foreach ($Pattern in @($Token.patterns)) {
        $Escaped = [regex]::Escape([string]$Pattern)
        $Regex = [regex]::new("(^|[^A-Za-z0-9])($Escaped)(?![A-Za-z0-9])", [Text.RegularExpressions.RegexOptions]::IgnoreCase)
        $Match = $Regex.Match($Text)
        if (-not $Match.Success) {
            continue
        }

        $Index = $Match.Index + $Match.Groups[1].Length
        if ($null -eq $Best -or $Index -lt $Best) {
            $Best = $Index
        }
    }
    return $Best
}

function Get-FundingOrder {
    param([string]$FundingText)
    $Hits = @()
    foreach ($Token in $RankTokens) {
        $Index = Find-TokenIndex -Text $FundingText -Token $Token
        if ($null -ne $Index) {
            $Hits += [pscustomobject]@{
                Id = [string]$Token.id
                Index = $Index
            }
        }
    }

    $Seen = @{}
    $Order = @()
    foreach ($Hit in ($Hits | Sort-Object Index)) {
        if ($Seen.ContainsKey($Hit.Id)) {
            continue
        }
        $Seen[$Hit.Id] = $true
        $Order += $Hit.Id
    }
    return @($Order)
}

function Add-Candidate {
    param(
        [System.Collections.Specialized.OrderedDictionary]$CandidatesByLabel,
        [string]$RelativePdfPath,
        [int]$ItemOrder
    )
    $Stem = [IO.Path]::GetFileNameWithoutExtension($RelativePdfPath)
    $Label = Get-DefaultLabelFromStem $Stem
    $Candidate = [pscustomobject]@{
        Label = $Label
        IsArxiv = ($Stem -match "_arXiv$")
        RelativePdfPath = $RelativePdfPath
        AbsolutePdfPath = Join-Path $Root $RelativePdfPath
        Order = $ItemOrder
    }

    if (-not $CandidatesByLabel.Contains($Label)) {
        $CandidatesByLabel[$Label] = @()
    }
    $CandidatesByLabel[$Label] = @($CandidatesByLabel[$Label]) + $Candidate
}

function Get-PublicationEntries {
    $Html = Get-Content -Raw -Encoding UTF8 $PublicationsPath
    $ItemRegex = [regex]::new('<!--\s*([\s\S]*?)\s*-->\s*<li class="pub-item[^"]*"[^>]*>([\s\S]*?)</li>', [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    $HrefRegex = [regex]::new('<a\b[^>]*href="([^"]+\.pdf)"', [Text.RegularExpressions.RegexOptions]::IgnoreCase)
    $CandidatesByLabel = [ordered]@{}
    $Order = 0

    foreach ($ItemMatch in $ItemRegex.Matches($Html)) {
        $Block = $ItemMatch.Groups[2].Value
        $Hrefs = @()
        foreach ($HrefMatch in $HrefRegex.Matches($Block)) {
            $Href = $HrefMatch.Groups[1].Value
            if (($Href -replace "^\./", "") -match "(^|/)publications/") {
                $Hrefs += $Href
            }
        }
        if ($Hrefs.Count -eq 0) {
            continue
        }

        $ItemOrder = $Order
        $Order++
        foreach ($Href in $Hrefs) {
            $RelativePdfPath = Normalize-Href $Href
            Add-Candidate -CandidatesByLabel $CandidatesByLabel -RelativePdfPath $RelativePdfPath -ItemOrder $ItemOrder

            if ($RelativePdfPath -match "_arXiv\.pdf$") {
                $FinalPdfPath = [regex]::Replace($RelativePdfPath, "_arXiv\.pdf$", ".pdf", [Text.RegularExpressions.RegexOptions]::IgnoreCase)
                if (Test-Path (Join-Path $Root $FinalPdfPath)) {
                    Add-Candidate -CandidatesByLabel $CandidatesByLabel -RelativePdfPath $FinalPdfPath -ItemOrder $ItemOrder
                }
            }
        }
    }

    $Entries = @()
    foreach ($Label in $CandidatesByLabel.Keys) {
        $Candidates = @($CandidatesByLabel[$Label]) | Sort-Object Order, IsArxiv
        $Chosen = @($Candidates | Where-Object { -not $_.IsArxiv } | Select-Object -First 1)
        if ($Chosen.Count -eq 0) {
            $Chosen = @($Candidates | Select-Object -First 1)
        }
        $Chosen = $Chosen[0]
        $Status = if ($Chosen.IsArxiv) { "coming" } else { "completed" }

        try {
            $PdfText = Invoke-PdfToText $Chosen.AbsolutePdfPath
            $FundingText = Get-FundingText $PdfText
            $FundingOrder = if ($FundingText) { Get-FundingOrder $FundingText } else { Get-FundingOrder $PdfText }
            if ($Verbose) {
                $TokenText = if ($FundingOrder.Count -gt 0) { $FundingOrder -join ", " } else { "no funding tokens found" }
                Write-Host "$Label [$Status, $($Chosen.RelativePdfPath)]: $TokenText"
            }
        }
        catch {
            Write-Warning "Failed to read $($Chosen.RelativePdfPath): $($_.Exception.Message)"
            $FundingOrder = @()
        }

        if ($FundingOrder.Count -gt 0) {
            $Entries += [pscustomobject]@{
                Label = $Label
                Status = $Status
                FundingOrder = @($FundingOrder)
                Source = $Chosen.RelativePdfPath
                Order = $Chosen.Order
            }
        }
    }

    return @($Entries)
}

function Get-ProjectRank {
    param(
        $Entry,
        $Project
    )
    foreach ($GrantId in @($Project.grantIds)) {
        for ($Index = 0; $Index -lt $Entry.FundingOrder.Count; $Index++) {
            if ($Entry.FundingOrder[$Index] -eq $GrantId) {
                return ($Index + 1)
            }
        }
    }
    return $null
}

function Sort-Entries {
    param($Entries)
    return @($Entries | Sort-Object Order)
}

function Get-ProjectResults {
    param($Entries)
    $Results = @{}
    foreach ($Project in $Projects) {
        $ProjectResults = @{}
        foreach ($Entry in $Entries) {
            $Rank = Get-ProjectRank -Entry $Entry -Project $Project
            if ($null -eq $Rank) {
                continue
            }

            if (-not $ProjectResults.ContainsKey($Entry.Status)) {
                $ProjectResults[$Entry.Status] = @()
            }
            $ProjectResults[$Entry.Status] = @($ProjectResults[$Entry.Status]) + [pscustomobject]@{
                Label = $Entry.Label
                Status = $Entry.Status
                FundingOrder = $Entry.FundingOrder
                Source = $Entry.Source
                Order = $Entry.Order
                Rank = $Rank
            }
        }

        foreach ($Status in @($ProjectResults.Keys)) {
            $ProjectResults[$Status] = Sort-Entries $ProjectResults[$Status]
        }
        $Results[[string]$Project.id] = $ProjectResults
    }
    return $Results
}

function Get-CountSummary {
    param($Entries)
    if ($Entries.Count -eq 0) {
        return "0"
    }

    $MaxRank = ($Entries | Measure-Object -Property Rank -Maximum).Maximum
    $Counts = @(for ($Index = 0; $Index -lt $MaxRank; $Index++) { 0 })
    foreach ($Entry in $Entries) {
        $Counts[$Entry.Rank - 1]++
    }
    return ($Counts -join "+")
}

function ConvertTo-Roman {
    param([int]$Value)
    $Map = @(
        @(1000, "m"), @(900, "cm"), @(500, "d"), @(400, "cd"),
        @(100, "c"), @(90, "xc"), @(50, "l"), @(40, "xl"),
        @(10, "x"), @(9, "ix"), @(5, "v"), @(4, "iv"), @(1, "i")
    )
    $Remaining = $Value
    $Result = ""
    foreach ($Pair in $Map) {
        while ($Remaining -ge $Pair[0]) {
            $Result += $Pair[1]
            $Remaining -= $Pair[0]
        }
    }
    return $Result
}

function Render-Block {
    param(
        [string]$ProjectId,
        [string]$Status,
        [string]$Indent,
        $Results,
        [string]$Eol
    )
    $Entries = @($Results[$ProjectId][$Status])
    $Label = if ($Config.statusLabels.$Status) { $Config.statusLabels.$Status } else { "($Status)" }
    $Lines = @("$Indent$Label $(Get-CountSummary $Entries) papers:$(if ($Entries.Count -gt 0) { '<br>' } else { '' })")
    if ($Entries.Count -gt 0) {
        $Total = $Entries.Count
        $Items = @()
        for ($Index = 0; $Index -lt $Entries.Count; $Index++) {
            $Entry = $Entries[$Index]
            $Items += "$(ConvertTo-Roman ($Total - $Index))) $($Entry.Label)_R$($Entry.Rank)"
        }
        $Lines += "$Indent$($Items -join '; ');"
    }
    return ($Lines -join $Eol)
}

function Update-Research {
    param($Results)
    $Html = Get-Content -Raw -Encoding UTF8 $ResearchPath
    $CrlfCount = ([regex]::Matches($Html, "`r`n")).Count
    $LfOnlyCount = ([regex]::Matches($Html, "(?<!`r)`n")).Count
    $Eol = if ($CrlfCount -ge $LfOnlyCount) { "`r`n" } else { "`n" }
    $Regex = [regex]::new('([ \t]*)<!-- funding:auto:start\s+project=([^\s]+)\s+status=([^\s]+)\s*-->\r?\n[\s\S]*?\r?\n[ \t]*<!-- funding:auto:end -->')

    $script:FundingReplacementCount = 0
    $Updated = $Regex.Replace($Html, {
        param($Match)
        $script:FundingReplacementCount++
        $Indent = $Match.Groups[1].Value
        $ProjectId = $Match.Groups[2].Value
        $Status = $Match.Groups[3].Value
        return @(
            "$Indent<!-- funding:auto:start project=$ProjectId status=$Status -->",
            (Render-Block -ProjectId $ProjectId -Status $Status -Indent $Indent -Results $Results -Eol $Eol),
            "$Indent<!-- funding:auto:end -->"
        ) -join $Eol
    })

    if ($script:FundingReplacementCount -eq 0) {
        throw "No funding:auto markers found in research.html."
    }

    if (-not $DryRun) {
        $Utf8NoBom = New-Object System.Text.UTF8Encoding($false)
        [IO.File]::WriteAllText($ResearchPath, $Updated, $Utf8NoBom)
    }
    return $script:FundingReplacementCount
}

function Write-Summary {
    param($Results)
    foreach ($Project in $Projects) {
        $ProjectId = [string]$Project.id
        foreach ($Status in @($Results[$ProjectId].Keys)) {
            $Entries = @($Results[$ProjectId][$Status])
            if ($Entries.Count -eq 0) {
                continue
            }
            $Detail = (($Entries | ForEach-Object { "$($_.Label)_R$($_.Rank)" }) -join ", ")
            Write-Host "$ProjectId ${Status}: $(Get-CountSummary $Entries) ($Detail)"
        }
    }
}

$Entries = Get-PublicationEntries
$Results = Get-ProjectResults $Entries
$ReplacementCount = Update-Research $Results
Write-Summary $Results
Write-Host "$(if ($DryRun) { 'Checked' } else { 'Updated' }) $ReplacementCount funding block(s) in research.html."
