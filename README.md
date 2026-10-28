## Thermal Soaring Animation

Place new animation files in the repository root using the name
`Thermal_Soaring_V<number>.html`, for example `Thermal_Soaring_V29.html`.
Use positive integer versions without leading zeros and keep these standalone
HTML files without YAML front matter. During each Jekyll build, `research.html`
automatically embeds the highest numeric version from `site.static_files`.
Older versions may be kept or removed; no Research page edit is needed.
Commit and publish the new file through the normal website workflow.

The local `_/preview-research.cjs` preview also discovers versions on each page
request, so refreshing the page picks up newly added files. Opening the raw
Research HTML directly with `file://` does not render Liquid templates.

## Funding Statistics

Update the project publication counts and lists in `research.html` from the local
PDFs linked in `publications.html`. Run these commands from the repository root.
Node.js and `pdftotext` must be available on `PATH`; `pdftotextCommand` in
`data/funding-projects.json` can also specify the executable's absolute path.

```powershell
# Check the results without changing research.html.
node scripts/update-funding.js --dry-run

# Show the funding IDs detected in each PDF.
node scripts/update-funding.js --dry-run --verbose

# Write the updated counts and publication lists.
node scripts/update-funding.js
```

The PowerShell entry point delegates to the same Node.js script and returns its
exit code. The existing parameters remain supported:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/update-funding.ps1 -DryRun
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/update-funding.ps1 -DryRun -Verbose
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/update-funding.ps1
```

Only the content between `funding:auto:start` and `funding:auto:end` comments is
regenerated. Project titles, funding amounts, targets and page layout are preserved.
An empty result is rendered as `0 papers`, without a placeholder publication.

- PDF filenames ending in `_arXiv.pdf` are classified as `coming`; other PDFs are
  classified as `completed`. When a matching final PDF exists beside an arXiv PDF,
  the final version takes precedence. Publication status is not checked online.
- Publication labels come from the PDF filename, without its numeric prefix or
  `_arXiv` suffix. Renaming a PDF changes its generated label on the next run.
- `R1`, `R2`, etc. indicate the order of recognized funding IDs in the PDF's funding
  or acknowledgements section, with a full-text fallback when no section is found.
  For example, `1+4+2` means one first-ranked, four second-ranked and two third-ranked
  acknowledgements. Funding IDs absent from `rankTokens` do not participate in ranking.
- Missing project configuration, missing ranking tokens, invalid markers or PDF
  extraction failures stop the update with a nonzero exit code before any write.

To add a project, add its recognized funding IDs to `rankTokens` and its project
entry to `projects` in `data/funding-projects.json`. Each project's `grantIds` must
refer to IDs in `rankTokens` (omitting `grantIds` uses the project ID). Add the
corresponding `completed` and/or `coming` blocks to its research-page entry:

```html
<!-- funding:auto:start project=GRANT_ID status=completed -->
(Completed) 0 papers:
<!-- funding:auto:end -->
<!-- funding:auto:start project=GRANT_ID status=coming -->
<br>(Coming..) 0 papers:
<!-- funding:auto:end -->
```

Run the focused regression checks with `node --test scripts/update-funding.test.js`.

## News Archives

News content is maintained in `_data/news.json` and rendered by `news.html`.
Each event has one story card with a `reports` array. Add reports there, then run
the command below. Each original report is downloaded and packaged into its own
verified ZIP under `news-backups`; only its lightweight cover remains outside
the package:

```powershell
node scripts/archive-news.js
```

Force-refresh every cold archive and its downloaded images:

```powershell
node scripts/archive-news.js --refresh
```

Cover selection follows the original webpage metadata (`og:image`), matching the cover chosen by WeChat. A substantial static body image is used only when the page does not provide cover metadata. For source pages with neither cover metadata nor article images, use a clearly editorial text cover rather than a site header or unrelated image. The HUST seminar announcements are examples.

Generate a labeled sheet for visual review with:

```powershell
powershell.exe -ExecutionPolicy Bypass -File scripts/build-news-cover-audit.ps1 -OutputPath news-cover-audit.png
```

Restore one story temporarily when its original link becomes unavailable:

```powershell
powershell.exe -ExecutionPolicy Bypass -File scripts/restore-news-archive.ps1 -Slug ARTICLE-SLUG
```

Running `node scripts/archive-news.js` again repackages any restored story. Package metadata and SHA-256 checksums are recorded in `news-backups/manifest.json`.

The `news-backups` directory stays in the website source repository but is excluded from the generated Jekyll site, so cold archives do not add weight to normal page delivery.

A story has `title`, `category`, `summary`, and one or more `reports`. Each report
retains its original `title`, `source`, `date`, ISO `datetime`, `url`, local
`image`, and unique `archive` slug. Optional `label` distinguishes reports from
the same publisher. `source_note` identifies a platform without counting it as
a separate source (for example, CSTAM's 澎湃号·政务 account).

The first report supplies the story's cover and main link. Each story is dated
by its earliest report, and stories are stored in reverse chronological order
using that date. Additional reports keep their own publication dates and links.
Stories about distinct events remain separate even when their research topics
overlap. Merging cards does not merge or delete source ZIPs.

The hero calculates Stories, Reports, Sources (unique source names across all
reports), and Archive span from this data. Archive and cover-audit scripts use
the same file. Report ZIP names and restore commands remain unchanged.

The Paper articles use the page's structured article body and configured `pic`
cover; its generic platform `og:image` icon is not an article cover. Run the
archive parser checks with `node --test scripts/archive-news.test.js`.
