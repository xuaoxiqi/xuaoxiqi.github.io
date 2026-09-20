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

Archive new stories added to `news.html`. Each story is downloaded, packaged into its own verified ZIP under `news-backups`, and removed from the live site except for its lightweight cover image:

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
