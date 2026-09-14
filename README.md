powershell.exe -ExecutionPolicy Bypass -File scripts/update-funding.ps1

Archive new stories added to `news.html`. Each story is downloaded, packaged into its own verified ZIP under `news-backups`, and removed from the live site except for its lightweight cover image:

```powershell
node scripts/archive-news.js
```

Force-refresh every cold archive and its downloaded images:

```powershell
node scripts/archive-news.js --refresh
```

Restore one story temporarily when its original link becomes unavailable:

```powershell
powershell.exe -ExecutionPolicy Bypass -File scripts/restore-news-archive.ps1 -Slug ARTICLE-SLUG
```

Running `node scripts/archive-news.js` again repackages any restored story. Package metadata and SHA-256 checksums are recorded in `news-backups/manifest.json`.

The `news-backups` directory stays in the website source repository but is excluded from the generated Jekyll site, so cold archives do not add weight to normal page delivery.
