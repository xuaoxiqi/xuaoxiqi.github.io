const fs = require('fs');
const path = require('path');
const childProcess = require('child_process');

const root = path.resolve(__dirname, '..');
const newsPath = path.join(root, 'news.html');
const archiveRoot = path.join(root, 'news-archive');
const imageRoot = path.join(root, 'images', 'news');
const backupRoot = path.join(root, 'news-backups');
const userAgent = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/140.0 Safari/537.36';

function decodeEntities(value) {
  return value
    .replace(/&amp;/g, '&')
    .replace(/&quot;/g, '"')
    .replace(/&#39;|&apos;/g, "'")
    .replace(/&lt;/g, '<')
    .replace(/&gt;/g, '>');
}

function escapeHtml(value) {
  return String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function parseQuoted(value) {
  const trimmed = value.trim();
  return trimmed.startsWith('"') && trimmed.endsWith('"') ? JSON.parse(trimmed) : trimmed;
}

function parseNews(source) {
  const frontMatter = source.split(/^---\s*$/m)[1];
  const items = [];
  let current = null;
  for (const line of frontMatter.split(/\r?\n/)) {
    const start = line.match(/^  - ([a-z_]+):\s*(.*)$/);
    const property = line.match(/^    ([a-z_]+):\s*(.*)$/);
    if (start) {
      current = {};
      current[start[1]] = parseQuoted(start[2]);
      items.push(current);
    } else if (property && current) {
      current[property[1]] = parseQuoted(property[2]);
    }
  }
  return items;
}

function articleSlug(url) {
  const parsed = new URL(url);
  const token = parsed.pathname.split('/').filter(Boolean).pop();
  if (!token) {
    throw new Error('Cannot derive a safe archive name from ' + url);
  }
  if (parsed.hostname === 'mp.weixin.qq.com' && /^[A-Za-z0-9_-]+$/.test(token)) {
    return token;
  }
  const stem = token.replace(/\.[^.]+$/, '').replace(/[^A-Za-z0-9_-]+/g, '-');
  const host = parsed.hostname.replace(/^www\./, '').replace(/[^A-Za-z0-9]+/g, '-');
  if (!stem || !host) throw new Error('Cannot derive a safe archive name from ' + url);
  return host + '-' + stem;
}

async function fetchResponse(url, retries) {
  let lastError;
  const retryCount = retries === undefined ? 2 : retries;
  for (let attempt = 0; attempt <= retryCount; attempt += 1) {
    const controller = new AbortController();
    const timeout = setTimeout(function () { controller.abort(); }, 45000);
    try {
      const response = await fetch(url, {
        headers: {
          'User-Agent': userAgent,
          Referer: new URL(url).origin + '/'
        },
        redirect: 'follow',
        signal: controller.signal
      });
      if (!response.ok) throw new Error('HTTP ' + response.status);
      return response;
    } catch (error) {
      lastError = error;
      if (attempt < retryCount) {
        await new Promise(function (resolve) {
          setTimeout(resolve, 700 * (attempt + 1));
        });
      }
    } finally {
      clearTimeout(timeout);
    }
  }
  throw lastError;
}

async function fetchText(url) {
  return (await fetchResponse(url)).text();
}

function metaImage(html) {
  const match = html.match(/<meta\s+property=["']og:image["']\s+content=["']([^"']*)["']/i);
  return match ? decodeEntities(match[1]) : '';
}

function extractDivContent(html, markerPattern) {
  const marker = html.search(markerPattern);
  if (marker < 0) return '';
  const start = html.lastIndexOf('<div', marker);
  if (start < 0) return '';
  const tagPattern = /<div\b[^>]*>|<\/div\s*>/gi;
  tagPattern.lastIndex = start;
  let depth = 0;
  let openingEnd = -1;
  let match;
  while ((match = tagPattern.exec(html))) {
    if (/^<div\b/i.test(match[0])) {
      depth += 1;
      if (openingEnd < 0) openingEnd = tagPattern.lastIndex;
    } else {
      depth -= 1;
      if (depth === 0) return html.slice(openingEnd, match.index);
    }
  }
  return '';
}

function extractMainContent(html) {
  const selectors = [
    /id=["']js_content["']/i,
    /id=["']vsb_content_2["']/i,
    /class=["'][^"']*\bv_news_content\b[^"']*["']/i
  ];
  for (const selector of selectors) {
    const content = extractDivContent(html, selector);
    if (content) return content;
  }
  throw new Error('A supported article body container was not found');
}

function attributeValue(tag, name) {
  const expression = new RegExp('\\s' + name + '\\s*=\\s*(?:"([^"]*)"|\\x27([^\\x27]*)\\x27|([^\\s>]+))', 'i');
  const match = tag.match(expression);
  return match ? decodeEntities(match[1] || match[2] || match[3] || '') : '';
}

function normalizedRemoteUrl(value, baseUrl) {
  if (!value || value.startsWith('data:')) return '';
  let candidate = value.trim();
  if (candidate.startsWith('//')) candidate = 'https:' + candidate;
  if (candidate.startsWith('http://')) candidate = 'https://' + candidate.slice(7);
  try {
    const parsed = new URL(candidate, baseUrl);
    return parsed.protocol === 'https:' ? parsed.href : '';
  } catch (error) {
    return '';
  }
}

function extensionFor(url, contentType) {
  const type = (contentType || '').toLowerCase();
  if (type.includes('png')) return 'png';
  if (type.includes('gif')) return 'gif';
  if (type.includes('webp')) return 'webp';
  if (type.includes('svg')) return 'svg';
  if (type.includes('jpeg') || type.includes('jpg')) return 'jpg';
  try {
    const parsed = new URL(url);
    const wxFormat = parsed.searchParams.get('wx_fmt');
    if (wxFormat && /^(png|gif|webp|jpeg|jpg)$/i.test(wxFormat)) {
      return wxFormat.toLowerCase().replace('jpeg', 'jpg');
    }
  } catch (error) {}
  return 'jpg';
}

function imageDimensions(bytes, extension) {
  try {
    if (extension === 'png' && bytes.length >= 24 && bytes.toString('ascii', 1, 4) === 'PNG') {
      return { width: bytes.readUInt32BE(16), height: bytes.readUInt32BE(20) };
    }
    if (extension === 'gif' && bytes.length >= 10 && bytes.toString('ascii', 0, 3) === 'GIF') {
      return { width: bytes.readUInt16LE(6), height: bytes.readUInt16LE(8) };
    }
    if (extension === 'jpg' && bytes.length >= 12 && bytes[0] === 0xff && bytes[1] === 0xd8) {
      let offset = 2;
      while (offset + 9 < bytes.length) {
        if (bytes[offset] !== 0xff) { offset += 1; continue; }
        const marker = bytes[offset + 1];
        if ([0xc0, 0xc1, 0xc2, 0xc3, 0xc5, 0xc6, 0xc7, 0xc9, 0xca, 0xcb, 0xcd, 0xce, 0xcf].includes(marker)) {
          return { width: bytes.readUInt16BE(offset + 7), height: bytes.readUInt16BE(offset + 5) };
        }
        if (marker === 0xd8 || marker === 0xd9) { offset += 2; continue; }
        const segmentLength = bytes.readUInt16BE(offset + 2);
        if (segmentLength < 2) break;
        offset += segmentLength + 2;
      }
    }
    if (extension === 'webp' && bytes.length >= 30 && bytes.toString('ascii', 0, 4) === 'RIFF' && bytes.toString('ascii', 8, 12) === 'WEBP') {
      const format = bytes.toString('ascii', 12, 16);
      if (format === 'VP8X') {
        return {
          width: 1 + bytes.readUIntLE(24, 3),
          height: 1 + bytes.readUIntLE(27, 3)
        };
      }
      if (format === 'VP8 ' && bytes.length >= 30) {
        return { width: bytes.readUInt16LE(26) & 0x3fff, height: bytes.readUInt16LE(28) & 0x3fff };
      }
      if (format === 'VP8L' && bytes.length >= 25) {
        const b1 = bytes[21];
        const b2 = bytes[22];
        const b3 = bytes[23];
        const b4 = bytes[24];
        return {
          width: 1 + (((b2 & 0x3f) << 8) | b1),
          height: 1 + (((b4 & 0x0f) << 10) | (b3 << 2) | ((b2 & 0xc0) >> 6))
        };
      }
    }
  } catch (error) {}
  return { width: 0, height: 0 };
}

async function downloadImage(url, directory, stem) {
  const response = await fetchResponse(url);
  const bytes = Buffer.from(await response.arrayBuffer());
  const extension = extensionFor(url, response.headers.get('content-type'));
  const dimensions = imageDimensions(bytes, extension);
  const fileName = stem + '.' + extension;
  fs.mkdirSync(directory, { recursive: true });
  const finalPath = path.join(directory, fileName);
  const temporaryPath = path.join(directory, '.' + fileName + '.download-' + process.pid + '-' + Date.now());
  fs.writeFileSync(temporaryPath, bytes);
  let saved = false;
  try {
    for (let attempt = 0; attempt < 6; attempt += 1) {
      try {
        if (fs.existsSync(finalPath)) fs.unlinkSync(finalPath);
        fs.renameSync(temporaryPath, finalPath);
        saved = true;
        break;
      } catch (error) {
        if (attempt === 5) throw error;
        await new Promise(function (resolve) { setTimeout(resolve, 250 * (attempt + 1)); });
      }
    }
  } finally {
    if (!saved && fs.existsSync(temporaryPath)) fs.unlinkSync(temporaryPath);
  }
  return { fileName: fileName, extension: extension, bytes: bytes.length, width: dimensions.width, height: dimensions.height };
}

async function mapWithConcurrency(values, limit, worker) {
  const results = new Array(values.length);
  let next = 0;
  async function run() {
    while (next < values.length) {
      const index = next;
      next += 1;
      results[index] = await worker(values[index], index);
    }
  }
  await Promise.all(Array.from({ length: Math.min(limit, values.length) }, run));
  return results;
}

async function localizeImages(content, slug, directory, pageUrl, preferredCoverIndex) {
  const tags = content.match(/<img\b[^>]*>/gi) || [];
  const urls = [];
  for (const tag of tags) {
    const source = normalizedRemoteUrl(attributeValue(tag, 'data-src') || attributeValue(tag, 'src'), pageUrl);
    if (source && !urls.includes(source)) urls.push(source);
  }
  const localMap = new Map();
  const downloads = new Array(urls.length);
  let downloadedBytes = 0;
  let failedImages = 0;
  await mapWithConcurrency(urls, 5, async function (url, index) {
    try {
      const saved = await downloadImage(url, directory, String(index + 1).padStart(3, '0'));
      localMap.set(url, './images/news/' + slug + '/' + saved.fileName);
      downloads[index] = {
        url: url,
        extension: saved.extension,
        width: saved.width,
        height: saved.height,
        bytes: saved.bytes
      };
      downloadedBytes += saved.bytes;
    } catch (error) {
      failedImages += 1;
      process.stderr.write('  image skipped (' + error.message + '): ' + url + '\n');
    }
  });
  const rewritten = content.replace(/<img\b[^>]*>/gi, function (tag) {
    const source = normalizedRemoteUrl(attributeValue(tag, 'data-src') || attributeValue(tag, 'src'), pageUrl);
    const local = localMap.get(source);
    if (!local) return '<span class="news-archive-image-missing">[Image unavailable in archive]</span>';
    const alt = attributeValue(tag, 'alt') || 'Archived article image';
    return '<img src="' + escapeHtml(local) + '" alt="' + escapeHtml(alt) + '" loading="lazy" decoding="async">';
  });
  const preferredCover = /^\d+$/.test(String(preferredCoverIndex || '')) ? downloads[Number(preferredCoverIndex) - 1] : null;
  const coverCandidate = preferredCover || downloads.filter(Boolean).find(function (image) {
    if (image.extension === 'gif' || image.bytes < 5000) return false;
    const substantialLandscape = image.width >= 480 && image.height >= 180;
    const substantialPortrait = image.width >= 300 && image.height >= 320;
    if (!substantialLandscape && !substantialPortrait) return false;
    const ratio = image.width / image.height;
    return ratio >= 0.3 && ratio <= 5;
  });
  return {
    content: rewritten,
    imageCount: localMap.size,
    failedImages: failedImages,
    downloadedBytes: downloadedBytes,
    coverUrl: coverCandidate ? coverCandidate.url : ''
  };
}

function sanitizeArticleHtml(content, pageUrl) {
  let safe = content
    .replace(/<!--[\s\S]*?-->/g, '')
    .replace(/<(script|style|iframe|form|button|textarea|select|object|embed|canvas|video|audio)\b[^>]*>[\s\S]*?<\/\1\s*>/gi, '')
    .replace(/<(input|source|track|link|meta)\b[^>]*\/?\s*>/gi, '');
  const allowed = new Set([
    'a', 'b', 'blockquote', 'br', 'caption', 'code', 'div', 'em', 'figcaption', 'figure',
    'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'hr', 'i', 'img', 'li', 'ol', 'p', 'pre',
    'section', 'small', 'span', 'strong', 'sub', 'sup', 'table', 'tbody', 'td', 'th',
    'thead', 'tr', 'u', 'ul'
  ]);
  safe = safe.replace(/<[^>]+>/g, function (tag) {
    const match = tag.match(/^<\s*(\/?)\s*([a-z0-9-]+)/i);
    if (!match) return '';
    const closing = Boolean(match[1]);
    const name = match[2].toLowerCase();
    if (!allowed.has(name)) return '';
    if (closing) return ['br', 'hr', 'img'].includes(name) ? '' : '</' + name + '>';
    if (name === 'br' || name === 'hr') return '<' + name + '>';
    if (name === 'img') {
      const source = attributeValue(tag, 'src');
      if (!source.startsWith('./images/news/')) return '';
      const alt = attributeValue(tag, 'alt') || 'Archived article image';
      return '<img src="' + escapeHtml(source) + '" alt="' + escapeHtml(alt) + '" loading="lazy" decoding="async">';
    }
    if (name === 'a') {
      const href = normalizedRemoteUrl(attributeValue(tag, 'href'), pageUrl);
      return href ? '<a href="' + escapeHtml(href) + '" target="_blank" rel="noopener noreferrer">' : '<a>';
    }
    if (name === 'td' || name === 'th') {
      const colspan = attributeValue(tag, 'colspan');
      const rowspan = attributeValue(tag, 'rowspan');
      const attrs = (/^\d+$/.test(colspan) ? ' colspan="' + colspan + '"' : '') +
        (/^\d+$/.test(rowspan) ? ' rowspan="' + rowspan + '"' : '');
      return '<' + name + attrs + '>';
    }
    return '<' + name + '>';
  });
  return safe
    .replace(/\{\{/g, '&#123;&#123;')
    .replace(/\{%/g, '&#123;%')
    .replace(/<p>\s*<\/p>/g, '')
    .trim();
}

function archivePage(item, body, archivedAt) {
  return [
    '---',
    'layout: main',
    'title: News Archive',
    'description: Archived media coverage from ' + item.source,
    'base_root: true',
    'noindex: true',
    '---',
    '',
    '<div id="divPost">',
    '  <article id="articlePost" class="news-archive-page">',
    '    <header class="news-archive-header">',
    '      <a class="news-archive-back" href="./news.html">&larr; Back to News</a>',
    '      <p class="news-archive-eyebrow">Locally preserved media coverage</p>',
    '      <h1>' + escapeHtml(item.title) + '</h1>',
    '      <div class="news-archive-meta">',
    '        <time datetime="' + escapeHtml(item.datetime) + '">' + escapeHtml(item.date) + '</time>',
    '        <span aria-hidden="true">·</span>',
    '        <span>' + escapeHtml(item.source) + '</span>',
    '      </div>',
    '      <a class="news-archive-original" href="' + escapeHtml(item.url) + '" target="_blank" rel="noopener noreferrer">Visit original webpage <span aria-hidden="true">&nearr;</span></a>',
    '    </header>',
    '',
    '    <aside class="news-archive-notice">',
    '      This is a locally preserved copy for continuity if the original webpage becomes unavailable. Copyright remains with the original publisher. Archived ' + escapeHtml(archivedAt) + '.',
    '    </aside>',
    '',
    '    <div class="news-archive-content">',
    body,
    '    </div>',
    '',
    '    <footer class="news-archive-source">',
    '      <strong>Original source:</strong> ' + escapeHtml(item.source) + ' ·',
    '      <a href="' + escapeHtml(item.url) + '" target="_blank" rel="noopener noreferrer">' + escapeHtml(item.url) + '</a>',
    '    </footer>',
    '  </article>',
    '  <div class="div0"></div>',
    '</div>',
    ''
  ].join('\n');
}

function existingCover(directory, configuredImage) {
  if (!fs.existsSync(directory)) return '';
  const configured = configuredImage && configuredImage.startsWith('./images/news/') ? path.basename(configuredImage) : '';
  if (configured && /^cover\.(?:avif|gif|jpe?g|png|svg|webp)$/i.test(configured) && fs.existsSync(path.join(directory, configured))) {
    return configured;
  }
  const cover = fs.readdirSync(directory).find(function (fileName) {
    return /^cover\.(?:avif|gif|jpe?g|png|svg|webp)$/i.test(fileName);
  });
  return cover || '';
}

async function archiveItem(item, refresh) {
  const slug = articleSlug(item.url);
  const directory = path.join(imageRoot, slug);
  const archivePath = path.join(archiveRoot, slug + '.html');
  const backupPath = path.join(backupRoot, slug + '.zip');
  const savedCover = existingCover(directory, item.image);
  if (!refresh && (fs.existsSync(backupPath) || fs.existsSync(archivePath)) && savedCover) {
    process.stdout.write('[keep] ' + item.datetime + ' ' + item.title + '\n');
    return {
      oldImage: item.image,
      coverPath: './images/news/' + slug + '/' + savedCover,
      sourceUrl: item.url,
      archiveSlug: slug,
      imageCount: 0,
      failedImages: 0,
      bytes: 0,
      reused: true
    };
  }
  process.stdout.write('[archive] ' + item.datetime + ' ' + item.title + '\n');
  const html = await fetchText(item.url);
  const body = extractMainContent(html);
  const localized = await localizeImages(body, slug, directory, item.url, item.cover_index);
  const safeBody = sanitizeArticleHtml(localized.content, item.url);
  if (safeBody.replace(/<[^>]+>/g, '').trim().length < 80) {
    throw new Error('Extracted article body is unexpectedly short');
  }
  const metadataCover = normalizedRemoteUrl(metaImage(html), item.url);
  const coverUrl = metadataCover || localized.coverUrl;
  const coverMethod = metadataCover ? 'page metadata' : 'article image';
  let coverPath = item.image;
  let coverBytes = 0;
  if (coverUrl) {
    const cover = await downloadImage(coverUrl, directory, 'cover');
    coverPath = './images/news/' + slug + '/' + cover.fileName;
    coverBytes = cover.bytes;
    for (const fileName of fs.readdirSync(directory)) {
      if (/^cover\.(?:avif|gif|jpe?g|png|svg|webp)$/i.test(fileName) && fileName !== cover.fileName) {
        fs.unlinkSync(path.join(directory, fileName));
      }
    }
  }
  process.stdout.write('  cover selected from ' + coverMethod + '\n');
  const archivedAt = new Date().toISOString().slice(0, 10);
  fs.mkdirSync(archiveRoot, { recursive: true });
  fs.writeFileSync(archivePath, archivePage(item, safeBody, archivedAt), 'utf8');
  process.stdout.write('  saved ' + localized.imageCount + ' body images, ' + localized.failedImages + ' skipped\n');
  return {
    oldImage: item.image,
    coverPath: coverPath,
    sourceUrl: item.url,
    archiveSlug: slug,
    imageCount: localized.imageCount,
    failedImages: localized.failedImages,
    bytes: localized.downloadedBytes + coverBytes,
    reused: false
  };
}

async function main() {
  const source = fs.readFileSync(newsPath, 'utf8');
  const items = parseNews(source);
  const limitFlag = process.argv.find(function (arg) { return arg.startsWith('--limit='); });
  const limit = limitFlag ? Number(limitFlag.split('=')[1]) : items.length;
  const slugFlag = process.argv.find(function (arg) { return arg.startsWith('--slug='); });
  const requestedSlug = slugFlag ? slugFlag.slice('--slug='.length) : '';
  const refresh = process.argv.includes('--refresh');
  const selected = requestedSlug ? items.filter(function (item) { return articleSlug(item.url) === requestedSlug; }) : items.slice(0, limit);
  if (requestedSlug && selected.length === 0) throw new Error('No news item found for slug ' + requestedSlug);
  const results = [];
  for (const item of selected) {
    try {
      results.push(await archiveItem(item, refresh));
    } catch (error) {
      process.stderr.write('[failed] ' + item.url + ': ' + error.message + '\n');
      process.exitCode = 1;
    }
  }
  if (results.length) {
    let updated = source;
    for (const result of results) {
      const oldLine = '    image: "' + result.oldImage + '"';
      const newLine = '    image: "' + result.coverPath + '"';
      if (updated.includes(oldLine)) updated = updated.replace(oldLine, newLine);
      const defaultSlug = new URL(result.sourceUrl).pathname.split('/').filter(Boolean).pop();
      const archiveLine = '    archive: "' + result.archiveSlug + '"';
      const urlLine = '    url: "' + result.sourceUrl + '"';
      if (result.archiveSlug !== defaultSlug && !updated.includes(archiveLine) && updated.includes(urlLine)) {
        updated = updated.replace(urlLine, archiveLine + '\n' + urlLine);
      }
    }
    fs.writeFileSync(newsPath, updated, 'utf8');
  }
  const expandedArchives = fs.existsSync(archiveRoot) && fs.readdirSync(archiveRoot).some(function (fileName) {
    return fileName.toLowerCase().endsWith('.html');
  });
  if (expandedArchives) {
    const shell = process.platform === 'win32' ? 'powershell.exe' : 'pwsh';
    const packScript = path.join(root, 'scripts', 'pack-news-archives.ps1');
    const packArguments = ['-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', packScript];
    if (process.argv.includes('--keep-expanded')) packArguments.push('-KeepExpanded');
    const packed = childProcess.spawnSync(shell, packArguments, {
      cwd: root,
      stdio: 'inherit'
    });
    if (packed.error) throw packed.error;
    if (packed.status !== 0) throw new Error('Cold-archive packaging failed with exit code ' + packed.status);
  }
  const totalBytes = results.reduce(function (sum, result) { return sum + result.bytes; }, 0);
  const totalImages = results.reduce(function (sum, result) { return sum + result.imageCount + (result.reused ? 0 : 1); }, 0);
  const failures = results.reduce(function (sum, result) { return sum + result.failedImages; }, 0);
  const reused = results.filter(function (result) { return result.reused; }).length;
  const archived = results.length - reused;
  process.stdout.write('Archived ' + archived + ' new/updated articles; kept ' + reused + ' existing archives; saved ' +
    totalImages + ' images (' + (totalBytes / 1024 / 1024).toFixed(1) + ' MB); ' + failures + ' body images skipped.\n');
}

main().catch(function (error) {
  process.stderr.write((error.stack || error.message) + '\n');
  process.exitCode = 1;
});
