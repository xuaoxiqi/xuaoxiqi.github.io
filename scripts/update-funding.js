const fs = require("fs");
const path = require("path");
const { execFileSync } = require("child_process");

const ROOT = path.resolve(__dirname, "..");
const CONFIG_PATH = path.join(ROOT, "data", "funding-projects.json");
const PUBLICATIONS_PATH = path.join(ROOT, "publications.html");
const RESEARCH_PATH = path.join(ROOT, "research.html");

const args = new Set(process.argv.slice(2));
const dryRun = args.has("--dry-run");
const verbose = args.has("--verbose");

const config = JSON.parse(fs.readFileSync(CONFIG_PATH, "utf8"));
const rankTokens = config.rankTokens || [];
const projects = config.projects || [];

function escapeRegex(value) {
  return value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function defaultLabelFromStem(stem) {
  return stem.replace(/^\d+_/, "").replace(/_arXiv$/i, "");
}

function normalizeHref(href) {
  return href.replace(/^\.\//, "").replace(/\//g, path.sep);
}

function runPdftotext(pdfPath) {
  return execFileSync(config.pdftotextCommand || "pdftotext", ["-layout", pdfPath, "-"], {
    cwd: ROOT,
    encoding: "utf8",
    maxBuffer: 20 * 1024 * 1024,
    stdio: ["ignore", "pipe", "pipe"]
  });
}

function extractFundingText(text) {
  const markerPattern = /\b(?:Funding\s*[\.:]|Acknowledgements?\.?|Acknowledgments?\.?|ACKNOWLEDGMENTS)\b/gi;
  const markers = [...text.matchAll(markerPattern)];
  if (markers.length === 0) return "";

  const start = markers[markers.length - 1].index;
  let chunk = text.slice(start, start + 3500);
  const stopPattern =
    /\b(?:Declaration of interests|Conflict of interest|Author ORCIDs|Data availability|AUTHOR DECLARATIONS|DATA AVAILABILITY|APPENDIX|Appendix|REFERENCES|References|Supplementary materials)\b/i;
  const stop = chunk.slice(1).search(stopPattern);
  if (stop >= 0) {
    chunk = chunk.slice(0, stop + 1);
  }
  return chunk;
}

function findToken(text, token) {
  let best = null;
  for (const pattern of token.patterns || [token.id]) {
    const re = new RegExp(`(^|[^A-Za-z0-9])(${escapeRegex(pattern)})(?![A-Za-z0-9])`, "i");
    const match = re.exec(text);
    if (!match) continue;
    const index = match.index + match[1].length;
    if (best === null || index < best) best = index;
  }
  return best;
}

function extractFundingOrder(fundingText) {
  const hits = [];
  for (const token of rankTokens) {
    const index = findToken(fundingText, token);
    if (index !== null) hits.push({ id: token.id, index });
  }
  hits.sort((a, b) => a.index - b.index);

  const seen = new Set();
  return hits
    .filter((hit) => {
      if (seen.has(hit.id)) return false;
      seen.add(hit.id);
      return true;
    })
    .map((hit) => hit.id);
}

function parsePublications() {
  const html = fs.readFileSync(PUBLICATIONS_PATH, "utf8");
  const itemPattern = /<!--\s*([^]*?)\s*-->\s*<li class="pub-item[^"]*"[^>]*>([\s\S]*?)<\/li>/g;
  const candidatesByLabel = new Map();
  const readErrors = [];
  let match;
  let order = 0;

  function addCandidate(relativePdfPath, itemOrder) {
    const stem = path.basename(relativePdfPath, ".pdf");
    const label = defaultLabelFromStem(stem);
    const isArxiv = /_arXiv$/i.test(stem);
    const candidate = {
      label,
      isArxiv,
      relativePdfPath,
      absolutePdfPath: path.join(ROOT, relativePdfPath),
      order: itemOrder
    };
    const candidates = candidatesByLabel.get(label) || [];
    candidates.push(candidate);
    candidatesByLabel.set(label, candidates);
  }

  while ((match = itemPattern.exec(html)) !== null) {
    const block = match[2];
    const hrefs = [...block.matchAll(/<a\b[^>]*href="([^"]+\.pdf)"/gi)]
      .map((hrefMatch) => hrefMatch[1])
      .filter((href) => /(?:^|\/)publications\//.test(href.replace(/^\.\//, "")));
    if (hrefs.length === 0) continue;

    const itemOrder = order++;
    for (const href of hrefs) {
      const relativePdfPath = normalizeHref(href);
      addCandidate(relativePdfPath, itemOrder);

      if (/_arXiv\.pdf$/i.test(relativePdfPath)) {
        const finalPdfPath = relativePdfPath.replace(/_arXiv\.pdf$/i, ".pdf");
        if (fs.existsSync(path.join(ROOT, finalPdfPath))) {
          addCandidate(finalPdfPath, itemOrder);
        }
      }
    }
  }

  const entries = [];
  for (const [label, candidates] of candidatesByLabel.entries()) {
    candidates.sort((a, b) => a.order - b.order || Number(a.isArxiv) - Number(b.isArxiv));
    const chosen = candidates.find((candidate) => !candidate.isArxiv) || candidates[0];
    const status = chosen.isArxiv ? "coming" : "completed";

    let fundingOrder = [];
    try {
      const pdfText = runPdftotext(chosen.absolutePdfPath);
      const fundingText = extractFundingText(pdfText);
      fundingOrder = extractFundingOrder(fundingText || pdfText);
      if (verbose) {
        console.log(`${label} [${status}, ${chosen.relativePdfPath}]: ${fundingOrder.join(", ") || "no funding tokens found"}`);
      }
    } catch (error) {
      const message = `${chosen.relativePdfPath}: ${error.message}`;
      console.warn(`Warning: failed to read ${message}`);
      readErrors.push(message);
    }

    if (fundingOrder.length > 0) {
      entries.push({
        label,
        status,
        fundingOrder,
        source: chosen.relativePdfPath,
        order: chosen.order
      });
    }
  }

  if (readErrors.length > 0) {
    throw new Error(`Failed to read ${readErrors.length} PDF(s):\n${readErrors.join("\n")}`);
  }

  return entries;
}

function projectRank(entry, project) {
  for (const grantId of project.grantIds || [project.id]) {
    const index = entry.fundingOrder.indexOf(grantId);
    if (index >= 0) return index + 1;
  }
  return null;
}

function sortEntries(entries, status) {
  const preferred = config.statusOrder?.[status] || [];
  const position = (label) => {
    const index = preferred.indexOf(label);
    return index === -1 ? Number.POSITIVE_INFINITY : index;
  };
  return [...entries].sort((a, b) => {
    const preferredDiff = position(a.label) - position(b.label);
    if (preferredDiff !== 0) return preferredDiff;
    return a.order - b.order;
  });
}

function buildProjectResults(entries) {
  const results = {};
  for (const project of projects) {
    results[project.id] = {};
    for (const entry of entries) {
      const rank = projectRank(entry, project);
      if (rank === null) continue;
      const statusEntries = results[project.id][entry.status] || [];
      statusEntries.push({ ...entry, rank });
      results[project.id][entry.status] = statusEntries;
    }
    for (const status of Object.keys(results[project.id])) {
      results[project.id][status] = sortEntries(results[project.id][status], status);
    }
  }
  return results;
}

function countSummary(entries) {
  if (entries.length === 0) return "0";
  const maxRank = Math.max(...entries.map((entry) => entry.rank));
  const counts = Array.from({ length: maxRank }, () => 0);
  for (const entry of entries) {
    counts[entry.rank - 1] += 1;
  }
  return counts.join("+");
}

function roman(value) {
  const map = [
    [1000, "m"],
    [900, "cm"],
    [500, "d"],
    [400, "cd"],
    [100, "c"],
    [90, "xc"],
    [50, "l"],
    [40, "xl"],
    [10, "x"],
    [9, "ix"],
    [5, "v"],
    [4, "iv"],
    [1, "i"]
  ];
  let remaining = value;
  let result = "";
  for (const [amount, symbol] of map) {
    while (remaining >= amount) {
      result += symbol;
      remaining -= amount;
    }
  }
  return result;
}

function renderBlock(projectId, status, indent, results, eol = "\n") {
  const entries = results[projectId]?.[status] || [];
  const label = config.statusLabels?.[status] || `(${status})`;
  const labelPrefix = status === "coming" ? "<br>" : "";
  const lines = [`${indent}${labelPrefix}${label} ${countSummary(entries)} papers:${entries.length > 0 ? "<br>" : ""}`];
  if (entries.length > 0) {
    const total = entries.length;
    const items = entries
      .map((entry, index) => `${roman(total - index)}) ${entry.label}_R${entry.rank}`)
      .join("; ");
    lines.push(`${indent}${items};`);
  }
  return lines.join(eol);
}

function updateResearch(results) {
  const html = fs.readFileSync(RESEARCH_PATH, "utf8");
  const crlfCount = (html.match(/\r\n/g) || []).length;
  const lfOnlyCount = (html.match(/(?<!\r)\n/g) || []).length;
  const eol = crlfCount >= lfOnlyCount ? "\r\n" : "\n";
  const markerPattern =
    /([ \t]*)<!-- funding:auto:start\s+project=([^\s]+)\s+status=([^\s]+)\s*-->\r?\n[\s\S]*?\r?\n[ \t]*<!-- funding:auto:end -->/g;
  let replacements = 0;
  const updated = html.replace(markerPattern, (match, indent, projectId, status) => {
    replacements += 1;
    return [
      `${indent}<!-- funding:auto:start project=${projectId} status=${status} -->`,
      renderBlock(projectId, status, indent, results, eol),
      `${indent}<!-- funding:auto:end -->`
    ].join(eol);
  });

  if (replacements === 0) {
    throw new Error("No funding:auto markers found in research.html.");
  }

  if (!dryRun) {
    fs.writeFileSync(RESEARCH_PATH, updated, "utf8");
  }
  return replacements;
}

function printSummary(results) {
  for (const project of projects) {
    const statuses = results[project.id] || {};
    for (const status of Object.keys(statuses)) {
      const entries = statuses[status];
      if (entries.length === 0) continue;
      const detail = entries.map((entry) => `${entry.label}_R${entry.rank}`).join(", ");
      console.log(`${project.id} ${status}: ${countSummary(entries)} (${detail})`);
    }
  }
}

const entries = parsePublications();
const results = buildProjectResults(entries);
const replacementCount = updateResearch(results);
printSummary(results);
console.log(`${dryRun ? "Checked" : "Updated"} ${replacementCount} funding block(s) in research.html.`);
