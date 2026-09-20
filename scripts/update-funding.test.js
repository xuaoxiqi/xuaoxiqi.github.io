const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { spawnSync } = require("node:child_process");
const test = require("node:test");

const projectId = "12672289";
const defaultConfig = {
  pdftotextCommand: "missing-pdftotext-test-command",
  rankTokens: [{ id: projectId }],
  projects: [{ id: projectId, grantIds: [projectId] }],
  statusLabels: { completed: "(Completed)", coming: "(Coming..)" }
};
const researchHtml = [
  '<h4 class="project-title">Project title</h4>',
  '<div class="project-meta">RMB 530,000</div>',
  '<details class="project-publications">',
  '    <summary>Publications</summary>',
  '    <div class="project-target">(Target) 6-8 papers</div>',
  `    <!-- funding:auto:start project=${projectId} status=completed -->`,
  "    stale completed content",
  "    <!-- funding:auto:end -->",
  `    <!-- funding:auto:start project=${projectId} status=coming -->`,
  "    stale coming content",
  "    <!-- funding:auto:end -->",
  "</details>",
  ""
].join("\n");

function fixture(t, { config = defaultConfig, html = researchHtml } = {}) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), "update-funding-test-"));
  t.after(() => {
    const resolved = fs.realpathSync(root);
    assert.equal(path.dirname(resolved), fs.realpathSync(os.tmpdir()));
    assert.ok(path.basename(resolved).startsWith("update-funding-test-"));
    fs.rmSync(resolved, { recursive: true, force: true });
  });
  fs.mkdirSync(path.join(root, "scripts"));
  fs.mkdirSync(path.join(root, "data"));
  for (const name of ["update-funding.js", "update-funding.ps1"]) {
    fs.copyFileSync(path.join(__dirname, name), path.join(root, "scripts", name));
  }
  fs.writeFileSync(path.join(root, "data", "funding-projects.json"), JSON.stringify(config));
  fs.writeFileSync(path.join(root, "publications.html"), "<ol></ol>");
  const researchPath = path.join(root, "research.html");
  fs.writeFileSync(researchPath, html);
  return {
    root,
    read: () => fs.readFileSync(researchPath, "utf8"),
    node: (...args) => spawnSync(process.execPath, [path.join(root, "scripts", "update-funding.js"), ...args], {
      cwd: os.tmpdir(), encoding: "utf8"
    }),
    powershell: (...args) => spawnSync("powershell.exe", [
      "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
      path.join(root, "scripts", "update-funding.ps1"), ...args
    ], { cwd: os.tmpdir(), encoding: "utf8" })
  };
}

function expectedEmptyResult(html) {
  return html.replace("stale completed content", "(Completed) 0 papers:")
    .replace("stale coming content", "<br>(Coming..) 0 papers:");
}

for (const [name, eol] of [["LF", "\n"], ["CRLF", "\r\n"]]) {
  test(`empty results preserve layout and ${name} line endings`, t => {
    const html = researchHtml.replace(/\n/g, eol);
    const env = fixture(t, { html });
    const result = env.node();
    assert.equal(result.status, 0, result.stderr);
    assert.equal(env.read(), expectedEmptyResult(html));
    assert.doesNotMatch(env.read(), /i\) _R/);
    assert.match(result.stdout, /Updated 2 funding block/);
  });
}

test("dry-run reports blocks without changing research.html", t => {
  const env = fixture(t);
  const result = env.node("--dry-run", "--verbose");
  assert.equal(result.status, 0, result.stderr);
  assert.equal(env.read(), researchHtml);
  assert.match(result.stdout, /Checked 2 funding block/);
});

const invalidCases = [
  {
    name: "project missing from configuration",
    config: { ...defaultConfig, projects: [] },
    error: /12672289.*missing from projects/
  },
  {
    name: "grant missing from rankTokens",
    config: { ...defaultConfig, rankTokens: [] },
    error: /12672289.*missing from rankTokens/
  },
  {
    name: "project has an empty grantIds list",
    config: { ...defaultConfig, projects: [{ id: projectId, grantIds: [] }] },
    error: /12672289 has no grantIds/
  },
  {
    name: "unknown status",
    html: researchHtml.replace("status=coming", "status=pending"),
    error: /unsupported funding status "pending"/
  },
  {
    name: "unmatched block marker",
    html: researchHtml.replace("<!-- funding:auto:end -->", ""),
    error: /Malformed or unmatched/
  }
];

for (const scenario of invalidCases) {
  test(`${scenario.name} fails before PDF processing or writes`, t => {
    const env = fixture(t, scenario);
    // A missing PDF must not mask the earlier configuration error.
    fs.writeFileSync(path.join(env.root, "publications.html"),
      '<!-- Paper --><li class="pub-item"><a href="./publications/missing.pdf">PDF</a></li>');
    const before = env.read();
    const result = env.node();
    assert.equal(result.status, 1);
    assert.match(result.stderr, scenario.error);
    assert.doesNotMatch(result.stderr, /failed to read|ENOENT/);
    assert.equal(env.read(), before);
  });
}

test("unknown command options do not trigger an update", t => {
  const env = fixture(t);
  const result = env.node("--dryrun");
  assert.equal(result.status, 1);
  assert.match(result.stderr, /Unknown argument/);
  assert.equal(env.read(), researchHtml);
});

test("PowerShell forwards flags and uses the same empty-result behavior", { skip: process.platform !== "win32" }, t => {
  const env = fixture(t);
  const dry = env.powershell("-DryRun", "-Verbose");
  assert.equal(dry.status, 0, dry.stderr);
  assert.match(dry.stdout, /Checked 2 funding block/);
  assert.equal(env.read(), researchHtml);
  const write = env.powershell();
  assert.equal(write.status, 0, write.stderr);
  assert.equal(write.stderr, "");
  assert.equal(env.read(), expectedEmptyResult(researchHtml));
});

test("PowerShell returns Node's failure exit code without writing", { skip: process.platform !== "win32" }, t => {
  const env = fixture(t, { config: { ...defaultConfig, projects: [] } });
  const result = env.powershell();
  assert.equal(result.status, 1);
  assert.match(result.stderr, /12672289.*missing from projects/);
  assert.equal(env.read(), researchHtml);
});
