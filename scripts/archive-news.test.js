const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { parseNews, archiveItems, articleSlug, metaImage, extractMainContent } = require('./archive-news');

function fixtureReport(name, datetime, source) {
  return { title: name, datetime, date: datetime, source, image: './images/news/' + name + '/cover.jpg',
    archive: 'example-org-' + name, url: 'https://example.org/' + name + '.htm' };
}

test('a merged event keeps three independently restorable reports with their own metadata', function () {
  const reports = [fixtureReport('college', '2024-01-03', 'College'),
    fixtureReport('society', '2023-12-25', 'Society'), fixtureReport('university', '2023-12-29', 'University')];
  const news = parseNews(JSON.stringify([{ title: 'One award', summary: 'Award announcement', category: 'Recognition', reports }]));
  const articles = archiveItems(news);
  assert.equal(news.length, 1);
  assert.equal(articles.length, 3);
  for (let i = 0; i < reports.length; i++) {
    for (const field of Object.keys(reports[i])) assert.equal(articles[i][field], reports[i][field]);
    assert.equal(articles[i].storyIndex, 0);
    assert.equal(articles[i].reportIndex, i);
    assert.equal(articleSlug(articles[i].url), reports[i].archive);
  }
  assert.throws(() => parseNews(JSON.stringify([{ ...news[0], reports: [] }])), /no reports/);
  assert.throws(() => parseNews(JSON.stringify([{ ...news[0], reports: [reports[0], reports[0]] }])), /Duplicate/);
  const older = { ...news[0], reports: [fixtureReport('older', '2020-01-01', 'Old source')] };
  assert.throws(() => parseNews(JSON.stringify([older, news[0]])), /earliest report date/);
});

test('all site reports remain mapped to their original archive metadata and available cover', function () {
  const root = path.resolve(__dirname, '..');
  const stories = parseNews(fs.readFileSync(path.join(root, '_data/news.json'), 'utf8'));
  const articles = archiveItems(stories);
  const manifest = JSON.parse(fs.readFileSync(path.join(root, 'news-backups/manifest.json'), 'utf8'));
  const byUrl = new Map(manifest.archives.map(entry => [entry.originalUrl, entry]));
  assert.equal(articles.length, byUrl.size);
  for (const article of articles) {
    const archive = byUrl.get(article.url);
    assert.ok(archive, 'Missing archive for ' + article.url);
    assert.equal(article.archive, archive.slug);
    assert.equal(article.title, archive.title);
    assert.equal(article.source, archive.source);
    assert.equal(article.datetime, archive.publishedAt);
    assert.equal(article.image, archive.cover);
    assert.ok(fs.existsSync(path.join(root, article.image)));
    assert.ok(fs.existsSync(path.join(root, archive.package)));
  }
});

test('The Paper archives the full article and uses its configured cover instead of its platform icon', function () {
  const content = '<p>Society report</p><img src=\\"https://example.org/body.jpg\\">';
  const data = { props: { pageProps: { detailData: { contentDetail: {
    content, pic: 'https://example.org/article-cover.jpg'
  } } } } };
  const html = '<meta property="og:image" content="https://file.thepaper.cn/clt/img/gov_small.png">' +
    '<div>Platform navigation</div><script id="__NEXT_DATA__" type="application/json">' + JSON.stringify(data) + '</script>';
  assert.equal(metaImage(html), 'https://example.org/article-cover.jpg');
  assert.equal(extractMainContent(html), '<p>Society report</p><img src="https://example.org/body.jpg">');
  assert.equal(metaImage('<meta property="og:image" content="https://file.thepaper.cn/clt/img/gov_small.png">'), '');
});

test('university articles with different numbered body containers exclude surrounding navigation', function () {
  for (const id of ['vsb_content', 'vsb_content_4', 'vsb_content_1109']) {
    assert.equal(extractMainContent('<div>Navigation</div><div id="' + id + '"><p>Report</p><div>Image caption</div></div><div>Footer</div>'),
      '<p>Report</p><div>Image caption</div>');
  }
});

test('existing WeChat metadata cover and article-body extraction remain supported', function () {
  const html = '<meta property="og:image" content="https://example.org/wechat.jpg?a=1&amp;b=2">' +
    '<div id="js_content"><p>WeChat article</p></div>';
  assert.equal(metaImage(html), 'https://example.org/wechat.jpg?a=1&b=2');
  assert.equal(extractMainContent(html), '<p>WeChat article</p>');
});
