# Builds docs/outputs/final_report/report.html: every new-data analysis step as one page.
# Each chapter is the published report embedded as a self-contained, auto-height frame, so the
# chapters keep their own styles, scripts and equations and none of the source pages is modified.
# Rerun after any chapter's report is regenerated:   julia build_final_report.jl
# Every build that changes the page is also archived, never overwritten, under
# docs/outputs/final_report/history/<version>/ with its modelling decisions (report_history.jl); the page opens
# with those decisions against the previous version, and history/index.html compares every version.

using Dates, TOML
include(joinpath(@__DIR__, "report_history.jl"))
# length of the long dose + fasting run (control_config.toml), for chapter 7
const LONG_DAYS = round(Int, TOML.parsefile(joinpath(@__DIR__, "Modeling_Approaches", "02_mechanical_automatic_package", "control_config.toml"))["schedule"]["long_days"])

const DOCS = joinpath(@__DIR__, "docs")
const OUT  = joinpath(DOCS, "outputs", "final_report")
const PKG  = "Modeling_Approaches/02_mechanical_automatic_package/outputs"

# (id, title, what the step does, path under docs/)
chapters = [
    ("data", "1 · Data", "Every CSV and XLSX sheet in <code>New Datasets</code>, plotted as mean cells per day on shared axes.",
     "outputs/new_data_report/report.html"),
    ("models", "2 · Compartment models (the models used)", "Stage 1 untreated → 2a low resource → 2b treated → 3a Ce0 + Ce1 with a shared carrying capacity → 3b low-resource coculture → 4 low resource + cisplatin coculture (both inputs at once). Every model is autonomous: delays and accumulation come from drug-effect, stress or damaged-cell compartments driven by the dose, never from time, with at most 3 compartments per mechanism. The model carried forward at each stage has the fewest free parameters within 2 BIC of the lowest, and both lines share one cisplatin mechanism.",
     "$PKG/zoo_fits/report.html"),
    ("bootstrap", "3 · Bootstrapping", "200-replicate moving-block residual bootstrap with every candidate refitted and the selection rule replayed, compared with the package's within-trajectory bootstrap: how often each model is chosen, 95% parameter intervals and bands.",
     "$PKG/zoo_analyses/bootstrap/report.html"),
    ("sensitivity", "4 · Sensitivity analysis", "Local sensitivity index of every fitted parameter, profile-likelihood identifiability and the bound check.",
     "$PKG/zoo_analyses/sensitivity/report.html"),
    ("validation", "5 · Validation", "Held-out treated cocultures predicted from the untreated cocultures and the monoculture drug fits, blocked leave-one-dose / ratio / run-out validation, and forecasts of the last 30% of days.",
     "$PKG/zoo_analyses/validation/report.html"),
    ("control", "6 · Optimal control", "Sequential (successive-linearisation) model-predictive control of the cisplatin dose at each media change, on the carried Stage 3 model (normal medium) and the same model with the fitted low-resource modules (low resource; equal to the Stage 4 fit at 1 µM), with one stated goal (the longest time to progression: N + C held under B_prog = K/2, with the resistant cells left at the end priced), a literature toxicity proxy, a burden penalty and a toxicity barrier, baseline schedules, bootstrap robustness, identifiability warnings and suggested experiments.",
     "$PKG/control/report.html"),
    ("control_30d", "7a · Optimal control over 30 days (cisplatin; beyond the data)", "The same controller and comparators as chapter 6, run for 30 days with the look-ahead ending at the run's end (user, 2026-10-27). Days 15–30 follow the fitted equations beyond the data.",
     "$PKG/control_30d/report.html"),
    ("control_60d", "7b · Optimal control over 60 days (cisplatin; beyond the data)", "The same controller and comparators, run for 60 days with the look-ahead ending at the run's end (user, 2026-10-27). Days 15–60 are extrapolation.",
     "$PKG/control_60d/report.html"),
    ("sweep", "7c · Control sweep: barriers, thresholds, run lengths", "The control study rerun under four burden-barrier shapes (relaxed log, quadratic, hard, none), thresholds from 0.5 to 0.8 K, 14/30/60-day runs, fixed vs run-ending look-ahead, other objectives, and with and without the A2780cis 1 µM monoculture: which control results come from the biology and which from modelling decisions.",
     "$PKG/control_sweep/report.html"),
    ("control30", "8 · Optimal control over $LONG_DAYS days: cisplatin and fasting", "Model-predictive control choosing both the cisplatin dose and fasting (low resource) at every media change for $LONG_DAYS days, read off at days 14, 30, 60 and $LONG_DAYS of the same run, with a body-weight limit (cachexia threshold) and fasting protection of healthy tissue from the literature. Days 14–$LONG_DAYS are beyond the data, and the fasting effect is the in-vitro effect at full strength.",
     "$PKG/control_$(LONG_DAYS)d/report.html"),
    ("findings", "9 · Key findings, assumptions and limitations", "What did not carry over between conditions, the estimated starting values, validation, the controller against the best fixed policies, linear vs Hill drug response, the mechanistic vs neural-network comparison, the assumptions register and the limitations with the experiments that would resolve them. Generated from the outputs of the same run.",
     "$PKG/findings/report.html"),
    ("nn", "10 · Mechanistic vs neural-network models", "The neural-network and universal-differential-equation track and its matched comparison with the mechanistic models on the same held-out data.",
     "outputs/nn_track/comparison/report.html"),
]

# Turns each chapter into dropdowns so nothing has to be scrolled past: every top-level (h2)
# section and every card-style section (a section/card whose first element is an h3) becomes a
# closed <details>, and any table longer than 8 rows is folded behind "Show table". The content is
# unchanged; only its visibility is.
const COLLAPSE_JS = raw"""<style>
details.fr-sec>summary,details.fr-tab>summary{cursor:pointer;list-style:none;display:flex;align-items:center;gap:8px}
details.fr-sec>summary::-webkit-details-marker,details.fr-tab>summary::-webkit-details-marker{display:none}
details.fr-sec>summary::before,details.fr-tab>summary::before{content:'▸';color:#176b87;font-size:15px;transition:transform .15s;flex:none}
details.fr-sec[open]>summary::before,details.fr-tab[open]>summary::before{transform:rotate(90deg)}
details.fr-sec>summary>h2,details.fr-sec>summary>h3{display:inline;margin:0}
details.fr-h2{margin:14px 0;border-bottom:1px solid #d9e2e6;padding-bottom:6px}details.fr-h2>summary{padding:6px 0}
details.fr-box>summary{padding:2px 0}
details.fr-tab{margin:6px 0}details.fr-tab>summary{font:600 13px system-ui,sans-serif;color:#176b87;padding:4px 0}
.fr-tools{display:flex;gap:8px;margin:0 0 10px}.fr-tools button{border:1px solid #c5d1d8;background:#f7fbfc;color:#176b87;border-radius:999px;padding:4px 12px;font:600 13px system-ui,sans-serif;cursor:pointer}
</style><script>
document.addEventListener('DOMContentLoaded', function () {
  var mk = function (cls, head) {
    var d = document.createElement('details'), s = document.createElement('summary');
    d.className = 'fr-sec ' + cls; s.appendChild(head); d.appendChild(s); return d;
  };
  // 1. card-style sections: the h3 becomes the summary of a details inside the card
  document.querySelectorAll('section, .panel, .card, .workbook').forEach(function (box) {
    var h = box.firstElementChild;
    if (!h || h.tagName !== 'H3' || box.closest('summary')) return;
    var d = mk('fr-box', h);
    while (box.firstChild) d.appendChild(box.firstChild);
    box.appendChild(d);
  });
  // 2. long tables fold behind "Show table"
  document.querySelectorAll('table').forEach(function (t) {
    var n = t.querySelectorAll('tr').length;
    if (n <= 9 || t.closest('summary') || t.closest('details.fr-tab')) return;
    var host = t.parentElement && t.parentElement.classList.contains('tscroll') ? t.parentElement : t;
    var d = document.createElement('details'), s = document.createElement('summary');
    d.className = 'fr-tab'; s.textContent = 'Show table (' + (n - 1) + ' rows)'; d.appendChild(s);
    host.parentNode.insertBefore(d, host); d.appendChild(host);
  });
  // 3. each h2 and everything up to the next h2 becomes one dropdown
  var parents = [];
  document.querySelectorAll('h2').forEach(function (h) { if (!h.closest('summary') && parents.indexOf(h.parentElement) < 0) parents.push(h.parentElement); });
  parents.forEach(function (par) {
    var cur = null;
    Array.prototype.slice.call(par.children).forEach(function (el) {
      if (el.tagName === 'H2') { var d = mk('fr-h2', el.cloneNode(true)); par.replaceChild(d, el); cur = d; return; }
      if (cur && el.tagName !== 'SCRIPT' && el.tagName !== 'STYLE') cur.appendChild(el);
    });
  });
  // open the dropdowns around any in-page link target
  var reveal = function (id) { var t = id && document.getElementById(id); for (var e = t; e; e = e.parentElement) if (e.tagName === 'DETAILS') e.open = true; };
  reveal(location.hash.slice(1));
  document.addEventListener('click', function (ev) { var a = ev.target.closest && ev.target.closest('a[href^="#"]'); if (a) reveal(a.getAttribute('href').slice(1)); });
  // expand / collapse everything in this chapter
  var first = document.querySelector('details.fr-h2, details.fr-box');
  if (first) {
    var bar = document.createElement('div'); bar.className = 'fr-tools';
    bar.innerHTML = '<button type="button" data-v="1">Expand all</button><button type="button" data-v="0">Collapse all</button>';
    bar.addEventListener('click', function (ev) { var v = ev.target.dataset && ev.target.dataset.v; if (v === undefined) return;
      document.querySelectorAll('details.fr-sec, details.fr-tab').forEach(function (d) { d.open = v === '1'; }); });
    first.parentNode.insertBefore(bar, first);
  }
  if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise();
});
</script>"""

# hide each chapter's own tab bar; links resolve (and open) relative to the original page
const CHAPTER_HEAD = "<style>nav.topnav, nav.tab-nav, nav[aria-label='Low-resource report navigation'], a.back, a.back-home{display:none !important}</style>" * COLLAPSE_JS

function json_str(s::AbstractString)
    io = IOBuffer(); print(io, '"')
    for c in s
        c == '"'  ? print(io, "\\\"") :
        c == '\\' ? print(io, "\\\\") :
        c == '\n' ? print(io, "\\n") :
        c == '\r' ? print(io, "\\r") :
        c == '\t' ? print(io, "\\t") :
        c == ' ' ? print(io, "\\u2028") :
        c == ' ' ? print(io, "\\u2029") :
        c < ' '   ? print(io, "\\u", string(UInt16(c); base=16, pad=4)) : print(io, c)
    end
    print(io, '"')
    replace(String(take!(io)), "</" => "<\\/")      # never close the host <script> early
end

parts = String[]
chapter_stamps = Tuple{String,String}[]               # (title, report time) of every chapter, for the version record
toc = String[]
missing_parts = String[]
for (id, title, blurb, rel) in chapters
    path = joinpath(DOCS, rel)
    if !isfile(path)
        push!(missing_parts, rel); continue
    end
    html = read(path, String)
    base = "../../" * rel                                          # relative to docs/outputs/final_report/
    inject = "<base href=\"$base\" target=\"_top\">$CHAPTER_HEAD"
    html = occursin(r"<head[^>]*>"i, html) ? replace(html, r"<head[^>]*>"i => s -> s * inject; count=1) : inject * html
    stamp = Dates.format(unix2datetime(mtime(path)), "yyyy-mm-dd HH:MM")
    push!(parts, "{id:$(json_str(id)),title:$(json_str(title)),blurb:$(json_str(blurb)),src:$(json_str(base)),stamp:$(json_str(stamp)),html:$(json_str(html))}")
    push!(toc, "<a href='#$id'>$title</a>"); push!(chapter_stamps, (title, stamp))
end

# ── version and modelling decisions ─────────────────────────────────────────────
settings = report_settings(); results = report_results()
content_sha = bytes2hex(sha256(join(parts) * hjson([collect(x) for x in settings]) * hjson([collect(x) for x in results]) * hjson(decision_log())))
current = joinpath(OUT, "report.html")
archive_unrecorded(current)                          # a report.html from before the history is kept, not lost
prev_all = load_versions()
same = !isempty(prev_all) && get(prev_all[end], "content_sha", "") == content_sha
version_id = same ? prev_all[end]["id"] : Dates.format(now(), "yyyy-mm-dd_HHMM")
while !same && isdir(joinpath(HIST, version_id))     # two builds in one minute
    global version_id = version_id * "b"
end
prev = same ? (length(prev_all) > 1 ? prev_all[end-1] : nothing) : (isempty(prev_all) ? nothing : prev_all[end])
DECISIONS = decisions_panel(version_id, settings, results, prev; chapters=chapter_stamps)
HISTORY = history_panel(version_id, settings)

page = """<!doctype html>
<html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>Final report — A2780 new data</title>
<style>
:root{--ink:#17212b;--muted:#52636b;--accent:#176b87;--line:#d9e2e6;--page:#f4f7f8}
*{box-sizing:border-box}body{margin:0;background:var(--page);color:var(--ink);font:15px/1.5 system-ui,sans-serif}
header{max-width:1400px;margin:auto;padding:24px 24px 8px}h1{color:#124b6e;margin:8px 0 6px;font-size:32px}
.lead{color:var(--muted);max-width:980px;margin:0 0 12px}
.home{display:inline-flex;align-items:center;padding:6px 9px;border-radius:8px;background:#eef5f6;color:var(--accent);text-decoration:none;font-weight:650}
nav.toc{position:sticky;top:0;z-index:5;background:var(--page);border-bottom:1px solid var(--line)}
nav.toc div{max-width:1400px;margin:auto;padding:10px 24px;display:flex;gap:8px;flex-wrap:wrap}
nav.toc a{padding:6px 11px;border-radius:8px;background:#eef5f6;color:var(--accent);text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}
nav.toc a.on{background:var(--accent);color:#fff}
section.chapter{max-width:1440px;margin:14px auto 0;padding:0 12px}
section.chapter>details{background:#fff;border:1px solid var(--line);border-radius:12px}
section.chapter>details>summary{list-style:none;cursor:pointer;padding:14px 16px;display:flex;gap:12px;align-items:flex-start}
section.chapter>details>summary::-webkit-details-marker{display:none}
section.chapter>details>summary::before{content:'▸';color:var(--accent);font-size:20px;line-height:1.2;transition:transform .15s}
section.chapter>details[open]>summary::before{transform:rotate(90deg)}
section.chapter>details[open]>summary{border-bottom:1px solid var(--line)}
.chead h2{margin:0;color:#124b6e;font-size:22px}.chead p{margin:4px 0 0;color:var(--muted)}.chead small{color:#7b8a91}
.tools{display:flex;gap:8px;margin:0 0 6px}.tools button{border:1px solid var(--line);background:#eef5f6;color:var(--accent);border-radius:999px;padding:5px 12px;font:650 13px system-ui,sans-serif;cursor:pointer}
iframe{display:block;width:100%;border:0;min-height:400px;background:var(--page)}
.warn{background:#fff6d9;border-radius:8px;padding:10px 14px}
section.decisions{max-width:1400px;margin:8px auto 0;padding:14px 18px;background:#fff;border:1px solid var(--line);border-left:6px solid var(--accent);border-radius:12px}
section.decisions h2{margin:0 0 4px;color:#124b6e;font-size:22px}section.decisions .tscroll{overflow-x:auto}
section.decisions table{border-collapse:collapse;font-size:13px;width:100%;margin:6px 0}section.decisions th,section.decisions td{border-bottom:1px solid #e5eaec;padding:5px 8px;text-align:left;vertical-align:top}
section.decisions th{background:#eef3f5}section.decisions tr.chg td{background:#fff3c4}section.decisions .tag{color:#9a6700;font-size:12px}section.decisions .mut{color:var(--muted)}
section.decisions details{margin:6px 0}section.decisions summary{cursor:pointer;font-weight:650;color:var(--accent)}
.version{color:var(--muted);font-size:13px}
</style></head><body>
<header>
<a class='home' href='../../index.html' aria-label='Reports home'>⌂ Reports home</a>
<h1>Final report — A2780 / A2780cis new data</h1>
<p class='lead'>The analysis on one page: the data, the compartment models used at every stage, their bootstrapping, sensitivity analysis and validation, and optimal control of the cisplatin schedule built on them, over 14 days (and, beyond the data, 30 and 60 days), a sweep over the control's modelling decisions, and, with fasting as a second control, over $LONG_DAYS days (read off at days 14, 30, 60 and $LONG_DAYS). The modelling decisions of this version come first, next to the previous version. Earlier time-based models are kept only in the <a href='../archive_time/index.html'>time-based archive</a>. Each chapter is the full report for that step (assembled $(Dates.format(now(), "yyyy-mm-dd HH:MM"))); the date under each heading is when that report was last generated. Rebuild with <code>julia build_final_report.jl</code>.</p>
<p class='version'>Version $(version_id) · <a href='#history'>report history</a> · <a href='history/index.html'>all versions and what changed</a></p>
$(isempty(missing_parts) ? "" : "<p class='warn'>Not found, so left out: " * join(missing_parts, ", ") * "</p>")
</header>
<div style='padding:0 12px'>$DECISIONS</div>
<div style='padding:0 12px'>$HISTORY</div>
<nav class='toc' aria-label='Final report chapters'><div>$(join(toc))</div></nav>
<div style='max-width:1400px;margin:10px auto 0;padding:0 24px' class='tools'><button type='button' id='openall'>Open all chapters</button><button type='button' id='closeall'>Close all</button><span style='color:var(--muted);font-size:13px;align-self:center'>Every chapter and section is a dropdown: click a heading to show it.</span></div>
<main id='chapters'></main>
<script>
const PARTS = [$(join(parts, ",\n"))];
const main = document.getElementById('chapters');
for (const p of PARTS) {
  const sec = document.createElement('section'); sec.className = 'chapter'; sec.id = p.id;
  sec.innerHTML = '<details><summary><div class="chead"><h2></h2><p></p><small></small></div></summary></details>';
  sec.querySelector('h2').textContent = p.title;
  sec.querySelector('p').innerHTML = p.blurb;
  sec.querySelector('small').innerHTML = 'Report generated ' + p.stamp + ' · <a href="' + p.src + '">open on its own</a>';
  const det = sec.querySelector('details');
  // the chapter loads the first time it is opened, so the page itself stays light
  det.addEventListener('toggle', () => {
    if (!det.open || det.querySelector('iframe')) return;
    const fr = document.createElement('iframe'); fr.title = p.title; fr.srcdoc = p.html;
    fr.addEventListener('load', () => {
      const doc = fr.contentDocument; if (!doc) return;
      const fit = () => { fr.style.height = Math.max(doc.documentElement.scrollHeight, doc.body ? doc.body.scrollHeight : 0) + 'px'; };
      fit(); new ResizeObserver(fit).observe(doc.documentElement);
      if (doc.body) new ResizeObserver(fit).observe(doc.body);
      doc.addEventListener('toggle', () => setTimeout(fit, 30), true);
      setTimeout(fit, 1500); setTimeout(fit, 5000);
    });
    det.appendChild(fr);
  });
  main.appendChild(sec);
}
const openChapter = id => { const d = document.querySelector('#' + id + ' > details'); if (d) d.open = true; };
document.querySelectorAll('nav.toc a').forEach(a => a.addEventListener('click', () => openChapter(a.getAttribute('href').slice(1))));
if (location.hash) openChapter(location.hash.slice(1));
document.getElementById('openall').onclick = () => document.querySelectorAll('section.chapter > details').forEach(d => d.open = true);
document.getElementById('closeall').onclick = () => document.querySelectorAll('section.chapter > details').forEach(d => d.open = false);
// highlight the chapter in view
const links = [...document.querySelectorAll('nav.toc a')];
const obs = new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) links.forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + e.target.id)); }), {rootMargin: '-40% 0px -55% 0px'});
document.querySelectorAll('section.chapter').forEach(s => obs.observe(s));
</script>
</body></html>
"""

mkpath(OUT)
include(joinpath(@__DIR__, "Modeling_Approaches", "02_mechanical_automatic_package", "stage_names.jl"))   # thesis stage names (1, 2a, 2b, 3a, 3b, 4)
final = stage_names(page)
write(current, final)
println("Wrote $current ($(round(filesize(current) / 1e6; digits=2)) MB, $(length(parts)) chapters), version $version_id")
if same
    # same content, rebuilt page (e.g. a new page layout): remember its hash so it is not archived as unrecorded
    let f = joinpath(HIST, version_id, "version.toml"), r = TOML.parsefile(f), sha = bytes2hex(sha256(final))
        sha == get(r, "report_sha", "") || sha in get(r, "report_shas", String[]) ||
            (r["report_shas"] = vcat(get(r, "report_shas", String[]), sha); open(io -> TOML.print(io, r), f, "w"))
    end
    println("Unchanged since version $version_id: not archived again")
else
    archive_version(version_id, final; settings=settings, results=results, chapters=chapter_stamps, source="build_final_report.jl", content_sha=content_sha)
    println("Archived as docs/outputs/final_report/history/$version_id/")
end
write_history_index()
println("Wrote docs/outputs/final_report/history/index.html ($(length(load_versions())) versions)")

# the same top bar on every main page (docs/assets/site-nav.js; user, 2026-11-01). Regenerated reports lose the
# script tag, so it is put back here, as docs/assets/add_site_nav.py does; pages that have it are left alone.
const NAV_PAGES = ["index.html", "papers/index.html", "papers/from-culture-to-control/index.html",
    "outputs/results/summary.html", "outputs/results/presentation.html", "outputs/results/explorer.html",
    "outputs/final_report/report.html", "outputs/final_report/history/index.html", "outputs/time_keeping/index.html",
    (joinpath(PKG, d, "report.html") for d in ("control", "control_30d", "control_60d", "control_120d", "control_sweep", "findings", "zoo_fits",
        "zoo_analyses/bootstrap", "zoo_analyses/sensitivity", "zoo_analyses/validation", "zoo_analyses/resistant_kill"))...]
function add_site_nav(rel)
    f = joinpath(DOCS, rel); isfile(f) || return false
    html = read(f, String)
    occursin("site-nav.js", first(html, 20000)) && return false
    src = replace(relpath(joinpath(DOCS, "assets", "site-nav.js"), dirname(f)), '\\' => '/')
    tag = "<script defer src=\"$src\"></script>\n"
    m = findfirst(r"</head\s*>"i, html)
    write(f, m === nothing ? tag * html : html[1:prevind(html, first(m))] * tag * html[first(m):end])
    true
end
println("Site navigation bar added to $(count(add_site_nav, NAV_PAGES)) page(s)")
isempty(missing_parts) || println("Missing: ", join(missing_parts, ", "))
