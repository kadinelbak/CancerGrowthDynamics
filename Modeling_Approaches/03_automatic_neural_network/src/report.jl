# Minimal static HTML report writer (same spirit as the 02 reports: tables + figures).
# Every table and figure carries an "i" tip that says how to read it and which direction is better
# (the same component as 02_mechanical_automatic_package/info_tips.jl). `html_table` and `figure` require
# `tip = (text, better)`, where better is :lower, :higher, :none (descriptive) or a phrase.

const INFO_TIPS_HEAD = raw"""
<style>
.ii{display:inline-flex;align-items:center;justify-content:center;width:17px;height:17px;margin:0 4px;padding:0;border-radius:50%;border:1.5px solid #176b87;background:#fff;color:#176b87;font:italic 700 11px/1 Georgia,serif;cursor:help;vertical-align:middle;flex:none}
.ii:hover,.ii:focus-visible,.ii[aria-expanded=true]{background:#176b87;color:#fff;outline:none}
#ii-pop{position:fixed;z-index:9999;max-width:min(400px,calc(100vw - 24px));background:#17212b;color:#f4f7f8;font:13px/1.45 system-ui,sans-serif;font-style:normal;font-weight:400;text-align:left;padding:10px 12px;border-radius:8px;box-shadow:0 6px 20px rgba(0,0,0,.25);display:none;white-space:normal}
</style>
<script>
(function () {
  if (window.__iiInit) return; window.__iiInit = true;
  var pop = null, cur = null, pinned = false;
  function place(b) {
    if (!pop) { pop = document.createElement('div'); pop.id = 'ii-pop'; pop.setAttribute('role', 'tooltip'); document.body.appendChild(pop); }
    pop.textContent = b.getAttribute('data-tip');
    var bt = b.getAttribute('data-better');
    if (bt) { var d = document.createElement('div'); d.style.cssText = 'margin-top:6px;font-weight:700;color:#ffd479'; d.textContent = 'What is better: ' + bt; pop.appendChild(d); }
    pop.style.display = 'block';
    var r = b.getBoundingClientRect(), w = pop.offsetWidth, h = pop.offsetHeight;
    var x = Math.min(Math.max(8, r.left + r.width / 2 - w / 2), window.innerWidth - w - 8);
    var y = r.bottom + 8; if (y + h > window.innerHeight - 8) y = Math.max(8, r.top - h - 8);
    pop.style.left = x + 'px'; pop.style.top = y + 'px';
    if (cur && cur !== b) cur.setAttribute('aria-expanded', 'false');
    cur = b; b.setAttribute('aria-expanded', 'true');
  }
  function hide() { if (pop) pop.style.display = 'none'; if (cur) cur.setAttribute('aria-expanded', 'false'); cur = null; pinned = false; }
  document.addEventListener('click', function (e) {
    var b = e.target.closest ? e.target.closest('.ii') : null;
    if (b) { e.preventDefault(); e.stopPropagation(); if (cur === b && pinned) hide(); else { place(b); pinned = true; } }
    else if (!(pop && pop.contains(e.target))) hide();
  }, true);
  document.addEventListener('mouseover', function (e) { var b = e.target.closest ? e.target.closest('.ii') : null; if (b && !pinned) place(b); });
  document.addEventListener('mouseout', function (e) { var b = e.target.closest ? e.target.closest('.ii') : null; if (b && !pinned && !b.contains(e.relatedTarget)) hide(); });
  document.addEventListener('focusin', function (e) { var b = e.target.closest ? e.target.closest('.ii') : null; if (b) place(b); });
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') hide(); });
  window.addEventListener('scroll', function () { if (!pinned) hide(); else if (cur) place(cur); }, true);
})();
</script>
"""

tip_esc(s) = replace(string(s), "&" => "&amp;", "\"" => "&quot;", "<" => "&lt;", ">" => "&gt;")
const BETTER_TEXT = Dict(:lower => "Lower is better.", :higher => "Higher is better.",
                         :none => "Neither: this describes the model or data; it is not a score.")
better_text(b::Symbol) = BETTER_TEXT[b]
better_text(b::AbstractString) = b
"The \"i\" button: `text` explains the item, `better` says which direction is better."
tip(text, better) = "<button type=\"button\" class=\"ii\" aria-label=\"How to read this\" data-tip=\"$(tip_esc(text))\" data-better=\"$(tip_esc(better_text(better)))\">i</button>"

const OUT_ROOT = get(ENV, "NN_OUT", normpath(joinpath(@__DIR__, "..", "outputs")))
outdir(stage) = mkpath(joinpath(OUT_ROOT, stage))

_esc(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
_fmt(x::AbstractFloat) = isfinite(x) ? @sprintf("%.4g", x) : string(x)
_fmt(x) = _esc(x)

"Table with a mandatory info tip `tip = (text, better)`."
function html_table(df::DataFrame; tip::Tuple)
    head = join(("<th>$(_esc(n))</th>" for n in names(df)))
    body = join(("<tr>" * join(("<td>$(_fmt(v))</td>" for v in row)) * "</tr>" for row in eachrow(df)))
    return "<div class=\"tl\">$(NNTrack.tip(tip...))<span class=\"note\">How to read this table</span></div>" *
           "<div class=\"tw\"><table><thead><tr>$head</tr></thead><tbody>$body</tbody></table></div>"
end

"""
write_report(path, title, sections) where each section is (heading, html_body).
Figures are referenced by relative path with `figure(src, caption)`.
"""
function write_report(path, title, sections)
    css = """
    :root{--bg:#fbfbfa;--fg:#1f2328;--mut:#5b6470;--line:#d9dde3;--acc:#2f6f9f}
    @media (prefers-color-scheme: dark){:root{--bg:#16181c;--fg:#e6e8eb;--mut:#9aa3ad;--line:#343a42;--acc:#7fb4dc}}
    body{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,sans-serif;margin:0 auto;max-width:1100px;padding:24px 16px}
    h1{font-size:1.6rem}h2{font-size:1.2rem;margin-top:2rem;border-bottom:1px solid var(--line);padding-bottom:4px}
    p,li{color:var(--fg)} .note{color:var(--mut)} .tw{overflow-x:auto}
    table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{border:1px solid var(--line);padding:3px 8px;text-align:right}
    th{background:color-mix(in srgb,var(--acc) 12%,transparent)} td:first-child,th:first-child{text-align:left}
    figure{margin:12px 0}img{max-width:100%;border:1px solid var(--line);background:#fff}figcaption{color:var(--mut);font-size:13px}
    code{font-size:13px} .tl{display:flex;align-items:center;gap:2px;margin-top:10px;font-size:13px}
    """
    parts = ["<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
             "<title>$(_esc(title))</title><style>$css</style>$INFO_TIPS_HEAD</head><body><h1>$(_esc(title))</h1>"]
    for (h, body) in sections
        push!(parts, "<h2>$(_esc(h))</h2>$body")
    end
    push!(parts, "</body></html>")
    write(path, join(parts))
    return path
end

"Figure with a mandatory info tip `tip = (text, better)`."
figure(src, caption = ""; tip::Tuple) = "<figure><img src=\"$(src)\" alt=\"$(_esc(caption))\"><figcaption>$(NNTrack.tip(tip...))$(_esc(caption))</figcaption></figure>"
para(s) = "<p>$s</p>"
