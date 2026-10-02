# Minimal static HTML report writer (same spirit as the 02 reports: tables + figures).

const OUT_ROOT = get(ENV, "NN_OUT", normpath(joinpath(@__DIR__, "..", "outputs")))
outdir(stage) = mkpath(joinpath(OUT_ROOT, stage))

_esc(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
_fmt(x::AbstractFloat) = isfinite(x) ? @sprintf("%.4g", x) : string(x)
_fmt(x) = _esc(x)

function html_table(df::DataFrame)
    head = join(("<th>$(_esc(n))</th>" for n in names(df)))
    body = join(("<tr>" * join(("<td>$(_fmt(v))</td>" for v in row)) * "</tr>" for row in eachrow(df)))
    return "<div class=\"tw\"><table><thead><tr>$head</tr></thead><tbody>$body</tbody></table></div>"
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
    code{font-size:13px}
    """
    parts = ["<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
             "<title>$(_esc(title))</title><style>$css</style></head><body><h1>$(_esc(title))</h1>"]
    for (h, body) in sections
        push!(parts, "<h2>$(_esc(h))</h2>$body")
    end
    push!(parts, "</body></html>")
    write(path, join(parts))
    return path
end

figure(src, caption = "") = "<figure><img src=\"$(src)\" alt=\"$(_esc(caption))\"><figcaption>$(_esc(caption))</figcaption></figure>"
para(s) = "<p>$s</p>"
