using CSV
using XLSX
using DataFrames
using Statistics
using Printf
using Dates

const ROOT = @__DIR__
const DATA_DIR = joinpath(ROOT, "New Datasets")
const OUT_DIR = joinpath(ROOT, "outputs", "new_data_report")
const DOCS_OUT_DIR = joinpath(ROOT, "docs", "outputs", "new_data_report")

html_escape(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
slug(x) = replace(lowercase(string(x)), r"[^a-z0-9]+" => "-")

function num(x)
    x isa Missing && return missing
    x isa Number && return Float64(x)
    s = strip(string(x))
    isempty(s) && return missing
    occursin("#", s) && return missing
    try
        return parse(Float64, replace(s, "," => ""))
    catch
        return missing
    end
end

function tidy(df)
    rename!(df, Symbol.(strip.(string.(names(df)))))
    required = [:Day, Symbol("Mean Cells")]
    available = propertynames(df)
    all(c -> c in available, required) || return nothing
    for c in (:Day, Symbol("Mean Cells"), Symbol("SD Cells"), Symbol("SEM Cells"), Symbol("N Samples"))
        c in available && (df[!, c] = [num(v) for v in df[!, c]])
    end
    filter!(r -> !ismissing(r.Day) && !ismissing(r[Symbol("Mean Cells")]), df)
    isempty(df) ? nothing : df
end

function matrix_df(matrix)
    isempty(matrix) && return nothing
    raw_header = strip.(string.(matrix[1, :]))
    keep = [!isempty(h) && lowercase(h) != "missing" for h in raw_header]
    header = Symbol.(raw_header[keep])
    seen = Dict{Symbol,Int}()
    header = [begin
        k = get!(seen, h, 0) + 1
        seen[h] = k
        k == 1 ? h : Symbol("$(h)_$(k)")
    end for h in header]
    rows = Vector{Vector{Any}}()
    for i in 2:size(matrix, 1)
        vals = Any[matrix[i, j] for j in findall(keep)]
        all(v -> v === missing || isempty(strip(string(v))), vals) && continue
        push!(rows, vals)
    end
    isempty(rows) && return nothing
    df = DataFrame([getindex.(rows, j) for j in eachindex(header)], header)
    tidy(df)
end

function kind_and_notes(file, sheet)
    f = lowercase(file)
    s = lowercase(sheet)
    if startswith(f, "ce0_")
        return ("Co-culture untreated (Ce0)", "Low-resource co-culture; Ce0 = untreated.")
    elseif startswith(f, "ce1_")
        return ("Co-culture treated (Ce1, 1 uM IC50)", "Low-resource co-culture; Ce1 = 1 uM cisplatin treatment at the IC50.")
    elseif startswith(f, "all ") && (occursin("1-1", s) || occursin("1-3", s) || occursin("3-1", s))
        return ("Co-culture low-resource", "Low-resource co-culture; ratio is sensitive:resistant.")
    elseif occursin("untreated", f)
        return ("Mono-culture untreated", "Mono-culture at 30,000 cells/mL; untreated.")
    else
        return ("Mono-culture treated", "Mono-culture at 30,000 cells/mL; treatment is stated in the filename.")
    end
end

function record(section, title, notes, sources)
    (section=section, title=title, notes=notes, sources=sources)
end

function series(file, sheet, df, label)
    (file=file, sheet=sheet, df=df, label=label)
end

function ratio_label(sheet)
    m = match(r"(\d-\d)", sheet)
    m === nothing ? "matched coculture" : "sensitive:resistant $(m.captures[1])"
end

function counterpart(sheet)
    lower = lowercase(sheet)
    if occursin("a2780cis", lower)
        return replace(sheet, "A2780cis" => "A2780")
    elseif occursin("a2780", lower)
        return replace(sheet, "A2780" => "A2780cis")
    elseif occursin("tyknu-cpr", lower)
        return replace(sheet, "Tyknu-cpr" => "Tyknu")
    elseif occursin("tyknu", lower)
        return replace(sheet, "Tyknu" => "Tyknu-cpr")
    end
    nothing
end

function read_sources()
    records = NamedTuple[]
    files = sort(filter(f -> endswith(lowercase(f), ".csv") || endswith(lowercase(f), ".xlsx"), readdir(DATA_DIR)))
    for file in files
        path = joinpath(DATA_DIR, file)
        if endswith(lowercase(file), ".csv")
            df = tidy(CSV.read(path, DataFrame; silencewarnings=true))
            df === nothing && continue
            untreated = occursin("untreated", lowercase(file))
            section = untreated ? "untreated" : "treated"
            kind = untreated ? "Mono-culture untreated" : "Mono-culture treated"
            notes = untreated ? "Mono-culture at 30,000 cells/mL; untreated." : "Mono-culture at 30,000 cells/mL; treatment is stated in the filename."
            push!(records, record(section, "$(kind): $(file) / CSV", notes, [series(file, "CSV", df, "$(file) / CSV")]))
        else
            book = XLSX.readxlsx(path)
            sheets = XLSX.sheetnames(book)
            parsed = Dict{String,DataFrame}()
            for sheet in sheets
                df = matrix_df(book[sheet][:])
                df !== nothing && (parsed[sheet] = df)
            end
            if startswith(lowercase(file), "ce0_") || startswith(lowercase(file), "ce1_")
                isempty(parsed) && continue
                state = startswith(lowercase(file), "ce0_") ? "untreated (Ce0)" : "treated (Ce1, 1 uM IC50)"
                notes = "Low-resource co-culture; Ce0 = untreated and Ce1 = 1 uM cisplatin at the IC50."
                src = [series(file, s, d, s) for (s, d) in sort(collect(parsed); by=first)]
                push!(records, record("coculture", "Co-culture $(state): $(file) / sensitive + resistant sheets", notes, src))
            else
                used = Set{String}()
                for sheet in sheets
                    haskey(parsed, sheet) && !(sheet in used) || continue
                    df = parsed[sheet]
                    if occursin("mono", lowercase(sheet))
                        push!(records, record("low_resource", "Low-resource mono-culture: $(file) / $(sheet)", "Low-resource mono-culture source.", [series(file, sheet, df, sheet)]))
                        push!(used, sheet)
                    else
                        other = counterpart(sheet)
                        paired = other !== nothing && haskey(parsed, other) && !(other in used)
                        src = paired ? [series(file, sheet, df, sheet), series(file, other, parsed[other], other)] : [series(file, sheet, df, sheet)]
                        label = paired ? "$(ratio_label(sheet)); sensitive/naive + resistant/cis" : "$(ratio_label(sheet))"
                        push!(records, record("low_resource", "Low-resource co-culture: $(file) / $(label)", "Low-resource co-culture; ratio is sensitive:resistant. The sensitive/naive and resistant/cis sheets are plotted together.", src))
                        push!(used, sheet)
                        paired && push!(used, other)
                    end
                end
            end
        end
    end
    order = Dict("untreated" => 1, "treated" => 2, "coculture" => 3, "low_resource" => 4)
    sort(records; by=r -> (order[r.section], r.title))
end

function svg_plot(sources, title)
    curves = NamedTuple[]
    for (i, src) in enumerate(sources)
        x = Float64.(collect(skipmissing(src.df.Day)))
        y = Float64.(collect(skipmissing(src.df[!, Symbol("Mean Cells")])))
        n = min(length(x), length(y)); n == 0 && continue
        x = x[1:n]; y = y[1:n]
        errcol = Symbol("SEM Cells") in propertynames(src.df) ? Symbol("SEM Cells") : Symbol("SD Cells")
        e = errcol in propertynames(src.df) ? [ismissing(v) ? 0.0 : Float64(v) for v in src.df[1:n, errcol]] : zeros(n)
        push!(curves, (x=x, y=y, e=e, label=src.label, color=["#176b87", "#d97706", "#8b5cf6", "#dc2626"][mod1(i, 4)]))
    end
    isempty(curves) && return "<div class=\"empty\">No numeric trajectory available.</div>"
    xmin = minimum(minimum(c.x) for c in curves); xmax = maximum(maximum(c.x) for c in curves); xmax == xmin && (xmax = xmin + 1)
    ymin = 0.0; ymax = max(maximum(maximum(c.y .+ c.e) for c in curves) * 1.12, 1.0)
    W, H, L, R, T, B = 760, 300, 72, 22, 28, 48
    sx(v) = L + (v - xmin) / (xmax - xmin) * (W - L - R)
    sy(v) = H - B - (v - ymin) / (ymax - ymin) * (H - T - B)
    lines = String[]
    bars = String[]
    dots = String[]
    legend = String[]
    for c in curves
        points = join([@sprintf("%.1f,%.1f", sx(c.x[i]), sy(c.y[i])) for i in eachindex(c.x)], " ")
        push!(lines, "<polyline fill=\"none\" stroke=\"$(c.color)\" stroke-width=\"2.5\" points=\"$points\"/>")
        push!(bars, join([@sprintf("<line stroke=\"%s\" stroke-width=\"1.2\" x1=\"%.1f\" x2=\"%.1f\" y1=\"%.1f\" y2=\"%.1f\"/><line stroke=\"%s\" stroke-width=\"1.2\" x1=\"%.1f\" x2=\"%.1f\" y1=\"%.1f\" y2=\"%.1f\"/>", c.color, sx(c.x[i]), sx(c.x[i]), sy(c.y[i]-c.e[i]), sy(c.y[i]+c.e[i]), c.color, sx(c.x[i])-3, sx(c.x[i])+3, sy(c.y[i]-c.e[i]), sy(c.y[i]-c.e[i])) for i in eachindex(c.x)], ""))
        push!(dots, join([@sprintf("<circle fill=\"%s\" cx=\"%.1f\" cy=\"%.1f\" r=\"3\"/>", c.color, sx(c.x[i]), sy(c.y[i])) for i in eachindex(c.x)], ""))
        push!(legend, "<span style=\"color:$(c.color)\">&#9679;</span> $(html_escape(c.label))")
    end
    xticks = join([@sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\">%.0f</text>", sx(v), H-12, v) for v in range(xmin, xmax; length=5)], "")
    yticks = join([@sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\">%.0f</text><line x1=\"%d\" x2=\"%d\" y1=\"%.1f\" y2=\"%.1f\"/>", L-8, sy(v)+4, v, L, W-R, sy(v), sy(v)) for v in range(ymin, ymax; length=5)], "")
    "<svg viewBox=\"0 0 $W $H\" role=\"img\" aria-label=\"$(html_escape(title))\"><style>text{font:12px sans-serif;fill:#344054}.axis{stroke:#667085;stroke-width:1}.legend{font:11px sans-serif}.legend span{font-size:15px}</style><text x=\"$(W/2)\" y=\"16\" text-anchor=\"middle\" font-weight=\"700\">$(html_escape(title))</text><text x=\"$(W/2)\" y=\"$H\" text-anchor=\"middle\">Day</text><text transform=\"translate(14 $(H/2)) rotate(-90)\" text-anchor=\"middle\">Mean Cells</text>$yticks<line class=\"axis\" x1=\"$L\" x2=\"$(W-R)\" y1=\"$(H-B)\" y2=\"$(H-B)\"/><line class=\"axis\" x1=\"$L\" x2=\"$L\" y1=\"$T\" y2=\"$(H-B)\"/>$xticks$(join(lines, ""))$(join(bars, ""))$(join(dots, ""))<foreignObject x=\"$(W-330)\" y=\"25\" width=\"310\" height=\"60\"><div xmlns=\"http://www.w3.org/1999/xhtml\" class=\"legend\">$(join(legend, "<br>"))</div></foreignObject></svg>"
end

function preview(df)
    cols = names(df)[1:min(5, ncol(df))]
    head = first(df, min(3, nrow(df)))
    th = join(["<th>$(html_escape(c))</th>" for c in cols], "")
    row_html(row) = join(["<td>$(html_escape(row[c]))</td>" for c in cols], "")
    rows = join(["<tr>$(row_html(row))</tr>" for row in eachrow(head)], "")
    "<table class=\"preview\"><thead><tr>$th</tr></thead><tbody>$rows</tbody></table>"
end

function preview_toggle(df, label)
    """
    <details class="preview-toggle">
      <summary>Expand preview for $(html_escape(label))</summary>
      <div class="preview-shell">$(preview(df))</div>
    </details>
    """
end

function ensure_back_buttons(dir)
    isdir(dir) || return
    link = "<p><a class=\"back-home\" href=\"../../../../index.html\" style=\"display:inline-block;margin:0 0 16px;color:#176b87;font-weight:700;text-decoration:none\">&larr; Back to reports home</a></p>"
    for path in readdir(dir; join=true)
        endswith(lowercase(path), ".html") || continue
        html = read(path, String)
        occursin("reports home", lowercase(html)) && continue
        body_tag = match(r"<body[^>]*>", html)
        body_tag === nothing && continue
        pos = body_tag.offset + ncodeunits(body_tag.match)
        html = html[1:pos-1] * link * html[pos:end]
        write(path, html)
    end
end

function main()
    mkpath(OUT_DIR)
    records = read_sources()
    cards_by_section = Dict{String,Vector{String}}()
    for (i, r) in enumerate(records)
        previews = join(["<h4>Sheet preview: $(html_escape(s.label))</h4>$(preview_toggle(s.df, s.label))" for s in r.sources], "")
        push!(get!(cards_by_section, r.section, String[]), "<article class=\"card section-$(r.section)\"><h3>$(html_escape(r.title))</h3><p class=\"meta\">$(html_escape(r.notes))</p>$(svg_plot(r.sources, r.title))$previews</article>")
    end
    files = Set(s.file for r in records for s in r.sources)
    manifest = "{\n  \"generated_by\": \"Julia $(VERSION)\",\n  \"generated_at\": \"$(Dates.now())\",\n  \"file_count\": $(length(files)),\n  \"plot_count\": $(length(records)),\n  \"sheet_count\": $(sum(length(r.sources) for r in records))\n}\n"
    write(joinpath(OUT_DIR, "manifest.json"), manifest)
    section_titles = Dict("untreated" => "1. Untreated Mono-culture", "treated" => "2. Treated Mono-culture", "coculture" => "3. Co-culture", "low_resource" => "4. Low-Resource Data")
    section_order = ["untreated", "treated", "coculture", "low_resource"]
    body = join(["<section class=\"report-section\"><h2>$(section_titles[s])</h2>$(join(get(cards_by_section, s, String[]), "\\n"))</section>" for s in section_order], "\n")
    html = """<!doctype html><html><head><meta charset=\"utf-8\"><title>New Data Report</title><style>
body{margin:0;background:#f4f7f8;color:#172b36;font:15px system-ui,sans-serif}main{max-width:1280px;margin:auto;padding:28px}h1{margin-bottom:8px}h2{margin-top:30px}.intro,.workbook{background:white;border:1px solid #d9e2e6;border-radius:14px;padding:20px;box-shadow:0 5px 18px #173b4d0d}.intro{border-left:6px solid #176b87}.card{margin-top:16px;border-top:1px solid #e6ecef;padding-top:16px}.meta{color:#52636b}.back{display:inline-block;margin:0 0 18px;color:#176b87;font-weight:700;text-decoration:none}.preview-toggle{margin:10px 0 0}.preview-toggle>summary{display:inline-flex;align-items:center;gap:8px;cursor:pointer;list-style:none;padding:9px 13px;border:1px solid #c5d1d8;border-radius:999px;background:#f7fbfc;color:#175e79;font-weight:700}.preview-toggle>summary::-webkit-details-marker{display:none}.preview-toggle>summary::marker{content:\"\"}.preview-toggle[open]>summary{background:#e8f3f6;border-color:#9fc4d0}.preview-shell{margin-top:10px}.preview{border-collapse:collapse;width:100%;font-size:12px}.preview th,.preview td{border:1px solid #d9e2e6;padding:5px;text-align:right}.preview th{background:#eef5f6;text-align:left}.empty{padding:30px;background:#fff7ed;border-radius:8px}svg{width:100%;max-height:300px;background:#fbfdfd;border:1px solid #e4e7ec;border-radius:8px}code{background:#eef5f6;padding:2px 5px;border-radius:4px}</style></head><body><main><a class=\"back\" href=\"../../index.html\">&larr; Back to reports home</a><h1>New Data Report</h1><div class=\"intro\"><p>This report was generated entirely in Julia $(VERSION) from every CSV and XLSX source in <code>New Datasets</code>.</p><p><b>How to read the labels:</b> LR means low resource. Ce0 is untreated co-culture. Ce1 is treated co-culture at the 1 uM IC50. Ratios 1-1, 3-1, and 1-3 are co-culture sensitive:resistant ratios. Mono-culture files are at 30,000 cells/mL, with treatment specified in the filename. Each chart includes its source file and sheet name, and each card includes a three-row preview.</p><p>Plotted values are Mean Cells over Day. Error bars use SEM when available, otherwise SD.</p></div><h2>All source sheets ($(length(records)))</h2>$body</main></body></html>"""
    write(joinpath(OUT_DIR, "report.html"), html)
    mkpath(DOCS_OUT_DIR)
    cp(joinpath(OUT_DIR, "report.html"), joinpath(DOCS_OUT_DIR, "report.html"); force=true)
    cp(joinpath(OUT_DIR, "manifest.json"), joinpath(DOCS_OUT_DIR, "manifest.json"); force=true)
    ensure_back_buttons(joinpath(ROOT, "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports"))
    ensure_back_buttons(joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports"))
    println("Generated $(length(records)) plots from $(length(files)) files")
    println(joinpath(OUT_DIR, "report.html"))
end

main()
