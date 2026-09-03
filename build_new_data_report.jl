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

function read_sources()
    records = NamedTuple[]
    files = sort(filter(f -> endswith(lowercase(f), ".csv") || endswith(lowercase(f), ".xlsx"), readdir(DATA_DIR)))
    for file in files
        path = joinpath(DATA_DIR, file)
        if endswith(lowercase(file), ".csv")
            df = tidy(CSV.read(path, DataFrame; silencewarnings=true))
            df === nothing && continue
            kind, notes = kind_and_notes(file, "CSV")
            push!(records, (file=file, sheet="CSV", kind=kind, notes=notes, df=df))
        else
            book = XLSX.readxlsx(path)
            for sheet in XLSX.sheetnames(book)
                df = matrix_df(book[sheet][:])
                df === nothing && continue
                kind, notes = kind_and_notes(file, sheet)
                push!(records, (file=file, sheet=sheet, kind=kind, notes=notes, df=df))
            end
        end
    end
    records
end

function svg_plot(df, title, subtitle)
    x = Float64.(collect(skipmissing(df.Day)))
    y = Float64.(collect(skipmissing(df[!, Symbol("Mean Cells")])))
    n = min(length(x), length(y)); x = x[1:n]; y = y[1:n]
    isempty(x) && return "<div class=\"empty\">No numeric trajectory available.</div>"
    errcol = Symbol("SEM Cells") in propertynames(df) ? Symbol("SEM Cells") : Symbol("SD Cells")
    e = errcol in propertynames(df) ? [ismissing(v) ? 0.0 : Float64(v) for v in df[1:n, errcol]] : zeros(n)
    W, H, L, R, T, B = 760, 300, 72, 22, 28, 48
    xmin, xmax = extrema(x); xmax == xmin && (xmax = xmin + 1)
    ymin = 0.0; ymax = max(maximum(y .+ e) * 1.12, 1.0)
    sx(v) = L + (v - xmin) / (xmax - xmin) * (W - L - R)
    sy(v) = H - B - (v - ymin) / (ymax - ymin) * (H - T - B)
    points = join([@sprintf("%.1f,%.1f", sx(x[i]), sy(y[i])) for i in eachindex(x)], " ")
    bars = join([@sprintf("<line x1=\"%.1f\" x2=\"%.1f\" y1=\"%.1f\" y2=\"%.1f\"/><line x1=\"%.1f\" x2=\"%.1f\" y1=\"%.1f\" y2=\"%.1f\"/>", sx(x[i]), sx(x[i]), sy(y[i]-e[i]), sy(y[i]+e[i]), sx(x[i])-3, sx(x[i])+3, sy(y[i]-e[i]), sy(y[i]-e[i])) for i in eachindex(x)], "")
    xticks = join([@sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\">%.0f</text>", sx(v), H-12, v) for v in range(xmin, xmax; length=5)], "")
    yticks = join([@sprintf("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\">%.0f</text><line x1=\"%d\" x2=\"%d\" y1=\"%.1f\" y2=\"%.1f\"/>", L-8, sy(v)+4, v, L, W-R, sy(v), sy(v)) for v in range(ymin, ymax; length=5)], "")
    "<svg viewBox=\"0 0 $W $H\" role=\"img\" aria-label=\"$(html_escape(title))\"><style>text{font:12px sans-serif;fill:#344054} .grid{stroke:#e4e7ec;stroke-width:1} .axis{stroke:#667085;stroke-width:1} .err{stroke:#d97706;stroke-width:1.2} .line{fill:none;stroke:#176b87;stroke-width:2.5} .dot{fill:#176b87}</style><text x=\"$(W/2)\" y=\"16\" text-anchor=\"middle\" font-weight=\"700\">$(html_escape(title))</text><text x=\"$(W/2)\" y=\"$H\" text-anchor=\"middle\">Day</text><text transform=\"translate(14 $(H/2)) rotate(-90)\" text-anchor=\"middle\">Mean Cells</text>$yticks<line class=\"axis\" x1=\"$L\" x2=\"$(W-R)\" y1=\"$(H-B)\" y2=\"$(H-B)\"/><line class=\"axis\" x1=\"$L\" x2=\"$L\" y1=\"$T\" y2=\"$(H-B)\"/>$xticks<polyline class=\"line\" points=\"$points\"/>$bars$(join([@sprintf("<circle class=\\\"dot\\\" cx=\\\"%.1f\\\" cy=\\\"%.1f\\\" r=\\\"3\\\"/>", sx(x[i]), sy(y[i])) for i in eachindex(x)], ""))</svg>"
end

function preview(df)
    cols = names(df)[1:min(5, ncol(df))]
    head = first(df, min(3, nrow(df)))
    th = join(["<th>$(html_escape(c))</th>" for c in cols], "")
    row_html(row) = join(["<td>$(html_escape(row[c]))</td>" for c in cols], "")
    rows = join(["<tr>$(row_html(row))</tr>" for row in eachrow(head)], "")
    "<table class=\"preview\"><thead><tr>$th</tr></thead><tbody>$rows</tbody></table>"
end

function main()
    mkpath(OUT_DIR)
    records = read_sources()
    cards = String[]
    for (i, r) in enumerate(records)
        title = "$(r.kind): $(r.file) / $(r.sheet)"
        push!(cards, "<article class=\"card\"><h3>$(html_escape(title))</h3><p class=\"meta\">$(html_escape(r.notes))</p>$(svg_plot(r.df, title, r.notes))<h4>Sheet preview</h4>$(preview(r.df))</article>")
    end
    groups = Dict{String,Vector{String}}()
    for r in records
        push!(get!(groups, r.file, String[]), r.sheet)
    end
    manifest = "{\n  \"generated_by\": \"Julia $(VERSION)\",\n  \"generated_at\": \"$(Dates.now())\",\n  \"file_count\": $(length(groups)),\n  \"sheet_count\": $(length(records))\n}\n"
    write(joinpath(OUT_DIR, "manifest.json"), manifest)
    body = join(cards, "\n")
    html = """<!doctype html><html><head><meta charset=\"utf-8\"><title>New Data Report</title><style>
body{margin:0;background:#f4f7f8;color:#172b36;font:15px system-ui,sans-serif}main{max-width:1280px;margin:auto;padding:28px}h1{margin-bottom:8px}h2{margin-top:30px}.intro,.workbook{background:white;border:1px solid #d9e2e6;border-radius:14px;padding:20px;box-shadow:0 5px 18px #173b4d0d}.intro{border-left:6px solid #176b87}.card{margin-top:16px;border-top:1px solid #e6ecef;padding-top:16px}.meta{color:#52636b}.preview{border-collapse:collapse;width:100%;font-size:12px}.preview th,.preview td{border:1px solid #d9e2e6;padding:5px;text-align:right}.preview th{background:#eef5f6;text-align:left}.empty{padding:30px;background:#fff7ed;border-radius:8px}svg{width:100%;max-height:300px;background:#fbfdfd;border:1px solid #e4e7ec;border-radius:8px}code{background:#eef5f6;padding:2px 5px;border-radius:4px}</style></head><body><main><h1>New Data Report</h1><div class=\"intro\"><p>This report was generated entirely in Julia $(VERSION) from every CSV and XLSX source in <code>New Datasets</code>.</p><p><b>How to read the labels:</b> LR means low resource. Ce0 is untreated co-culture. Ce1 is treated co-culture at the 1 uM IC50. Ratios 1-1, 3-1, and 1-3 are co-culture sensitive:resistant ratios. Mono-culture files are at 30,000 cells/mL, with treatment specified in the filename. Each chart includes its source file and sheet name, and each card includes a three-row preview.</p><p>Plotted values are Mean Cells over Day. Error bars use SEM when available, otherwise SD.</p></div><h2>All source sheets ($(length(records)))</h2>$body</main></body></html>"""
    write(joinpath(OUT_DIR, "report.html"), html)
    mkpath(DOCS_OUT_DIR)
    cp(joinpath(OUT_DIR, "report.html"), joinpath(DOCS_OUT_DIR, "report.html"); force=true)
    cp(joinpath(OUT_DIR, "manifest.json"), joinpath(DOCS_OUT_DIR, "manifest.json"); force=true)
    println("Generated $(length(records)) plots from $(length(groups)) files")
    println(joinpath(OUT_DIR, "report.html"))
end

main()
