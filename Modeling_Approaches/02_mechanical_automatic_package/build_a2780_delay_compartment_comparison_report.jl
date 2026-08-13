using CSV
using DataFrames
using Dates
using Printf
using Plots

const ROOT = @__DIR__
const OUT_CSV = joinpath(ROOT, "outputs", "csv", "coculture_treated")
const OUT_IMG = joinpath(ROOT, "outputs", "images", "coculture_treated", "figures", "delay_compartment_comparison")
const OUT_REPORT = joinpath(ROOT, "outputs", "reports")

const DATA_COLORS = Dict(
    "sensitive" => :crimson,
    "resistant" => :steelblue,
)

const DATA_LABELS = Dict(
    "sensitive" => "Naive / sensitive component",
    "resistant" => "cis / resistant component",
)

const MODEL_COLORS = Dict(
    "time_delay" => :darkorange,
    "compartment" => :purple,
)

const MODEL_SPECS = [
    (
        key = "time_delay",
        title = "Best Time-Delay-Only Candidate",
        model = "dual_delayed_ramp_kill",
        pooling = "partial_5pct",
        source = "coculture_treated_joint_overlays.csv",
        ranking = "coculture_treated_pooling_model_ranking.csv",
        context = "",
        mechanism = "Delayed ramp kill: drug pressure turns on over time, without adding a damaged-visible compartment.",
        equation = raw"\frac{dX_i}{dt}=G_i^{\mathrm{co}}-A_i(t)H_i(z)X_i",
    ),
    (
        key = "compartment",
        title = "Best Transit / Compartment Candidate",
        model = "dual_transit_load_scaled",
        pooling = "shared",
        source = "coculture_treated_joint_overlays.csv",
        ranking = "coculture_treated_pooling_model_ranking.csv",
        context = "",
        mechanism = "Transit damage with load scaling: cells move through a damaged-visible state before clearance, and drug damage is scaled by effective load.",
        equation = raw"\frac{dP_i}{dt}=G_i^{\mathrm{co}}-M_iA_i(t)H_i(z)P_i,\qquad \frac{dD_i}{dt}=M_iA_i(t)H_i(z)P_i-k_{\mathrm{clear},i}D_i",
    ),
]

html_escape(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
as_text(value) = ismissing(value) ? "" : string(value)
text_column(values) = [as_text(value) for value in values]
legend_panel(density, mix) = density == "30k" && mix == "50-50"

function finite_metric(value)
    try
        parsed = Float64(value)
        return isfinite(parsed) && abs(parsed) < 1e11
    catch
        return false
    end
end

function metric_for(spec)
    path = joinpath(OUT_CSV, spec.ranking)
    df = CSV.read(path, DataFrame)
    rows = df[text_column(df.model) .== spec.model, :]
    if :pooling_mode in names(rows)
        rows = rows[text_column(rows.pooling_mode) .== spec.pooling, :]
    elseif :pooling in names(rows)
        rows = rows[text_column(rows.pooling) .== spec.pooling, :]
    end
    nrow(rows) == 0 && return (bic = NaN, delta_bic = NaN, n_parameters = missing)
    valid = rows[[finite_metric(value) for value in rows.bic], :]
    row = nrow(valid) > 0 ? valid[argmin(Float64.(valid.bic)), :] : first(rows)
    delta = :delta_bic in names(rows) ? tryparse(Float64, string(row.delta_bic)) : nothing
    return (
        bic = tryparse(Float64, string(row.bic)),
        delta_bic = delta === nothing ? NaN : delta,
        n_parameters = :n_parameters in names(rows) ? row.n_parameters : missing,
    )
end

function source_overlay(spec)
    path = joinpath(OUT_CSV, "figures", spec.source)
    df = CSV.read(path, DataFrame)
    mask = (text_column(df.model) .== spec.model) .& (text_column(df.pooling_mode) .== spec.pooling)
    if spec.context != "" && :context in names(df)
        mask .&= text_column(df.context) .== spec.context
    end
    selected = df[mask, :]
    nrow(selected) > 0 || error("No overlay rows for $(spec.model) | $(spec.pooling)")
    return selected
end

function panel_data(df, density, mix, component)
    rows = df[
        (text_column(df.density) .== density) .&
        (text_column(df.mix) .== mix) .&
        (text_column(df.component) .== component),
        :,
    ]
    sort!(rows, :time)
    return rows
end

function render_single_model_grid(spec, overlay)
    panels = Any[]
    model_color = MODEL_COLORS[spec.key]
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        panel = plot(
            title = "$(density), mix $(mix)",
            titlefontsize = 9,
            xlabel = "Time (day)",
            ylabel = mix == "25-75" ? "Measured population" : "",
            legend = false,
        )
        for component in ("sensitive", "resistant")
            rows = panel_data(overlay, density, mix, component)
            scatter!(
                panel,
                rows.time,
                rows.observed;
                color = DATA_COLORS[component],
                ms = 3.0,
                markerstrokewidth = 0,
                alpha = 0.78,
                label = DATA_LABELS[component],
            )
            plot!(
                panel,
                rows.time,
                rows.predicted;
                color = model_color,
                lw = component == "sensitive" ? 2.8 : 2.8,
                linestyle = component == "sensitive" ? :solid : :dash,
                label = component == "sensitive" ? "$(spec.title) fit" : "$(spec.title) fit, cis/resistant",
            )
        end
        push!(panels, panel)
    end
    figure = plot(
        panels...;
        layout = (2, 3),
        size = (1600, 860),
        plot_title = "$(spec.title): $(spec.model) | $(spec.pooling)",
        margin = 4 * Plots.mm,
    )
    path = joinpath(OUT_IMG, "$(spec.key)_grid.png")
    mkpath(dirname(path))
    savefig(figure, path)
    return path
end

function render_all_model_grid(overlays)
    panels = Any[]
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        panel = plot(
            title = "$(density), mix $(mix)",
            titlefontsize = 9,
            xlabel = "Time (day)",
            ylabel = mix == "25-75" ? "Measured population" : "",
            legend = false,
        )
        reference = overlays["time_delay"]
        for component in ("sensitive", "resistant")
            observed = panel_data(reference, density, mix, component)
            scatter!(
                panel,
                observed.time,
                observed.observed;
                color = DATA_COLORS[component],
                ms = 3.1,
                markerstrokewidth = 0,
                alpha = 0.75,
                label = DATA_LABELS[component],
            )
        end
        for spec in MODEL_SPECS
            overlay = overlays[spec.key]
            for component in ("sensitive", "resistant")
                rows = panel_data(overlay, density, mix, component)
                plot!(
                    panel,
                    rows.time,
                    rows.predicted;
                    color = MODEL_COLORS[spec.key],
                    lw = 2.5,
                    linestyle = component == "sensitive" ? :solid : :dash,
                    label = component == "sensitive" ? spec.title : "",
                )
            end
        end
        push!(panels, panel)
    end
    figure = plot(
        panels...;
        layout = (2, 3),
        size = (1700, 900),
        plot_title = "A2780 Treated Coculture: Time Delay vs Transit/Compartment Model",
        margin = 4 * Plots.mm,
    )
    path = joinpath(OUT_IMG, "delay_vs_compartment_comparison_grid.png")
    mkpath(dirname(path))
    savefig(figure, path)
    return path
end

function rel_report_path(path)
    return replace(relpath(path, OUT_REPORT), "\\" => "/")
end

function model_table(specs, metrics)
    rows = String[]
    for spec in specs
        metric = metrics[spec.key]
        bic = metric.bic isa Number && isfinite(metric.bic) ? @sprintf("%.3f", metric.bic) : "NA"
        delta = metric.delta_bic isa Number && isfinite(metric.delta_bic) ? @sprintf("%.3f", metric.delta_bic) : "NA"
        push!(rows, """
        <tr>
          <td>$(html_escape(spec.title))</td>
          <td><code>$(html_escape(spec.model))</code></td>
          <td><code>$(html_escape(spec.pooling))</code></td>
          <td>$(bic)</td>
          <td>$(delta)</td>
          <td>$(html_escape(metric.n_parameters))</td>
        </tr>
        """)
    end
    return join(rows, "\n")
end

function stage_walkthrough_html()
    return """
        <section class="stage-card">
          <h2>Stage 1: Untreated Monoculture Growth</h2>
          <p><strong>Question:</strong> how does each cell line grow when it is alone and untreated?</p>
          <p><strong>Inherited object:</strong> each lineage gets its own growth law, written as <code>G_N</code> for A2780Naive and <code>G_C</code> for A2780cis.</p>
          <div class="math">\\[G_i(X)=r_iX\\left[1-\\left(\\frac{X}{K_i}\\right)^{\\theta_i}\\right]\\]</div>
          <p><strong>Meaning:</strong> <code>r_i</code> controls early growth, <code>K_i</code> is the carrying capacity, and <code>theta_i</code> controls how sharply growth slows as the population approaches carrying capacity.</p>
        </section>

        <section class="stage-card">
          <h2>Stage 2: Treated Monoculture Timing</h2>
          <p><strong>Question:</strong> does cisplatin killing act immediately, gradually, after an onset, or through a visible-damage delay?</p>
          <p><strong>Inherited object:</strong> the drug-response strength <code>H_i(z)</code> and time-activation term <code>A_i(t)</code>.</p>
          <div class="math">\\[\\frac{dX_i}{dt}=G_i(X_i)-A_i(t)H_i(z)X_i\\]</div>
          <div class="math">\\[H_i(z)=\\frac{E_{\\max,i}z^{h_i}}{EC_{50,i}^{h_i}+z^{h_i}}\\]</div>
          <div class="math">\\[A_i(t)=\\begin{cases}0, & t\\le t_{\\mathrm{on},i} \\\\ 1-\\exp[-\\lambda_i(t-t_{\\mathrm{on},i})], & t>t_{\\mathrm{on},i}\\end{cases}\\]</div>
          <p><strong>Meaning:</strong> <code>H_i(z)</code> is the dose effect, while <code>A_i(t)</code> lets the effect turn on over time. A rising <code>A_i(t)</code> can produce early growth followed by later decline without adding a new visible compartment.</p>
        </section>

        <section class="stage-card">
          <h2>Stage 3: Untreated Coculture Competition</h2>
          <p><strong>Question:</strong> how do Naive and cis cells alter each other's growth when mixed without drug?</p>
          <p><strong>Inherited object:</strong> effective-load variables and competition coefficients.</p>
          <div class="math">\\[L_N=N+\\alpha_{NC}C,\\qquad L_C=C+\\alpha_{CN}N\\]</div>
          <div class="math">\\[\\frac{dN}{dt}=G_N(N,L_N)-d_NN,\\qquad \\frac{dC}{dt}=G_C(C,L_C)-d_CC\\]</div>
          <p><strong>Meaning:</strong> each lineage grows according to its own Stage-1 growth law, but growth is limited by effective coculture load. If cis cells decline, <code>L_N</code> falls and Naive cells can experience less competition.</p>
        </section>

        <section class="stage-card">
          <h2>Stage 4A: Time-Delay-Only Treated Coculture Candidate</h2>
          <p><strong>Question:</strong> can inherited growth, inherited competition, and a time-delayed drug effect explain the treated coculture data without a damaged-visible compartment?</p>
          <p><strong>Best candidate here:</strong> <code>dual_delayed_ramp_kill | partial_5pct</code>.</p>
          <div class="math">\\[\\frac{dN}{dt}=G_N(N,L_N)-A_N(t)H_N(z)N\\]</div>
          <div class="math">\\[\\frac{dC}{dt}=G_C(C,L_C)-A_C(t)H_C(z)C\\]</div>
          <div class="math">\\[A_i(t)=\\begin{cases}0, & t\\le t_{\\mathrm{on},i} \\\\ 1-\\exp[-\\lambda_i(t-t_{\\mathrm{on},i})], & t>t_{\\mathrm{on},i}\\end{cases}\\]</div>
          <p><strong>Meaning:</strong> the only special treated behavior is that drug killing is time-activated. The visible population is still modeled directly, so the hump must come from growth winning early and drug effect winning later.</p>
        </section>

        <section class="stage-card">
          <h2>Stage 4B: Transit / Compartment Treated Coculture Candidate</h2>
          <p><strong>Question:</strong> does the data prefer a visible-damage compartment, where cells are injured before disappearing from the measured population?</p>
          <p><strong>Best candidate here:</strong> <code>dual_transit_load_scaled | shared</code>.</p>
          <div class="math">\\[\\frac{dP_i}{dt}=G_i^{\\mathrm{co}}-M_iA_i(t)H_i(z)P_i\\]</div>
          <div class="math">\\[\\frac{dD_i}{dt}=M_iA_i(t)H_i(z)P_i-k_{\\mathrm{clear},i}D_i\\]</div>
          <div class="math">\\[Y_i(t)=P_i(t)+D_i(t),\\qquad M_i=\\exp\\left(\\beta_i\\frac{L_i}{K_i}\\right)\\]</div>
          <div class="math">\\[A_i(t)=\\begin{cases}0, & t\\le t_{\\mathrm{on},i} \\\\ 1-\\exp[-\\lambda_i(t-t_{\\mathrm{on},i})], & t>t_{\\mathrm{on},i}\\end{cases}\\]</div>
          <p><strong>Meaning:</strong> cells can stop proliferating or become drug-damaged before they are cleared from the measurement. The observed hump can therefore reflect delayed disappearance, not just delayed death onset.</p>
        </section>

        <section class="stage-card">
          <h2>Identifiability And Global Search</h2>
          <p><strong>Identifiability matrix idea:</strong> build a parameter-by-parameter view of whether two parameters can trade off while producing nearly the same fitted curves.</p>
          <div class="math">\\[J_{k,j}=\\frac{\\partial \\hat y_k}{\\partial p_j},\\qquad F=J^TWJ\\]</div>
          <div class="math">\\[\\operatorname{corr}(p_a,p_b)=\\frac{(F^{-1})_{ab}}{\\sqrt{(F^{-1})_{aa}(F^{-1})_{bb}}}\\]</div>
          <p><strong>Meaning:</strong> correlations near +1 or -1 mean two parameters are hard to distinguish. A near-singular matrix or broad profile likelihood means the data do not uniquely identify that mechanism.</p>
          <p><strong>Global search idea:</strong> yes, it makes sense to run broader multistart or global searches for these candidate models. Delay onset, ramp speed, kill strength, clearance, and load scaling can form ridges where several parameter combinations produce similar hump-shaped trajectories.</p>
        </section>
    """
end

function write_report(image_paths, metrics)
    mkpath(OUT_REPORT)
    generated = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    cards = String[]
    for spec in MODEL_SPECS
        push!(cards, """
        <section class="model-card">
          <h2>$(html_escape(spec.title))</h2>
          <p><strong>Model:</strong> <code>$(html_escape(spec.model))</code> | <strong>Pooling:</strong> <code>$(html_escape(spec.pooling))</code></p>
          <p><strong>Mechanism:</strong> $(html_escape(spec.mechanism))</p>
          <div class="math">\\[$(spec.equation)\\]</div>
          <figure>
            <img src="$(rel_report_path(image_paths[spec.key]))" alt="$(html_escape(spec.title)) overlay grid">
            <figcaption>Observed data are red for Naive/sensitive component and blue for cis/resistant component. The model fit is $(html_escape(string(MODEL_COLORS[spec.key]))); solid line is Naive/sensitive, dashed line is cis/resistant.</figcaption>
          </figure>
        </section>
        """)
    end
    html = """
    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <title>A2780 Delay and Compartment Model Comparison</title>
      <script>
        window.MathJax = {
          tex: { inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']] },
          svg: { fontCache: 'global' }
        };
      </script>
      <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
      <script>
        function showComparisonTab(tabName) {
          document.querySelectorAll('[data-comparison-panel]').forEach(function(panel) {
            panel.hidden = panel.dataset.comparisonPanel !== tabName;
          });
          document.querySelectorAll('[data-comparison-tab]').forEach(function(button) {
            const selected = button.dataset.comparisonTab === tabName;
            button.classList.toggle('active', selected);
            button.setAttribute('aria-selected', selected ? 'true' : 'false');
          });
        }
        window.addEventListener('DOMContentLoaded', function() { showComparisonTab('both'); });
      </script>
      <style>
        :root { --ink: #1f2933; --muted: #52606d; --line: #cbd2d9; --bg: #f7f9fb; --panel: #ffffff; }
        body { margin: 0; font-family: Arial, Helvetica, sans-serif; color: var(--ink); background: var(--bg); line-height: 1.45; }
        header, main { max-width: 1180px; margin: 0 auto; padding: 24px; }
        header { padding-top: 32px; }
        h1 { margin: 0 0 8px; font-size: 30px; }
        h2 { margin: 0 0 10px; font-size: 22px; }
        p { margin: 8px 0; }
        code { background: #eef2f6; padding: 1px 4px; border-radius: 4px; }
        .math { overflow-x: auto; padding: 8px 0; }
        .tabs { display: flex; gap: 8px; flex-wrap: wrap; margin: 14px 0 12px; }
        .tabs button { border: 1px solid var(--line); background: #fff; color: var(--ink); padding: 8px 12px; border-radius: 6px; cursor: pointer; font-weight: 600; }
        .tabs button.active { background: #1f2933; color: #fff; border-color: #1f2933; }
        .legend-strip { display: flex; justify-content: center; gap: 22px; flex-wrap: wrap; padding: 10px 12px; background: #ffffff; border: 1px solid var(--line); border-top: 0; font-size: 14px; }
        .legend-item { display: inline-flex; align-items: center; gap: 7px; white-space: nowrap; }
        .dot { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
        .line { width: 28px; height: 0; border-top: 3px solid currentColor; display: inline-block; }
        .dash { border-top-style: dashed; }
        table { width: 100%; border-collapse: collapse; background: var(--panel); margin: 18px 0 24px; }
        th, td { border: 1px solid var(--line); padding: 9px 10px; text-align: left; vertical-align: top; }
        th { background: #e8edf2; }
        .model-card { background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 18px; margin: 22px 0; }
        .stage-card { background: var(--panel); border-left: 4px solid #7b8794; padding: 16px 18px; margin: 18px 0; }
        figure { margin: 14px 0 0; }
        img { max-width: 100%; height: auto; border: 1px solid var(--line); background: white; }
        figcaption { color: var(--muted); font-size: 14px; margin-top: 6px; }
        .note { color: var(--muted); }
      </style>
    </head>
    <body>
      <header>
        <h1>A2780 Treated Coculture: Time Delay vs Transit / Compartment</h1>
        <p class="note">Generated $(html_escape(generated)). This report uses existing fitted prediction overlays; it does not refit the models.</p>
      </header>
      <main>
        <section>
          <h2>Direct Overlay Comparison</h2>
          <p>The best time-delay-only fit and the best transit/compartment fit are overlaid on the same treated-coculture data panels. Data are red/blue by lineage component; model colors identify the candidate mechanism.</p>
          <div class="tabs" role="tablist" aria-label="Model overlay view">
            <button type="button" data-comparison-tab="both" role="tab" onclick="showComparisonTab('both')">Both models</button>
            <button type="button" data-comparison-tab="time" role="tab" onclick="showComparisonTab('time')">Time-delay only</button>
            <button type="button" data-comparison-tab="compartment" role="tab" onclick="showComparisonTab('compartment')">Transit / compartment only</button>
          </div>
          <figure data-comparison-panel="both">
            <img src="$(rel_report_path(image_paths["all"]))" alt="Delay versus compartment model comparison grid">
            <figcaption>Both candidate mechanisms overlaid on the same observed treated-coculture data.</figcaption>
          </figure>
          <figure data-comparison-panel="time" hidden>
            <img src="$(rel_report_path(image_paths["time_delay"]))" alt="Time-delay-only model overlay grid">
            <figcaption>Only the best time-delay-only model is shown on the observed treated-coculture data.</figcaption>
          </figure>
          <figure data-comparison-panel="compartment" hidden>
            <img src="$(rel_report_path(image_paths["compartment"]))" alt="Transit compartment model overlay grid">
            <figcaption>Only the best transit/compartment model is shown on the observed treated-coculture data.</figcaption>
          </figure>
          <div class="legend-strip" aria-label="Overlay key">
            <span class="legend-item"><span class="dot" style="background:#dc143c"></span>Naive / sensitive observed</span>
            <span class="legend-item"><span class="dot" style="background:#4682b4"></span>cis / resistant observed</span>
            <span class="legend-item" style="color:#ff8c00"><span class="line"></span>Time-delay fit, Naive/sensitive</span>
            <span class="legend-item" style="color:#ff8c00"><span class="line dash"></span>Time-delay fit, cis/resistant</span>
            <span class="legend-item" style="color:#800080"><span class="line"></span>Transit/compartment fit, Naive/sensitive</span>
            <span class="legend-item" style="color:#800080"><span class="line dash"></span>Transit/compartment fit, cis/resistant</span>
          </div>
        </section>
        <section>
          <h2>How The Math Builds Up</h2>
          <p>These sections show the inherited equations stage by stage. The important idea is that the treated coculture candidates do not start from scratch: they inherit baseline growth, drug timing, and untreated competition structure, then test whether the treated hump is better explained by time activation alone or by a damaged-visible compartment.</p>
        </section>
        $(stage_walkthrough_html())
        <section>
          <h2>Model Summary</h2>
          <table>
            <thead>
              <tr><th>Comparison</th><th>Model</th><th>Pooling</th><th>BIC</th><th>Delta BIC</th><th>Parameters</th></tr>
            </thead>
            <tbody>
              $(model_table(MODEL_SPECS, metrics))
            </tbody>
          </table>
        </section>
        $(join(cards, "\n"))
      </main>
    </body>
    </html>
    """
    path = joinpath(OUT_REPORT, "a2780_delay_compartment_comparison.html")
    write(path, html)
    return path
end

function main()
    overlays = Dict{String,DataFrame}()
    image_paths = Dict{String,String}()
    metrics = Dict{String,Any}()
    for spec in MODEL_SPECS
        overlays[spec.key] = source_overlay(spec)
        metrics[spec.key] = metric_for(spec)
        image_paths[spec.key] = render_single_model_grid(spec, overlays[spec.key])
    end
    image_paths["all"] = render_all_model_grid(overlays)
    report = write_report(image_paths, metrics)
    println("report_path=$(report)")
end

main()
