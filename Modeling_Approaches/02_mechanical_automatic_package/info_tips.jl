# Small "i" info buttons that explain how to read a graph or table. Used by zoo_report.jl,
# zoo_analyses.jl and build_new_data_report.jl; docs/outputs/results/site.js has the same
# component for the Summary, Presentation and Explorer.
#   INFO_TIPS_HEAD  put once in the page <head> (styles + script)
#   tip(text)       the button; text is plain text shown in a pop-up on hover, focus or click
#   warn(text)      a small amber ⚠ button next to the "i": a caveat shown in the same pop-up on click,
#                   so a warning is written once (or kept in WARNINGS) instead of repeated in the prose;
#                   tip(text, better; warn=...) and withtip(html, text, better; warn=...) put it beside the i
# The pop-up is position:fixed and kept inside the window, so it also works inside the Final
# Report's chapter frames. Clicks on a button do not toggle an enclosing <details>/<summary>.

const INFO_TIPS_HEAD = raw"""
<style>
.ii{display:inline-flex;align-items:center;justify-content:center;width:17px;height:17px;margin:0 4px;padding:0;border-radius:50%;border:1.5px solid #176b87;background:#fff;color:#176b87;font:italic 700 11px/1 Georgia,serif;cursor:help;vertical-align:middle;flex:none}
.ii:hover,.ii:focus-visible,.ii[aria-expanded=true]{background:#176b87;color:#fff;outline:none}
#ii-pop.ww-pop{background:#3a2606;border-left:5px solid #f59e0b}
.ww{display:inline-flex;align-items:center;justify-content:center;width:19px;height:17px;margin:0 4px 0 0;padding:0;border:1.5px solid #b45309;border-radius:4px;background:#fff6d9;color:#b45309;font:700 11px/1 system-ui,sans-serif;cursor:help;vertical-align:middle;flex:none}
.ww:hover,.ww:focus-visible,.ww[aria-expanded=true]{background:#b45309;color:#fff;outline:none}
.tipbtns{display:inline-flex;align-items:center;vertical-align:middle}.tipwrap>.tipbtns{position:absolute;top:8px;right:8px;z-index:2}
#ii-pop{position:fixed;z-index:9999;max-width:min(400px,calc(100vw - 24px));background:#17212b;color:#f4f7f8;font:13px/1.45 system-ui,sans-serif;font-style:normal;font-weight:400;text-align:left;padding:10px 12px;border-radius:8px;box-shadow:0 6px 20px rgba(0,0,0,.25);display:none;white-space:normal}
</style>
<script>
(function () {
  if (window.__iiInit) return; window.__iiInit = true;
  var pop = null, cur = null, pinned = false;
  function place(b) {
    if (!pop) { pop = document.createElement('div'); pop.id = 'ii-pop'; pop.setAttribute('role', 'tooltip'); document.body.appendChild(pop); }
    var wt = b.getAttribute('data-warn');
    pop.className = wt ? 'ww-pop' : '';
    if (wt) { pop.textContent = ''; var hd = document.createElement('div'); hd.style.cssText = 'font-weight:700;color:#ffd479;margin-bottom:4px'; hd.textContent = '⚠ Caution'; pop.appendChild(hd); pop.appendChild(document.createTextNode(wt)); }
    else pop.textContent = b.getAttribute('data-tip');
    var bt = wt ? null : b.getAttribute('data-better');
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
    var b = e.target.closest ? e.target.closest('.ii, .ww') : null;
    if (b) { e.preventDefault(); e.stopPropagation(); if (cur === b && pinned) hide(); else { place(b); pinned = true; } }
    else if (!(pop && pop.contains(e.target))) hide();
  }, true);
  document.addEventListener('mouseover', function (e) { var b = e.target.closest ? e.target.closest('.ii, .ww') : null; if (b && !pinned) place(b); });
  document.addEventListener('mouseout', function (e) { var b = e.target.closest ? e.target.closest('.ii, .ww') : null; if (b && !pinned && !b.contains(e.relatedTarget)) hide(); });
  document.addEventListener('focusin', function (e) { var b = e.target.closest ? e.target.closest('.ii, .ww') : null; if (b) place(b); });
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') hide(); });
  window.addEventListener('scroll', function () { if (!pinned) hide(); else if (cur) place(cur); }, true);
})();
</script>
"""

tip_esc(s) = replace(string(s), "&" => "&amp;", "\"" => "&quot;", "<" => "&lt;", ">" => "&gt;")

# Every tip must say which direction is better; it is shown as the pop-up's last line.
#   :lower / :higher / :none (descriptive, not a score), or a phrase such as "Closer to 1 is better."
const BETTER_TEXT = Dict(:lower => "Lower is better.", :higher => "Higher is better.",
                         :none => "Neither: this describes the model or data; it is not a score.")
better_text(b::Symbol) = BETTER_TEXT[b]
better_text(b::AbstractString) = b
tip_i(text, better) = "<button type=\"button\" class=\"ii\" aria-label=\"How to read this\" data-tip=\"$(tip_esc(text))\" data-better=\"$(tip_esc(better_text(better)))\">i</button>"
# warn(:key) uses the shared wording in WARNINGS; warn(text) any other caveat; a vector gives one ⚠ per item
warn_btn(text::AbstractString) = "<button type=\"button\" class=\"ww\" aria-label=\"Caution\" data-warn=\"$(tip_esc(text))\">&#9888;</button>"
warn_btn(key::Symbol) = warn_btn(WARNINGS[key])
warn_btn(ws::AbstractVector) = join(warn_btn.(ws))
warn_btn(::Nothing) = ""
warn(x) = warn_btn(x)
tip(text, better; warn=nothing) = warn === nothing ? tip_i(text, better) : "<span class='tipbtns'>$(tip_i(text, better))$(warn_btn(warn))</span>"
# a figure or table with its tip (and caveats) in the top-right corner
withtip(html, text, better; warn=nothing) = "<div class='tipwrap'><span class='tipbtns'>$(tip_i(text, better))$(warn_btn(warn))</span>$html</div>"

# Caveats that recur across the Final Report chapters, written once so every ⚠ says the same thing.
const WARNINGS = Dict(
    :ce1_selection => "Not a fully independent test: the Ce1 cocultures were looked at when the modelling setup was chosen (the residual weighting and the drug model were compared partly on how well they predicted Ce1), so scores on Ce1 are somewhat optimistic. The intermittent-dosing data are the untouched validation set.",
    :extrapolation => "Beyond the data: the experiments ran 14 days, so every day after 14 is a model extrapolation (no passaging, no re-plating, the fitted parameters held constant).",
    :at_bound => "A parameter at its search bound is not pinned down by the data: its interval and sensitivity are not meaningful, and the profile-likelihood test does not apply.",
    :tox_proxy => "Toxicity was never measured: T is a literature proxy (de Pillis & Radunskaya 2001) with an assumed limit, so toxicity numbers are illustrative; with the measured doses the limit never binds.",
    :soft_penalty => "B_cancer is a soft penalty, not a hard cap: the controller pays a steep cost above it but may cross it when the alternative costs more.",
    :refeeding => "Recovery after refeeding was never measured: how fast the low-resource drug interaction fades after a fast is an untested assumption.",
    :assumption => "An assumption, not fitted to data: the results change if this value changes (see the sensitivity tables).",
    :not_patient_days => "Days are culture days in a well, not patient days: the host terms (toxicity, weight, protection) are assumed from the literature, not fitted.",
    :bootstrap_narrow => "The within-trajectory bootstrap treats errors as independent, so its intervals are usually too narrow (over-confident); trust the block-bootstrap interval.",
    :hypothesis_test => "A hypothesis test, not a validation: it asks whether an extra term is supported, not whether the model predicts new data.",
)

# Table headers: a 3-tuple (label, text, better) is used as is; a 2-tuple looks its direction up by
# label, and an unknown label stops the build so no tip goes out without a direction.
function header_better(dict, h)
    length(h) >= 3 && h[3] !== nothing && return h[3]     # (label, text, better[, warning]); better = nothing looks it up
    haskey(dict, h[1]) || error("info tip for column \"$(h[1])\" has no better-direction; add it to the BETTER table")
    dict[h[1]]
end
