# Small "i" info buttons that explain how to read a graph or table. Used by zoo_report.jl,
# zoo_analyses.jl and build_new_data_report.jl; docs/outputs/results/site.js has the same
# component for the Summary, Presentation and Explorer.
#   INFO_TIPS_HEAD  put once in the page <head> (styles + script)
#   tip(text)       the button; text is plain text shown in a pop-up on hover, focus or click
# The pop-up is position:fixed and kept inside the window, so it also works inside the Final
# Report's chapter frames. Clicks on a button do not toggle an enclosing <details>/<summary>.

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

# Every tip must say which direction is better; it is shown as the pop-up's last line.
#   :lower / :higher / :none (descriptive, not a score), or a phrase such as "Closer to 1 is better."
const BETTER_TEXT = Dict(:lower => "Lower is better.", :higher => "Higher is better.",
                         :none => "Neither: this describes the model or data; it is not a score.")
better_text(b::Symbol) = BETTER_TEXT[b]
better_text(b::AbstractString) = b
tip(text, better) = "<button type=\"button\" class=\"ii\" aria-label=\"How to read this\" data-tip=\"$(tip_esc(text))\" data-better=\"$(tip_esc(better_text(better)))\">i</button>"

# Table headers: a 3-tuple (label, text, better) is used as is; a 2-tuple looks its direction up by
# label, and an unknown label stops the build so no tip goes out without a direction.
function header_better(dict, h)
    length(h) == 3 && return h[3]
    haskey(dict, h[1]) || error("info tip for column \"$(h[1])\" has no better-direction; add it to the BETTER table")
    dict[h[1]]
end
