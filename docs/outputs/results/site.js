// Shared helpers for the Summary, Presentation and Explorer pages. Everything is drawn from
// window.RESULTS (bundle.js), which each full pipeline run rewrites.
'use strict';

const RES = window.RESULTS || null;

// ── formatting ───────────────────────────────────────────────────────────────
function fmt(v) {
  if (v === null || v === undefined || !isFinite(v)) return '—';
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1000) return v.toFixed(0);
  if (a >= 100) return v.toFixed(1);
  if (a < 1e-3) return a < 1e-8 ? '≈ 0' : v.toExponential(1);
  return String(+v.toPrecision(3));
}
const pctf = v => (v === null || v === undefined || !isFinite(v)) ? "—" : Math.round(v) + "%";
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const texBlock = t => `<div class="eq">\\[${t}\\]</div>`;
const texInline = t => `\\(${t}\\)`;

// ── "i" info buttons: infoTip(text) returns a button; hovering, focusing or clicking it shows how
// to read the graph or column next to it (same component as info_tips.jl in the Julia reports).
// Every tip states what is better (same rule as info_tips.jl): 'lower', 'higher', 'none' or a phrase.
const BETTER_TEXT = { lower: 'Lower is better.', higher: 'Higher is better.', none: 'Neither: this describes the model or data; it is not a score.' };
const infoTip = (t, b) => {
  if (b === undefined) console.warn('info tip without a better-direction:', t);
  const bt = b === undefined ? '' : (BETTER_TEXT[b] || b);
  return `<button type="button" class="ii" aria-label="How to read this" data-tip="${esc(t)}" data-better="${esc(bt)}">i</button>`;
};
// the shared tips below, with their direction
const tipK = k => infoTip(TIPS[k], TIPS_BETTER[k]);
const TIPS = {
  plot: "Points: measured mean cells (open points and dashed lines: treated Ce1 cocultures; grey × = excluded outlier). Line: the optimal model fitted to the data. Shaded band: where 95% of the curves refitted to the bootstrap data sets lie, i.e. the uncertainty of the fitted curve, not of single measurements.",
  param: "Parameters fitted at this stage. Rows under 'Carried over / fixed' come from earlier stages and are held at those values.",
  value: "Value fitted to the data.",
  boot: "95% of the values refitted to the bootstrap data sets (the data resampled in 3-day blocks) fall in this range. A wide range = the data determine this parameter poorly.",
  flag: "'at bound' = the value sits at its search limit; 'interval reaches 0' = the data are consistent with this term being absent.",
  model: "The optimal model: fewest free parameters within 3 BIC of the lowest, with its equations. N, C: sensitive (A2780) and resistant (A2780cis) counts; D: cisplatin dose (µM); u_R: low-resource level; A, S: delay compartments driven by D or u_R; K: carrying capacity; α: crowding coefficients.",
  chosen: "Share of the bootstrap data sets in which the same selection rule picks this model. Above about 70% = a robust choice; low = another model is about as good.",
  bic: "Bayesian information criterion, n·ln(SSE/n) + k·ln(n): the fit error penalised for the number of parameters. When the fits are weighted, SSE uses each residual divided by that point's replicate SD, so BICs are on that scale. Lower is better; differences under about 2 are negligible, above 10 decisive. Only comparable within one fit.",
  runner: "The next candidate and its ΔBIC relative to the optimal model. Negative = it fits slightly better but needs more parameters, so the rule keeps the smaller model.",
  k: "Number of free parameters fitted at this stage (values carried over from earlier stages are not counted).",
  family: "The delay mechanism: immediate, effect compartment, transit compartments, second order, accumulation, or damaged cells.",
  dbic: "BIC minus the optimal model's BIC. Negative = slightly better fit but more parameters; positive = worse fit.",
  fit: "The stage and data set.",
  stats: "k = free parameters; BIC = fit quality penalised for parameters (lower is better); chosen in = how often the selection rule picks this model when the data are resampled; reaches zero = whether the population dies out in a long-term simulation.",
  eq: "The model's equations. N, C: sensitive and resistant counts; V: viable (non-damaged) cells; A or S: delay compartments driven by the dose D or the low-resource level u_R. Values are in the table below.",
  note: "How the shared-mechanism rule was applied across the two cell lines.",
  why: "Why the rule kept this model rather than the other candidates within 3 BIC of the lowest: fewest free parameters first, then the lower BIC.",
};
const TIPS_BETTER = {
  plot: 'Better: the line runs through the points and the shaded band is narrow.',
  param: 'none', value: 'none', model: 'none', family: 'none', fit: 'none', eq: 'none', note: 'none', why: 'none',
  boot: 'Narrower is better: the parameter is more precisely determined.',
  flag: 'Better: no flag.',
  chosen: 'higher', bic: 'lower',
  runner: 'Neither: it shows how close the next model is. A larger positive ΔBIC means a clearer choice.',
  k: 'Lower is better (fewer parameters), as long as the fit stays within 3 BIC.',
  dbic: 'Lower means a better fit, but within 3 of 0 the carried model wins by having fewer parameters; above +3 is clearly worse.',
  stats: 'Lower BIC and fewer parameters are better; a higher "chosen in" is better; "reaches zero" describes the prediction and is neither.',
};
(function () {
  if (typeof document === 'undefined' || typeof window === 'undefined' || window.__iiInit) return; window.__iiInit = true;
  let pop = null, cur = null, pinned = false;
  function place(b) {
    if (!pop) { pop = document.createElement('div'); pop.id = 'ii-pop'; pop.setAttribute('role', 'tooltip'); document.body.appendChild(pop); }
    pop.textContent = b.getAttribute('data-tip');
    const bt = b.getAttribute('data-better');
    if (bt) { const d = document.createElement('div'); d.style.cssText = 'margin-top:6px;font-weight:700;color:#ffd479'; d.textContent = 'What is better: ' + bt; pop.appendChild(d); }
    pop.style.display = 'block';
    const r = b.getBoundingClientRect(), w = pop.offsetWidth, h = pop.offsetHeight;
    const x = Math.min(Math.max(8, r.left + r.width / 2 - w / 2), innerWidth - w - 8);
    let y = r.bottom + 8; if (y + h > innerHeight - 8) y = Math.max(8, r.top - h - 8);
    pop.style.left = x + 'px'; pop.style.top = y + 'px';
    if (cur && cur !== b) cur.setAttribute('aria-expanded', 'false');
    cur = b; b.setAttribute('aria-expanded', 'true');
  }
  function hide() { if (pop) pop.style.display = 'none'; if (cur) cur.setAttribute('aria-expanded', 'false'); cur = null; pinned = false; }
  document.addEventListener('click', e => {
    const b = e.target.closest && e.target.closest('.ii');
    if (b) { e.preventDefault(); e.stopPropagation(); if (cur === b && pinned) hide(); else { place(b); pinned = true; } }
    else if (!(pop && pop.contains(e.target))) hide();
  }, true);
  document.addEventListener('mouseover', e => { const b = e.target.closest && e.target.closest('.ii'); if (b && !pinned) place(b); });
  document.addEventListener('mouseout', e => { const b = e.target.closest && e.target.closest('.ii'); if (b && !pinned && !b.contains(e.relatedTarget)) hide(); });
  document.addEventListener('focusin', e => { const b = e.target.closest && e.target.closest('.ii'); if (b) place(b); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') hide(); });
  addEventListener('scroll', () => { if (!pinned) hide(); else if (cur) place(cur); }, true);
})();
function typeset(el) { if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise(el ? [el] : undefined); }

// ── navigation shared by the three pages ─────────────────────────────────────
// the site-wide bar now comes from docs/assets/site-nav.js (same tabs on every main page); this returns
// nothing unless that script is missing, so the page still has a way home
function navHTML(active) {
  if (document.querySelector('script[src*="site-nav.js"]')) return '';
  const l = (id, href, label) => `<a href="${href}" class="${id === active ? 'on' : ''}">${label}</a>`;
  return `<nav class="topnav" aria-label="Results navigation">${l('home', '../../index.html', '⌂ Home')}` +
    l('summary', 'summary.html', 'Summary') + l('presentation', 'presentation.html', 'Presentation') +
    l('explorer', 'explorer.html', 'Explorer') + '<span class="sep"></span>' +
    l('final', '../final_report/report.html', 'Final report (all details)') + '</nav>';
}
function runInfo() {
  return RES ? `<span class="runinfo">Latest full run: ${esc(RES.generated)} · ${RES.replicates} bootstrap replicates</span>` : '';
}
function missingBundle(el) {
  el.innerHTML = '<div class="card note">No results yet: <code>bundle.js</code> is written at the end of <code>run_zoo_fits.jl</code>. Run it once to fill these pages.</div>';
}

// ── SVG plot ─────────────────────────────────────────────────────────────────
// o: {w,h,xmax,ymax,title,points:[{x,y,color,filled}],excluded:[{x,y}],lines:[{t,y,color,dash,width,opacity}],
//     bands:[{t,lo,hi,color}], legend:[{label,color,dash,kind:'line'|'dot'|'open'}]}
function plotSVG(o) {
  const W = o.w || 760, H = o.h || 330, L = 58, R = 16, T = o.title ? 30 : 14;
  const legendRows = o.legend ? Math.ceil(o.legend.length / 4) : 0;
  const B = 38 + legendRows * 18;
  const sx = v => L + (v / o.xmax) * (W - L - R);
  const sy = v => H - B - (Math.max(0, Math.min(v, o.ymax)) / o.ymax) * (H - T - B);
  const p = [`<svg class="plot" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.title || 'plot')}"><rect width="100%" height="100%" fill="#fff"/>`];
  if (o.title) p.push(`<text x="${L}" y="19" font-size="14" font-weight="700" fill="#17212b">${esc(o.title)}</text>`);
  for (let i = 0; i <= 4; i++) {
    const v = o.ymax * i / 4;
    p.push(`<line x1="${L}" x2="${W - R}" y1="${sy(v)}" y2="${sy(v)}" stroke="#e8edf0"/><text x="${L - 6}" y="${sy(v) + 3}" font-size="10" text-anchor="end" fill="#52636b">${Math.round(v)}</text>`);
  }
  const step = o.xmax > 16 ? 5 : 2;
  for (let v = 0; v <= o.xmax + 1e-9; v += step) p.push(`<text x="${sx(v)}" y="${H - B + 14}" font-size="10" text-anchor="middle" fill="#52636b">${v}</text>`);
  p.push(`<line x1="${L}" y1="${H - B}" x2="${W - R}" y2="${H - B}" stroke="#667"/><line x1="${L}" y1="${T}" x2="${L}" y2="${H - B}" stroke="#667"/>`);
  p.push(`<text x="${(L + W - R) / 2}" y="${H - B + 28}" font-size="11" text-anchor="middle" fill="#52636b">Day</text>`);
  p.push(`<text transform="translate(13 ${(T + H - B) / 2}) rotate(-90)" font-size="11" text-anchor="middle" fill="#52636b">Mean cells</text>`);
  for (const b of o.bands || []) {
    const up = [], dn = [];
    b.t.forEach((t, i) => { if (b.hi[i] != null && isFinite(b.hi[i])) up.push(`${sx(t).toFixed(1)},${sy(b.hi[i]).toFixed(1)}`); });
    for (let i = b.t.length - 1; i >= 0; i--) if (b.lo[i] != null && isFinite(b.lo[i])) dn.push(`${sx(b.t[i]).toFixed(1)},${sy(b.lo[i]).toFixed(1)}`);
    if (up.length) p.push(`<polygon points="${up.concat(dn).join(' ')}" fill="${b.color}" opacity="0.16"/>`);
  }
  for (const l of o.lines || []) {
    const pts = l.t.map((t, i) => `${sx(t).toFixed(1)},${sy(l.y[i]).toFixed(1)}`).join(' ');
    p.push(`<polyline points="${pts}" fill="none" stroke="${l.color}" stroke-width="${l.width || 2.2}"${l.dash ? ` stroke-dasharray="${l.dash}"` : ''}${l.opacity ? ` opacity="${l.opacity}"` : ''}/>`);
  }
  for (const q of o.points || []) p.push(`<circle cx="${sx(q.x).toFixed(1)}" cy="${sy(q.y).toFixed(1)}" r="3.4" fill="${q.filled ? q.color : '#fff'}" stroke="${q.color}" stroke-width="1.5"/>`);
  for (const q of o.excluded || []) {
    const x = sx(q.x), y = sy(q.y);
    p.push(`<path d="M${x - 4} ${y - 4}L${x + 4} ${y + 4}M${x - 4} ${y + 4}L${x + 4} ${y - 4}" stroke="#9ca3af" stroke-width="1.8"/>`);
  }
  (o.legend || []).forEach((g, i) => {
    const lx = L + (i % 4) * 176, ly = H - B + 46 + Math.floor(i / 4) * 18;
    if (g.kind === 'line') p.push(`<line x1="${lx}" x2="${lx + 24}" y1="${ly - 4}" y2="${ly - 4}" stroke="${g.color}" stroke-width="2.4"${g.dash ? ` stroke-dasharray="${g.dash}"` : ''}/>`);
    else p.push(`<circle cx="${lx + 12}" cy="${ly - 4}" r="4" fill="${g.kind === 'open' ? '#fff' : g.color}" stroke="${g.color}" stroke-width="1.5"/>`);
    p.push(`<text x="${lx + 30}" y="${ly}" font-size="11" fill="#17212b">${esc(g.label)}</text>`);
  });
  p.push('</svg>');
  return p.join('');
}

// data points and the optimal fit (solid, with 95% band) for one panel of a group
function groupPlot(g, panel, opts = {}) {
  const ser = g.series.filter(s => s.panel === panel);
  const points = [], excluded = [], lines = [], bands = [];
  for (const s of ser) {
    s.x.forEach((x, i) => points.push({ x, y: s.y[i], color: s.color, filled: s.filled }));
    (s.xbad || []).forEach((x, i) => excluded.push({ x, y: s.ybad[i] }));
    const key = String(s.state);
    const opt = g.optimal;
    if (opt.band && opt.band.lo[key]) bands.push({ t: opt.band.t, lo: opt.band.lo[key], hi: opt.band.hi[key], color: s.color });
    lines.push({ t: opt.curve.t, y: opt.curve.y[key], color: s.color, dash: s.dashed ? '8,4' : '', width: 2.4 });
  }
  const legend = [];
  const seen = new Set();
  for (const s of ser) if (!seen.has(s.name)) { seen.add(s.name); legend.push({ label: s.name, color: s.color, kind: s.filled ? 'dot' : 'open' }); }
  legend.push({ label: 'optimal (95% band)', color: '#374151', kind: 'line' });
  if (excluded.length) legend.push({ label: 'excluded point', color: '#9ca3af', kind: 'open' });
  const title = opts.title !== undefined ? opts.title : (g.panels[panel - 1] || '');
  return `<div class="tipwrap">${tipK('plot')}` + plotSVG({ w: opts.w, h: opts.h, xmax: g.xmax, ymax: g.ymax, title, points, excluded, lines, bands, legend }) + '</div>';
}

// ── parameter tables ─────────────────────────────────────────────────────────
function paramTable(m, compact = false) {
  const rows = m.params.map(p => {
    const flag = p.atbound ? '<span class="flag">at bound</span>' : p.removable ? '<span class="flag">interval reaches 0</span>' : '';
    return `<tr><td>${texInline(p.tex)}</td><td class="num">${fmt(p.value)}</td><td class="num">${fmt(p.lo)} – ${fmt(p.hi)}</td>${compact ? '' : `<td>${flag}</td>`}</tr>`;
  }).join('');
  const inh = m.inherited.map(p => `<tr class="inh"><td>${texInline(p.sym)}</td><td class="num">${fmt(p.value)}</td><td class="muted">${esc(p.note)}</td>${compact ? '' : '<td></td>'}</tr>`).join('');
  const head = `<tr><th>Parameter${tipK('param')}</th><th>Value${tipK('value')}</th><th>95% bootstrap${tipK('boot')}</th>${compact ? '' : `<th>${tipK('flag')}</th>`}</tr>`;
  const fitted = m.params.length ? rows : `<tr><td colspan="${compact ? 3 : 4}" class="muted">No free parameters: every value is carried over.</td></tr>`;
  return `<div class="tscroll"><table>${head}${fitted}${inh ? `<tr><th colspan="${compact ? 3 : 4}">Carried over / fixed</th></tr>${inh}` : ''}</table></div>`;
}

// ── in-browser ODE simulation (fixed-step RK4, states clamped at 0) ─────────
const sig = z => 1 / (1 + Math.exp(-z));
function erlang(x, n) { let s = 0, term = 1; for (let k = 0; k < n; k++) { if (k > 0) term *= x / k; s += term; } return 1 - Math.exp(-x) * s; }

const MODELS = {
  // Stage 1 — p: fitted params in Julia order
  logistic: (c, p) => (t, u) => [p[0] * u[0] * (1 - u[0] / p[1])],
  theta_logistic: (c, p) => (t, u) => [p[0] * u[0] * (1 - Math.pow(Math.max(u[0], 0) / p[1], p[2]))],
  gompertz: (c, p) => (t, u) => [p[0] * u[0] * Math.log(p[1] / Math.max(u[0], 1e-8))],
  // Stage 2a — c: {r, K, uR}
  lr_subtractive_loss: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * c.uR * u[0]],
  lr_rho_modifiers: (c, p) => (t, u) => [p[0] * c.r * u[0] * (1 - u[0] / (p[1] * c.K))],
  lr_separate_r: (c, p) => (t, u) => [p[0] * u[0] * (1 - u[0] / c.K)],
  delayed_loss_exp: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * (1 - Math.exp(-p[1] * t)) * u[0]],
  delayed_loss_sigmoid: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * sig(t - p[1]) * u[0]],
  delayed_loss_sigmoid_w: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * sig((t - p[1]) / p[2]) * u[0]],
  decaying_K: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / Math.max(c.K * Math.exp(-p[0] * t), 1e-6))],
  decaying_r_loss: (c, p) => (t, u) => [Math.exp(-p[0] * t) * c.r * u[0] * (1 - u[0] / c.K) - p[1] * u[0]],
  resource_depletion: (c, p) => (t, u) => { const R = Math.max(u[1], 0); return [R * c.r * u[0] * (1 - u[0] / c.K) - p[1] * (1 - R) * u[0], -p[0] * R * u[0] / 1000]; },
  toxin: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * u[1] * u[0], u[0] / 1000],
  toxin_clearance: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * u[1] * u[0], u[0] / 1000 - p[1] * u[1]],
  quiescence_death: (c, p) => (t, u) => { const P = Math.max(u[0] - u[1], 0); return [c.r * P * (1 - u[0] / c.K) - p[1] * u[1], p[0] * P - p[1] * u[1]]; },
  transit_damage: (c, p) => (t, u) => { const P = Math.max(u[0] - u[1] - u[2], 0); return [c.r * P * (1 - u[0] / c.K) - p[1] * u[2], p[0] * P - p[1] * u[1], p[1] * (u[1] - u[2])]; },
  lr_accumulated_loss4: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - p[0] * c.uR * erlang(p[1] * t, 4) * u[0]],
  lr_loss_plus_accum4: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - (p[0] + p[1] * erlang(p[2] * t, 4)) * c.uR * u[0]],
  lr_loss_plus_accum8: (c, p) => (t, u) => [c.r * u[0] * (1 - u[0] / c.K) - (p[0] + p[1] * erlang(p[2] * t, 8)) * c.uR * u[0]],
  // Stage 2b — c: {r, K, EC50, doses}; one state per dose
  delayed_hill_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * (1 - Math.exp(-p[1] * t)) * hillf(c.doses[i], p[2], c.EC50) * x),
  delayed_linear_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * (1 - Math.exp(-p[1] * t)) * c.doses[i] * x),
  immediate_linear_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * c.doses[i] * x),
  immediate_hill_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * hillf(c.doses[i], p[1], c.EC50) * x),
  // Stage 3a — six environments [N, C]; unfitted terms: α = 1, d = 0, β = inherited value (0 if none)
  coculture: (c, p) => {
    const v = n => { const i = c.pnames.indexOf(n); return i >= 0 ? p[i] : null; };
    const K = v('K'), aNC = v('α_NC') ?? 1, aCN = v('α_CN') ?? 1, dN = v('d_N') ?? 0, dC = v('d_C') ?? 0;
    const bN = v('β_N') ?? c.betaN, bC = v('β_C') ?? c.betaC;
    return (t, u) => {
      const du = new Array(12), AN = 1 - Math.exp(-c.lN * t), AC = 1 - Math.exp(-c.lC * t);
      for (let i = 0; i < 6; i++) {
        const N = u[2 * i], C = u[2 * i + 1], on = c.envDose[i] > 0 ? 1 : 0;
        du[2 * i] = c.rN * N * (1 - (N + aNC * C) / K) - dN * N - on * bN * AN * c.HN * N;
        du[2 * i + 1] = c.rC * C * (1 - (C + aCN * N) / K) - dC * C - on * bC * AC * c.HC * C;
      }
      return du;
    };
  },
  // Stage 3b — six low-resource environments [N, C]
  lr_coculture: (c, p) => {
    const map = { 'K': 'K', 'α_NC': 'aNC', 'α_CN': 'aCN', 'd_N^LR': 'dN', 'd_C^LR': 'dC' };
    const v = Object.assign({}, c.defaults);
    c.pnames.forEach((n, i) => { if (map[n]) v[map[n]] = p[i]; });
    return (t, u) => {
      const du = new Array(12);
      for (let i = 0; i < 6; i++) {
        const N = u[2 * i], C = u[2 * i + 1];
        du[2 * i] = c.rN * N * (1 - (N + v.aNC * C) / v.K) - v.dN * c.uR * N;
        du[2 * i + 1] = c.rC * C * (1 - (C + v.aCN * N) / v.K) - v.dC * c.uR * C;
      }
      return du;
    };
  },
};
function hillf(D, h, EC50) { const dh = Math.pow(D, h); return dh / (Math.pow(EC50, h) + dh); }

// ── autonomous model zoo (mirrors zoo_models.jl: no term uses t) ─────────────
const ZOO_NCELL = { damage: 1, transit: 2 }, ZOO_NKER = { comp1: 1, auc: 1, chain2: 2, chain3: 3, chain8: 8, second: 2 }, ZOO_NAUX = { resource: 1, toxin: 1 };
function zooSystem(sim, p) {
  const S = sim.sys, uR = sim.uR, EC50 = sim.EC50;
  const sv = s => (s[0] > 0 ? p[s[0] - 1] : s[1]);
  const growth = (g, reff, V, x, th) => {
    x = Math.max(x, 0);
    if (g === 'logistic') return reff * V * (1 - x);
    if (g === 'theta') return reff * V * (1 - Math.pow(x, th));
    return reff * V * Math.log(1 / Math.max(x, 1e-9));
  };
  // One mechanism module (mirrors zoo_models.jl: module_effect!): delay kernel + action.
  // o: one before the module's first hidden state (hidden state j is u[o + j]).
  const width = M => (ZOO_NCELL[M.action] || 0) + (ZOO_NKER[M.kernel] || 0) + (ZOO_NAUX[M.action] || 0);
  // ks: Norton–Simon factor on the damaged-cell transfer β·E·V (1 = no coupling), as in zoo_models.jl
  function mod(du, u, o, M, w, X, V, r, K, ks = 1) {
    const s = M.s, nc = ZOO_NCELL[M.action] || 0, nk = ZOO_NKER[M.kernel] || 0;
    w *= sv(s['ε']);
    const ko = o + nc;
    let E = w;
    if (M.kernel === 'comp1') { du[ko + 1] = sv(s.k) * (w - u[ko + 1]); E = Math.max(u[ko + 1], 0); }
    else if (M.kernel.startsWith('chain')) {
      const k = sv(s.k);
      du[ko + 1] = k * (w - u[ko + 1]);
      for (let j = 2; j <= nk; j++) du[ko + j] = k * (u[ko + j - 1] - u[ko + j]);
      E = Math.max(u[ko + nk], 0);
    } else if (M.kernel === 'second') {
      const om = sv(s['ω']), ze = sv(s['ζ']);
      du[ko + 1] = u[ko + 2]; du[ko + 2] = om * om * (w - u[ko + 1]) - 2 * ze * om * u[ko + 2];
      E = Math.max(u[ko + 1], 0);
    } else if (M.kernel === 'auc') { du[ko + 1] = w; E = Math.max(u[ko + 1], 0); }
    const ao = ko + nk;
    let reff = r, Keff = K, loss = 0, out = 0;
    const a = M.action;
    if (a === 'loss') loss += sv(s.d0) * w + sv(s['β']) * E;
    else if (a === 'hill') { const h = sv(s.h), ec = s.EC50 ? sv(s.EC50) : EC50; loss += E > 0 ? sv(s['β']) * Math.pow(E, h) / (Math.pow(ec, h) + Math.pow(E, h)) : 0; }
    else if (a === 'damage') { const d = sv(s['δ']), Dm = Math.max(u[o + 1], 0); du[o + 1] = ks * sv(s['β']) * E * V - d * Dm; out = d * Dm; }
    else if (a === 'transit') { const tau = sv(s['τ']), D1 = Math.max(u[o + 1], 0), D2 = Math.max(u[o + 2], 0); du[o + 1] = ks * sv(s['β']) * E * V - tau * D1; du[o + 2] = tau * (D1 - D2); out = tau * D2; }
    else if (a === 'capacity') Keff = K / (1 + sv(s.c) * E);
    else if (a === 'fade') { reff = r / (1 + sv(s.c) * E); loss += sv(s.d0) * w; }
    else if (a === 'rho') { const f = w / uR; reff = r * (1 + (sv(s['ρr']) - 1) * f); Keff = K * (1 + (sv(s['ρK']) - 1) * f); }
    else if (a === 'sep_r') reff = r + (sv(s.rLR) - r) * w / uR;
    else if (a === 'resource') { const R = Math.max(u[ao + 1], 0); du[ao + 1] = -sv(s.c) * E * R * X / 1000; reff = r * R; loss += sv(s['β']) * (1 - R); }
    else if (a === 'toxin') { const W = Math.max(u[ao + 1], 0); du[ao + 1] = E * X / 1000 - sv(s['κ']) * W; loss += sv(s['γ']) * W; }
    return [reff, Keff, loss, out];
  }
  // A line with one or more modules (mirrors line_rhs! / line_rhs_multi!): growth-rate and
  // capacity effects multiply, losses add; drug losses scale by (1 + χ·u_R/0.955) if L.chi is set.
  // xi: observed-count index; o: one before the line's first hidden state.
  function line(du, u, xi, o, L, env, crowd) {
    const mods = [L].concat(L.extra || []);
    const X = Math.max(u[xi], 0);
    let cells = 0, off = o;
    for (const M of mods) { const nc = ZOO_NCELL[M.action] || 0; for (let j = 1; j <= nc; j++) cells += Math.max(u[off + j], 0); off += width(M); }
    const V = Math.max(X - cells, 0);
    const r = sv(L.s.r), K = sv(L.s.K);
    let reff = r, Keff = K, loss = sv(L.s.d), dloss = 0, out = 0;
    const chi = L.chi ? sv(L.chi) : 0, th = sv(L.s['θ']);
    off = o;
    for (const M of mods) {
      const w = M.input === 'D' ? env.D : env.R;
      const gx = Math.max(crowd / K, 0), ks = L.ns && M.input === 'D' && (M.action === 'damage' || M.action === 'transit')
        ? (L.growth === 'theta' ? Math.max(1 - Math.pow(gx, th), 0) : Math.max(1 - gx, 0)) : 1;
      const [re, Ke, lo, ou] = mod(du, u, off, M, w, X, V, r, K, ks);
      reff *= re / r; Keff *= Ke / K;
      if (M.input === 'D') dloss += lo * (1 + chi * env.R / uR); else loss += lo;
      out += ou; off += width(M);
    }
    // Norton–Simon (L.ns, mirrors growth_factor in zoo_models.jl): drug losses × the line's growth factor
    if (L.ns) { const x = Math.max(crowd / Keff, 0); dloss *= L.growth === 'theta' ? Math.max(1 - Math.pow(x, th), 0) : Math.max(1 - x, 0); }
    loss += dloss;
    du[xi] = growth(L.growth, reff, V, crowd / Keff, th) - loss * V - out;
  }
  return (t, u) => {
    const du = new Array(u.length).fill(0);
    S.envs.forEach((env, e) => env.lines.forEach((l, j) => {
      const xi = S.xi[e][j], L = S.lines[l];
      let crowd = Math.max(u[xi], 0);
      env.lines.forEach((lm, m) => { if (m !== j) crowd += sv(S.alpha[l][lm]) * Math.max(u[S.xi[e][m]], 0); });
      line(du, u, xi, S.hb[e][j] - 1, L, env, crowd);
    }));
    return du;
  };
}
MODELS.zoo = (sim, p) => zooSystem(sim, p);

// ── model builder: any kernel × effect on a monoculture stage (mirrors line_spec in zoo_models.jl) ──
const ZOO_KEYS = ['r', 'K', 'θ', 'k', 'ω', 'ζ', 'β', 'h', 'δ', 'd0', 'c', 'ρr', 'ρK', 'rLR', 'γ', 'κ', 'τ', 'ε', 'd'];
const ZOO_DEFAULT = { θ: 1, k: 1, ω: 1, ζ: 1, β: 0, h: 1, δ: 1, d0: 0, c: 0, ρr: 1, ρK: 1, rLR: NaN, γ: 0, κ: 0, τ: 1, ε: 1, d: 0 };
function builderSlots(b, kernel, action, immediate) {
  const ks = b.kernels.find(k => k.id === kernel).slots, as = b.actions.find(a => a.id === action).slots;
  const extra = immediate && action === 'loss' && kernel !== 'none' ? ['d0'] : [];
  return [...new Set([...extra, ...as, ...ks])];
}
function zooBuild(b, kernel, action, free) {
  const s = {};
  for (const k of ZOO_KEYS) s[k] = [0, k === 'r' ? b.r : k === 'K' ? b.K : k === 'θ' ? b['θ'] : ZOO_DEFAULT[k]];
  free.forEach((k, i) => { s[k] = [i + 1, 0]; });
  const nc = ZOO_NCELL[action] || 0, nk = ZOO_NKER[kernel] || 0, na = ZOO_NAUX[action] || 0, w = nc + nk + na;
  const nE = b.envs.length, xi = [], hb = [], u0 = new Array(nE + nE * w).fill(0), clampMask = new Array(nE + nE * w).fill(true);
  b.envs.forEach((e, i) => {
    xi.push([i]); hb.push([nE + i * w]); u0[i] = b.X0[i];
    if (action === 'resource') u0[nE + i * w + nc + nk] = 1;
    if (kernel === 'second') clampMask[nE + i * w + nc + 1] = false;
  });
  return { id: 'zoo', uR: b.uR, EC50: b.EC50, u0, clamp: clampMask,
    sys: { lines: [{ growth: b.growth, kernel, action, input: b.input, s, ns: !!b.ns && b.input === 'D' }], envs: b.envs.map(e => ({ D: e.D, R: e.R, lines: [0] })),
           alpha: [[[0, 1]]], xi, hb, n: u0.length } };
}
// symbolic equations for one monoculture line (a port of line_tex in zoo_models.jl)
function zooTex(b, kernel, action, free) {
  const X = b.X || 'X', nm = k => (b.slots[k] ? b.slots[k].tex : k);
  const w = b.input === 'D' ? 'D' : 'u_R', A = b.input === 'D' ? `A_{${X}}` : `S_{${X}}`, lines = [];
  let E = w;
  const n = { chain2: 2, chain3: 3, chain8: 8 }[kernel];
  if (kernel === 'comp1') { lines.push(`\\frac{d${A}}{dt}=${nm('k')}\\left(${w}-${A}\\right)`); E = A; }
  else if (n) { const a = i => (b.input === 'D' ? `A_{${X},${i}}` : `S_{${X},${i}}`); lines.push(`\\frac{d${a(1)}}{dt}=${nm('k')}\\left(${w}-${a(1)}\\right),\\ \\frac{d${a('j')}}{dt}=${nm('k')}\\left(${a('j-1')}-${a('j')}\\right),\\ j=2,\\dots,${n}`); E = a(n); }
  else if (kernel === 'second') { lines.push(`\\frac{d^2${A}}{dt^2}={${nm('ω')}}^2\\left(${w}-${A}\\right)-2${nm('ζ')}${nm('ω')}\\frac{d${A}}{dt}`); E = A; }
  else if (kernel === 'auc') { lines.push(`\\frac{d${A}}{dt}=${w}`); E = A; }
  const V = action === 'damage' || action === 'transit' ? `V_{${X}}` : X;
  let r = `r_{${X}}`, K = `K_{${X}}`;
  if (action === 'capacity') K = `\\frac{K_{${X}}}{1+${nm('c')}${E}}`;
  if (action === 'fade') r = `\\frac{r_{${X}}}{1+${nm('c')}${E}}`;
  if (action === 'rho') { r = `\\rho_r r_{${X}}`; K = `\\rho_K K_{${X}}`; }
  if (action === 'sep_r') r = nm('rLR');
  if (action === 'resource') r = `r_{${X}} R_{${X}}`;
  const fr = `\\frac{${X}}{${K}}`;
  let rhs = b.growth === 'logistic' ? `${r}\\,${V}\\left(1-${fr}\\right)` : b.growth === 'theta' ? `${r}\\,${V}\\left[1-\\left(${fr}\\right)^{\\theta_{${X}}}\\right]` : `${r}\\,${V}\\ln\\left(\\frac{${K}}{${X}}\\right)`;
  if (action === 'loss') rhs += free.includes('d0') ? `-\\left(${nm('d0')}${w}+${nm('β')}${E}\\right)${V}` : `-${nm('β')}\\,${E}\\,${V}`;
  else if (action === 'hill') rhs += `-${nm('β')}\\frac{${E}^{${nm('h')}}}{EC_{50}^{${nm('h')}}+${E}^{${nm('h')}}}\\,${V}`;
  else if (action === 'damage') { rhs += `-${nm('δ')}\\,${X}_D`; lines.push(`\\frac{d${X}_D}{dt}=${nm('β')}\\,${E}\\,${V}-${nm('δ')}\\,${X}_D,\\quad ${V}=${X}-${X}_D`); }
  else if (action === 'transit') { rhs += `-${nm('τ')}\\,${X}_{D2}`; lines.push(`\\frac{d${X}_{D1}}{dt}=${nm('β')}${E}${V}-${nm('τ')}${X}_{D1},\\ \\frac{d${X}_{D2}}{dt}=${nm('τ')}\\left(${X}_{D1}-${X}_{D2}\\right),\\ ${V}=${X}-${X}_{D1}-${X}_{D2}`); }
  else if (action === 'fade') rhs += `-${nm('d0')}\\,${w}\\,${X}`;
  else if (action === 'resource') { rhs += `-${nm('β')}\\left(1-R_{${X}}\\right)${X}`; lines.push(`\\frac{dR_{${X}}}{dt}=-${nm('c')}\\,${E}\\,R_{${X}}\\frac{${X}}{1000},\\ R_{${X}}(0)=1`); }
  else if (action === 'toxin') { rhs += `-${nm('γ')}\\,W_{${X}}\\,${X}`; lines.push(`\\frac{dW_{${X}}}{dt}=${E}\\frac{${X}}{1000}-${nm('κ')}W_{${X}}`); }
  return [`\\frac{d${X}}{dt}=${rhs}`, ...lines];
}

// ── in-browser least-squares fit (bounded Nelder–Mead on the fraction of each bound range) ──
function sseOf(g, sim, p, fast = false) {
  const res = fast ? simulate(sim, p, g.xmax, 0.25, 4) : simulate(sim, p, g.xmax);
  if (!res) return Infinity;
  let s = 0;
  // residuals divided by each point's scale when the fits are weighted (ser.w, from the bundle), as in the fits
  for (const ser of g.series) ser.x.forEach((x, i) => { const d = (at(res, ser.state, x) - ser.y[i]) / (ser.w ? ser.w[i] : 1); s += d * d; });
  return isFinite(s) ? s : Infinity;
}
const dsign = d => (d >= 0 ? '+' : '') + fmt(d);
function nPoints(g) { return g.series.reduce((a, s) => a + s.x.length, 0); }
function bicOf(g, sse, k) { const n = nPoints(g); return n * Math.log(sse / n) + k * Math.log(n); }
function fitNM(g, sim, p0, bounds, iters = 300) {
  const k = p0.length;
  if (!k) return p0;
  const toP = z => z.map((v, i) => bounds[i][0] + (bounds[i][1] - bounds[i][0]) * Math.min(1, Math.max(0, v)));
  const toZ = p => p.map((v, i) => (v - bounds[i][0]) / (bounds[i][1] - bounds[i][0]));
  const f = z => sseOf(g, sim, toP(z), true);
  let simplex = [toZ(p0)];
  for (let i = 0; i < k; i++) { const z = simplex[0].slice(); z[i] = z[i] + (z[i] < 0.5 ? 0.05 : -0.05); simplex.push(z); }
  let vals = simplex.map(f);
  for (let it = 0; it < iters; it++) {
    const ord = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
    simplex = ord.map(i => simplex[i]); vals = ord.map(i => vals[i]);
    const c = simplex[0].map((_, j) => simplex.slice(0, k).reduce((a, z) => a + z[j], 0) / k);
    const pt = t => c.map((v, j) => v + t * (simplex[k][j] - v));
    const xr = pt(-1), fr = f(xr);
    if (fr < vals[0]) { const xe = pt(-2), fe = f(xe); if (fe < fr) { simplex[k] = xe; vals[k] = fe; } else { simplex[k] = xr; vals[k] = fr; } }
    else if (fr < vals[k - 1]) { simplex[k] = xr; vals[k] = fr; }
    else {
      const xc = pt(0.5), fc = f(xc);
      if (fc < vals[k]) { simplex[k] = xc; vals[k] = fc; }
      else { for (let i = 1; i <= k; i++) { simplex[i] = simplex[i].map((v, j) => simplex[0][j] + 0.5 * (v - simplex[0][j])); vals[i] = f(simplex[i]); } }
    }
  }
  const best = vals.indexOf(Math.min(...vals));
  return toP(simplex[best]);
}

function simulate(sim, params, tmax, dtOut = 0.1, sub = 20) {
  const make = MODELS[sim.id];
  if (!make) return null;
  const f = make(sim.id === 'zoo' ? sim : sim.consts, params);
  const h = dtOut / sub;
  let u = sim.u0.slice(), t = 0;
  const ts = [0], us = [u.slice()];
  const nOut = Math.round(tmax / dtOut);
  const mask = sim.clamp;
  const clamp = a => a.map((x, i) => (!isFinite(x) ? 0 : (mask && !mask[i]) ? x : Math.max(x, 0)));
  const add = (a, b, s) => a.map((x, i) => x + s * b[i]);
  for (let k = 0; k < nOut; k++) {
    for (let j = 0; j < sub; j++) {
      const k1 = f(t, u), k2 = f(t + h / 2, add(u, k1, h / 2)), k3 = f(t + h / 2, add(u, k2, h / 2)), k4 = f(t + h, add(u, k3, h));
      u = clamp(u.map((x, i) => x + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])));
      t += h;
    }
    ts.push(+(t.toFixed(6))); us.push(u.slice());
  }
  return { t: ts, u: us };
}
// value of state s (1-based) at time x, by linear interpolation of a simulation
function at(simres, s, x) {
  const dt = simres.t[1] - simres.t[0], i = Math.min(Math.floor(x / dt), simres.t.length - 2), f = (x - simres.t[i]) / dt;
  return simres.u[i][s - 1] * (1 - f) + simres.u[i + 1][s - 1] * f;
}
