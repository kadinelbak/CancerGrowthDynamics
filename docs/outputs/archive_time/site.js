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
function typeset(el) { if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise(el ? [el] : undefined); }

// ── navigation shared by the three pages ─────────────────────────────────────
function navHTML(active) {
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
  el.innerHTML = '<div class="card note">No results yet: <code>bundle.js</code> is written at the end of <code>run_bootstrap_fits.jl</code>. Run it once to fill these pages.</div>';
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
  return plotSVG({ w: opts.w, h: opts.h, xmax: g.xmax, ymax: g.ymax, title, points, excluded, lines, bands, legend });
}

// ── parameter tables ─────────────────────────────────────────────────────────
function paramTable(m, compact = false) {
  const rows = m.params.map(p => {
    const flag = p.atbound ? '<span class="flag">at bound</span>' : p.removable ? '<span class="flag">interval reaches 0</span>' : '';
    return `<tr><td>${texInline(p.tex)}</td><td class="num">${fmt(p.value)}</td><td class="num">${fmt(p.lo)} – ${fmt(p.hi)}</td>${compact ? '' : `<td>${flag}</td>`}</tr>`;
  }).join('');
  const inh = m.inherited.map(p => `<tr class="inh"><td>${texInline(p.sym)}</td><td class="num">${fmt(p.value)}</td><td class="muted">${esc(p.note)}</td>${compact ? '' : '<td></td>'}</tr>`).join('');
  const head = `<tr><th>Parameter</th><th>Value</th><th>95% bootstrap</th>${compact ? '' : '<th></th>'}</tr>`;
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
  // Stage 1b — c: {r, K, uR}
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
  // Stage 2 — c: {r, K, EC50, doses}; one state per dose
  delayed_hill_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * (1 - Math.exp(-p[1] * t)) * hillf(c.doses[i], p[2], c.EC50) * x),
  delayed_linear_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * (1 - Math.exp(-p[1] * t)) * c.doses[i] * x),
  immediate_linear_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * c.doses[i] * x),
  immediate_hill_kill: (c, p) => (t, u) => u.map((x, i) => c.r * x * (1 - x / c.K) - p[0] * hillf(c.doses[i], p[1], c.EC50) * x),
  // Stage 3 — six environments [N, C]; unfitted terms: α = 1, d = 0, β = inherited value (0 if none)
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

function simulate(sim, params, tmax, dtOut = 0.1) {
  const make = MODELS[sim.id];
  if (!make) return null;
  const f = make(sim.consts, params);
  const sub = 20, h = dtOut / sub;
  let u = sim.u0.slice(), t = 0;
  const ts = [0], us = [u.slice()];
  const nOut = Math.round(tmax / dtOut);
  const clamp = a => a.map(x => (isFinite(x) && x > 0 ? x : 0));
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
