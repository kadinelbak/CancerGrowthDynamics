// Small interactive plots used inside equation explanations. Each widget gets an
// empty container and the widget options from equation-notes.js. Values are
// either quoted from the draft or marked as illustrative in the widget itself.
(() => {
  const COLORS = { sens: '#2563eb', res: '#dc2626', total: '#1f2933', dose: '#9aa5b1', a: '#2563eb', b: '#d97706', c: '#059669', d: '#7c3aed' };

  function niceStep(span, target = 5) {
    const raw = span / target;
    const mag = 10 ** Math.floor(Math.log10(raw));
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
  }
  const fmt = (v) => (Math.abs(v) >= 1000 ? Math.round(v).toLocaleString() : +v.toFixed(Math.abs(v) < 1 ? 2 : 1));

  // plot(canvas, {x:[min,max], y:[min,max], xLabel, yLabel, series:[{pts, color, dash, width}], hlines:[{y,color,label}], vbands:[{from,to,color}]})
  function plot(canvas, spec) {
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth || 600;
    const h = canvas.clientHeight || 220;
    canvas.width = w * dpr;
    canvas.height = h * dpr;
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const m = { l: 52, r: 12, t: 10, b: 34 };
    const [x0, x1] = spec.x;
    const [y0, y1] = spec.y;
    const X = (x) => m.l + ((x - x0) / (x1 - x0)) * (w - m.l - m.r);
    const Y = (y) => h - m.b - ((y - y0) / (y1 - y0)) * (h - m.t - m.b);
    ctx.font = '11px Arial, Helvetica, sans-serif';
    ctx.fillStyle = '#52606d';
    ctx.strokeStyle = '#e4e7eb';
    ctx.lineWidth = 1;
    (spec.vbands || []).forEach((b) => {
      ctx.fillStyle = b.color;
      ctx.fillRect(X(b.from), m.t, X(b.to) - X(b.from), h - m.t - m.b);
    });
    ctx.fillStyle = '#52606d';
    const ys = niceStep(y1 - y0);
    for (let v = Math.ceil(y0 / ys) * ys; v <= y1 + 1e-9; v += ys) {
      ctx.beginPath(); ctx.moveTo(m.l, Y(v)); ctx.lineTo(w - m.r, Y(v)); ctx.stroke();
      ctx.textAlign = 'right'; ctx.fillText(fmt(v), m.l - 6, Y(v) + 4);
    }
    const xs = niceStep(x1 - x0, 7);
    const xTicks = spec.xTicks || [];
    if (!spec.xTicks) for (let v = Math.ceil(x0 / xs) * xs; v <= x1 + 1e-9; v += xs) xTicks.push([v, fmt(v)]);
    xTicks.forEach(([v, label]) => { ctx.textAlign = 'center'; ctx.fillText(label, X(v), h - m.b + 14); });
    ctx.strokeStyle = '#9aa5b1';
    ctx.beginPath(); ctx.moveTo(m.l, m.t); ctx.lineTo(m.l, h - m.b); ctx.lineTo(w - m.r, h - m.b); ctx.stroke();
    if (spec.xLabel) { ctx.textAlign = 'center'; ctx.fillText(spec.xLabel, (m.l + w - m.r) / 2, h - 4); }
    if (spec.yLabel) {
      ctx.save(); ctx.translate(12, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center'; ctx.fillText(spec.yLabel, 0, 0); ctx.restore();
    }
    ctx.save();
    ctx.beginPath(); ctx.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); ctx.clip();
    (spec.hlines || []).forEach((l) => {
      ctx.strokeStyle = l.color || '#9aa5b1'; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(m.l, Y(l.y)); ctx.lineTo(w - m.r, Y(l.y)); ctx.stroke();
      ctx.setLineDash([]);
      if (l.label) { ctx.fillStyle = l.color || '#52606d'; ctx.textAlign = 'right'; ctx.fillText(l.label, w - m.r - 4, Y(l.y) - 4); }
    });
    spec.series.forEach((s) => {
      ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []);
      ctx.beginPath();
      s.pts.forEach(([x, y], i) => { if (i) ctx.lineTo(X(x), Y(y)); else ctx.moveTo(X(x), Y(y)); });
      ctx.stroke();
      (s.dots || []).forEach(([x, y]) => { ctx.fillStyle = s.color; ctx.beginPath(); ctx.arc(X(x), Y(y), 3.5, 0, 7); ctx.fill(); });
    });
    ctx.restore();
    ctx.setLineDash([]);
  }

  // Builds the canvas, legend, sliders and readout; calls draw(values) on every change.
  function scaffold(box, { controls, legend = [], note, draw }) {
    box.innerHTML = '';
    const canvas = document.createElement('canvas');
    box.append(canvas);
    if (legend.length) {
      const lg = document.createElement('div');
      lg.className = 'legend';
      lg.innerHTML = legend.map(([label, color, dashed]) => `<span><i style="background:${dashed ? `repeating-linear-gradient(90deg,${color} 0 4px,transparent 4px 7px)` : color}"></i>${label}</span>`).join('');
      box.append(lg);
    }
    const readout = document.createElement('div');
    readout.className = 'readout';
    const wrap = document.createElement('div');
    wrap.className = 'controls';
    const values = {};
    const outs = {};
    controls.forEach((c) => {
      const label = document.createElement('label');
      if (c.options) {
        label.innerHTML = `<span>${c.label}</span>`;
        const sel = document.createElement('select');
        c.options.forEach(([v, text]) => sel.add(new Option(text, v)));
        sel.value = c.value;
        values[c.id] = c.value;
        sel.addEventListener('input', () => { values[c.id] = sel.value; redraw(); });
        label.append(sel);
      } else {
        label.innerHTML = `<span>${c.label}: <b></b></span>`;
        const input = document.createElement('input');
        Object.assign(input, { type: 'range', min: c.min, max: c.max, step: c.step, value: c.value });
        values[c.id] = Number(c.value);
        outs[c.id] = label.querySelector('b');
        input.addEventListener('input', () => { values[c.id] = Number(input.value); redraw(); });
        label.append(input);
      }
      wrap.append(label);
    });
    box.append(readout, wrap);
    if (note) {
      const n = document.createElement('div');
      n.className = 'small';
      n.innerHTML = note;
      box.append(n);
    }
    function redraw() {
      controls.forEach((c) => { if (outs[c.id]) outs[c.id].textContent = c.show ? c.show(values[c.id]) : values[c.id]; });
      readout.innerHTML = draw(canvas, values) || '';
    }
    redraw();
    new ResizeObserver(() => redraw()).observe(canvas);
  }

  // Fixed-step RK4 used by every simulation below.
  function simulate(f, y0, t1, h, input) {
    let y = y0.slice();
    const out = [[0, y.slice()]];
    const steps = Math.round(t1 / h);
    const add = (a, b, s) => a.map((v, i) => v + s * b[i]);
    for (let i = 0; i < steps; i += 1) {
      const t = i * h;
      const w = input ? input(t) : null;
      const k1 = f(y, w, t);
      const k2 = f(add(y, k1, h / 2), w, t + h / 2);
      const k3 = f(add(y, k2, h / 2), w, t + h / 2);
      const k4 = f(add(y, k3, h), w, t + h);
      y = y.map((v, j) => v + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
      out.push([t + h, y.slice()]);
    }
    return out;
  }
  const col = (sol, j) => sol.map(([t, y]) => [t, y[j]]);
  const schedules = {
    constant: [(t) => 1, 'Constant 1 µM'],
    max: [(t) => 1.47, 'Constant 1.47 µM (highest measured)'],
    pulse: [(t) => (Math.floor(t / 2) % 2 === 0 ? 1 : 0), '2 days on / 2 days off at 1 µM'],
    stop: [(t) => (t < 6 ? 1 : 0), '1 µM for days 0–6, then none'],
    delayed: [(t) => (t >= 6 && t < 12 ? 1 : 0), '1 µM on days 6–12 only'],
    none: [() => 0, 'No drug'],
  };
  const scheduleOptions = (keys) => keys.map((k) => [k, schedules[k][1]]);

  const W = {};

  W.growth = (box, o) => scaffold(box, {
    legend: [['Logistic (2.6)', COLORS.a], ['Theta-logistic (2.7)', COLORS.b], ['Gompertz (2.8)', COLORS.c]],
    controls: [
      { id: 'r', label: 'Growth rate r (per day)', min: 0.1, max: 1.5, step: 0.01, value: o.r ?? 0.47 },
      { id: 'theta', label: 'θ (theta-logistic only)', min: 0.2, max: 8, step: 0.1, value: o.theta ?? 3.3 },
      { id: 'K', label: 'Carrying capacity K', min: 2000, max: 6000, step: 50, value: o.K ?? 3841 },
    ],
    note: 'Start: 100 cells per imaged area (the 30,000-cell seeding). Defaults are the A2780cis values quoted in Section 3.2.',
    draw(canvas, v) {
      const N0 = 100;
      const T = 14;
      const lg = simulate((y) => [v.r * y[0] * (1 - y[0] / v.K)], [N0], T, 0.05);
      const th = simulate((y) => [v.r * y[0] * (1 - (y[0] / v.K) ** v.theta)], [N0], T, 0.05);
      const gz = simulate((y) => [v.r * y[0] * Math.log(v.K / y[0])], [N0], T, 0.05);
      const bold = (k) => (o.show === k ? 3.2 : 1.6);
      plot(canvas, {
        x: [0, T], y: [0, 6200], xLabel: 'day', yLabel: 'cells per image',
        hlines: [{ y: v.K, label: 'K' }],
        series: [
          { pts: col(lg, 0), color: COLORS.a, width: bold('logistic') },
          { pts: col(th, 0), color: COLORS.b, width: bold('theta') },
          { pts: col(gz, 0), color: COLORS.c, width: bold('gompertz') },
        ],
      });
      const half = (sol) => (sol.find(([, y]) => y[0] >= v.K / 2) || [NaN])[0];
      return `Day the population reaches half of K: logistic <b>${half(lg).toFixed(1)}</b>, theta-logistic <b>${half(th).toFixed(1)}</b>, Gompertz <b>${half(gz).toFixed(1)}</b>.`;
    },
  });

  W.effect = (box, o) => scaffold(box, {
    legend: [['Dose u(t)', COLORS.dose, true], ['Drug effect E (2.3)', COLORS.a], ['Time ramp 1 − e^(−γt) from (2.2)', COLORS.b, true]],
    controls: [
      { id: 'k', label: 'Rate k_on = k_off = γ (per day)', min: 0.1, max: 3, step: 0.05, value: 0.6 },
      { id: 'sched', label: 'Dose schedule', value: o.schedule || 'stop', options: scheduleOptions(['constant', 'stop', 'pulse']) },
    ],
    note: 'With k_on = k_off, E settles at the dose itself. The dashed ramp is what the rejected time-based model (2.2) does: it keeps rising after the drug is gone.',
    draw(canvas, v) {
      const u = schedules[v.sched][0];
      const sol = simulate((y, w) => [v.k * w - v.k * y[0]], [0], 14, 0.05, u);
      const dose = sol.map(([t]) => [t, u(t)]);
      const ramp = sol.map(([t]) => [t, 1 - Math.exp(-v.k * t)]);
      plot(canvas, { x: [0, 14], y: [0, 1.2], xLabel: 'day', yLabel: 'µM equivalent', series: [
        { pts: dose, color: COLORS.dose, dash: [5, 4] }, { pts: ramp, color: COLORS.b, dash: [6, 4] }, { pts: col(sol, 0), color: COLORS.a, width: 2.6 },
      ] });
      return `Half-life of the effect after the drug stops: <b>${(Math.log(2) / v.k).toFixed(1)} days</b>.`;
    },
  });

  W.hill = (box) => scaffold(box, {
    legend: [['𝓗(E) = Eʰ / (EC50ʰ + Eʰ)', COLORS.a], ['measured doses', COLORS.res]],
    controls: [
      { id: 'h', label: 'Hill slope h', min: 0.5, max: 6, step: 0.1, value: 2 },
      { id: 'ec', label: 'EC50 (µM, fixed at 1 in the draft)', min: 0.3, max: 3, step: 0.05, value: 1 },
    ],
    draw(canvas, v) {
      const H = (e) => e ** v.h / (v.ec ** v.h + e ** v.h);
      const pts = Array.from({ length: 151 }, (_, i) => [i * 0.02, H(i * 0.02)]);
      const doses = [0.62, 1, 1.47];
      plot(canvas, { x: [0, 3], y: [0, 1], xLabel: 'drug effect E (µM equivalent at steady state)', yLabel: 'fraction of max kill',
        hlines: [{ y: 0.5, label: 'half-maximal' }],
        series: [{ pts, color: COLORS.a, width: 2.6 }, { pts: [], color: COLORS.res, dots: doses.map((d) => [d, H(d)]) }] });
      return `At the measured doses the kill runs at ${doses.map((d) => `<b>${Math.round(100 * H(d))}%</b> (${d} µM)`).join(', ')} of its maximum β.`;
    },
  });

  W.chain = (box, o) => scaffold(box, {
    legend: [['Dose u', COLORS.dose, true], ['A₁', COLORS.c], ['A₂', COLORS.b], ['A₃ = E (drives the kill)', COLORS.a]],
    controls: [
      { id: 'k', label: 'Transit rate k (per day)', min: 0.1, max: 3, step: 0.05, value: o.k ?? 0.7 },
      { id: 'sched', label: 'Dose schedule', value: 'stop', options: scheduleOptions(['constant', 'stop', 'pulse', 'delayed']) },
    ],
    note: 'Each compartment smooths and delays the one before it, so the kill arrives late and outlasts the dose. Mean delay through three compartments = 3/k days.',
    draw(canvas, v) {
      const u = schedules[v.sched][0];
      const sol = simulate((y, w) => [v.k * (w - y[0]), v.k * (y[0] - y[1]), v.k * (y[1] - y[2])], [0, 0, 0], 14, 0.05, u);
      plot(canvas, { x: [0, 14], y: [0, 1.2], xLabel: 'day', yLabel: 'µM equivalent', series: [
        { pts: sol.map(([t]) => [t, u(t)]), color: COLORS.dose, dash: [5, 4] },
        { pts: col(sol, 0), color: COLORS.c, width: 1.5 }, { pts: col(sol, 1), color: COLORS.b, width: 1.5 }, { pts: col(sol, 2), color: COLORS.a, width: 2.8 },
      ] });
      const peak = col(sol, 2).reduce((a, b) => (b[1] > a[1] ? b : a));
      return `Mean delay 3/k = <b>${(3 / v.k).toFixed(1)} days</b>; E peaks at <b>${peak[1].toFixed(2)}</b> on day <b>${peak[0].toFixed(1)}</b>.`;
    },
  });

  // Coculture, optionally with cisplatin: equations (2.10) and (2.24)-(2.27).
  function coculture(box, o) {
    const treated = !!o.treated;
    const controls = [
      { id: 'aNC', label: 'α_NC (A2780cis → A2780)', min: 0, max: 3, step: 0.01, value: 0.61 },
      { id: 'aCN', label: 'α_CN (A2780 → A2780cis)', min: 0, max: 3, step: 0.01, value: 1.51 },
      { id: 'f', label: 'Starting share of A2780', min: 0.1, max: 0.9, step: 0.05, value: 0.5, show: (x) => `${Math.round(x * 100)}%` },
    ];
    if (treated) {
      controls.push(
        { id: 'sched', label: 'Dose schedule', value: 'constant', options: scheduleOptions(['none', 'constant', 'max', 'pulse', 'delayed']) },
        { id: 'bN', label: 'β_N max kill of A2780', min: 0, max: 1.5, step: 0.01, value: 0.97 },
        { id: 'bC', label: 'β_C max kill of A2780cis (0 = scenario A)', min: 0, max: 1, step: 0.01, value: 0 },
        { id: 'dN', label: 'd_N extra A2780 death (illustrative)', min: 0, max: 0.4, step: 0.01, value: 0.24 },
        { id: 'k', label: 'transit rate k (illustrative)', min: 0.1, max: 2, step: 0.05, value: 0.3 },
        { id: 'h', label: 'Hill slope h (illustrative)', min: 0.5, max: 5, step: 0.1, value: 4 },
      );
    }
    scaffold(box, {
      legend: [['A2780 (N)', COLORS.sens], ['A2780cis (C)', COLORS.res], ['Total', COLORS.total, true]],
      controls,
      note: treated
        ? 'r_N = 0.82, r_C = 0.47, θ_C = 3.3, α values and β_N from Sections 3.2–3.5; K = 4,794 (Stage 3a, from the cost-function table). k, h and d_N are not reported in the draft yet (Appendix A), so they are set by hand (d_N = 0.24, k = 0.3, h = 4) so that constant 1 µM lands near the 3,893 cells reported in Section 3.9. Both lines share k and h here. For intuition only, not the fitted model.'
        : 'r_N = 0.82, r_C = 0.47, θ_C = 3.3 (Section 3.2), K = 4,794 (Stage 3a). Start: 100 cells per image. For intuition only.',
      draw(canvas, v) {
        const p = { rN: 0.82, rC: 0.47, th: 3.3, K: 4794 };
        const u = treated ? schedules[v.sched][0] : () => 0;
        const H = (e) => (e <= 0 ? 0 : e ** v.h / (1 + e ** v.h));
        const T = treated ? 14 : 30;
        const sol = simulate((y, w) => {
          const [N, C, a1, a2, a3] = y;
          const dN = p.rN * N * (1 - (N + v.aNC * C) / p.K) - (treated ? v.dN * N + v.bN * H(a3) * N : 0);
          const crowd = Math.max(0, (C + v.aCN * N) / p.K);
          const dC = p.rC * C * (1 - crowd ** p.th) - (treated ? v.bC * H(a3) * C : 0);
          const k = treated ? v.k : 0;
          return [dN, dC, k * (w - a1), k * (a1 - a2), k * (a2 - a3)];
        }, [100 * v.f, 100 * (1 - v.f), 0, 0, 0], T, 0.05, u);
        const N = col(sol, 0);
        const C = col(sol, 1);
        const tot = sol.map(([t, y]) => [t, y[0] + y[1]]);
        plot(canvas, { x: [0, T], y: [0, 5200], xLabel: 'day', yLabel: 'cells per image',
          hlines: treated ? [{ y: 0.8 * p.K, label: 'B_cancer = 0.8 K', color: '#b45309' }] : [{ y: p.K, label: 'K' }],
          series: [{ pts: tot, color: COLORS.total, dash: [6, 4], width: 1.5 }, { pts: N, color: COLORS.sens, width: 2.4 }, { pts: C, color: COLORS.res, width: 2.4 }] });
        const [, last] = sol[sol.length - 1];
        const total = last[0] + last[1];
        let dose = 0;
        if (treated) for (let t = 0; t < T; t += 0.05) dose += u(t) * 0.05;
        return `Day ${T}: <b>${Math.round(total).toLocaleString()}</b> cells, resistant fraction <b>${(last[1] / total).toFixed(2)}</b>${treated ? `, cumulative dose <b>${dose.toFixed(1)} µM·day</b>` : ''}.`;
      },
    });
  }
  W.coculture = coculture;

  W.eps = (box) => scaffold(box, {
    legend: [['ε_N(u), A2780', COLORS.sens], ['ε_C(u), A2780cis', COLORS.res]],
    controls: [
      { id: 'eN', label: 'ε_N at 1 µM (Stage 4)', min: 0, max: 10, step: 0.1, value: 7.7 },
      { id: 'eC', label: 'ε_C at 1 µM (Stage 4)', min: 0, max: 2, step: 0.01, value: 0.38 },
    ],
    draw(canvas, v) {
      const e = (s4) => Array.from({ length: 61 }, (_, i) => { const u = i * 0.025; return [u, 1 + (s4 - 1) * Math.min(u, 1)]; });
      plot(canvas, { x: [0, 1.5], y: [0, Math.max(2, v.eN + 0.5)], xLabel: 'cisplatin u (µM)', yLabel: 'low-resource strength ε',
        hlines: [{ y: 1, label: 'inherited unchanged' }],
        series: [{ pts: e(v.eN), color: COLORS.sens, width: 2.4 }, { pts: e(v.eC), color: COLORS.res, width: 2.4 }] });
      return `At 0.5 µM the low-resource effect is scaled by <b>${(1 + (v.eN - 1) * 0.5).toFixed(2)}</b> on A2780 and <b>${(1 + (v.eC - 1) * 0.5).toFixed(2)}</b> on A2780cis; above 1 µM it stays at the Stage 4 value.`;
    },
  });

  W.barrier = (box) => scaffold(box, {
    legend: [['Relaxed barrier b_μ(z) (2.45)', COLORS.a], ['Pure log barrier −μ ln z', COLORS.b, true]],
    controls: [
      { id: 'mu', label: 'μ', min: 0.01, max: 0.2, step: 0.005, value: 0.05 },
      { id: 'd', label: 'δ', min: 0.005, max: 0.3, step: 0.005, value: 0.02 },
    ],
    note: 'z is the remaining room before a limit: z = 1 means far from it, z = 0 means at it, z < 0 means past it.',
    draw(canvas, v) {
      const b = (z) => (z >= v.d ? -v.mu * Math.log(z) : v.mu * (0.5 * ((z - 2 * v.d) / v.d) ** 2 - 0.5 - Math.log(v.d)));
      const relaxed = Array.from({ length: 281 }, (_, i) => { const z = -0.2 + i * 0.005; return [z, b(z)]; });
      const pure = relaxed.filter(([z]) => z > 0.001).map(([z]) => [z, -v.mu * Math.log(z)]);
      plot(canvas, { x: [-0.2, 1.2], y: [-0.1, 1], xLabel: 'z (room left before the limit)', yLabel: 'cost',
        series: [{ pts: relaxed, color: COLORS.a, width: 2.6 }, { pts: pure, color: COLORS.b, dash: [6, 4], width: 1.6 }] });
      return `At the limit (z = 0) the relaxed barrier costs <b>${b(0).toFixed(2)}</b> instead of infinity; at z = δ both equal <b>${(-v.mu * Math.log(v.d)).toFixed(3)}</b>.`;
    },
  });

  W.toxicity = (box) => scaffold(box, {
    legend: [['T* without fasting', COLORS.res], ['T* with protection P', COLORS.a]],
    controls: [
      { id: 'u', label: 'Dose u (µM)', min: 0, max: 3, step: 0.01, value: 1.47 },
      { id: 'P', label: 'Fasting protection P', min: 0, max: 1, step: 0.01, value: 0 },
    ],
    note: 'a_H = 0.1, r_H = 1, s_B = 1, d_B = 1, p_H = 0.5 as listed under (2.33)–(2.36).',
    draw(canvas, v) {
      const T = (u, P) => (0.1 * (1 - 0.5 * P) * (1 - Math.exp(-u))) / 1;
      const xs = Array.from({ length: 121 }, (_, i) => i * 0.025);
      plot(canvas, { x: [0, 3], y: [0, 0.3], xLabel: 'dose u (µM)', yLabel: 'steady toxicity T*',
        hlines: [{ y: 0.25, label: 'B_tox = 0.25', color: '#b45309' }, { y: 0.1, label: 'ceiling a_H / r_H = 0.1' }],
        series: [{ pts: xs.map((u) => [u, T(u, 0)]), color: COLORS.res }, { pts: xs.map((u) => [u, T(u, v.P)]), color: COLORS.a, width: 2.4, dots: [[v.u, T(v.u, v.P)]] }] });
      return `T* = <b>${T(v.u, v.P).toFixed(3)}</b> (healthy tissue at ${(100 * (1 - T(v.u, v.P))).toFixed(1)}% of capacity). The limit of 0.25 can never be reached with these values.`;
    },
  });

  W.fasting = (box) => scaffold(box, {
    legend: [['Fasting v', COLORS.dose, true], ['Protection P (2.36)', COLORS.a], ['Body weight W (2.37)', COLORS.res]],
    controls: [
      { id: 'pat', label: 'Fasting pattern', value: 'one', options: [['one', 'One 48 h fast (days 4–6)'], ['two', 'Two 48 h fasts (days 4–6, 14–16)'], ['long', 'One 4-day fast (days 4–8)']] },
      { id: 'kP', label: 'k_P protection build-up (per day)', min: 0.2, max: 3, step: 0.05, value: 1 },
      { id: 'kW', label: 'k_W weight loss per fasting day', min: 0.002, max: 0.03, step: 0.001, value: 0.012 },
      { id: 'rW', label: 'r_W recovery (per day)', min: 0.02, max: 0.5, step: 0.01, value: 0.15 },
    ],
    note: 'The second axis is shared: P runs from 0 to 1, and W is shown as weight relative to baseline (1 = no loss).',
    draw(canvas, v) {
      const windows = { one: [[4, 6]], two: [[4, 6], [14, 16]], long: [[4, 8]] }[v.pat];
      const fast = (t) => (windows.some(([a, b]) => t >= a && t < b) ? 1 : 0);
      const sol = simulate((y, w) => [v.kP * (w - y[0]), -v.kW * w + v.rW * (1 - y[1]) * (1 - w)], [0, 1], 30, 0.05, fast);
      const minW = Math.min(...sol.map(([, y]) => y[1]));
      const maxP = Math.max(...sol.map(([, y]) => y[0]));
      plot(canvas, { x: [0, 30], y: [0, 1.05], xLabel: 'day', yLabel: 'P   /   W',
        vbands: windows.map(([a, b]) => ({ from: a, to: b, color: 'rgba(154,165,177,.18)' })),
        hlines: [{ y: 0.95, label: 'W_min = 0.95', color: '#b45309' }],
        series: [{ pts: col(sol, 0), color: COLORS.a, width: 2.4 }, { pts: col(sol, 1), color: COLORS.res, width: 2.4 }] });
      const flag = minW < 0.95 ? ' — <b style="color:#b91c1c">below the cachexia limit</b>' : '';
      return `Peak protection <b>${maxP.toFixed(2)}</b>; lowest weight <b>${(100 * minW).toFixed(1)}%</b> of baseline${flag}.`;
    },
  });

  W.bic = (box) => scaffold(box, {
    controls: [
      { id: 'n', label: 'Observations n', min: 10, max: 200, step: 1, value: 84 },
      { id: 'k1', label: 'Model A parameters k', min: 1, max: 10, step: 1, value: 2 },
      { id: 'rss1', label: 'Model A RSS', min: 1e5, max: 5e6, step: 1e4, value: 1.2e6, show: (x) => x.toLocaleString() },
      { id: 'k2', label: 'Model B parameters k', min: 1, max: 10, step: 1, value: 4 },
      { id: 'rss2', label: 'Model B RSS', min: 1e5, max: 5e6, step: 1e4, value: 1.0e6, show: (x) => x.toLocaleString() },
    ],
    note: 'Numbers are made up to show the trade-off; drag them to see when an extra parameter pays for itself.',
    draw(canvas, v) {
      const bic = (rss, k) => v.n * Math.log(rss / v.n) + k * Math.log(v.n);
      const a = bic(v.rss1, v.k1);
      const b = bic(v.rss2, v.k2);
      const fitA = v.n * Math.log(v.rss1 / v.n);
      const fitB = v.n * Math.log(v.rss2 / v.n);
      const lo = Math.min(fitA, fitB) - 20;
      const hi = Math.max(a, b) + 20;
      plot(canvas, { x: [0, 2.6], y: [lo, hi], yLabel: 'BIC (lower is better)', xTicks: [[0.6, 'Model A'], [2.0, 'Model B']], series: [
        { pts: [[0.6, lo], [0.6, fitA]], color: COLORS.a, width: 26 }, { pts: [[0.6, fitA], [0.6, a]], color: '#93c5fd', width: 26 },
        { pts: [[2.0, lo], [2.0, fitB]], color: COLORS.b, width: 26 }, { pts: [[2.0, fitB], [2.0, b]], color: '#fcd34d', width: 26 },
      ] });
      const best = Math.min(a, b);
      const within = [['A', a, v.k1], ['B', b, v.k2]].filter(([, x]) => x - best <= 3);
      const pick = within.sort((x, y) => x[2] - y[2] || x[1] - y[1])[0][0];
      return `Model A: BIC <b>${a.toFixed(1)}</b> (fit ${fitA.toFixed(1)} + penalty ${(v.k1 * Math.log(v.n)).toFixed(1)}); Model B: BIC <b>${b.toFixed(1)}</b> (fit ${fitB.toFixed(1)} + penalty ${(v.k2 * Math.log(v.n)).toFixed(1)}). ΔBIC = <b>${Math.abs(a - b).toFixed(1)}</b>. The draft's rule picks <b>Model ${pick}</b> (fewest parameters within 3 BIC of the lowest). Dark bar = fit term, light bar = parameter penalty.`;
    },
  });

  W.rk4 = (box) => scaffold(box, {
    legend: [['Exact logistic curve', COLORS.total], ['RK4 with step h', COLORS.a]],
    controls: [{ id: 'h', label: 'Step size h (days)', min: 0.05, max: 2, step: 0.05, value: 1 }],
    note: 'Logistic growth with r = 0.82, K = 3,950 (A2780, Section 3.2) has an exact solution, so the RK4 error can be measured directly. The draft uses h = 0.05.',
    draw(canvas, v) {
      const r = 0.82;
      const K = 3950;
      const N0 = 100;
      const exact = (t) => K / (1 + ((K - N0) / N0) * Math.exp(-r * t));
      const sol = simulate((y) => [r * y[0] * (1 - y[0] / K)], [N0], 14, v.h);
      const err = Math.max(...sol.map(([t, y]) => Math.abs(y[0] - exact(t))));
      plot(canvas, { x: [0, 14], y: [0, 4500], xLabel: 'day', yLabel: 'cells per image', series: [
        { pts: Array.from({ length: 141 }, (_, i) => [i * 0.1, exact(i * 0.1)]), color: COLORS.total, width: 1.5 },
        { pts: col(sol, 0), color: COLORS.a, width: 1.5, dots: col(sol, 0) },
      ] });
      return `Largest error over 14 days: <b>${err < 0.01 ? err.toExponential(1) : err.toFixed(2)}</b> cells (${(100 * err / K).toPrecision(2)}% of K) using <b>${Math.round(14 / v.h)}</b> steps.`;
    },
  });

  window.PAPER_WIDGETS = W;
})();
