(() => {
  const plots = [...document.querySelectorAll('svg[role="img"]')];
  const axisMaximum = (plot) => Number(plot.querySelector('text[x="8"]')?.textContent.replace(/,/g, '')) || 0;
  const sharedMaximum = Math.max(...plots.map(axisMaximum));
  if (!sharedMaximum) return;

  plots.forEach((plot) => {
    const previousMaximum = axisMaximum(plot);
    if (!previousMaximum || previousMaximum === sharedMaximum) return;
    const rescale = (oldY) => 265 - ((265 - oldY) * previousMaximum / sharedMaximum);
    plot.querySelectorAll('polyline').forEach((line) => {
      const points = line.getAttribute('points').trim().split(/\s+/).map((point) => {
        const [x, y] = point.split(',').map(Number);
        return `${x.toFixed(1)},${rescale(y).toFixed(1)}`;
      });
      line.setAttribute('points', points.join(' '));
    });
    plot.querySelectorAll('circle[cy]').forEach((dot) => dot.setAttribute('cy', rescale(Number(dot.getAttribute('cy'))).toFixed(1)));
    const label = plot.querySelector('text[x="8"]');
    label.textContent = sharedMaximum.toLocaleString(undefined, { maximumFractionDigits: 0 });
  });

  const ledger = document.querySelector('h2');
  if (ledger && !document.querySelector('#joint-model-equation')) {
    const detail = document.createElement('section');
    detail.id = 'joint-model-equation';
    detail.className = 'equation-card';
    detail.innerHTML = `
      <h2>Top joint co-culture contestant: shared competition coefficients</h2>
      <p><strong>Dataset pair:</strong> ce0_75S25R.xlsx + ce1_75S25R.xlsx. <strong>Model:</strong> joint competition logistic + constant treatment loss. <strong>BIC:</strong> 693.28 (the lowest LR-4 BIC among the three matched Ce0/Ce1 pairs).</p>
      <p>This is a single simultaneous fit of the untreated and treated trajectories. The growth, carrying-capacity, and competition coefficients below are <strong>identical</strong> in both conditions; only the loss terms apply to treated cells.</p>
      <table><tr><th>Parameter</th><th>Estimate</th><th>Meaning</th></tr>
      <tr><td>r<sub>S</sub></td><td>0.68067 day<sup>-1</sup></td><td>Sensitive intrinsic growth rate</td></tr>
      <tr><td>K<sub>S</sub></td><td>3,558 cells</td><td>Sensitive carrying capacity</td></tr>
      <tr><td>&alpha;<sub>SR</sub></td><td>0.99475</td><td>Resistant-cell competitive load on sensitive cells</td></tr>
      <tr><td>r<sub>R</sub></td><td>0.64428 day<sup>-1</sup></td><td>Resistant intrinsic growth rate</td></tr>
      <tr><td>K<sub>R</sub></td><td>13,602 cells</td><td>Resistant carrying capacity</td></tr>
      <tr><td>&alpha;<sub>RS</sub></td><td>5.0000</td><td>Sensitive-cell competitive load on resistant cells</td></tr>
      <tr><td>d<sub>S</sub></td><td>0.065911 day<sup>-1</sup></td><td>Treated-only sensitive loss</td></tr>
      <tr><td>d<sub>R</sub></td><td>2.8073e-13 day<sup>-1</sup> (~0)</td><td>Treated-only resistant loss</td></tr></table>
      <p><strong>Full fitted equation</strong> (S = sensitive, R = resistant; u = untreated, t = treated):</p>
      <pre>dS_u/dt = 0.68067 S_u [1 - (S_u + 0.99475 R_u)/3558]
dR_u/dt = 0.64428 R_u [1 - (R_u + 5.0000 S_u)/13602]

dS_t/dt = 0.68067 S_t [1 - (S_t + 0.99475 R_t)/3558] - 0.065911 S_t
dR_t/dt = 0.64428 R_t [1 - (R_t + 5.0000 S_t)/13602] - (2.8073e-13) R_t</pre>
      <p class="note"><strong>Boundary warning:</strong> &alpha;<sub>RS</sub> reached its imposed upper bound of 5 and d<sub>R</sub> is effectively zero. The shared-coefficient model wins this candidate set, but those two estimates are not precise biological quantities.</p>`;
    ledger.after(detail);
  }

  // The excluded Run 1 overlays remain in the document so they can be reviewed
  // one at a time without changing the joint-fit card or its BIC comparison.
  const excludedRunOneCards = new Set(['A2780 30k 1-1', 'A2780 30k 3-1']);
  document.querySelectorAll('#joint section').forEach((card) => {
    const heading = card.querySelector('h3');
    const key = heading?.textContent.split(' — ')[0];
    if (!excludedRunOneCards.has(key) || card.querySelector('.run-one-toggle')) return;

    const chart = card.querySelector('svg[role="img"]');
    const runOneLabel = [...chart.querySelectorAll('text')].find((text) => text.textContent.includes('Run 1'));
    if (!runOneLabel) return;
    const runOneColour = runOneLabel.getAttribute('fill');
    const runOneMarks = [...chart.querySelectorAll(`polyline[stroke="${runOneColour}"], circle[fill="${runOneColour}"]`)];
    const chartTitle = chart.querySelector('text[font-weight="bold"]');
    const originalTitle = chartTitle?.textContent;
    const control = document.createElement('button');
    control.className = 'run-one-toggle';

    const setRunOneVisible = (visible) => {
      runOneMarks.forEach((mark) => { mark.style.display = visible ? '' : 'none'; });
      runOneLabel.style.display = visible ? '' : 'none';
      if (chartTitle) chartTitle.textContent = visible ? originalTitle : "Joint fitted shape restored to Run 2's observed N0";
      control.textContent = visible ? 'Hide excluded Run 1' : 'Show excluded Run 1';
      control.setAttribute('aria-pressed', String(visible));
    };
    control.addEventListener('click', () => setRunOneVisible(control.getAttribute('aria-pressed') !== 'true'));
    heading.after(control);
    setRunOneVisible(false);
  });

  const style = document.createElement('style');
  style.textContent = `.equation-card{border:1px solid #ccd6dd;border-radius:10px;padding:16px;margin:20px 0;background:#fff}.equation-card h2{margin-top:0}.equation-card pre{overflow:auto;padding:12px;background:#f7fbfd;border:1px solid #ccd6dd;line-height:1.55}.run-one-toggle{border:1px solid #176b87;background:#f7fbfd;color:#124b6e;border-radius:5px;padding:6px 10px;margin:0 0 10px;font-weight:700;cursor:pointer}.run-one-toggle:hover{background:#176b87;color:#fff}.tab-panel.active{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.tab-panel.active section{margin:0;border:1px solid #ccd6dd;border-radius:10px;padding:14px;background:#fff}.tab-panel.active section h3{margin-top:0;font-size:16px}.tab-panel.active section svg{width:100%}@media(max-width:560px){.tab-panel.active{grid-template-columns:1fr}}`;
  document.head.append(style);
})();
