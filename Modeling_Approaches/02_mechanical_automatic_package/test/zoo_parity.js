// JS ↔ Julia parity: every stage's optimal model simulated by the site's own code (site.js, the Explorer and
// Summary) against the Julia solution stored in the bundle (bundle.js "curve").
//   node test/zoo_parity.js            → prints the largest relative difference per stage, then PARITY OK / FAILED
// Run after every bundle export; the site must draw the same model the fits used (e.g. Norton–Simon, χ, ε).
const fs = require('fs'), vm = require('vm'), path = require('path');
const DOCS = path.join(__dirname, '..', '..', '..', 'docs', 'outputs', 'results');
global.window = global; global.addEventListener = () => {};
global.document = { write() {}, addEventListener() {}, querySelectorAll() { return []; }, getElementById() { return null; }, createElement() { return { style: {}, setAttribute() {}, appendChild() {} }; }, body: { appendChild() {} } };
vm.runInThisContext(fs.readFileSync(path.join(DOCS, 'bundle.js'), 'utf8'));
vm.runInThisContext(fs.readFileSync(path.join(DOCS, 'site.js'), 'utf8'));
const TOL = 2e-3;                         // relative, on max(|y|, 50 cells): RK4 in the browser vs the Julia solver
let worst = 0, failed = [];
for (const g of RESULTS.groups) {
  const o = g.optimal, c = o.curve;
  if (!c) continue;
  const res = simulate(o.sim, o.params.map(p => p.value), Math.max(...c.t), 0.05, 20);   // params are {name, value, ...}
  if (!res) { failed.push(g.id + ' (no simulation)'); continue; }
  let m = 0;
  for (const [s, ys] of Object.entries(c.y)) c.t.forEach((t, i) => {
    const d = Math.abs(at(res, +s, t) - ys[i]) / Math.max(Math.abs(ys[i]), 50); if (d > m) m = d;
  });
  console.log(`${g.id.padEnd(14)} max rel diff ${m.toExponential(2)}${o.sim.sys.lines.some(L => L.ns) ? '  (Norton–Simon)' : ''}`);
  worst = Math.max(worst, m); if (!(m < TOL)) failed.push(g.id);
}
console.log(failed.length ? `PARITY FAILED: ${failed.join(', ')}` : `PARITY OK (worst ${worst.toExponential(2)})`);
process.exit(failed.length ? 1 : 0);
