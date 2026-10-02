(() => {
  const treatedHeading = [...document.querySelectorAll('.report-section > h2')]
    .find((heading) => heading.textContent.includes('Treated Mono-culture'));
  if (treatedHeading) {
    const grid = document.createElement('div');
    grid.className = 'treated-grid';
    [...treatedHeading.parentElement.querySelectorAll(':scope > article')].forEach((card) => grid.append(card));
    treatedHeading.after(grid);
  }

  const plots = [...document.querySelectorAll('svg[role="img"]')];
  const axisMaximum = (plot) => {
    const labels = [...plot.querySelectorAll('text[x="64.0"]')];
    return Number(labels.at(-1)?.textContent.replace(/,/g, '')) || 0;
  };
  const sharedMaximum = Math.max(...plots.map(axisMaximum));
  if (!sharedMaximum) return;

  plots.forEach((plot) => {
    const previousMaximum = axisMaximum(plot);
    if (!previousMaximum || previousMaximum === sharedMaximum) return;
    const rescale = (oldY) => 252 - ((252 - oldY) * previousMaximum / sharedMaximum);

    plot.querySelectorAll('polyline').forEach((line) => {
      const points = line.getAttribute('points').trim().split(/\s+/).map((point) => {
        const [x, y] = point.split(',').map(Number);
        return `${x.toFixed(1)},${rescale(y).toFixed(1)}`;
      });
      line.setAttribute('points', points.join(' '));
    });
    plot.querySelectorAll('circle[cy]').forEach((dot) => dot.setAttribute('cy', rescale(Number(dot.getAttribute('cy'))).toFixed(1)));
    plot.querySelectorAll('line[stroke]').forEach((bar) => {
      ['y1', 'y2'].forEach((attribute) => bar.setAttribute(attribute, rescale(Number(bar.getAttribute(attribute))).toFixed(1)));
    });

    [...plot.querySelectorAll('text[x="64.0"]')].forEach((label, index) => {
      const value = sharedMaximum * index / 4;
      label.textContent = Math.round(value).toLocaleString();
      label.setAttribute('y', (256 - index * 56).toFixed(1));
      const gridline = label.nextElementSibling;
      if (gridline?.tagName.toLowerCase() === 'line') {
        const y = (252 - index * 56).toFixed(1);
        gridline.setAttribute('y1', y);
        gridline.setAttribute('y2', y);
      }
    });
  });

  const style = document.createElement('style');
  style.textContent = `.treated-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.treated-grid .card{margin-top:0;background:#fff;border:1px solid #d9e2e6;border-radius:10px;padding:14px}.treated-grid h3{font-size:16px;margin:0}.treated-grid .meta{font-size:13px;min-height:38px}.treated-grid svg{max-height:210px}@media(max-width:560px){.treated-grid{grid-template-columns:1fr}}`;
  document.head.append(style);
})();
