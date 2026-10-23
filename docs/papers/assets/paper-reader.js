// Interactive reader for manuscript drafts: navigation, progress, equation notes,
// glossary tooltips and citation previews. Content comes from equation-notes.js
// (per paper) and paper-glossary.js (shared).
(() => {
  const doc = document.querySelector('.doc');
  if (!doc) return;
  const NOTES = window.EQUATION_NOTES || {};
  const GLOSSARY = window.PAPER_GLOSSARY || [];
  const WIDGETS = window.PAPER_WIDGETS || {};
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
  const store = {
    get(key) { try { return localStorage.getItem(`paper:${key}`); } catch { return null; } },
    set(key, value) { try { localStorage.setItem(`paper:${key}`, value); } catch { /* private mode */ } },
  };
  // Shrink any formula wider than the space it sits in, so every symbol stays visible without scrollbars.
  function fitMath(root = document) {
    root.querySelectorAll('mjx-container').forEach((c) => {
      c.style.fontSize = '';
      const math = c.firstElementChild;
      if (!math) return;
      const display = c.getAttribute('display') === 'true';
      const host = display ? c.parentElement : c.closest('p, li, td, th, dd, dt, h5, .tip, .readout, figcaption') || c.parentElement;
      let avail = host.clientWidth - (display ? 0 : 4);
      if (host.tagName === 'DT') avail = Math.max(120, (host.closest('.note')?.clientWidth || 600) * 0.45);
      const width = math.getBoundingClientRect().width;
      if (avail > 0 && width > avail) c.style.fontSize = `${Math.max(30, Math.floor((avail / width) * 98))}%`;
    });
  }
  const typeset = (nodes) => {
    if (window.MathJax?.typesetPromise) {
      return MathJax.typesetPromise(nodes).then(() => nodes.forEach((n) => fitMath(n))).catch(() => {});
    }
    return Promise.resolve();
  };
  const escapeHtml = (s) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  // Keep each equation's TeX source before MathJax replaces it.
  $$('.eq-body').forEach((body) => { body.dataset.tex = body.textContent.trim().replace(/^\\\[|\\\]$/g, ''); });

  // ---------- Sections, contents and status ----------
  const headings = $$('h1, h2, h3, h4', doc);
  const level = (h) => ({ H1: 1, H2: 1, H3: 2, H4: 3 }[h.tagName]);
  const statusColor = { Drafted: 'var(--drafted)', Now: 'var(--now)', Update: 'var(--update)', Pending: 'var(--pending)', 'Write last': 'var(--pending)' };
  const toc = $('.toc');
  const tocLinks = new Map();
  let lastLevel2 = null;
  headings.forEach((h) => {
    const li = document.createElement('li');
    li.className = `lvl-${level(h)}`;
    const a = document.createElement('a');
    a.href = `#${h.id}`;
    const label = document.createElement('span');
    label.textContent = h.firstChild?.textContent.trim() || h.textContent.trim();
    a.append(label);
    if (h.dataset.status) {
      const dot = document.createElement('span');
      dot.className = 'dot';
      dot.style.background = statusColor[h.dataset.status];
      dot.title = h.dataset.status;
      a.append(dot);
    }
    li.append(a);
    toc.append(li);
    tocLinks.set(h, a);
    if (level(h) === 2) lastLevel2 = li;
    if (level(h) === 3 && lastLevel2) li.dataset.parent = lastLevel2.querySelector('a').getAttribute('href');
    if (level(h) === 1) lastLevel2 = null;
  });

  // Draft status overview in the header.
  const statuses = $$('[data-status]', doc).map((h) => h.dataset.status);
  if (statuses.length) {
    const counts = {};
    statuses.forEach((s) => { counts[s] = (counts[s] || 0) + 1; });
    const bar = $('.status-bar');
    const legend = document.createElement('div');
    legend.className = 'status-legend';
    Object.entries(counts).forEach(([s, n]) => {
      const seg = document.createElement('span');
      seg.style.width = `${(100 * n) / statuses.length}%`;
      seg.style.background = statusColor[s];
      seg.title = `${s}: ${n} sections`;
      bar.append(seg);
      legend.insertAdjacentHTML('beforeend', `<span><i style="background:${statusColor[s]}"></i>${s} ${n}</span>`);
    });
    bar.after(legend);
  }

  // Words per section, for the time-left estimate.
  const WPM = 220;
  const wordsIn = (text) => (text.match(/\S+/g) || []).length;
  const sectionWords = headings.map((h, i) => {
    const range = document.createRange();
    range.setStartBefore(h);
    if (headings[i + 1]) range.setEndBefore(headings[i + 1]); else range.setEndAfter(doc.lastChild);
    return wordsIn(range.toString());
  });
  const totalWords = sectionWords.reduce((a, b) => a + b, 0);
  $('.read-time').textContent = `About ${Math.round(totalWords / WPM)} min read · ${$$('.eq-num').length} numbered equations, each with an explanation`;

  // ---------- Scroll: progress, active section ----------
  const fill = $('.progress-fill');
  const whereSection = $('.where-section');
  const whereTime = $('.where-time');
  let active = -1;
  const onScroll = () => {
    const max = document.documentElement.scrollHeight - innerHeight;
    fill.style.width = `${Math.min(100, (100 * scrollY) / Math.max(1, max))}%`;
    let idx = -1;
    for (let i = 0; i < headings.length; i += 1) {
      if (headings[i].getBoundingClientRect().top < 110) idx = i; else break;
    }
    if (idx === active) return;
    active = idx;
    headings.forEach((h, i) => {
      const a = tocLinks.get(h);
      a.classList.toggle('active', i === idx);
      a.classList.toggle('read', i < idx);
    });
    // Show the sub-subsections only under the section being read.
    let parentHref = null;
    for (let i = idx; i >= 0; i -= 1) {
      if (level(headings[i]) === 2) { parentHref = `#${headings[i].id}`; break; }
      if (level(headings[i]) === 1) break;
    }
    $$('.toc .lvl-3').forEach((li) => li.classList.toggle('open', li.dataset.parent === parentHref));
    const current = headings[idx];
    whereSection.textContent = current ? current.firstChild.textContent.trim() : 'Start';
    const left = sectionWords.slice(Math.max(0, idx)).reduce((a, b) => a + b, 0);
    whereTime.textContent = `${Math.max(1, Math.round(left / WPM))} min left`;
    const link = current && tocLinks.get(current);
    if (link && !$('#pane-toc').hidden) {
      const pane = $('#pane-toc');
      const r = link.getBoundingClientRect();
      const pr = pane.getBoundingClientRect();
      if (r.top < pr.top || r.bottom > pr.bottom) pane.scrollTop += r.top - pr.top - pr.height / 3;
    }
  };
  let ticking = false;
  addEventListener('scroll', () => {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(() => { ticking = false; onScroll(); });
  }, { passive: true });

  // ---------- Sidebar tabs and mobile drawer ----------
  const sidebar = $('#sidebar');
  const menuButton = $('.tb-menu');
  const closeDrawer = () => { sidebar.classList.remove('open'); menuButton.setAttribute('aria-expanded', 'false'); };
  $$('.side-tabs button').forEach((tab) => tab.addEventListener('click', () => {
    $$('.side-tabs button').forEach((t) => t.setAttribute('aria-selected', String(t === tab)));
    $$('.side-pane').forEach((p) => { p.hidden = p.id !== `pane-${tab.dataset.pane}`; });
  }));
  const showPane = (name) => $(`.side-tabs [data-pane="${name}"]`).click();
  menuButton.addEventListener('click', () => {
    const open = sidebar.classList.toggle('open');
    menuButton.setAttribute('aria-expanded', String(open));
  });
  sidebar.addEventListener('click', (e) => { if (e.target.closest('a')) closeDrawer(); });

  // ---------- Equation notes ----------
  const eqTargets = [
    ...$$('.eq').map((el) => ({ key: el.dataset.eq, el, numbered: !!$('.eq-num', el) })),
    ...$$('[data-eq-row]').map((el) => ({ key: el.dataset.eqRow, el, row: true })),
  ];
  const byKey = new Map(eqTargets.map((t) => [t.key, t]));
  const labelFor = (key) => {
    const t = byKey.get(key);
    if (!t) return key;
    if (t.row) return `Table ${key.split('-')[0]}, row ${key.split('-')[1]}`;
    return t.numbered ? `Eq. ${key}` : 'Unnumbered equation';
  };

  function renderNote(key) {
    const n = NOTES[key];
    const note = document.createElement('section');
    note.className = 'note';
    note.setAttribute('aria-label', `Explanation of ${labelFor(key)}`);
    if (!n) {
      note.innerHTML = `<div class="note-head"><h5>${labelFor(key)}</h5><button class="close" aria-label="Close">×</button></div><p>No explanation has been written for this equation yet.</p>`;
      return note;
    }
    const parts = [`<div class="note-head"><h5>${labelFor(key)} · ${n.title}</h5><button class="close" aria-label="Close">×</button></div>`];
    if (n.plain) parts.push(`<h6>In plain words</h6><p>${n.plain}</p>`);
    if (n.terms?.length) {
      parts.push('<h6>Term by term</h6><dl class="terms">');
      n.terms.forEach(([sym, meaning]) => parts.push(`<dt>\\(${sym}\\)</dt><dd>${meaning}</dd>`));
      parts.push('</dl>');
    }
    if (n.model) parts.push(`<h6>Where it sits in the model</h6><p class="connect">${n.model}</p>`);
    if (n.widget && WIDGETS[n.widget.type]) parts.push(`<h6>Try it</h6><div class="widget" data-widget="${n.widget.type}"></div>`);
    if (n.caveat) parts.push(`<h6>Reading carefully</h6><p class="caveat">${n.caveat}</p>`);
    if (n.related?.length) {
      parts.push('<h6>Related</h6><div class="related">');
      n.related.forEach((r) => {
        if (byKey.has(r)) parts.push(`<a href="#${byKey.get(r).el.id}" data-jump="${r}">${labelFor(r)}${NOTES[r] ? ` · ${NOTES[r].title}` : ''}</a>`);
      });
      parts.push('</div>');
    }
    const t = byKey.get(key);
    if (t && !t.row) parts.push('<div class="foot"><button type="button" class="copy">Copy LaTeX</button></div>');
    note.innerHTML = parts.join('');
    const widgetBox = $('.widget', note);
    if (widgetBox) {
      try { WIDGETS[n.widget.type](widgetBox, n.widget); } catch (err) { widgetBox.textContent = 'This plot could not be drawn.'; }
    }
    return note;
  }

  function texOf(eqEl) {
    const body = $('.eq-body', eqEl);
    return body.dataset.tex || '';
  }

  function setOpen(key, open, { scroll = false } = {}) {
    const t = byKey.get(key);
    if (!t) return;
    const button = t.row ? $(`.eq-info[data-eq="${key}"]`) : $('.eq-info', t.el);
    const isOpen = button.getAttribute('aria-expanded') === 'true';
    if (open === isOpen) {
      if (scroll) t.el.scrollIntoView({ behavior: 'smooth', block: 'center' });
      return;
    }
    button.setAttribute('aria-expanded', String(open));
    if (open) {
      const note = renderNote(key);
      note.dataset.for = key;
      if (t.row) {
        const tr = document.createElement('tr');
        tr.className = 'note-row';
        const td = document.createElement('td');
        td.colSpan = t.el.children.length;
        td.append(note);
        tr.append(td);
        t.el.after(tr);
        t.el.classList.add('row-open');
      } else {
        t.el.after(note);
        t.el.classList.add('open');
      }
      $('.close', note).addEventListener('click', () => { setOpen(key, false); button.focus(); });
      $('.copy', note)?.addEventListener('click', (e) => {
        navigator.clipboard?.writeText(texOf(t.el)).then(() => { e.target.textContent = 'Copied'; });
      });
      typeset([note]);
      if (scroll) t.el.scrollIntoView({ behavior: 'smooth', block: 'start' });
    } else {
      const note = $(`.note[data-for="${CSS.escape(key)}"]`);
      (t.row ? note?.closest('tr') : note)?.remove();
      t.el.classList.remove('open', 'row-open');
    }
  }

  doc.addEventListener('click', (e) => {
    const button = e.target.closest('.eq-info');
    if (button) {
      const key = button.dataset.eq || button.closest('.eq').dataset.eq;
      setOpen(key, button.getAttribute('aria-expanded') !== 'true');
      return;
    }
    const jump = e.target.closest('[data-jump]');
    if (jump) {
      e.preventDefault();
      goToEquation(jump.dataset.jump);
    }
  });

  function goToEquation(key) {
    const t = byKey.get(key);
    if (!t) return;
    setOpen(key, true);
    t.el.scrollIntoView({ behavior: 'smooth', block: 'start' });
    t.el.classList.remove('flash');
    void t.el.offsetWidth;
    t.el.classList.add('flash');
    history.replaceState(null, '', `#${t.el.id}`);
  }

  // Equations pane, grouped by the subsection each equation sits in.
  const eqList = $('.eq-list');
  let lastGroup = null;
  eqTargets
    .slice()
    .sort((a, b) => (a.el.compareDocumentPosition(b.el) & Node.DOCUMENT_POSITION_FOLLOWING ? -1 : 1))
    .forEach((t) => {
      const heading = headings.filter((h) => level(h) <= 2 && (h.compareDocumentPosition(t.el) & Node.DOCUMENT_POSITION_FOLLOWING)).pop();
      const group = heading ? heading.firstChild.textContent.trim() : '';
      if (group !== lastGroup) {
        eqList.insertAdjacentHTML('beforeend', `<li class="grp">${escapeHtml(group)}</li>`);
        lastGroup = group;
      }
      const num = t.row ? `T${t.key.split('-')[0]}·${t.key.split('-')[1]}` : (t.numbered ? `(${t.key})` : '—');
      const title = NOTES[t.key]?.title || '';
      eqList.insertAdjacentHTML('beforeend', `<li><a href="#${t.el.id}" data-jump="${t.key}"><span class="n">${num}</span><span class="t">${escapeHtml(title)}</span></a></li>`);
    });
  eqList.addEventListener('click', (e) => {
    const a = e.target.closest('[data-jump]');
    if (!a) return;
    e.preventDefault();
    closeDrawer();
    goToEquation(a.dataset.jump);
  });

  // ---------- Glossary ----------
  const entries = GLOSSARY.map((g, i) => ({ ...g, i }));
  const aliasMap = new Map();
  const patterns = [];
  entries.forEach((g) => [g.term, ...(g.aliases || [])].forEach((alias) => {
    const exact = alias !== alias.toLowerCase();
    aliasMap.set(exact ? alias : alias.toLowerCase(), g);
    patterns.push(alias);
  }));
  patterns.sort((a, b) => b.length - a.length);
  const esc = (s) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const termRe = patterns.length ? new RegExp(`(?<![\\w-])(${patterns.map(esc).join('|')})(?![\\w-])`, 'gi') : null;
  const firstUse = new Map();

  function wrapTerms() {
    if (!termRe) return;
    const seenInSection = new Set();
    const walker = document.createTreeWalker(doc, NodeFilter.SHOW_ELEMENT | NodeFilter.SHOW_TEXT, {
      acceptNode(node) {
        if (node.nodeType === 1) {
          if (node.matches('h1, h2, h3, h4')) { seenInSection.clear(); return NodeFilter.FILTER_SKIP; }
          if (node.matches('.math, .eq, a, script, style, .reference, figcaption, button, .note')) return NodeFilter.FILTER_REJECT;
          return NodeFilter.FILTER_SKIP;
        }
        return NodeFilter.FILTER_ACCEPT;
      },
    });
    const jobs = [];
    let node;
    // Section boundaries are tracked while walking, so collect matches in order.
    while ((node = walker.nextNode())) {
      if (node.nodeType !== 3) continue;
      const text = node.nodeValue;
      termRe.lastIndex = 0;
      let m;
      const hits = [];
      while ((m = termRe.exec(text))) {
        const g = aliasMap.get(m[1]) || aliasMap.get(m[1].toLowerCase());
        if (!g) continue;
        if (seenInSection.has(g.i)) continue;
        seenInSection.add(g.i);
        hits.push({ start: m.index, end: m.index + m[1].length, g });
      }
      if (hits.length) jobs.push({ node, hits });
    }
    jobs.forEach(({ node: textNode, hits }) => {
      const frag = document.createDocumentFragment();
      let last = 0;
      const text = textNode.nodeValue;
      hits.forEach(({ start, end, g }) => {
        frag.append(text.slice(last, start));
        const span = document.createElement('span');
        span.className = 'term';
        span.tabIndex = 0;
        span.dataset.term = g.i;
        span.textContent = text.slice(start, end);
        frag.append(span);
        if (!firstUse.has(g.i)) firstUse.set(g.i, span);
        last = end;
      });
      frag.append(text.slice(last));
      textNode.replaceWith(frag);
    });
  }
  wrapTerms();

  const termList = $('.term-list');
  entries.slice().sort((a, b) => a.term.localeCompare(b.term, undefined, { sensitivity: 'base' })).forEach((g) => {
    const sym = g.symbol ? ` <span class="sym">\\(${g.symbol}\\)</span>` : '';
    const jump = firstUse.has(g.i) ? ` <a href="#" data-first="${g.i}">first use</a>` : '';
    termList.insertAdjacentHTML('beforeend', `<div class="term-item" data-search="${escapeHtml([g.term, ...(g.aliases || []), g.def].join(' ').toLowerCase())}"><dt>${escapeHtml(g.term)}${sym}</dt><dd>${g.def}${jump}</dd></div>`);
  });
  termList.addEventListener('click', (e) => {
    const a = e.target.closest('[data-first]');
    if (!a) return;
    e.preventDefault();
    closeDrawer();
    const span = firstUse.get(Number(a.dataset.first));
    span.scrollIntoView({ behavior: 'smooth', block: 'center' });
    setTimeout(() => { span.focus(); showTip(span); }, 450);
  });
  const search = $('.term-search');
  search.addEventListener('input', () => {
    const q = search.value.trim().toLowerCase();
    $$('.term-item', termList).forEach((item) => { item.hidden = q && !item.dataset.search.includes(q); });
  });
  let glossaryTyped = false;
  $('.side-tabs [data-pane="terms"]').addEventListener('click', () => {
    if (!glossaryTyped) { glossaryTyped = true; typeset([termList]); }
  });

  // ---------- Tooltips for terms and citations ----------
  const tip = $('.tip');
  let tipFor = null;
  function showTip(el) {
    let content = '';
    if (el.classList.contains('term')) {
      const g = entries[Number(el.dataset.term)];
      content = `<b>${escapeHtml(g.term)}${g.symbol ? ` \\(${g.symbol}\\)` : ''}</b>${g.def}`;
    } else {
      const ref = document.getElementById(el.getAttribute('href').slice(1));
      if (!ref) return;
      content = `${ref.innerHTML}<span class="more">Click to jump to the reference list</span>`;
    }
    tip.innerHTML = content;
    tip.hidden = false;
    tipFor = el;
    const r = el.getBoundingClientRect();
    const w = Math.min(360, innerWidth - 24);
    tip.style.maxWidth = `${w}px`;
    const left = Math.max(12, Math.min(scrollX + r.left, scrollX + innerWidth - w - 12));
    tip.style.left = `${left}px`;
    tip.style.top = `${scrollY + r.bottom + 8}px`;
    typeset([tip]).then(() => {
      if (tipFor !== el) return;
      const tr = tip.getBoundingClientRect();
      if (tr.bottom > innerHeight - 70 && r.top > tr.height + 16) tip.style.top = `${scrollY + r.top - tr.height - 8}px`;
    });
  }
  function hideTip() { tip.hidden = true; tipFor = null; }
  const tipTarget = (e) => e.target.closest?.('.term, a.cite');
  document.addEventListener('mouseover', (e) => { const el = tipTarget(e); if (el && el !== tipFor) showTip(el); });
  document.addEventListener('mouseout', (e) => {
    const el = tipTarget(e);
    if (el && !el.contains(e.relatedTarget)) hideTip();
  });
  document.addEventListener('focusin', (e) => { const el = tipTarget(e); if (el) showTip(el); });
  document.addEventListener('focusout', (e) => { if (tipTarget(e)) hideTip(); });
  doc.addEventListener('click', (e) => {
    const term = e.target.closest('.term');
    if (term) { if (tipFor === term) hideTip(); else showTip(term); return; }
    const cite = e.target.closest('a.cite');
    if (cite) {
      const ref = document.getElementById(cite.getAttribute('href').slice(1));
      if (ref) {
        ref.classList.remove('flash');
        void ref.offsetWidth;
        ref.classList.add('flash');
      }
    }
  });
  addEventListener('scroll', () => { if (tipFor && !tipFor.matches(':focus')) hideTip(); }, { passive: true });

  // ---------- Toolbar and keys ----------
  const sectionHeads = headings.filter((h) => level(h) <= 2);
  const goSection = (dir) => {
    const y = 120;
    const tops = sectionHeads.map((h) => h.getBoundingClientRect().top);
    let target;
    if (dir > 0) target = sectionHeads[tops.findIndex((t) => t > y + 5)];
    else {
      const idx = tops.findLastIndex((t) => t < y - 5);
      target = sectionHeads[Math.max(0, idx - (tops[idx] > -40 ? 1 : 0))];
    }
    target?.scrollIntoView({ behavior: 'smooth' });
  };
  $('.tb-prev').addEventListener('click', () => goSection(-1));
  $('.tb-next').addEventListener('click', () => goSection(1));
  $('.tb-top').addEventListener('click', () => scrollTo({ top: 0, behavior: 'smooth' }));
  let size = Number(store.get('size')) || 17;
  const applySize = () => { document.documentElement.style.setProperty('--text-size', `${size}px`); store.set('size', size); };
  applySize();
  $('.tb-smaller').addEventListener('click', () => { size = Math.max(14, size - 1); applySize(); });
  $('.tb-larger').addEventListener('click', () => { size = Math.min(23, size + 1); applySize(); });
  const allNotes = $('.tb-notes');
  allNotes.addEventListener('click', () => {
    const on = allNotes.getAttribute('aria-pressed') !== 'true';
    allNotes.setAttribute('aria-pressed', String(on));
    const anchor = headings[Math.max(0, active)];
    const before = anchor?.getBoundingClientRect().top;
    eqTargets.forEach((t) => setOpen(t.key, on));
    if (anchor) scrollBy(0, anchor.getBoundingClientRect().top - before);
  });

  document.addEventListener('keydown', (e) => {
    if (e.target.matches('input, textarea, select') || e.metaKey || e.ctrlKey || e.altKey) {
      if (e.key === 'Escape') e.target.blur();
      return;
    }
    if (e.key === ']') goSection(1);
    else if (e.key === '[') goSection(-1);
    else if (e.key === '/') {
      e.preventDefault();
      showPane('terms');
      sidebar.classList.add('open');
      search.focus();
    } else if (e.key === 'Escape') {
      hideTip();
      closeDrawer();
      lightbox.hidden = true;
    }
  });

  // ---------- Figure zoom ----------
  const lightbox = $('.lightbox');
  doc.addEventListener('click', (e) => {
    const img = e.target.closest('img.zoomable');
    if (!img) return;
    $('img', lightbox).src = img.src;
    $('img', lightbox).alt = img.alt;
    lightbox.hidden = false;
  });
  lightbox.addEventListener('click', () => { lightbox.hidden = true; });

  // ---------- Start ----------
  onScroll();
  const whenTypeset = () => (window.MathJax?.startup?.promise
    ? MathJax.startup.promise.then(() => fitMath(document))
    : setTimeout(whenTypeset, 200));
  whenTypeset();
  let fitTimer;
  addEventListener('resize', () => { clearTimeout(fitTimer); fitTimer = setTimeout(() => fitMath(document), 150); });
  [$('.tb-smaller'), $('.tb-larger')].forEach((b) => b.addEventListener('click', () => setTimeout(() => fitMath(document), 50)));
  const hashKey = location.hash.startsWith('#eq-') && eqTargets.find((t) => `#${t.el.id}` === location.hash);
  if (hashKey) setTimeout(() => goToEquation(hashKey.key), 300);
})();
