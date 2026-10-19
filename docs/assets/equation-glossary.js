// Click any equation to see what its variables mean.
// Every MathJax equation (and the few plain-text equations marked with class "eqg-text")
// gets a hover highlight; a click opens a small box under the equation listing each
// variable with its meaning. Meanings come from equation-glossary-data.js (loaded next to
// this file): a project-wide table plus per-report overrides where a symbol means something
// else. The script also adds itself to same-origin iframes, such as the Final Report's
// chapters. Add it to a page with <script defer src="<relative path>/assets/equation-glossary.js"></script>;
// docs/assets/add_equation_glossary.py does that for every report page.
(function () {
  if (window.__eqgInit) return; window.__eqgInit = true;

  var me = document.currentScript;
  var base = me && me.src ? me.src.replace(/[^/]*$/, '') : (window.__eqgBase || '');
  window.__eqgBase = base;

  // ---------- TeX -> variable keys ----------
  var GREEK = ('alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi pi ' +
    'varpi rho varrho sigma varsigma tau upsilon phi varphi chi psi omega Gamma Delta Theta Lambda Xi Pi Sigma ' +
    'Upsilon Phi Psi Omega ell').split(' ');
  var TEXTCMD = { text: 1, mathrm: 1, textrm: 1, operatorname: 1, mathit: 1, textit: 1, mbox: 1, textbf: 1, mathsf: 1, textsf: 1 };
  var FONTCMD = { mathbb: 1, mathbf: 1, mathcal: 1, boldsymbol: 1, bm: 1 };
  var ACCENT = { hat: 1, widehat: 1, bar: 1, overline: 1, tilde: 1, widetilde: 1, dot: 1, ddot: 1, vec: 1 };
  var ENVCMD = { begin: 1, end: 1 };
  var WORDFN = { logit: 1, softplus: 1, tanh: 1, sigmoid: 1, logistic: 1 };
  var UNITS = { d: 1, k: 1, M: 1, uM: 1, nM: 1, mL: 1, h: 1, day: 1, days: 1, mg: 1, kg: 1, cells: 1 };
  var ISGREEK = {}; GREEK.forEach(function (g) { ISGREEK[g] = 1; });

  function readGroup(s, i) {            // s[i] === '{' -> [content, next index]
    var depth = 0, j = i;
    for (; j < s.length; j++) {
      if (s[j] === '\\') { j++; continue; }
      if (s[j] === '{') depth++;
      else if (s[j] === '}') { depth--; if (depth === 0) return [s.slice(i + 1, j), j + 1]; }
    }
    return [s.slice(i + 1), s.length];
  }
  function readArg(s, i) {               // one TeX argument: {group}, \command or one character
    while (s[i] === ' ') i++;
    if (s[i] === '{') return readGroup(s, i);
    if (s[i] === '\\') { var m = /^\\([a-zA-Z]+|.)/.exec(s.slice(i)); return [m[0], i + m[0].length]; }
    return [s[i] || '', i + 1];
  }
  function unwrap(t) {                   // tidy a sub/superscript: drop \text{} wrappers, braces and spaces
    var prev;
    do { prev = t; t = t.replace(/\\(?:text|mathrm|textrm|operatorname|mathit|textit|mbox|textbf|mathsf|textsf)\s*\{([^{}]*)\}/g, '$1'); } while (t !== prev);
    t = t.replace(/\\[,;:! ]|~/g, '').replace(/\s+/g, '');
    while (/^\{[^{}]*\}$/.test(t)) t = t.slice(1, -1);
    return t;
  }
  function isLabel(sup) {                // a superscript that names a variant (N^{R}, d^{acc}) rather than a power
    return /^[A-Za-z]{1,6}$/.test(sup) && !/^(n|m|k|p|q|h|x|y|T|top)$/.test(sup);
  }

  // Returns [{key, sub, sup}] in reading order, one per distinct variable.
  function texKeys(tex) {
    var out = [], seen = {};
    function push(baseKey, sub, sup) {
      var key = baseKey + (sub ? '_{' + sub + '}' : '') + (sup ? '^{' + sup + '}' : '');
      if (!seen[key]) { seen[key] = 1; out.push(key); }
    }
    function scripts(s, i) {             // read any _ and ^ after a base
      var sub = '', sup = '', r;
      for (var n = 0; n < 3; n++) {
        while (s[i] === ' ') i++;
        if (s[i] === "'") { i++; continue; }
        if (s[i] === '_' && !sub) { r = readArg(s, i + 1); sub = unwrap(r[0]); i = r[1]; continue; }
        if (s[i] === '^' && !sup) { r = readArg(s, i + 1); sup = r[0]; i = r[1]; continue; }
        break;
      }
      return [sub, sup, i];
    }
    function walk(s) {
      var i = 0;
      while (i < s.length) {
        var c = s[i], r, sc, b = null;
        if (c === '\\') {
          var m = /^\\([a-zA-Z]+)/.exec(s.slice(i));
          if (!m) { i += 2; continue; }
          var name = m[1]; i += m[0].length;
          if (name === 'mu' && /^\s*(\\(text|mathrm)\s*\{\s*M\s*\}|M\b)/.test(s.slice(i))) {   // the unit µM
            i += /^\s*(\\(text|mathrm)\s*\{\s*M\s*\}|M)/.exec(s.slice(i))[0].length; continue;
          }
          if (ISGREEK[name]) b = '\\' + name;
          else if (TEXTCMD[name]) {
            r = readArg(s, i); i = r[1];
            var w = unwrap(r[0]);
            if (/^[A-Za-z]{2,}$/.test(w) && /^\s*[_^]/.test(s.slice(i))) b = w;     // \mathrm{EC}_{50}
            else if (/^[A-Za-z][A-Za-z0-9]*$/.test(w) && !UNITS[w] && w.length <= 4 && !/\s/.test(r[0])) b = w; // \mathrm{IC}, \text{BIC}
            else continue;                                                              // prose such as \text{ if } or \mathrm{d}
          } else if (FONTCMD[name]) { r = readArg(s, i); i = r[1]; b = '\\' + name + '{' + unwrap(r[0]) + '}'; }
          else if (ACCENT[name]) {
            r = readArg(s, i); i = r[1];
            var inner = texKeys(r[0]);
            if (inner.length === 1) b = '\\' + name + '{' + inner[0] + '}'; else { walk(r[0]); continue; }
          } else if (ENVCMD[name]) { r = readArg(s, i); i = r[1]; continue; }
          else continue;                 // \frac, \left, \exp, \max, spacing ... their arguments are walked as normal text
        } else if (/[A-Za-z]/.test(c)) {
          var ab = /[A-Za-z]$/.test(s.slice(0, i)) ? null : (/^((?:EC|IC|LD|ED)(?:\d+)?)(?![A-Za-z])/.exec(s.slice(i)) || /^(SSE|RMSE|BIC|AIC|SI)(?![A-Za-z])/.exec(s.slice(i)) || /^([a-z]{4,})(?![A-Za-z])/.exec(s.slice(i))); // EC_{50}, SSE_w are one symbol
          if (ab) { b = ab[1]; i += ab[1].length; }
          else if (c === 'e' && /[0-9.]$/.test(s.slice(0, i)) && /^e[-+]?\d/.test(s.slice(i))) { i++; continue; } // 2.9e-07
          else { b = c; i++; }
        }
        else { i++; continue; }
        sc = scripts(s, i); i = sc[2];
        var sub = sc[0], sup = sc[1], label = '';
        if (sup) { var u = unwrap(sup); if (isLabel(u)) label = u; else walk(sup); }
        if (sub && /\\/.test(sub) && !/^[A-Za-z0-9,\\{}_]+$/.test(sub)) { walk(sub); sub = ''; }
        if (b === 'd' && !sub && !label) continue;
        if (WORDFN[b]) continue;          // the d of dN/dt
        if (b === 'e' && !sub && sup && !label) continue;   // e^{...} is Euler's number
        push(b, sub, label);
      }
    }
    walk(tex);
    return out;
  }

  // Plain-text equations: identifiers such as r_i, k_clear_i, alpha_NC, kill_S, rho
  var TEXT_SKIP = { exp: 1, ln: 1, log: 1, max: 1, min: 1, sum: 1, softplus: 1, NN: 1, tanh: 1, if: 1, else: 1, same: 1, for: 1, and: 1, or: 1, with: 1, where: 1, only: 1, wells: 1, low: 1, resource: 1, left: 1, line: 1 };
  function textKeys(t) {
    var out = [], seen = {}, re = /[A-Za-zα-ωΑ-Ω][A-Za-z0-9α-ωΑ-Ω_]*/g, m;
    function add(w) { if (!seen[w]) { seen[w] = 1; out.push(w); } }
    while ((m = re.exec(t))) {
      var w = m[0], after = t.slice(re.lastIndex);
      if (w === 'dt' || TEXT_SKIP[w]) continue;
      if (/^\s*\//.test(after) && /^d[A-Za-zα-ωΑ-Ω]/.test(w) && /^\s*\/\s*dt/.test(after)) { add(w.slice(1)); continue; } // dN/dt
      add(w);
    }
    return out;
  }

  // ---------- glossary lookup ----------
  var G = null;
  function pageId() {
    if (window.__eqgPage) return window.__eqgPage;
    var p = decodeURIComponent(location.pathname || '');
    if (/srcdoc/.test(location.href) && G && G.titles && G.titles[document.title]) return G.titles[document.title];
    return p;
  }
  // G.scopes: [{match: [path fragments], symbols: {key: meaning}}], most specific first
  function tables() {
    var id = pageId(), list = [];
    ((G && G.scopes) || []).forEach(function (sc) {
      if (sc.match.some(function (m) { return id.indexOf(m) !== -1; })) list.push(sc.symbols);
    });
    if (G && G.common) list.push(G.common);
    return list;
  }
  function find(key, list) {
    for (var n = 0; n < list.length; n++) if (Object.prototype.hasOwnProperty.call(list[n], key)) return list[n][key] === null ? false : list[n][key];
    return null;
  }
  function variants(key) {               // A_{N,3} -> A_{N,j}; S_{N}^{R} -> S_{N}; x^{H} -> x
    var v = [key];
    var m = /^(.*?)(?:_\{([^}]*(?:\{[^}]*\}[^}]*)*)\})?(?:\^\{([^}]*)\})?$/.exec(key);
    if (!m) return v;
    var b = m[1], sub = m[2], sup = m[3];
    function mk(s, p) { return b + (s ? '_{' + s + '}' : '') + (p ? '^{' + p + '}' : ''); }
    if (sub) {
      var parts = sub.split(',');
      var gen = parts.map(function (p, idx) { return idx > 0 || parts.length === 1 ? p.replace(/^([1-9]|j-1|j\+1|i|k|n)$/, 'j') : p; }).join(',');
      if (gen !== sub) v.push(mk(gen, sup));
    }
    if (sup) { v.push(mk(sub, '')); }
    return v;
  }
  function meaning(key) {
    var list = tables(), vs = variants(key);
    for (var n = 0; n < vs.length; n++) {
      var d = find(vs[n], list);
      if (d === false) return { skip: true };     // listed as null: not a variable on this page
      if (d) return { text: d, generic: n > 0 && /_\{[^}]*j[^}]*\}/.test(vs[n]) };
    }
    return null;
  }

  // ---------- UI ----------
  var css = '' +
    '.eqg{cursor:pointer;border-radius:5px;transition:background-color .12s,box-shadow .12s}' +
    '.eqg:hover,.eqg:focus-visible{background:rgba(23,107,135,.10);box-shadow:0 0 0 2px rgba(23,107,135,.28);outline:none}' +
    '.eqg.eqg-open{background:rgba(23,107,135,.16);box-shadow:0 0 0 2px rgba(23,107,135,.55)}' +
    'mjx-container[display="true"].eqg{display:block;width:fit-content;max-width:100%;margin-left:auto!important;margin-right:auto!important;padding:2px 6px}' +
    '#eqg-pop{position:fixed;z-index:10000;display:none;max-width:min(460px,calc(100vw - 16px));max-height:min(360px,60vh);overflow:auto;' +
    'background:#fff;color:#17212b;border:1px solid #9fc3cf;border-top:3px solid #176b87;border-radius:8px;box-shadow:0 8px 24px rgba(0,0,0,.18);' +
    'font:13px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif;text-align:left;padding:8px 10px}' +
    '#eqg-pop .eqg-h{font-weight:700;color:#124b6e;margin:0 0 4px;font-size:12px;letter-spacing:.02em;text-transform:uppercase}' +
    '#eqg-pop table{border-collapse:collapse;width:100%;margin:0;background:none}' +
    '#eqg-pop td{border:0;border-top:1px solid #e3ecef;padding:4px 6px;vertical-align:top;background:none;font-size:13px;white-space:normal;text-align:left}' +
    '#eqg-pop tr:first-child td{border-top:0}' +
    '#eqg-pop td.eqg-s{white-space:nowrap;width:1%;padding-left:0}' +
    '#eqg-pop .eqg-none{color:#7b8794;font-style:italic}' +
    '#eqg-pop .eqg-gen{color:#52606d}';
  var style = document.createElement('style'); style.textContent = css;
  (document.head || document.documentElement).appendChild(style);

  var pop = null, cur = null;
  function texOf(el) {
    if (el.__eqgTex !== undefined) return el.__eqgTex;
    var t = null;
    try {
      var doc = window.MathJax && MathJax.startup && MathJax.startup.document;
      if (doc && doc.math) for (var it of doc.math) { if (it.typesetRoot === el) { t = it.math; break; } }
    } catch (e) { }
    return t;
  }
  function keysOf(el) {
    if (el.classList.contains('eqg-text')) {
      if (el.getAttribute('data-vars')) return el.getAttribute('data-vars').split(',');
      var plain = '';                    // r<sub>S</sub> reads as r_S
      (function walk(n) {
        n.childNodes.forEach(function (c) {
          if (c.nodeType === 3) plain += c.nodeValue;
          else if (c.tagName === 'SUB') plain += '_' + c.textContent + ' ';
          else if (c.tagName === 'BR') plain += '; ';
          else walk(c);
        });
      })(el);
      return textKeys(plain);
    }
    var tex = texOf(el);
    if (tex) return texKeys(tex);
    // fallback: read the MathML MathJax keeps for screen readers
    var out = [], seen = {};
    el.querySelectorAll('mjx-assistive-mml mi, math mi').forEach(function (n) {
      var k = n.textContent.trim();
      if (k && k !== 'd' && !seen[k]) { seen[k] = 1; out.push(k); }
    });
    return out;
  }
  function esc(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }
  function textSym(k) {                  // r_i -> r<sub>i</sub> for plain-text equations
    var m = /^([^_]+)_(.+)$/.exec(k);
    return m ? '<i>' + esc(m[1]) + '</i><sub>' + esc(m[2].replace(/_/g, ',')) + '</sub>' : '<i>' + esc(k) + '</i>';
  }
  function build(el) {
    var keys = keysOf(el), isText = el.classList.contains('eqg-text');
    var rows = [];
    keys.forEach(function (k) {
      var m = meaning(k);
      if (m && m.skip) return;
      var sym = isText ? textSym(k) : '\\(' + k + '\\)';
      var txt = m ? esc(m.text) : '<span class="eqg-none">not defined in this report</span>';
      rows.push('<tr><td class="eqg-s">' + (isText ? sym : esc(sym)) + '</td><td>' + txt + '</td></tr>');
    });
    if (!rows.length) rows.push('<tr><td class="eqg-none">This expression has no variables (only numbers or operators).</td></tr>');
    return '<div class="eqg-h">Variables</div><table>' + rows.join('') + '</table>';
  }
  function place() {
    if (!cur || !pop) return;
    var target = cur.querySelector('mjx-math, svg') || cur;
    var r = target.getBoundingClientRect();
    if (r.bottom < 0 || r.top > window.innerHeight) { hide(); return; }
    var w = pop.offsetWidth, h = pop.offsetHeight;
    var x = Math.min(Math.max(8, r.left), window.innerWidth - w - 8);
    var y = r.bottom + 6;
    if (y + h > window.innerHeight - 8 && r.top - h - 6 > 8) y = r.top - h - 6;
    pop.style.left = x + 'px'; pop.style.top = Math.max(8, y) + 'px';
  }
  function show(el) {
    if (!pop) {
      pop = document.createElement('div'); pop.id = 'eqg-pop'; pop.setAttribute('role', 'dialog'); pop.setAttribute('aria-label', 'Variables in this equation');
      document.body.appendChild(pop);
    }
    if (cur) { cur.classList.remove('eqg-open'); cur.setAttribute('aria-expanded', 'false'); }
    cur = el; el.classList.add('eqg-open'); el.setAttribute('aria-expanded', 'true');
    pop.innerHTML = build(el);
    pop.style.display = 'block'; place();
    if (!el.classList.contains('eqg-text') && window.MathJax && MathJax.typesetPromise) {
      MathJax.typesetPromise([pop]).then(place).catch(function () { });
    }
  }
  function hide() {
    if (pop) pop.style.display = 'none';
    if (cur) { cur.classList.remove('eqg-open'); cur.setAttribute('aria-expanded', 'false'); }
    cur = null;
  }

  function mark(root) {
    (root || document).querySelectorAll('mjx-container, .eqg-text').forEach(function (el) {
      if (el.__eqg || (pop && pop.contains(el)) || el.closest('#eqg-pop')) return;
      if (el.parentElement && el.parentElement.closest('mjx-container')) return;
      el.__eqg = true; el.classList.add('eqg');
      el.setAttribute('tabindex', '0'); el.setAttribute('role', 'button');
      el.setAttribute('aria-expanded', 'false');
      el.setAttribute('title', 'Click to see what the variables mean');
    });
  }

  document.addEventListener('click', function (e) {
    var el = e.target.closest ? e.target.closest('.eqg') : null;
    if (el && !(pop && pop.contains(el))) {
      e.preventDefault(); e.stopPropagation();
      if (cur === el) hide(); else show(el);
    } else if (!(pop && pop.contains(e.target))) hide();
  }, true);
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') { hide(); return; }
    var el = e.target.closest ? e.target.closest('.eqg') : null;
    if (el && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); if (cur === el) hide(); else show(el); }
  });
  window.addEventListener('scroll', place, true);
  window.addEventListener('resize', place);

  // equations appear when MathJax typesets (also later, e.g. when a dropdown or a slide is rendered)
  var pending = false;
  function soon() { if (pending) return; pending = true; setTimeout(function () { pending = false; mark(); }, 60); }
  function start() {
    mark();
    new MutationObserver(function (ms) {
      for (var n = 0; n < ms.length; n++) if (ms[n].addedNodes.length) { soon(); break; }
    }).observe(document.body, { childList: true, subtree: true });
    if (window.MathJax && MathJax.startup && MathJax.startup.promise) MathJax.startup.promise.then(soon);
  }

  // same-origin iframes (the Final Report loads each chapter in one)
  function adopt(fr) {
    try {
      var w = fr.contentWindow, d = fr.contentDocument;
      if (!w || !d || w.__eqgInit || !d.body) return;
      w.__eqgBase = base;
      var s = d.createElement('script'); s.src = base + 'equation-glossary.js'; d.body.appendChild(s);
    } catch (e) { }
  }
  document.addEventListener('load', function (e) { if (e.target && e.target.tagName === 'IFRAME') adopt(e.target); }, true);

  function load(cb) {
    if (window.EQ_GLOSSARY) { G = window.EQ_GLOSSARY; cb(); return; }
    if (window.parent !== window) { try { if (window.parent.EQ_GLOSSARY) { G = window.EQ_GLOSSARY = window.parent.EQ_GLOSSARY; cb(); return; } } catch (e) { } }
    var s = document.createElement('script'); s.src = base + 'equation-glossary-data.js';
    s.onload = function () { G = window.EQ_GLOSSARY || {}; cb(); };
    s.onerror = function () { G = {}; cb(); };
    (document.head || document.documentElement).appendChild(s);
  }
  function boot() { load(function () { start(); document.querySelectorAll('iframe').forEach(adopt); }); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot); else boot();

  window.__eqg = { texKeys: texKeys, textKeys: textKeys, keysOf: keysOf, meaning: meaning, pageId: pageId, mark: mark };
})();
