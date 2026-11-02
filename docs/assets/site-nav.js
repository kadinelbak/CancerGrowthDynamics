/* One top navigation bar for the main pages of the site (user, 2026-11-01: "make the top nav bar as synchronous
   as possible across summary, presentation, explorer, final report and the optimal control section").
   Loaded with <script defer src=".../assets/site-nav.js"> by docs/assets/add_site_nav.py. It replaces each page's
   own bar (nav.site-tabs, nav.topnav and the "Reports home" button) with the same tabs everywhere, marks the
   current one, and on the optimal-control pages adds a second row for the 14/30/60/120-day runs and the sweep.
   Pages shown inside another page (Final Report chapters, the version viewer) keep no bar. */
(function () {
  // a page shown inside one of our own pages (Final Report chapters, the version viewer) gets no bar; a page
  // shown inside someone else's frame (a hosting viewer, cross-origin) still does
  let inOurFrame = false; try { inOurFrame = !!window.frameElement; } catch (e) {}
  if (inOurFrame) return;
  const me = document.currentScript || [...document.scripts].find(s => /site-nav\.js/.test(s.src));
  const ROOT = new URL('..', me.src);                           // docs/
  const here = location.pathname;
  const R = p => new URL(p, ROOT).href;
  const PKG = 'Modeling_Approaches/02_mechanical_automatic_package/outputs/';
  const TABS = [
    ['home', 'Reports', 'index.html', p => p === ROOT.pathname || p === new URL(R('index.html')).pathname],
    ['summary', 'Summary', 'outputs/results/summary.html', p => /results\/summary\.html$/.test(p)],
    ['presentation', 'Presentation', 'outputs/results/presentation.html', p => /results\/(presentation|staged)\.html$/.test(p)],
    ['explorer', 'Explorer', 'outputs/results/explorer.html', p => /results\/explorer\.html$/.test(p)],
    ['control', 'Optimal control', PKG + 'control/report.html', p => /\/outputs\/control(_\w+)?\/report\.html$/.test(p)],
    ['final', 'Final report', 'outputs/final_report/history/index.html', p => /\/final_report\//.test(p)],
    ['papers', 'Papers', 'papers/index.html', p => /\/papers\//.test(p)],
    ['education', 'Education', 'education/index.html', p => /\/education\//.test(p)],
  ];
  const CONTROL = [
    ['14 days', PKG + 'control/report.html', /\/outputs\/control\/report\.html$/],
    ['30 days', PKG + 'control_30d/report.html', /control_30d\//],
    ['60 days', PKG + 'control_60d/report.html', /control_60d\//],
    ['120 days, dose + fasting', PKG + 'control_120d/report.html', /control_120d\//],
    ['Control sweep', PKG + 'control_sweep/report.html', /control_sweep\//],
  ];
  const css = `
nav.cgd-nav{position:sticky;top:0;z-index:50;display:flex;gap:2px;padding:0 16px;min-height:50px;align-items:stretch;background:#fff;border-bottom:1px solid #d9e2e6;font:700 14px/1.2 system-ui,-apple-system,"Segoe UI",sans-serif;overflow-x:auto;scrollbar-width:none;margin:0}
nav.cgd-nav::-webkit-scrollbar{display:none}
nav.cgd-nav a{display:flex;align-items:center;padding:0 11px;color:#52636b;text-decoration:none;border-bottom:3px solid transparent;white-space:nowrap}
nav.cgd-nav a:hover{color:#17212b}
nav.cgd-nav a[aria-current="page"]{color:#176b87;border-bottom-color:#176b87}
nav.cgd-nav .cgd-name{color:#17212b;padding-left:0;margin-right:10px;font-size:15px}
nav.cgd-sub{display:flex;gap:6px;flex-wrap:wrap;padding:8px 16px;background:#f4f7f8;border-bottom:1px solid #d9e2e6;font:650 13px system-ui,sans-serif}
nav.cgd-sub a{padding:4px 10px;border-radius:999px;background:#eef5f6;color:#176b87;text-decoration:none}
nav.cgd-sub a[aria-current="page"]{background:#176b87;color:#fff}
@media (max-width:760px){nav.cgd-nav .cgd-name{display:none}nav.cgd-nav{padding:0 8px}}
@media print{nav.cgd-nav,nav.cgd-sub{display:none}}`;
  function build() {
    const st = document.createElement('style'); st.textContent = css; document.head.appendChild(st);
    const nav = document.createElement('nav');
    nav.className = 'cgd-nav site-tabs'; nav.setAttribute('aria-label', 'Site');
    const cur = TABS.find(t => t[3](here));
    nav.innerHTML = `<a class="cgd-name" href="${R('index.html')}">Cancer Growth Dynamics</a>` +
      TABS.map(t => `<a href="${R(t[2])}"${cur === t ? ' aria-current="page"' : ''}>${t[1]}</a>`).join('');
    // drop the page's own bars
    document.querySelectorAll('nav.site-tabs:not(.cgd-nav), nav.topnav, a.home[href$="index.html"], a.back-home, a.nav-home').forEach(e => {
      const nav2 = e.closest('nav.topnav'); (nav2 || e).remove();
    });
    document.body.insertBefore(nav, document.body.firstChild);
    if (cur && cur[0] === 'control') {
      const sub = document.createElement('nav');
      sub.className = 'cgd-sub'; sub.setAttribute('aria-label', 'Optimal control runs');
      sub.innerHTML = CONTROL.map(c => `<a href="${R(c[1])}"${c[2].test(here) ? ' aria-current="page"' : ''}>${c[0]}</a>`).join('');
      nav.after(sub);
    }
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', build); else build();
})();
