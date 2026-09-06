param(
    [string]$ReportPath = "outputs\low_resource_staged_fits\report.html"
)

$html = [IO.File]::ReadAllText((Resolve-Path -LiteralPath $ReportPath))
$panels = [regex]::Matches($html, '<section>.*?</section>', [Text.RegularExpressions.RegexOptions]::Singleline) | ForEach-Object Value
$run1 = [Collections.Generic.List[string]]::new()
$run2 = [Collections.Generic.List[string]]::new()
$joint = [Collections.Generic.List[string]]::new()
$coculture = [Collections.Generic.List[string]]::new()

foreach ($panel in $panels) {
    if ($panel -match 'coculture') { $coculture.Add($panel) }
    elseif ($panel -match 'joint Run 1 \+ Run 2') { $joint.Add($panel) }
    elseif ($panel -match '— Run 2') { $run2.Add($panel) }
    else { $run1.Add($panel) }
}

$tabs = @"
<h2>Fits and BIC model comparisons</h2>
<div class='tabs' role='tablist' aria-label='Fit comparison groups'>
<button class='tab active' role='tab' aria-selected='true' data-tab='run1'>Run 1</button>
<button class='tab' role='tab' aria-selected='false' data-tab='run2'>Run 2</button>
<button class='tab' role='tab' aria-selected='false' data-tab='joint'>Joint</button>
<button class='tab' role='tab' aria-selected='false' data-tab='coculture'>Coculture</button>
</div>
<div class='tab-panel active' id='run1' role='tabpanel'>$($run1 -join "`n")</div>
<div class='tab-panel' id='run2' role='tabpanel'>$($run2 -join "`n")</div>
<div class='tab-panel' id='joint' role='tabpanel'>$($joint -join "`n")</div>
<div class='tab-panel' id='coculture' role='tabpanel'>$($coculture -join "`n")</div>
<script>for(const b of document.querySelectorAll('.tab'))b.addEventListener('click',()=>{document.querySelectorAll('.tab,.tab-panel').forEach(x=>{x.classList.remove('active');if(x.classList.contains('tab'))x.setAttribute('aria-selected','false')});b.classList.add('active');b.setAttribute('aria-selected','true');document.getElementById(b.dataset.tab).classList.add('active')})</script>
"@

$tabCss = '.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:12px 0}.tab{border:1px solid #176b87;background:#f7fbfd;color:#124b6e;border-radius:5px;padding:8px 14px;font-weight:700;cursor:pointer}.tab.active,.tab:hover{background:#176b87;color:#fff}.tab-panel{display:none}.tab-panel.active{display:block}'
$html = $html.Replace('</style>', "$tabCss</style>")
$html = [regex]::Replace($html, '<h2>Fits and BIC model comparisons</h2>.*?</body>', "$tabs</body>", [Text.RegularExpressions.RegexOptions]::Singleline)
[IO.File]::WriteAllText((Resolve-Path -LiteralPath $ReportPath), $html, [Text.UTF8Encoding]::new($false))
Write-Output "Organized $($panels.Count) fit panels: Run 1=$($run1.Count), Run 2=$($run2.Count), Joint=$($joint.Count), Coculture=$($coculture.Count)."
