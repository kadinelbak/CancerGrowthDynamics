$files = @(
  "C:\Users\elbak\OneDrive\Desktop\Research\Results\CancerGrowthDynamics\outputs\new_data_report\report.html",
  "C:\Users\elbak\OneDrive\Desktop\Research\Results\CancerGrowthDynamics\docs\outputs\new_data_report\report.html"
)

$previewPattern = '<h4>Sheet preview: (.*?)</h4><table class="preview">(.*?)</table>'
$previewReplacement = '<h4>Sheet preview: $1</h4><details class="preview-toggle"><summary>Expand preview for $1</summary><div class="preview-shell"><table class="preview">$2</table></div></details>'
$stylePattern = '\.back\{display:inline-block;margin:0 0 18px;color:#176b87;font-weight:700;text-decoration:none\}\.preview\{border-collapse:collapse;width:100%;font-size:12px\}'
$styleReplacement = '.back{display:inline-block;margin:0 0 18px;color:#176b87;font-weight:700;text-decoration:none}.preview-toggle{margin:10px 0 0}.preview-toggle>summary{display:inline-flex;align-items:center;gap:8px;cursor:pointer;list-style:none;padding:9px 13px;border:1px solid #c5d1d8;border-radius:999px;background:#f7fbfc;color:#175e79;font-weight:700}.preview-toggle>summary::-webkit-details-marker{display:none}.preview-toggle>summary::marker{content:""}.preview-toggle[open]>summary{background:#e8f3f6;border-color:#9fc4d0}.preview-shell{margin-top:10px}.preview{border-collapse:collapse;width:100%;font-size:12px}'

foreach ($f in $files) {
  $html = Get-Content -Raw $f
  $html = [regex]::Replace($html, $previewPattern, $previewReplacement, [System.Text.RegularExpressions.RegexOptions]::Singleline)
  $html = [regex]::Replace($html, $stylePattern, $styleReplacement)
  [System.IO.File]::WriteAllText($f, $html)
}

Write-Host "Updated preview toggles in $($files.Count) report files."
