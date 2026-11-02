$ErrorActionPreference='Stop'
$previewRoot=$PSScriptRoot
$staticRoot=Join-Path (Split-Path $previewRoot) 'KO_CHECKER\static'
$html=[IO.File]::ReadAllText((Join-Path $staticRoot 'index.html'))
$css=[IO.File]::ReadAllText((Join-Path $staticRoot 'style.css'))
$app=[IO.File]::ReadAllText((Join-Path $staticRoot 'app.js'))
$app=[regex]::Replace($app,'init\(\);\s*$','')
$structure=[IO.File]::ReadAllText((Join-Path $staticRoot 'structure.js'))
$sample=[IO.File]::ReadAllText((Join-Path $previewRoot 'sample.js'))
$html=[regex]::Replace($html,'<script src="[^"]+" defer></script>','')
$html=$html.Replace('<link rel="stylesheet" href="/style.css">',"<style>$css</style>")
$logo=[Convert]::ToBase64String([IO.File]::ReadAllBytes((Join-Path $staticRoot 'stem-cell-boys.png')))
$html=$html.Replace('src="/stem-cell-boys.png"',('src="data:image/png;base64,'+$logo+'"'))
$html=$html.Replace('</body>',"<script>$app`n$structure`n$sample</script></body>")
[IO.File]::WriteAllText((Join-Path $previewRoot 'completed.html'),$html,(New-Object Text.UTF8Encoding($false)))
