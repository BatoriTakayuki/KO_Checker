$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

function Test-DockerReady {
    # Windows PowerShell can treat native stderr as an error even when redirected.
    $ErrorActionPreference = 'Continue'
    & docker info *> $null
    return ($LASTEXITCODE -eq 0)
}

try {
    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        throw 'Docker was not found. Install Docker Desktop, then run start.bat again.'
    }

    if (-not (Test-DockerReady)) {
        $desktopPath = Join-Path $env:ProgramFiles 'Docker\Docker\Docker Desktop.exe'
        if (-not (Test-Path -LiteralPath $desktopPath)) {
            throw 'Docker is not running. Start Docker Desktop in Linux containers mode, then try again.'
        }
        Write-Host 'Starting Docker Desktop...'
        Start-Process -FilePath $desktopPath -WindowStyle Hidden
        $dockerDeadline = (Get-Date).AddSeconds(120)
        do {
            Start-Sleep -Seconds 3
            if (Test-DockerReady) { break }
        } while ((Get-Date) -lt $dockerDeadline)
        if (-not (Test-DockerReady)) {
            throw 'Docker did not become ready within 120 seconds. Check Docker Desktop, then try again.'
        }
    }

    & docker compose version
    if ($LASTEXITCODE -ne 0) { throw 'Docker Compose is required. Check your Docker Desktop installation.' }

    Write-Host 'Starting KO_Checker (the first build may take several minutes)...'
    & docker compose up -d --build
    if ($LASTEXITCODE -ne 0) { throw 'KO_Checker could not start. See the Docker error above.' }

    Write-Host 'Waiting for KO_Checker...'
    $appDeadline = (Get-Date).AddSeconds(60)
    $ready = $false
    do {
        try {
            $response = Invoke-WebRequest -Uri 'http://127.0.0.1:8080/api/health' -UseBasicParsing -TimeoutSec 3
            if ($response.StatusCode -eq 200) { $ready = $true }
        } catch {
            # The container may still be starting.
        }
        if (-not $ready) { Start-Sleep -Seconds 2 }
    } while (-not $ready -and (Get-Date) -lt $appDeadline)

    if (-not $ready) {
        throw 'KO_Checker did not respond within 60 seconds. Run: docker compose logs --tail 100 rnako'
    }
    Start-Process 'http://localhost:8080'
    Write-Host 'KO_Checker is ready: http://localhost:8080'
    exit 0
} catch {
    Write-Host ("Error: " + $_.Exception.Message) -ForegroundColor Red
    exit 1
}
