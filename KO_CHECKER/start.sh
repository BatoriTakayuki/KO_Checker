#!/usr/bin/env bash
# Linux server launcher: builds and starts KO_Checker, reachable from the same network.
# Usage: ./start.sh            (port 8083)
#        KO_PORT=9000 ./start.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

export KO_BIND="${KO_BIND:-0.0.0.0}"
export KO_PORT="${KO_PORT:-8083}"
# Analyses run at the same time and the default CPU threads of each, e.g. KO_PARALLEL=2 KO_THREADS=32 ./start.sh
export KO_PARALLEL="${KO_PARALLEL:-4}"
export KO_THREADS="${KO_THREADS:-16}"
# RAM work area (tmpfs): 50% of physical memory unless set, e.g. KO_RAM_SIZE=200g ./start.sh
# (TE_Checker on the same server uses its own tmpfs, so both are kept below the physical memory together).
export KO_RAM_SIZE="${KO_RAM_SIZE:-$(awk '/MemTotal/ {printf "%dg", $2 * 0.5 / 1024 / 1024}' /proc/meminfo)}"

die() { echo -e "\033[31mError: $*\033[0m" >&2; exit 1; }

command -v docker >/dev/null || die 'Docker was not found. Install Docker Engine, then run start.sh again.'

# Use sudo when the current user is not in the docker group.
DOCKER=(docker)
if ! docker info >/dev/null 2>&1; then
    if sudo docker info >/dev/null 2>&1; then
        DOCKER=(sudo --preserve-env=KO_BIND,KO_PORT,KO_RAM_SIZE,KO_PARALLEL,KO_THREADS,KO_REFERENCES_DIR,KO_DATA_DIR,KO_RESULTS_DIR,RNAKO_DATA_DIR docker)
    else
        die 'Docker is not running or not accessible. Try: sudo systemctl start docker'
    fi
fi

"${DOCKER[@]}" compose version >/dev/null || die 'Docker Compose is required (docker-compose-plugin).'

# Folders may be moved to another disk in .env (KO_REFERENCES_DIR, KO_DATA_DIR, KO_RESULTS_DIR).
if [ -f .env ]; then set -a; . ./.env; set +a; fi
FOLDERS=("${KO_REFERENCES_DIR:-references}" "${KO_RESULTS_DIR:-results}" "${KO_DATA_DIR:-workdata}" inputs)
# Folders on a root-owned disk (e.g. /mnt/...) are created with sudo and handed to the current user.
if ! mkdir -p "${FOLDERS[@]}" 2>/dev/null; then
    sudo mkdir -p "${FOLDERS[@]}" || die 'Could not create the data folders set in .env.'
    sudo chown "$(id -u):$(id -g)" "${FOLDERS[@]}" || true
fi

echo 'Starting KO_Checker (the first build may take several minutes)...'
"${DOCKER[@]}" compose up -d --build || die 'KO_Checker could not start. See the Docker error above.'

echo 'Waiting for KO_Checker...'
for _ in $(seq 30); do
    if curl -fsS -m 3 "http://127.0.0.1:${KO_PORT}/api/health" >/dev/null 2>&1; then
        echo
        echo "KO_Checker is ready:"
        echo "  This server:    http://localhost:${KO_PORT}"
        if [ "$KO_BIND" = 0.0.0.0 ]; then
            for ip in $(hostname -I); do
                case "$ip" in 172.1[7-9].*|172.2[0-9].*|172.3[01].*|*:*) continue ;; esac  # skip Docker bridges / IPv6
                echo "  Same network:   http://${ip}:${KO_PORT}"
            done
        fi
        echo
        echo "RAM work area: ${KO_RAM_SIZE} (used only by jobs with the RAM option)"
        echo "Analyses at the same time: ${KO_PARALLEL}, default CPU threads per analysis: ${KO_THREADS}"
        echo "Stop: ${DOCKER[*]} compose stop    Logs: ${DOCKER[*]} compose logs --tail 100 rnako"
        exit 0
    fi
    sleep 2
done
die "KO_Checker did not respond within 60 seconds. Run: ${DOCKER[*]} compose logs --tail 100 rnako"
