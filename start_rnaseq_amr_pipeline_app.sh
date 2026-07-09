#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
host="${HOST:-127.0.0.1}"
port="${PORT:-8791}"
exec python3 rnaseq_amr_pipeline_app.py serve --host "$host" --port "$port"
