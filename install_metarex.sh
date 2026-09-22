#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -gt 0 && "${1:-}" != -* ]]; then
    PREFIX="$1"
    shift
else
    DEFAULT_PREFIX="${HOME}/metarex"
    read -r -p "Enter the installation directory [${DEFAULT_PREFIX}]: " PREFIX
    PREFIX="${PREFIX:-${DEFAULT_PREFIX}}"
fi

if [[ -z "${PREFIX}" ]]; then
    echo "Error: installation directory cannot be empty." >&2
    exit 1
fi

echo "MetaREX will be installed into: ${PREFIX}"
echo "The installation includes Conda environments, package cache, runtime files and databases."

exec python3 "${SCRIPT_DIR}/install_metarex.py" install \
    --prefix "${PREFIX}" \
    "$@"
