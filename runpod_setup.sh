#!/bin/bash
set -e

echo "=== Installing uv ==="
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"

echo "=== Setting up project ==="
cd /workspace/repurposing

echo "=== Installing dependencies ==="
uv sync

echo "=== Running 06_improve.py ==="
uv run python 06_improve.py

echo "=== Done! ==="
