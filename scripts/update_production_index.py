#!/usr/bin/env python
"""Regenerate the campaign-level result index from case-local status files."""

from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.production_index import build_global_index  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("campaign_root", type=Path)
args = parser.parse_args()
print(json.dumps(build_global_index(args.campaign_root), indent=2, sort_keys=True))
