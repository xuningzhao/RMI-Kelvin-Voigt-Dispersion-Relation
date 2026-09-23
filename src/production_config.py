"""Versioned configuration for production parameter sweeps."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ProductionConfig:
    path: Path
    data: dict[str, Any]
    canonical_json: str
    sha256: str

    @property
    def levels(self) -> tuple[float, ...]:
        return tuple(float(value) for value in self.data["contrast_levels"])

    @property
    def solver_options(self) -> dict[str, Any]:
        return dict(self.data["solver"])

    @property
    def save_all_candidates(self) -> bool:
        return bool(self.data["storage"]["save_all_candidates"])


def _validate(data: dict[str, Any]) -> None:
    if data.get("schema_version") != 1:
        raise ValueError("unsupported production configuration schema")
    levels = data.get("contrast_levels")
    if not isinstance(levels, list) or len(levels) != 7:
        raise ValueError("contrast_levels must contain exactly seven values")
    numeric = [float(value) for value in levels]
    if numeric != sorted(numeric) or len(set(numeric)) != 7:
        raise ValueError("contrast_levels must be unique and increasing")
    for name in ("Ek", "contrast"):
        scan = data.get("scans", {}).get(name, {})
        if int(scan.get("samples", 0)) < 2:
            raise ValueError(f"{name} scan requires at least two samples")
        if float(scan["minimum"]) >= float(scan["maximum"]):
            raise ValueError(f"{name} scan bounds are invalid")
    execution = data.get("execution", {})
    for name in (
        "case_workers",
        "group_a_cases_per_worker",
        "point_workers",
    ):
        if int(execution.get(name, 0)) < 1:
            raise ValueError(f"{name} must be positive")
    if int(execution["point_workers"]) != 1:
        raise ValueError("point_workers must be one; nested point pools are unsupported")
    shape = execution.get("group_b_tile_shape")
    if not isinstance(shape, list) or len(shape) != 2 or any(int(x) < 1 for x in shape):
        raise ValueError("group_b_tile_shape must contain two positive integers")


def load_production_config(path: str | Path) -> ProductionConfig:
    source = Path(path)
    data = json.loads(source.read_text(encoding="utf-8"))
    _validate(data)
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":"))
    return ProductionConfig(
        path=source,
        data=data,
        canonical_json=canonical,
        sha256=hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
    )
