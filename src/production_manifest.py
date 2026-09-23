"""Deterministic Group A and Group B production manifests."""

from __future__ import annotations

from dataclasses import dataclass
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .production_config import ProductionConfig


GROUPS = ("A", "B_rho", "B_mu", "B_G")


@dataclass(frozen=True)
class ProductionCase:
    case_index: int
    case_id: str
    group: str
    fixed: dict[str, float]
    swept: tuple[str, ...]
    axes: dict[str, tuple[float, ...]]
    output_path: str
    config_sha256: str

    def to_json(self) -> dict[str, Any]:
        return {
            "case_index": self.case_index,
            "case_id": self.case_id,
            "group": self.group,
            "fixed": self.fixed,
            "swept": list(self.swept),
            "axes": {name: list(values) for name, values in self.axes.items()},
            "shape": [len(self.axes[name]) for name in self.swept],
            "point_count": int(np.prod([len(self.axes[name]) for name in self.swept])),
            "output_path": self.output_path,
            "config_sha256": self.config_sha256,
        }


def _axis(config: ProductionConfig, name: str) -> tuple[float, ...]:
    scan_name = "Ek" if name == "Ek" else "contrast"
    scan = config.data["scans"][scan_name]
    return tuple(
        float(value)
        for value in np.linspace(scan["minimum"], scan["maximum"], scan["samples"])
    )


def build_cases(config: ProductionConfig, group: str) -> tuple[ProductionCase, ...]:
    if group not in GROUPS:
        raise ValueError(f"unsupported production group: {group}")
    levels = config.levels
    ek = _axis(config, "Ek")
    cases: list[ProductionCase] = []
    if group == "A":
        for index, (ri, mi, gi) in enumerate(itertools.product(range(7), repeat=3)):
            case_id = f"a_r{ri:02d}_m{mi:02d}_g{gi:02d}"
            cases.append(ProductionCase(index, case_id, group,
                {"Arho": levels[ri], "Amu": levels[mi], "AG": levels[gi]},
                ("Ek",), {"Ek": ek}, f"group_a/{case_id}", config.sha256))
        return tuple(cases)

    definitions = {
        "B_rho": ("brho", "Arho", ("Amu", "AG"), ("m", "g")),
        "B_mu": ("bmu", "Amu", ("Arho", "AG"), ("r", "g")),
        "B_G": ("bg", "AG", ("Arho", "Amu"), ("r", "m")),
    }
    prefix, scanned, fixed_names, letters = definitions[group]
    contrast = _axis(config, scanned)
    for index, (first, second) in enumerate(itertools.product(range(7), repeat=2)):
        case_id = f"{prefix}_{letters[0]}{first:02d}_{letters[1]}{second:02d}"
        fixed = {fixed_names[0]: levels[first], fixed_names[1]: levels[second]}
        cases.append(ProductionCase(index, case_id, group, fixed,
            ("Ek", scanned), {"Ek": ek, scanned: contrast},
            f"group_b_{group[2:].lower()}/{case_id}", config.sha256))
    return tuple(cases)


def write_manifest(path: str | Path, cases: Iterable[ProductionCase]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    content = "".join(json.dumps(case.to_json(), sort_keys=True) + "\n" for case in cases)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(target)


def read_manifest(path: str | Path) -> tuple[dict[str, Any], ...]:
    with Path(path).open(encoding="utf-8") as handle:
        return tuple(json.loads(line) for line in handle if line.strip())


def bundle_records(records: tuple[dict[str, Any], ...], bundle_index: int, bundle_size: int):
    if bundle_index < 0 or bundle_size < 1:
        raise ValueError("invalid bundle selection")
    start = bundle_index * bundle_size
    return records[start : start + bundle_size]
