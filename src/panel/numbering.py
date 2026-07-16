"""Mature <-> precursor residue numbering for SOD1 variants.

Variant names (e.g. ``A4V``) use MATURE (153-residue) numbering. The precursor has
an extra N-terminal Met, so ``precursor_pos = mature_pos + mature_offset``. The
offset is a config value (``config/pipeline.yaml`` -> ``project.mature_offset``) and
must never be hardcoded elsewhere -- callers pass it in from the loaded config.

``tests/test_numbering.py`` guards this module (CLAUDE.md rule 3).
"""

from __future__ import annotations

import re

# Standard one-letter amino-acid codes.
_AA1 = frozenset("ACDEFGHIKLMNPQRSTVWY")
_VARIANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_variant(name: str) -> tuple[str, int, str]:
    """Parse a mature-numbered variant like ``A4V`` into ``(wt, pos, mut)``.

    Args:
        name: Variant string in ``<WT><pos><MUT>`` form (surrounding whitespace ok).

    Returns:
        Tuple of wild-type residue, mature position (1-based), mutant residue.

    Raises:
        ValueError: If the name is malformed or uses a non-standard amino acid.
    """
    match = _VARIANT_RE.match(name.strip())
    if match is None:
        raise ValueError(f"Malformed variant name: {name!r} (expected e.g. 'A4V')")
    wt, pos_str, mut = match.group(1), match.group(2), match.group(3)
    for aa in (wt, mut):
        if aa not in _AA1:
            raise ValueError(f"Non-standard amino acid {aa!r} in variant {name!r}")
    pos = int(pos_str)
    if pos < 1:
        raise ValueError(f"Residue position must be >= 1 in variant {name!r}")
    return wt, pos, mut


def mature_to_precursor(mature_pos: int, mature_offset: int) -> int:
    """Map a mature-numbered residue position to precursor numbering."""
    return mature_pos + mature_offset


def precursor_to_mature(precursor_pos: int, mature_offset: int) -> int:
    """Inverse of :func:`mature_to_precursor`."""
    return precursor_pos - mature_offset
