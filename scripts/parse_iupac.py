#!/usr/bin/env python3
"""
Parse IUPAC Dissociation-Constants CSV into Compound objects.

Reads the IUPAC high-confidence pKa dataset (v2.3) and extracts high-quality
aqueous acid dissociation constants suitable for the CChO acid ranking
question generator.

Filters applied:
  - assessment in (Reliable, Approximate, Uncertain)
  - cosolvent empty (water-only measurements)
  - T in (25, 20, '') — 25degC preferred
  - pka_type in acid dissociation types (pKa1, pKa, pKaH1, pKa2, ...)
  - acidity_label in (AH, A, '') — acid forms only
  - pka_value parseable as float in [-20, 60]

Deduplication per unique SMILES: prefers T=25 > T=20 > empty;
Reliable > Approximate > Uncertain; pKa1 > pKa > pKaH1 > pKa2 > ...;
and lowest pKa value among ties.

Usage:
    from parse_iupac import get_iupac_compounds
    compounds = get_iupac_compounds()
"""

import csv
import os
import sys


class Compound:
    """A single compound entry from the pKa database.

    Matches the Compound class in generate.py.
    """

    def __init__(self, smiles, name, pka_h2o=None, pka_dmso=None, category=None, pka_type=None):
        self.smiles = smiles
        self.name = name
        self.pka_h2o = pka_h2o
        self.pka_dmso = pka_dmso
        self.category = category
        self.pka_type = pka_type

    def best_pka(self):
        """Return the best available pKa value (prefer H2O)."""
        if self.pka_h2o is not None:
            return self.pka_h2o
        return self.pka_dmso

    def __repr__(self):
        return f"Compound({self.name}, pKa={self.best_pka()})"


# ── Constants ──────────────────────────────────────────────────────

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CSV_PATH = os.path.join(SKILL_DIR, "iupac_high-confidence_v2_3.csv")

# Column indices in the IUPAC CSV (v2.3)
COL_UNIQUE_ID = 0
COL_SMILES = 1
COL_INCHI = 2
COL_PKA_TYPE = 3
COL_PKA_VALUE = 4
COL_T = 5
COL_REMARKS = 6
COL_METHOD = 7
COL_ASSESSMENT = 8
COL_REF = 9
COL_REF_REMARKS = 10
COL_ENTRY_REMARKS = 11
COL_ORIGINAL_NAMES = 12
COL_NAME_CONTRIBUTORS = 13
COL_NUM_CONTRIBUTORS = 14
COL_NICKNAMES = 15
COL_SOURCE = 16
COL_PRESSURE = 17
COL_ACIDITY_LABEL = 18
COL_ORIGINAL_T = 19
COL_COSOLVENT = 20
MIN_COLS = 21

# pKa type precedence (lower rank = preferred)
PKA_TYPE_RANK = {
    'pKa1': 0,
    'pKa': 1,
    'pKaH1': 2,
    'pKa2': 3,
    'pKa3': 4,
    'pKa4': 5,
    'pKa5': 6,
}

# Assessment precedence (lower rank = preferred)
ASSESSMENT_RANK = {
    'Reliable': 0,
    'Approximate': 1,
    'Uncertain': 2,
}

# Temperature precedence (lower rank = preferred)
TEMP_RANK = {
    '25': 0,
    '20': 1,
    '': 2,
}

# Valid acidity labels for acid forms
VALID_ACIDITY_LABELS = {'AH', 'A', ''}

# Reasonable pKa bounds
PKA_MIN = -20
PKA_MAX = 60


# ── pKaH1 SMILES protonation helper ────────────────────────────────

def _is_azulene_core(smiles):
    """Check whether the SMILES contains an azulene core (fused 5+7 ring).
    Azulene and its derivatives undergo C-protonation on the ring,
    never on substituent heteroatoms.
    """
    # Azulene pattern: bridged fused 5+7 aromatic rings
    # Typical SMILES: c1ccc2cccc-2cc1 (unsubstituted)
    # Substituted: R-c1ccc2cccccc1-2 or similar
    return '-2' in smiles and ('cccccc' in smiles or 'cccc-2' in smiles)


def _protonate_smiles_for_pkah1(smiles):
    """If SMILES represents a neutral base whose conjugate-acid pKa
    is stored (pKaH1), attempt to add H+ to the most basic atom.

    Protonation order: N > carbonyl-O > alcohol/phenol-O > ether-O > S > P.
    Azulene cores are left untouched (C-protonation cannot be auto-detected).

    Returns (protonated_smiles, protonated) where protonated
    indicates whether the SMILES was changed.
    """
    if '+' in smiles:
        return smiles, False

    # Azulene core → C-protonation, cannot auto-protonate
    if _is_azulene_core(smiles):
        return smiles, False

    try:
        from rdkit import Chem
    except ImportError:
        return smiles, False

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, False

    # Build candidate list: (elem_priority, score, idx, symbol)
    # elem_priority ensures N always beats O/S/P
    max_valence = {7: 4, 8: 3, 16: 3, 15: 4}

    candidates = []
    for atom in mol.GetAtoms():
        anum = atom.GetAtomicNum()
        if anum not in max_valence:
            continue
        if atom.GetFormalCharge() != 0:
            continue

        heavy_bonds = sum(1 for bond in atom.GetBonds()
                          if bond.GetOtherAtom(atom).GetAtomicNum() != 1)
        if heavy_bonds >= max_valence[anum]:
            continue

        n_h = atom.GetTotalNumHs()
        score = -heavy_bonds + n_h * 2

        # Carbonyl oxygen (C=O) is more basic than alcohol/phenol OH
        if anum == 8:
            is_carbonyl = any(
                bond.GetBondTypeAsDouble() == 2.0
                and bond.GetOtherAtom(atom).GetAtomicNum() == 6
                for bond in atom.GetBonds()
            )
            if is_carbonyl:
                score += 3  # bonus over alcohol/phenol OH within same element

        elem_prio = {7: 0, 8: 1, 16: 1, 15: 2}[anum]
        candidates.append((-elem_prio, score, atom.GetIdx(), atom.GetSymbol()))

    if not candidates:
        return smiles, False

    # Sort: element priority first (N > O/S > P), then score within element
    candidates.sort(reverse=True)

    for _, _, best_idx, sym in candidates:
        try:
            mol_h = Chem.RWMol(mol)
            atom = mol_h.GetAtomWithIdx(best_idx)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            atom.SetFormalCharge(1)
            mol_h = mol_h.GetMol()
            Chem.SanitizeMol(mol_h)
            new_smiles = Chem.MolToSmiles(mol_h)
            if '+' in new_smiles:
                return new_smiles, True
        except Exception:
            continue

    return smiles, False


# ── Public API ─────────────────────────────────────────────────────

def get_iupac_compounds(csv_path=None):
    """Parse the IUPAC CSV and return deduplicated Compound objects.

    Parameters
    ----------
    csv_path : str or None
        Path to the IUPAC CSV file. Uses the default Downloads location if None.

    Returns
    -------
    list of Compound
        Filtered, deduplicated compounds with aqueous pKa values.
    """
    path = _resolve_path(csv_path)
    if path is None:
        return []

    entries_by_smiles = {}
    stats = _init_stats()

    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header

        for row in reader:
            stats['total_rows'] += 1
            entry = _process_row(row, stats)
            if entry is None:
                continue

            smiles = entry['smiles']
            if smiles not in entries_by_smiles:
                entries_by_smiles[smiles] = []
            entries_by_smiles[smiles].append(entry)

    compounds = _deduplicate(entries_by_smiles)
    _print_stats(stats, len(compounds))
    return compounds


# ── Internal helpers ───────────────────────────────────────────────

def _resolve_path(csv_path):
    """Find the CSV file, trying the supplied path, then skill directory."""
    if csv_path is not None and os.path.exists(csv_path):
        return csv_path

    if os.path.exists(CSV_PATH):
        return CSV_PATH

    print(f"ERROR: IUPAC CSV not found at {CSV_PATH}", file=sys.stderr)
    return None


def _init_stats():
    """Return a fresh statistics dict for tracking filter counts."""
    return {
        'total_rows': 0,
        'valid_assessment': 0,
        'water_only': 0,
        'valid_temp': 0,
        'acid_type': 0,
        'valid_acidity_label': 0,
        'valid_pka_value': 0,
        'valid_smiles': 0,
    }


def _process_row(row, stats):
    """Filter and parse a single CSV row.  Returns an entry dict or None."""
    if len(row) < MIN_COLS:
        return None

    smiles = row[COL_SMILES].strip()
    pka_type = row[COL_PKA_TYPE].strip()
    pka_value_str = row[COL_PKA_VALUE].strip()
    T = row[COL_T].strip()
    assessment = row[COL_ASSESSMENT].strip()
    acidity_label = row[COL_ACIDITY_LABEL].strip()
    cosolvent = row[COL_COSOLVENT].strip()
    original_names = row[COL_ORIGINAL_NAMES].strip()

    # 1) assessment filter
    if assessment not in ASSESSMENT_RANK:
        return None
    stats['valid_assessment'] += 1

    # 2) cosolvent filter (water only)
    if cosolvent:
        return None
    stats['water_only'] += 1

    # 3) temperature filter
    if T not in TEMP_RANK:
        return None
    stats['valid_temp'] += 1

    # 4) pKa type filter (acid dissociation constants)
    if pka_type not in PKA_TYPE_RANK:
        return None
    stats['acid_type'] += 1

    # 5) acidity label filter (acid forms)
    if acidity_label not in VALID_ACIDITY_LABELS:
        return None
    stats['valid_acidity_label'] += 1

    # 6) parse and validate pKa value
    try:
        pka_float = float(pka_value_str)
    except ValueError:
        return None
    if pka_float < PKA_MIN or pka_float > PKA_MAX:
        return None
    stats['valid_pka_value'] += 1

    # 7) non-empty SMILES
    if not smiles:
        return None
    stats['valid_smiles'] += 1

    # Build sort key for later deduplication
    sort_key = (
        TEMP_RANK[T],
        ASSESSMENT_RANK[assessment],
        PKA_TYPE_RANK[pka_type],
        pka_float,
    )

    return {
        'smiles': smiles,
        'pka': pka_float,
        'pka_type': pka_type,
        'name': original_names,
        'T': T,
        'assessment': assessment,
        'sort_key': sort_key,
    }


def _deduplicate(entries_by_smiles):
    """Per SMILES, keep the single best entry and build Compound objects."""
    compounds = []
    protonated_count = 0
    for smiles, entries in entries_by_smiles.items():
        entries.sort(key=lambda e: e['sort_key'])
        best = entries[0]
        pka_type = best['pka_type']

        # Protonate neutral SMILES for pKaH1 entries (conjugate acids)
        final_smiles = smiles
        if pka_type and pka_type.startswith('pKaH'):
            final_smiles, changed = _protonate_smiles_for_pkah1(smiles)
            if changed:
                protonated_count += 1

        compounds.append(Compound(
            smiles=final_smiles,
            name=best['name'],
            pka_h2o=best['pka'],
            pka_dmso=None,
            category='IUPAC',
            pka_type=pka_type,
        ))

    if protonated_count:
        print(f"  [pKaH1 protonation] {protonated_count} SMILES protonated")
    return compounds


def _print_stats(stats, n_compounds):
    """Print filtering statistics to stdout."""
    rows = stats['total_rows']
    print(f"IUPAC CSV: {rows} total rows")
    print(f"  After assessment filter: {stats['valid_assessment']} "
          f"({_pct(stats['valid_assessment'], rows)})")
    print(f"  After cosolvent filter:  {stats['water_only']} "
          f"({_pct(stats['water_only'], rows)})")
    print(f"  After temperature filter: {stats['valid_temp']} "
          f"({_pct(stats['valid_temp'], rows)})")
    print(f"  After pKa type filter:    {stats['acid_type']} "
          f"({_pct(stats['acid_type'], rows)})")
    print(f"  After acidity_label:      {stats['valid_acidity_label']} "
          f"({_pct(stats['valid_acidity_label'], rows)})")
    print(f"  After pKa value filter:   {stats['valid_pka_value']} "
          f"({_pct(stats['valid_pka_value'], rows)})")
    print(f"  With valid SMILES:        {stats['valid_smiles']} "
          f"({_pct(stats['valid_smiles'], rows)})")
    print(f"  Unique compounds (after dedup): {n_compounds}")


def _pct(part, whole):
    """Return a formatted percentage string."""
    if whole == 0:
        return "0.0%"
    return f"{100.0 * part / whole:.1f}%"


# ── CLI entry point ────────────────────────────────────────────────

if __name__ == '__main__':
    compounds = get_iupac_compounds()
    print()
    print("Sample compounds:")
    for c in compounds[:10]:
        print(f"  {c.smiles:30s} {c.name:45s} pKa={c.pka_h2o:6.1f}")
    print(f"\nTotal: {len(compounds)} IUPAC compounds loaded.")
