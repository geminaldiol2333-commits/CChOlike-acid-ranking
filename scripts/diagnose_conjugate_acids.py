#!/usr/bin/env python3
"""
Diagnose pKaH1 protonation failures.

Runs the same filtering + deduplication as parse_iupac.py, then for every
pKaH1 entry whose SMILES still lacks a '+' charge, diagnoses why the
_protonate_smiles_for_pkah1 function could not protonate it.

Output: a table of all failure cases with atom composition and a summary
breakdown by root cause.
"""

import csv
import os
import sys

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: rdkit not available. Will use string-based fallback.")


# ── Constants (mirrored from parse_iupac.py) ───────────────────────

CSV_PATH = os.path.join(
    os.path.expanduser("~"), "Downloads",
    "Dissociation-Constants", "iupac_high-confidence_v2_3.csv"
)

PKA_TYPE_RANK = {
    'pKa1': 0, 'pKa': 1, 'pKaH1': 2,
    'pKa2': 3, 'pKa3': 4, 'pKa4': 5, 'pKa5': 6,
}
ASSESSMENT_RANK = {
    'Reliable': 0, 'Approximate': 1, 'Uncertain': 2,
}
TEMP_RANK = {'25': 0, '20': 1, '': 2}
VALID_ACIDITY_LABELS = {'AH', 'A', ''}
PKA_MIN, PKA_MAX = -20, 60


def _resolve_csv():
    if os.path.exists(CSV_PATH):
        return CSV_PATH
    alt = os.path.normpath(os.path.join(
        os.path.dirname(__file__), "..", "..", "..",
        "Downloads", "Dissociation-Constants", "iupac_high-confidence_v2_3.csv"
    ))
    if os.path.exists(alt):
        return alt
    print(f"ERROR: CSV not found at {CSV_PATH} or {alt}", file=sys.stderr)
    return None


# ── Filtering (same logic as parse_iupac.py) ───────────────────────

def filter_and_collect(path):
    entries_by_smiles = {}
    stats = {'total_rows': 0, 'passed': 0}
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            stats['total_rows'] += 1
            if len(row) < 21:
                continue
            smiles = row[1].strip()
            pka_type = row[3].strip()
            pka_str = row[4].strip()
            T = row[5].strip()
            assessment = row[8].strip()
            acidity_label = row[18].strip()
            cosolvent = row[20].strip()
            name = row[12].strip()

            if assessment not in ASSESSMENT_RANK:
                continue
            if cosolvent:
                continue
            if T not in TEMP_RANK:
                continue
            if pka_type not in PKA_TYPE_RANK:
                continue
            if acidity_label not in VALID_ACIDITY_LABELS:
                continue
            try:
                pka = float(pka_str)
            except ValueError:
                continue
            if pka < PKA_MIN or pka > PKA_MAX:
                continue
            if not smiles:
                continue
            stats['passed'] += 1

            sort_key = (TEMP_RANK[T], ASSESSMENT_RANK[assessment],
                        PKA_TYPE_RANK[pka_type], pka)

            if smiles not in entries_by_smiles:
                entries_by_smiles[smiles] = []
            entries_by_smiles[smiles].append({
                'smiles': smiles, 'pka': pka, 'pka_type': pka_type,
                'name': name, 'T': T, 'assessment': assessment,
                'sort_key': sort_key,
            })
    return entries_by_smiles, stats


def deduplicate(entries_by_smiles):
    compounds = []
    for smiles, entries in entries_by_smiles.items():
        entries.sort(key=lambda e: e['sort_key'])
        best = entries[0]
        compounds.append(best)
    return compounds


# ── RDKit-based diagnostics ────────────────────────────────────────

def diagnose_rdkit(smiles):
    result = {
        'smiles': smiles,
        'parse_ok': False,
        'has_N': False, 'has_O': False, 'has_S': False,
        'has_P': False, 'has_Se': False, 'has_halogen': False,
        'N_count': 0, 'O_count': 0, 'S_count': 0, 'P_count': 0,
        'N_protonatable': [],   # list of (idx, symbol, bonds, h_count, charge)
        'O_protonatable': [],
        'S_protonatable': [],
        'P_protonatable': [],
        'sanitize_ok': False,
        'sanitize_error': '',
        'root_cause': '',
    }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result['root_cause'] = 'RDKit parse failed'
        return result
    result['parse_ok'] = True

    for atom in mol.GetAtoms():
        anum = atom.GetAtomicNum()
        sym = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        heavy_bonds = sum(1 for b in atom.GetBonds()
                          if b.GetOtherAtom(atom).GetAtomicNum() != 1)
        n_h = atom.GetTotalNumHs()

        if anum == 7:
            result['has_N'] = True
            result['N_count'] += 1
            if charge == 0 and heavy_bonds < 4:
                result['N_protonatable'].append(
                    (atom.GetIdx(), sym, heavy_bonds, n_h, charge))
        elif anum == 8:
            result['has_O'] = True
            result['O_count'] += 1
            if charge == 0 and heavy_bonds < 3:
                result['O_protonatable'].append(
                    (atom.GetIdx(), sym, heavy_bonds, n_h, charge))
        elif anum == 16:
            result['has_S'] = True
            result['S_count'] += 1
            if charge == 0 and heavy_bonds < 3:
                result['S_protonatable'].append(
                    (atom.GetIdx(), sym, heavy_bonds, n_h, charge))
        elif anum == 15:
            result['has_P'] = True
            result['P_count'] += 1
            if charge == 0 and heavy_bonds < 4:
                result['P_protonatable'].append(
                    (atom.GetIdx(), sym, heavy_bonds, n_h, charge))
        elif anum == 34:
            result['has_Se'] = True
        elif anum in (9, 17, 35, 53):
            result['has_halogen'] = True

    # Try the actual protonation logic from parse_iupac.py and see if
    # sanitization succeeds
    if result['N_protonatable']:
        # Score N atoms the same way parse_iupac does
        best_idx = None
        best_score = -999
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 7:
                continue
            if atom.GetFormalCharge() != 0:
                continue
            heavy_bonds = sum(1 for b in atom.GetBonds()
                              if b.GetOtherAtom(atom).GetAtomicNum() != 1)
            n_h = atom.GetTotalNumHs()
            score = -heavy_bonds + n_h * 2
            if score > best_score:
                best_score = score
                best_idx = atom.GetIdx()

        if best_idx is not None:
            try:
                mol_h = Chem.RWMol(mol)
                atom = mol_h.GetAtomWithIdx(best_idx)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(1)
                mol_h = mol_h.GetMol()
                Chem.SanitizeMol(mol_h)
                new_smiles = Chem.MolToSmiles(mol_h)
                if '+' in new_smiles:
                    result['sanitize_ok'] = True
                    result['root_cause'] = 'N protonation DID work on retry — check dedup logic'
                else:
                    result['sanitize_ok'] = False
                    result['sanitize_error'] = 'Sanitized OK but no + in result SMILES'
                    result['root_cause'] = 'Sanitize ok but SMILES lacks +'
            except Exception as e:
                result['sanitize_ok'] = False
                result['sanitize_error'] = str(e)[:200]
                result['root_cause'] = 'Sanitize failed after N protonation'
        else:
            result['root_cause'] = 'N atoms present but none protonatable (all N+ or N with 4 bonds)'
    elif result['O_protonatable'] or result['S_protonatable'] or result['P_protonatable']:
        result['root_cause'] = 'No nitrogen, but O/S/P sites available (protonation not handled)'
    else:
        result['root_cause'] = 'No protonatable atoms (all atoms saturated or charged)'

    return result


# ── String-based fallback (no RDKit) ────────────────────────────────

def diagnose_string(smiles):
    result = {
        'smiles': smiles,
        'parse_ok': False,
        'has_N': 'N' in smiles,
        'has_O': 'O' in smiles,
        'has_S': 'S' in smiles,
        'has_P': 'P' in smiles,
        'has_Se': 'Se' in smiles,
        'has_halogen': any(h in smiles for h in ('F', 'Cl', 'Br', 'I')),
        'N_count': smiles.count('N'),
        'O_count': smiles.count('O'),
        'S_count': smiles.count('S'),
        'P_count': smiles.count('P'),
        'N_protonatable': [], 'O_protonatable': [], 'S_protonatable': [],
        'P_protonatable': [],
        'sanitize_ok': False, 'sanitize_error': '',
        'root_cause': '(string-based — install rdkit for full diagnosis)',
    }
    if not result['has_N'] and not result['has_O'] and not result['has_S']:
        result['root_cause'] = 'No N/O/S atoms — possibly inorganic or exotic'
    elif not result['has_N']:
        result['root_cause'] = 'No nitrogen — protonation of O/S not handled'
    elif result['has_N']:
        result['root_cause'] = 'Has N — needs RDKit to check valence/sanitization'
    return result


# ── Main diagnosis ─────────────────────────────────────────────────

def main():
    path = _resolve_csv()
    if not path:
        sys.exit(1)

    print("=" * 100)
    print("DIAGNOSIS: pKaH1 protonation failures in parse_iupac.py")
    print("=" * 100)
    print()

    # Step 1: load and filter
    print("[1/3] Loading and filtering IUPAC CSV...")
    entries_by_smiles, stats = filter_and_collect(path)
    print(f"  Total rows: {stats['total_rows']}, passed filters: {stats['passed']}")

    # Step 2: deduplicate
    print("[2/3] Deduplicating per unique SMILES...")
    compounds = deduplicate(entries_by_smiles)
    print(f"  Unique compounds: {len(compounds)}")

    # Step 3: find pKaH1 entries and diagnose failures
    print("[3/3] Diagnosing pKaH1 entries...")
    diagnose = diagnose_rdkit if HAS_RDKIT else diagnose_string

    pkah_entries = [c for c in compounds
                    if c['pka_type'] and c['pka_type'].startswith('pKaH')]
    all_charged = [c for c in pkah_entries if '+' in c['smiles']]
    all_neutral = [c for c in pkah_entries if '+' not in c['smiles']]

    print(f"  Total pKaH1 entries (after dedup): {len(pkah_entries)}")
    print(f"    Already charged (+): {len(all_charged)}")
    print(f"    Neutral (can't protonate): {len(all_neutral)} ← THESE ARE THE PROBLEM")
    print()

    # Detailed diagnosis on neutral entries
    results = []
    for c in all_neutral:
        d = diagnose(c['smiles'])
        d['pka'] = c['pka']
        d['name'] = c['name']
        d['pka_type'] = c['pka_type']
        d['assessment'] = c['assessment']
        results.append(d)

    # ── Output: detailed table ──
    print(f"{'#':>3}  {'pKa':>6}  {'Assessment':>12}  {'N':>3} {'O':>3} {'S':>3} {'P':>3}  {'Root Cause'}")
    print("-" * 100)
    for i, r in enumerate(results, 1):
        n = r['N_count'] if HAS_RDKIT else int(r['has_N'])
        o = r['O_count'] if HAS_RDKIT else int(r['has_O'])
        s = r['S_count'] if HAS_RDKIT else int(r['has_S'])
        p = r['P_count'] if HAS_RDKIT else int(r['has_P'])
        root = r['root_cause'][:60] if r['root_cause'] else '???'
        print(f"{i:>3}  {r['pka']:>6.1f}  {r['assessment']:>12}  {n:>3} {o:>3} {s:>3} {p:>3}  {root}")
    print()

    # ── Output: per-compound detail ──
    print("=" * 100)
    print("PER-COMPOUND DETAIL")
    print("=" * 100)
    for i, r in enumerate(results, 1):
        print(f"\n--- [{i}] {r['name'][:80]} ---")
        print(f"  SMILES:     {r['smiles']}")
        print(f"  pKa:        {r['pka']:.1f}")
        print(f"  Type:       {r['pka_type']}")
        print(f"  Assessment: {r['assessment']}")
        print(f"  Root cause: {r['root_cause']}")
        if HAS_RDKIT:
            print(f"  Parse OK:   {r['parse_ok']}")
            if r['N_protonatable']:
                print(f"  N sites ({len(r['N_protonatable'])}): {r['N_protonatable']}")
            if r['O_protonatable']:
                print(f"  O sites ({len(r['O_protonatable'])}): {r['O_protonatable']}")
            if r['S_protonatable']:
                print(f"  S sites ({len(r['S_protonatable'])}): {r['S_protonatable']}")
            if r['P_protonatable']:
                print(f"  P sites ({len(r['P_protonatable'])}): {r['P_protonatable']}")
            if r['sanitize_error']:
                print(f"  Sanitize error: {r['sanitize_error']}")

    # ── Summary breakdown ──
    print("\n" + "=" * 100)
    print("SUMMARY BY ROOT CAUSE")
    print("=" * 100)
    from collections import Counter
    cause_counts = Counter(r['root_cause'] for r in results)
    for cause, count in cause_counts.most_common():
        print(f"  {count:>3}  {cause}")

    # Category breakdown
    atom_cats = Counter()
    for r in results:
        has_n = r['N_count'] > 0 if HAS_RDKIT else r['has_N']
        has_o = r['O_count'] > 0 if HAS_RDKIT else r['has_O']
        has_s = r['S_count'] > 0 if HAS_RDKIT else r['has_S']
        has_p = r['P_count'] > 0 if HAS_RDKIT else r['has_P']
        if not has_n and not has_o and not has_s:
            atom_cats['No N/O/S'] += 1
        elif not has_n:
            if has_o and not has_s:
                atom_cats['O only'] += 1
            elif has_s and not has_o:
                atom_cats['S only'] += 1
            elif has_o and has_s:
                atom_cats['O+S only'] += 1
            else:
                atom_cats['Other (no N)'] += 1
        elif has_n:
            atom_cats['Has N'] += 1

    print(f"\n  Atom-level categories:")
    for cat, count in atom_cats.most_common():
        print(f"    {count:>3}  {cat}")

    print(f"\n  Total neutral pKaH1 that need fixing: {len(results)}")


if __name__ == '__main__':
    main()
