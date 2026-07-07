import csv, os, sys
from rdkit import Chem

path = os.path.join(os.path.expanduser("~"), "Downloads",
                    "Dissociation-Constants", "iupac_high-confidence_v2_3.csv")
PKA_TYPE_RANK = {"pKa1":0,"pKa":1,"pKaH1":2,"pKa2":3,"pKa3":4,"pKa4":5,"pKa5":6}
ASSESSMENT_RANK = {"Reliable":0,"Approximate":1,"Uncertain":2}
TEMP_RANK = {"25":0,"20":1,"":2}
VALID = {"AH","A",""}

ebs = {}
with open(path, "r", encoding="utf-8", errors="replace") as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) < 21: continue
        smi = row[1].strip(); pt = row[3].strip(); ps = row[4].strip()
        T = row[5].strip(); a = row[8].strip(); al = row[18].strip()
        co = row[20].strip(); nm = row[12].strip()
        if a not in ASSESSMENT_RANK: continue
        if co: continue
        if T not in TEMP_RANK: continue
        if pt not in PKA_TYPE_RANK: continue
        if al not in VALID: continue
        try: pv = float(ps)
        except: continue
        if pv < -20 or pv > 60: continue
        if not smi: continue
        if not pt.startswith("pKaH"): continue
        if "+" in smi: continue
        sk = (TEMP_RANK[T], ASSESSMENT_RANK[a], PKA_TYPE_RANK[pt], pv)
        if smi not in ebs: ebs[smi] = []
        ebs[smi].append({"smiles":smi,"pka":pv,"pka_type":pt,"name":nm,"sort_key":sk})

def has_N(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return False
    for a in mol.GetAtoms():
        if a.GetAtomicNum()==7 and a.GetFormalCharge()==0:
            hb = sum(1 for b in a.GetBonds() if b.GetOtherAtom(a).GetAtomicNum()!=1)
            if hb < 4: return True
    return False

def protonate(smi):
    if "+" in smi: return smi, False
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return smi, False
    mv = {7:4,8:3,16:3,15:4}
    cands = []
    for a in mol.GetAtoms():
        an = a.GetAtomicNum()
        if an not in mv: continue
        if a.GetFormalCharge() != 0: continue
        hb = sum(1 for b in a.GetBonds() if b.GetOtherAtom(a).GetAtomicNum()!=1)
        if hb >= mv[an]: continue
        nh = a.GetTotalNumHs()
        score = -hb + nh*2
        ep = {7:0,8:1,16:1,15:2}[an]
        cands.append((score, -ep, a.GetIdx(), a.GetSymbol()))
    if not cands: return smi, False
    cands.sort(reverse=True)
    for _, _, idx, sym in cands:
        try:
            mh = Chem.RWMol(mol)
            a = mh.GetAtomWithIdx(idx)
            a.SetNumExplicitHs(a.GetNumExplicitHs()+1)
            a.SetFormalCharge(1)
            mh = mh.GetMol()
            Chem.SanitizeMol(mh)
            ps = Chem.MolToSmiles(mh)
            if "+" in ps:
                return ps, True
        except:
            continue
    return smi, False

print("Compounds failing all protonation attempts:")
print("=" * 80)
n = 0
for smi, entries in ebs.items():
    entries.sort(key=lambda e: e["sort_key"])
    best = entries[0]
    if not best["pka_type"].startswith("pKaH"): continue
    if has_N(smi): continue
    ps, ok = protonate(smi)
    if not ok:
        n += 1
        mol = Chem.MolFromSmiles(smi)
        atoms = set()
        if mol:
            for a in mol.GetAtoms(): atoms.add(a.GetSymbol())
        print(f"\n{n}. {best['name'][:90]}")
        print(f"   SMILES: {smi}")
        print(f"   pKaH1: {best['pka']:.1f}")
        print(f"   Atoms: {sorted(atoms)}")
        if mol:
            for a in mol.GetAtoms():
                an = a.GetAtomicNum()
                fc = a.GetFormalCharge()
                hb = sum(1 for b in a.GetBonds() if b.GetOtherAtom(a).GetAtomicNum()!=1)
                nh = a.GetTotalNumHs()
                print(f"     {a.GetSymbol():>3} idx={a.GetIdx():>3} hbonds={hb} h={nh} charge={fc} valence={hb+nh+fc}")

print(f"\nTotal still failing: {n}")
