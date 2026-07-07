import csv, os
from rdkit import Chem

path = os.path.join(os.path.expanduser("~"), "Downloads",
                    "Dissociation-Constants", "iupac_high-confidence_v2_3.csv")

keywords = ["formyl", "xanthen", "Xanthen", "azulene"]

with open(path, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        if len(row) < 21: continue
        name = row[12].strip()
        for kw in keywords:
            if kw.lower() in name.lower():
                smi = row[1].strip()
                pt = row[3].strip()
                pv = row[4].strip()
                assessment = row[8].strip()
                al = row[18].strip()
                print(f"pKa_type={pt}  pKa={pv}  label={al}  assessment={assessment}")
                print(f"  SMILES: {smi}")
                print(f"  Name: {name}")
                print()
                break
