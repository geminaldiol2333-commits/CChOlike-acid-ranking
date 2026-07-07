# CChO Acid Ranking Question Generator

Generate Chinese Chemistry Olympiad (CChO) style acidity ranking questions with
ACS 1996 standard chemical structures embedded in DOCX files.

## Trigger

Use when the user asks to:
- "生成酸性排序题"
- "generate acidity ranking questions"
- "出CChO酸性比较题"
- "make acid ranking problems"

## Material Library

The skill uses two data sources merged together:

1. **酸性数据.md** (~200 hand-curated compounds): Bordwell/Evans pKa tables with
   SMILES, names, and both H₂O and DMSO pKa values.

2. **IUPAC Dissociation-Constants** (~7,300 compounds after filtering):
   High-confidence aqueous pKa data from the IUPAC SC-Database. Automatically
   loaded from `~/Downloads/Dissociation-Constants/iupac_high-confidence_v2_3.csv`
   if the repository has been cloned.

Combined library: **~7,500 compounds** with aqueous pKa values.

## Output Format

- **5 questions**, each with **3 sub-questions** (4-5 compounds each)
- Font: **宋体 + Times New Roman** (五号 10.5pt), letter labels use **Arial**
- Numbering: **"1. "** for questions, **"1-1 "** for sub-questions (bold)
- Chemical structures in **ACS 1996 standard** (PNG embedded in DOCX)
- CDXML files exported alongside for ChemDraw editing
- Answer key at end of document

## Question Design Principles (CChO Style)

1. **取代基效应** (Substituent effects): same parent structure, different substituents
   (e.g. X-C₆H₄COOH, X-CH₂COOH, X-C₆H₄OH)
2. **官能团酸性差异** (Functional group comparisons): across carboxylic acids, phenols,
   alcohols, C-H acids, amides, heterocycles
3. **著名例外** (Famous exceptions): Meldrum's acid, dimedone, acetylacetone,
   nitromethane, fluorene, cyclopentadiene, etc.
4. **无机物** (Inorganic acids): freely inserted into any sub-question as distractors
   — students are expected to memorize their pKa values
5. **同小问内可跨主题但不交叉**: 允许混合取代基效应与官能团比较，但排序后
   同类化合物必须相邻（不可出现 A(取代基) > B(官能团) > C(取代基) 的穿插）
6. **由易到难**: 3 sub-questions per question progress from easy (level 1) to
   medium (level 2) to hard (level 3)

## Usage

```
python scripts/generate.py [--count N] [--output path.docx]
```

Options:
- `--count N`: Number of main questions (default: 5)
- `--output path.docx`: Output DOCX path (default: ~/Downloads/CChO酸性排序题.docx)
- `--library path.md`: Path to acidic data library (default: ~/Downloads/酸性数据.md)

## Dependencies

- Python 3.9+
- rdkit
- python-docx
- Pillow (PIL)
- cairosvg (for SVG rendering)
