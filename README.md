# CChO 酸性排序题生成器

生成**中国化学奥林匹克（CChO）**风格的酸性排序题目，附带 **ACS 1996 标准**化学结构式，输出为 DOCX 文件。

## 功能

- 自动生成 5 道大题，每道 3 小题（每题含 4-5 个化合物）
- 题型覆盖：取代基效应、官能团酸性差异、著名例外、无机酸穿插
- 化学结构式以 ACS 1996 标准渲染，嵌入 DOCX
- 同步导出 CDXML 文件，方便 ChemDraw 编辑
- 附带答案与解析

## 材料库

本 Skill 合并两个数据源：

1. **酸性数据.md**（约 200 个手工整理化合物）：来自 Bordwell/Evans pKa 表，包含 SMILES、名称、H₂O 和 DMSO pKa 值
2. **IUPAC Dissociation-Constants**（筛选后约 7,300 个化合物）：来自 IUPAC SC-Database 的高置信度水溶液 pKa 数据

> 数据文件已加入 `.gitignore`，需自行准备后放入仓库根目录。

## 安装

### 前置依赖

- Python 3.9+
- rdkit
- python-docx
- Pillow (PIL)
- cairosvg

```bash
pip install rdkit python-docx Pillow cairosvg
```

### 数据文件准备

将以下文件放入本 Skill 根目录：

| 文件 | 说明 |
|------|------|
| `酸性数据.md` | 手工整理的酸性数据 |
| `iupac_high-confidence_v2_3.csv` | IUPAC 高置信度 pKa 数据 |

## 使用方法

```
python scripts/generate.py [--count N] [--output path.docx]
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--count N` | 大题数量 | 5 |
| `--output path.docx` | 输出 DOCX 路径 | `~/Downloads/CChO酸性排序题.docx` |
| `--library path.md` | 酸性数据库路径 | `~/Downloads/酸性数据.md` |

### 示例

```bash
# 生成默认 5 道题
python scripts/generate.py

# 生成 10 道题，指定输出路径
python scripts/generate.py --count 10 --output ./output/酸性排序题.docx
```

## 出题设计原则

1. **取代基效应**：同一母体结构，不同取代基（如 X-C₆H₄COOH、X-CH₂COOH、X-C₆H₄OH）
2. **官能团酸性差异**：跨羧酸、酚、醇、C-H 酸、酰胺、杂环比较
3. **著名例外**：Meldrum's acid、dimedone、乙酰丙酮、硝基甲烷、芴、环戊二烯等
4. **无机物**：可自由插入任意小题中作为干扰项
5. **同小问内可跨主题但不交叉**：排序后同类化合物必须相邻
6. **由易到难**：每题 3 小题按难度递增

## 输出格式

- 字体：宋体 + Times New Roman（五号 10.5pt），选项标签使用 Arial
- 编号：大题用 **"1. "**，小题用 **"1-1 "**（加粗）
- 化学结构：ACS 1996 标准 PNG 嵌入 DOCX
- 文档末尾附带答案

## 目录结构

```
CChOlike-acid-ranking/
├── SKILL.md              # Skill 定义文件
├── README.md             # 本文件
├── extract.py            # 数据提取脚本
├── .gitignore            # 忽略数据文件
└── scripts/
    ├── generate.py       # 主生成脚本
    ├── parse_iupac.py    # IUPAC 数据解析
    ├── check_failures.py # 检查失败化合物
    ├── inspect_compounds.py # 化合物检查工具
    └── export_conjugate_acids.py # 共轭酸导出
```

## 许可

仅供学习与教学使用。
