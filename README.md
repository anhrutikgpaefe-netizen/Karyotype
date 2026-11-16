# 多物种染色体保守连锁与重排比较

根据基因在染色体上的归属和顺序，比较多个物种的染色体对应关系，寻找保守连锁和重排（融合/裂变）的候选。
当前使用的脚本只用 Python 标准库。

## 物种

物种登记在 [web/species.tsv](web/species.tsv)。

- **哺乳类（12 种）**：须鲸 `Balae`、河马 `Hippo`、大象 `elephant`、儒艮 `Dugong`、海牛 `inunguis`、猫 `cat`、狐獴 `meerkat`、
  马 `Equus`、黑猩猩 `hxx`、人 `homo`、穿山甲 `Manis`、假吸血蝠 `Megaderma`
- **鱼类及外群（8 种，`fish/`）**：矛尾鱼 `Latimeria`、棱皮龟 `coriacea`、雀鳝 `gar`，以及 `banmasha`、`jingsha`、`zhusha`、
  `xiaoyao`、`bsnake`（后五个的物种名待确认）

## 输入文件

| 文件 | 内容 |
|---|---|
| `<物种>.txt` | 两列：基因名、染色体 |
| `<物种>_positions.tsv` | Chromosome、Symbol、Start、End、Strand、GeneID、Gene_order |

两者都由 [fish/parse_ncbi_v2.py](fish/parse_ncbi_v2.py) 从 NCBI 下载的基因表（`*_dataset.tsv`）解析得到。
原始基因表体积大，没有放进仓库，见 `.gitignore`。

## 分析

### 1. 染色体对应关系表（`web/`）——当前主要结果

两个物种间，每对染色体共有多少个基因；三个物种时，可以查看两个物种对应到参照物种同一条染色体的那些染色体彼此是否同源。

```
cd web
python chrom_pairs.py species.tsv pairs.json                      # 计算所有物种的基因级数据
python build_chrom_page.py pairs.json chromosome_correspondence.html   # 生成网页
python chrom_pairs.py species.tsv --pair coriacea gar --min 20 --out coriacea_gar.tsv   # 单个物种对的 Excel 式表格
```

加物种：在 `species.tsv` 里加一行，再运行前两条命令。

### 2. 基因顺序比较（`compare_orderN.py`，结果在 `order5/`）

在对应的染色体片段内检验各物种的基因顺序是否共线，并结合物种树判断哪一侧是衍生排列。唯一使用坐标的分析。

```
python compare_orderN.py --species A=A_positions.tsv --species B=B_positions.tsv ... --tree "(elephant,(hxx,(meerkat,Equus,Balae)))" --out order5/mammal5
```

配套：`test_shared_breaks.py`（共享断点检验）、`report_multigroup.py`、`report_shared_events.py`。

### 3. 连锁单元（`units/`）

按基因在所有物种中的染色体组合划分互不重叠的单元，没有锚点物种。

```
python units/units.py --species NAME=FILE ... --min-genes 30 --min-present 5 --out mammal8
```

配套：`sensitivity.py`（阈值敏感性）、`unit_events.py`、`shared_pair.py`、`map_v4.py`、`plot_units.py`；带 `_e` 的是加入马之后的版本。

### 4. 连锁模块与事件候选（`compare8_v5.py`）——探索性

以各物种轮流为锚点找保守模块，再列出"两个模块在部分物种同染色体、在其他物种分开"的候选。
候选之间互相重叠，不等于独立的重排事件；结论以第 2、3 项为准。
用法见 [compare8_使用说明.txt](compare8_使用说明.txt)。

当前版本是 v5（输出 `mammal8_v5_*`）。根目录还留着第一版（`mammal8_*`）和 v4（`mammal8_v4_*`）的输出，
因为 `units/unit_events.py`、`units/map_v4.py` 和 `plot_submodule_events.py` 要读它们。

## archive/

已被取代、不再有脚本读取的旧版本，只为保留记录：

- `archive/compare_old/`：`compare4.py`、`compare8.py`、`compare8_v2.py`–`compare8_v4.py`，以及 4 物种版、v2、v3 的输出。
- `archive/web_xlsx/`：网页第一版用的 Excel 解析脚本 `parse_chrom_xlsx.py`（需要 `openpyxl`）和它的输出。

### 5. 鱼类及外群（`fish/`）

上述方法在 8 个鱼类及外群物种上的版本：`compare8_v2.py`、`compare3_v2.py`、`compare_order3.py`、
`report_one_differs.py`、`polarize_one_differs.py`。

## 局限

- 基因的对应只看基因名是否相同，不是经过验证的直系同源；物种越远，漏掉和错配的基因越多。
- 只使用在一个物种里只出现在一条染色体上的基因名。
- 不同物种的染色体编号互相独立，数字相同不代表同源。
- 模块和事件候选不给出融合或裂变的方向；只有基因顺序比较会借助物种树推断衍生的一侧。
