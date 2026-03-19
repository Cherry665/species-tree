### 列出所有层级

There are no noteworthy classification ranks other than species.

```shell
nwr member Trichoderma |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    rgr md stdin --num

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    rgr md stdin

```

| rank     | count |
| -------- | ----: |
| genus    |     1 |
| species  |   501 |
| no rank  |     1 |
| varietas |     2 |
| strain   |    14 |
| forma    |     2 |

| #rank      | sci_name          | tax_id |
| ---------- | ----------------- | ------ |
| kingdom    | Fungi             | 4751   |
| subkingdom | Dikarya           | 451864 |
| phylum     | Ascomycota        | 4890   |
| subphylum  | Pezizomycotina    | 147538 |
| class      | Sordariomycetes   | 147550 |
| subclass   | Hypocreomycetidae | 222543 |
| order      | Hypocreales       | 5125   |
| family     | Hypocreaceae      | 5129   |
| *genus*    | Trichoderma       | 5543   |

### Species with assemblies

The family Hypocreaceae as outgroups.

```shell
mkdir -p ~/Trichoderma/summary
cd ~/Trichoderma/summary

# 从 NCBI 分类数据库中找到肉座菌科（Hypocreaceae）下属的所有属
nwr member Hypocreaceae -r genus |
    grep -v " x " |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
#19 genus.list

# 针对每一个指定的属（Genus），从 NCBI 的本地数据库中筛选出具有高质量完整基因组的物种
# SELECT species_id, species, COUNT(*)：获取物种 ID、物种名称，并统计该物种下有多少个组装版本
# FROM ar：从 ar（Assembly Reports）表中提取数据
# WHERE genus_id = ${RANK_ID}：约束范围，只看当前循环到的这个“属”下面的物种
# AND species NOT LIKE '% x %'：生物学过滤，剔除名称中带有 " x " 的杂交种，保证分析的是纯种
# AND genome_rep IN ('Full')：质量控制，只拿 Complete Genome 或 Full 级别的结果，剔除碎片化的草图
# GROUP BY species_id：将同一个物种的不同组装版本（可能有多个版本）合并在一起
# HAVING count >= 1：确保该物种至少有一个满足上述条件的基因组
cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %' -- Crossbreeding of two species
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
#   10 RS1.tsv
#  91 GB1.tsv

# 计算每个文件中所有物种拥有的基因组组装总数
for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     10
#GB1     248

```

## Download all assemblies

### Create .assembly.tsv

This step is pretty important

* `nwr kb formats` will give the requirements for `.assembly.tsv`.

* The naming of assemblies has two aspects:
    * for program operation they are unique identifiers;
    * for researchers, they should provide taxonomic information.

If a RefSeq assembly is available, the corresponding GenBank one will not be listed

```shell
cd ~/Trichoderma/summary

# Reference genome
# 导出酿酒酵母的参考基因组
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND species IN ('Saccharomyces cerevisiae')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv

# RS
# 将所有 RS1.tsv 文件中的 TaxID 写在一行，并用逗号隔开，赋值给 SPECIES
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

# 根据筛选出的 RS1.tsv，去数据库里把这些物种所有符合高质量标准的“组装版本（Assemblies）”信息（物种名 + 株系名 + 组装号等）全部提取出来，并存入 raw.tsv
echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# 对那些木霉属中“未鉴定到种”的样本（即 sp. 样本）进行补充提取
echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Preference for refseq
# 提取已经选定的 RefSeq 基因组编号（组装号），防止在下一步加入 GenBank 数据时出现重复
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

# 重复上面操作
# gbrs_paired_asm 是 NCBI 内部的一个“配对指针”，如果一个 GenBank 组装（GCA）有一个对应的 RefSeq 版本（GCF），这个字段就会记录那个 GCF 的编号，如果它没有对应的 RefSeq 版本，这个字段通常为空或指向自己
# 利用 -e 排除模式，排除 gbrs_paired_asm 中和 rs.acc.tsv 中重复的内容
echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

echo "
    SELECT
        genus || ' sp. ' || infraspecific_name || ' ' || assembly_accession AS name,
        genus || ' sp.', genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species LIKE '% sp.%'
        AND species NOT LIKE '% x %'
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

# 深度去重并验证格式
cat raw.tsv |
    rgr dedup stdin |
    datamash check
#250 lines, 7 fields

# Create abbr.
# 从 raw.tsv 原始数据文件中，过滤、去重、标准化处理木霉属相关的基因组组装数据 → 生成缩写名 → 去重（保证 FTP 路径和缩写名唯一）→ 过滤有效 FTP/HTTP 链接 → 按物种 + 缩写名排序 → 最终输出结构化的 Trichoderma.assembly.tsv 文件
# kb abbr：nwr 的子命令，生成一个Perl 脚本（abbr.pl），用于对物种名等文本生成标准化缩写
# （）内用来重新添加表头
nwr kb abbr > abbr.pl
cat raw.tsv |
    grep -v '^#' |
    rgr dedup stdin |
    tsv-select -f 1-6 |
    perl abbr.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Trichoderma.assembly.tsv

datamash check < Trichoderma.assembly.tsv
#248 lines, 5 fields

# find potential duplicate strains or assemblies
# 找出第 1 列（缩写名）重复出现的所有行
cat Trichoderma.assembly.tsv |
    tsv-uniq -f 1 --repeated

# 找出第 2 列（ftp_path）不包含字符串「ftp」 的所有行
cat Trichoderma.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp

# Edit .assembly.tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Trichoderma.assembly.tsv
#
# Save the file to another directory to prevent accidentally changing it
# cp Trichoderma.assembly.tsv ~/Scripts/genomes/assembly

# Cleaning
rm raw*.*sv

```

### Count before download

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class

```shell
cd ~/Trichoderma

nwr template ./summary/Trichoderma.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
# 生成上面两个文件，分别是将从物种名追溯到它的属、科、目、纲和统计有多少个菌株数 (Strain)、种 (Species)、属 (Genus)、科 (Family)、目 (Order)、纲 (Class)（下表1）
bash Count/strains.sh

# 转换成 Markdown 表格格式
cat Count/taxa.tsv |
    rgr md stdin --fmt

# .lst and .count.tsv
# 生成上面两个文件，分别统计有多少属和每个属有多少种和菌株数（下表2）
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    rgr md stdin --num

```

| item    | count |
| ------- | ----: |
| strain  |   247 |
| species |    65 |
| genus   |     7 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus            | #species | #strains |
| ---------------- | -------: | -------: |
| Cladobotryum     |        3 |        4 |
| Escovopsis       |        2 |        7 |
| Hypomyces        |        4 |        4 |
| Mycogone         |        1 |        1 |
| Saccharomyces    |        1 |        1 |
| Sphaerostilbella |        1 |        1 |
| Trichoderma      |       53 |      229 |

### Download and check

* 当`rsync.sh`被中断时，在重新启动下载之前，请先运行`check.sh`（用于检查哪些文件已损坏或下载不完整）  

* 针对已完成下载但变更了品系名称（Strains）的项目：你可以运行`reorder.sh`来重新排序文件，以避免重复下载  
    相关的错误放置信息记录在`misplaced.tsv`中
    需要删除的文件列表记录在`remove.list`中
* `n50.sh`的参数设定：应根据描述性统计数据的分布情况来决定
* `collect.sh`的功能：它会生成一个`.tsv`格式的文件，旨在用电子表格软件（如 Excel）打开
    * 数据来源：组装信息是在下载完成后，从每个样本对应的`*_assembly_report.txt`文件中提取的
    * **Note**: `*_assembly_report.txt`文件的行尾使用的是 CRLF（Windows 换行符格式）
* `finish.sh`生成以下文件：
    * `omit.lst` - 遗漏清单 —— 记录了那些因为没有注释信息而被剔除的物种
    * `collect.pass.tsv` - 质控汇总表 —— 记录了通过了 N50 质量检测的详细信息
    * `pass.lst` - 通过名单 —— 记录了通过 N50 检测的物种名
    * `rep.lst` - 代表性名单 —— 记录了典型株系（Representative strains）或参考株系（Reference strains）
    * `counts.tsv`

```shell
cd ~/Trichoderma

nwr template ./summary/Trichoderma.assembly.tsv \
    --ass

# Run
# 下载基因组文件（基因组 DNA 序列（必选）、蛋白质序列（建树必选）、基因注释信息（必选）、组装报告（质控必选））
# 下载方式需修改
bash ASSEMBLY/rsync.sh

# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

## Put the misplaced directory into the right place
#bash ASSEMBLY/reorder.sh
#
## This operation will delete some files in the directory, so please be careful
#cat ASSEMBLY/remove.lst |
#    parallel --no-run-if-empty --linebuffer -k -j 1 '
#        if [[ -e "ASSEMBLY/{}" ]]; then
#            echo Remove {}
#            rm -fr "ASSEMBLY/{}"
#        fi
#    '

# N50 C S; create n50.tsv and n50.pass.tsv
# 只关注完整基因组 fasta 文件
# LEN_N50：N50值   N_CONTIG:contig数量(C)    LEN_SUM：基因组长度（S）
# 3个值有默认值100000 1000 1000000，也可以手动输入
bash ASSEMBLY/n50.sh 100000 1000 1000000

# Adjust parameters passed to `n50.sh`
#统计数据中最差指标：N50 的最小值、Contig 数量的最大值、总长度的最小值
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#579860  533     31700302

# 统计三个值的中位数、10%、90%阈值
cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#N50_pct10       139491
#N50_pct50       1289709
#C_pct50 147
#C_pct90 883.8
#S_pct10 32267925.4
#S_pct50 37251948

# After the above steps are completed, run the following commands.

# Collect收集信息; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```

| #item            | fields | lines |
| ---------------- | -----: | ----: |
| url.tsv          |      3 |   247 |
| check.lst        |      1 |   247 |
| collect.tsv      |     20 |   248 |
| n50.tsv          |      4 |   248 |
| n50.pass.tsv     |      4 |   223 |
| collect.pass.tsv |     23 |   223 |
| pass.lst         |      1 |   222 |
| omit.lst         |      1 |   176 |
| rep.lst          |      1 |    51 |
| sp.lst           |      1 |    29 |

### Rsync to hpcc

```shell
rsync -avP \
    ~/data/Trichoderma/ \
    wangq@202.119.37.251:data/Trichoderma

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/Trichoderma/ \
    wangq@58.213.64.36:data/Trichoderma

# rsync -avP wangq@202.119.37.251:data/Trichoderma/ ~/data/Trichoderma

# rsync -avP -e 'ssh -p 8804' wangq@58.213.64.36:data/Trichoderma/ ~/data/Trichoderma

```

## BioSample

ENA's BioSample missed many strains, so NCBI's was used.  
统计一些样本的信息

```shell
cd ~/Trichoderma

# 查一下系统最高允许开多少个文件，将当前权限提升到那个最高值
ulimit -n `ulimit -Hn`

nwr template ./summary/Trichoderma.assembly.tsv \
    --bs

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#245 lines, 42 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

菌株间核苷酸分化程度评估

* 异常菌株
    * [文献](https://doi.org/10.1038/s41467-018-07641-9) 指出：同物种内平均核苷酸一致性（ANI）＞95%，不同物种间ANI＜83%
    * 若某物种内部菌株间的最大 ANI 差异＞0.05，则需报告其中位数与最大值。无法通过中位数 ANI 关联的菌株（如在该物种中无相似菌株），将被判定为异常菌株
    * 异常通常分为两种情况：
        1. 物种鉴定错误
        2. 基因组组装质量较差

* 去冗余菌株
    * 若同一物种内两株菌株的 ANI 差异＜0.005，则视为冗余菌株
    * 需要文件：representative.lst and omit.lst

* MinHash 进化树
    * 通过 k‑均值聚类生成初步进化树

* T这些异常菌株需人工核查，以确定是否纳入后续分析步骤

```shell
cd ~/Trichoderma

nwr template ./summary/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
# 先筛选出需要分析的菌株，再并行检查并生成每个菌株的 MinHash 草图（.msh文件）
bash MinHash/compute.sh

# Non-redundant strains within species
# 生成每个物种的非冗余组装体 ID 列表 NR.lst，冗余组装体 ID 列表 redundant.lst
# 将 cluster 一行改为 hnsm clust cc stdin
bash MinHash/nr.sh

# 组合所有的 NR.lst 和 redundant.lst，合并后去重并排序
find MinHash -name "NR.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/NR.lst
find MinHash -name "redundant.lst" |
    xargs cat |
    sort |
    uniq \
    > summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
#  117 summary/NR.lst
#  68 summary/redundant.lst

# Abnormal strains
# 识别基因组相似度异常的菌株：筛选出菌株间最大距离超过阈值的物种，并定位其中的异常菌株
# 将 cluster 一行改为 hnsm clust cc stdin
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#22

# 计算草图间的距离，继而进行层次聚类分析
cd ~/Trichoderma/

nwr template ./summary/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/redundant.lst \
    --height 0.4

# 筛选非冗余的基因组 Mash 索引文件，计算所有基因组间的 Mash 距离矩阵,用 R 对基因组聚类，并按距离阈值划分聚类组，输出进化树（tree.nwk）和聚类结果（groups.tsv）
bash MinHash/dist.sh

```

### 合并minhash tree的分支

* 该系统发育树并非严格意义上规范 / 正确的版本，不应被用于解读系统发育关系
* 其仅用于筛选更多的异常菌株

```shell
mkdir -p ~/Trichoderma/tree
cd ~/Trichoderma/tree

# 将树的根设置在 Sa_cer_S288C 上
# 对进化树节点进行排序，--nd 按每个节点的 “后代数量” 升序排列（后代数量越少的分支越靠前），--an 按节点标签的字母数字顺序升序排列
nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |
    nwr ops order stdin --nd --an \
    > minhash.reroot.newick

# 将物种名映射到树上，并按物种层级合并树分支，清理树的注释信息
nwr pl-condense --map -r species \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr viz comment stdin -r "(S|member)=" |
    nwr viz comment stdin -r "^\d+$" |
    nwr ops order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

# 转换为 LaTeX 格式的可视化代码文件
nwr viz tex minhash.condensed.newick --bl -o Trichoderma.minhash.tex

# 编译LaTeX文件生成PDF
tectonic Trichoderma.minhash.tex

```

## 统计有效的物种和菌株数量

### For *genomic alignments*

```shell
cd ~/Trichoderma/

nwr template ./summary/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
# 筛选出 “质量合格 + 非异常” 的菌株,生成统计数据（下表1）
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
# 提取所有有效属名并生成列表，再对每个属批量统计其包含的唯一物种数和唯一菌株数（下表2）
bash Count/rank.sh

cat Count/genus.count.tsv |
    rgr md stdin --num

# Can accept N_COUNT
# 按「科→属→物种」层级统计菌株数量，筛选出数量达标（≥输入数量）的物种
bash Count/lineage.sh 1

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
| ------- | ----: |
| strain  |   200 |
| species |    58 |
| genus   |     5 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
| ------------- | -------: | -------: |
| Cladobotryum  |        3 |        4 |
| Escovopsis    |        2 |        7 |
| Hypomyces     |        4 |        4 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       48 |      184 |

| #family            | genus         | species                     | count |
| ------------------ | ------------- | --------------------------- | ----: |
| Hypocreaceae       | Cladobotryum  | Cladobotryum mycophilum     |     2 |
|                    |               | Cladobotryum protrusum      |     1 |
|                    |               | Cladobotryum sp.            |     1 |
|                    | Escovopsis    | Escovopsis sp.              |     5 |
|                    |               | Escovopsis weberi           |     2 |
|                    | Hypomyces     | Hypomyces aurantius         |     1 |
|                    |               | Hypomyces perniciosus       |     1 |
|                    |               | Hypomyces rosellus          |     1 |
|                    |               | Hypomyces semicircularis    |     1 |
|                    | Trichoderma   | Trichoderma aethiopicum     |     1 |
|                    |               | Trichoderma afarasin        |     1 |
|                    |               | Trichoderma afroharzianum   |     5 |
|                    |               | Trichoderma aggressivum     |     1 |
|                    |               | Trichoderma arundinaceum    |     4 |
|                    |               | Trichoderma asperelloides   |     4 |
|                    |               | Trichoderma asperellum      |    22 |
|                    |               | Trichoderma atrobrunneum    |     1 |
|                    |               | Trichoderma atroviride      |    12 |
|                    |               | Trichoderma austrokoningii  |     1 |
|                    |               | Trichoderma barbatum        |     1 |
|                    |               | Trichoderma breve           |     1 |
|                    |               | Trichoderma brevicrassum    |     1 |
|                    |               | Trichoderma camerunense     |     1 |
|                    |               | Trichoderma caribbaeum      |     1 |
|                    |               | Trichoderma ceciliae        |     1 |
|                    |               | Trichoderma cf. simile WF8  |     1 |
|                    |               | Trichoderma chlorosporum    |     1 |
|                    |               | Trichoderma citrinoviride   |     6 |
|                    |               | Trichoderma compactum       |     1 |
|                    |               | Trichoderma cornu-damae     |     1 |
|                    |               | Trichoderma endophyticum    |     4 |
|                    |               | Trichoderma erinaceum       |     2 |
|                    |               | Trichoderma evansii         |     1 |
|                    |               | Trichoderma gamsii          |     4 |
|                    |               | Trichoderma ghanense        |     1 |
|                    |               | Trichoderma gracile         |     2 |
|                    |               | Trichoderma guizhouense     |     1 |
|                    |               | Trichoderma hamatum         |     5 |
|                    |               | Trichoderma harzianum       |    16 |
|                    |               | Trichoderma koningii        |     1 |
|                    |               | Trichoderma koningiopsis    |     8 |
|                    |               | Trichoderma lentiforme      |     1 |
|                    |               | Trichoderma lixii           |     1 |
|                    |               | Trichoderma longibrachiatum |    11 |
|                    |               | Trichoderma novae-zelandiae |     1 |
|                    |               | Trichoderma orchidacearum   |     1 |
|                    |               | Trichoderma pleuroticola    |     1 |
|                    |               | Trichoderma polysporum      |     1 |
|                    |               | Trichoderma reesei          |    25 |
|                    |               | Trichoderma semiorbis       |     1 |
|                    |               | Trichoderma simmonsii       |     1 |
|                    |               | Trichoderma sp.             |    11 |
|                    |               | Trichoderma velutinum       |     1 |
|                    |               | Trichoderma virens          |     9 |
|                    |               | Trichoderma viride          |     4 |
|                    |               | Trichoderma virilente       |     1 |
|                    |               | Trichoderma yunnanense      |     1 |
| Saccharomycetaceae | Saccharomyces | Saccharomyces cerevisiae    |     1 |

### For *protein families*

```shell
cd ~/Trichoderma/

nwr template ./summary/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv（下表1）
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```

| item    | count |
| ------- | ----: |
| strain  |    59 |
| species |    36 |
| genus   |     4 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
| ------------- | -------: | -------: |
| Cladobotryum  |        1 |        1 |
| Escovopsis    |        1 |        1 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       33 |       56 |

## Collect proteins

```shell
cd ~/Trichoderma/

nwr template ./summary/Trichoderma.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# 为每个物种构建标准化的蛋白序列资源库：先根据输入参数筛选目标物种的菌株列表，再对每个物种批量提取其所有菌株的蛋白序列（*_protein.faa.gz），去重后生成该物种的非冗余蛋白库，同时整理蛋白注释和组装体关联信息并压缩保存
bash Protein/collect.sh

# clustering
# It may need to be run several times
# 先把物种的蛋白去重（菌株级去冗余→物种级代表序列）（95% 相似rep_seq.fa.gz），再按 80% 相似分功能家族（物种级代表序列→蛋白家族）（fam88_cluster.tsv），按 30% 相似分进化家族（远缘蛋白家族聚类）（fam38_cluster.tsv），最终得到不同层级的蛋白分类结果
bash Protein/cluster.sh

rm -fr Protein/tmp/

# info.tsv
# 先筛选目标物种的有效菌株列表，再将每个物种的蛋白序列、注释、聚类结果等多维度数据整合到 SQLite 数据库（seq.sqlite）中
# pro.fa.gz（非冗余蛋白库）、rep_cluster.tsv（95% 聚类结果）、anno.tsv.gz（蛋白注释）、asmseq.tsv.gz（蛋白 - 菌株关联）、fam88_cluster.tsv（80% 聚类）、fam38_cluster.tsv（30% 聚类）
bash Protein/info.sh

# counts
# 先筛选目标物种的有效菌株列表，再对每个物种的 seq.sqlite 数据库中提取统计指标（下表中的指标）
bash Protein/count.sh

cat Protein/counts.tsv |
    tsv-summarize -H --count --sum 2-7 |
    sed 's/^count/species/' |
    datamash transpose |
    (echo -e "#item\tcount" && cat) |
    rgr md stdin --fmt

```

| #item      |   count |
| ---------- | ------: |
| species    |      36 |
| strain_sum |      67 |
| total_sum  | 687,394 |
| dedup_sum  | 687,394 |
| rep_sum    | 529,592 |
| fam88_sum  | 466,019 |
| fam38_sum  | 390,874 |


## 用 BUSCO 进行系统发育分析

```shell
cd ~/Trichoderma/

rm -fr BUSCO

curl -L https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz |
    tar xvz
mv fungi_odb10/ BUSCO

#curl -L https://busco-data.ezlab.org/v5/data/lineages/ascomycota_odb10.2024-01-08.tar.gz |
#    tar xvz
#mv ascomycota_odb10/ BUSCO

```

### 通过`hmmsearch`筛选对应的代表性蛋白

```shell
cd ~/Trichoderma

# 只取 pass.lst 中的物种，排除 omit.lst 中没有注释的物种
cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

# 在每个物种的蛋白质组中，找到与 BUSCO（真菌通用单拷贝直系同源基因库）匹配的序列，格式化输出 BUSCO marker - 蛋白 ID 映射表
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ -s Protein/"${SPECIES}"/busco.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/rep_seq.fa.gz ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    cat BUSCO/scores_cutoff |
        parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 4 "
            gzip -dcf Protein/${SPECIES}/rep_seq.fa.gz |
                hmmsearch -T {2} --domT {2} --noali --notextw BUSCO/hmms/{1}.hmm - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({1}), \$1; '
        " \
        > Protein/${SPECIES}/busco.tsv
done

# 统计每个 BUSCO marker 出现的次数，计算下四分位数、中位数和上四分位数
fd --full-path "Protein/.+/busco.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-summarize --quantile 2:0.25,0.5,0.75
#40      42      45

# There are 36 species and 67 strains
# 需根据实际修改，这里是保留出现次数在 20 到 30 次之间的 Marker，丢弃少于20次和多于30次的（在 marker.omit.lst 中）
fd --full-path "Protein/.+/busco.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-filter --invert --ge 2:40 --le 2:75 |
    cut -f 1 \
    > Protein/marker.omit.lst

# 提取所有 BUSCO marker 列表
cat BUSCO/scores_cutoff |
    parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 1 "
        echo {1}
    " \
    > Protein/marker.lst

wc -l Protein/marker.lst Protein/marker.omit.lst
# 758 Protein/marker.lst
#   186 Protein/marker.omit.lst

# 去掉需要去除的基因，生成单拷贝基因列表
# 在本地 SQLite 数据库中，将 BUSCO Marker 的身份（ID）与实际的蛋白质序列建立“强关联”索引
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -s Protein/"${SPECIES}"/busco.tsv ]]; then
        continue
    fi
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    # single copy
    cat Protein/"${SPECIES}"/busco.tsv |
        grep -v -Fw -f Protein/marker.omit.lst \
        > Protein/"${SPECIES}"/busco.sc.tsv

    nwr seqdb -d Protein/${SPECIES} --rep f3=Protein/${SPECIES}/busco.sc.tsv

done

```

### 结构域相关的蛋白质序列

```shell
cd ~/Trichoderma

mkdir -p Domain

# each assembly
# 从 36 个物种的本地数据库中，通过 SQL 查询精准地把那些“单拷贝 BUSCO 基因”对应的蛋白质序列提取出来
cat Protein/species-f.tsv |
    tsv-select -f 2 |
    rgr dedup stdin |
while read SPECIES; do
    if [[ ! -f Protein/"${SPECIES}"/seq.sqlite ]]; then
        continue
    fi

    echo >&2 "${SPECIES}"

    echo "
        SELECT
            seq.name,
            asm.name,
            rep.f3
        FROM asm_seq
        JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
        JOIN seq ON asm_seq.seq_id = seq.id
        JOIN rep ON rep_seq.rep_id = rep.id
        JOIN asm ON asm_seq.asm_id = asm.id
        WHERE 1=1
            AND rep.f3 IS NOT NULL
        ORDER BY
            asm.name,
            rep.f3
        " |
        sqlite3 -tabs Protein/${SPECIES}/seq.sqlite \
        > Protein/${SPECIES}/seq_asm_f3.tsv

    hnsm some Protein/"${SPECIES}"/pro.fa.gz <(
            tsv-select -f 1 Protein/"${SPECIES}"/seq_asm_f3.tsv |
                rgr dedup stdin
        )
done |
    hnsm dedup stdin |
    hnsm gz stdin -o Domain/busco.fa

fd --full-path "Protein/.+/seq_asm_f3.tsv" -X cat \
    > Domain/seq_asm_f3.tsv

# 去冗余
cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### “比对并串联标记基因以构建物种树

```shell
cd ~/Trichoderma

# Extract proteins
# 对每一个 BUSCO Marker，提取对应的蛋白序列，放在一个文件夹中
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

        mkdir -p Domain/{}

        hnsm some Domain/busco.fa.gz <(
            cat Domain/seq_asm_f3.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) \
            > Domain/{}/{}.pro.fa
    '

# 利用 mafft 对每一个 BUSCO Marker 进行序列比对，对蛋白质序列进行“对齐”
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Domain/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Domain/{}/{}.aln.fa ]; then
            exit
        fi
#        muscle -quiet -in Domain/{}/{}.pro.fa -out Domain/{}/{}.aln.fa
        mafft --auto Domain/{}/{}.pro.fa > Domain/{}/{}.aln.fa
    '

# 将比对文件中的蛋白名改为菌株名
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
while read marker; do
    echo >&2 "==> marker [${marker}]"
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # Only NR strains
    # 1 name to many names
    cat Domain/seq_asm_f3.NR.tsv |
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

# Concat marker genes
# 将几百个 Marker 的比对结果全部合并进一个文件（Domain/busco.aln.fas）
cat Protein/marker.lst |
    grep -v -Fw -f Protein/marker.omit.lst |
while read marker; do
    if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    cat Domain/${marker}/${marker}.replace.fa

    # empty line for .fas
    echo
done \
    > Domain/busco.aln.fas

# 将基因片段的集合按菌株名串联起来
cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/busco.aln.fas stdin -o Domain/busco.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/busco.aln.fa -out Domain/busco.trim.fa -automated1

# 统计长度（上：原始串联长度，下：修剪后长度）
hnsm size Domain/busco.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#556558
#292922

# To make it faster
# 正式建树去掉 -fastest -noml，建更正式的 ML 树
FastTree -fastest -noml Domain/busco.trim.fa > Domain/busco.trim.newick

```

### The protein tree

```shell
cd ~/Trichoderma/tree

# （与 MinHash 中类似）
nw_reroot  ../Domain/busco.trim.newick Sa_cer_S288C |
    nwr ops order stdin --nd --an \
    > busco.reroot.newick

nwr pl-condense --map -r species \
    busco.reroot.newick ../Count/species.tsv |
    nwr viz comment stdin -r "(S|member)=" |
    nwr viz comment stdin -r "^\d+$" |
    nwr ops order stdin --nd --an \
    > busco.condensed.newick

mv condensed.tsv busco.condense.tsv

nwr viz tex minhash.condensed.newick --bl -o Trichoderma.busco.tex

tectonic Trichoderma.busco.tex

```
