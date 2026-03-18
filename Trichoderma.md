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

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --fmt

# .lst and .count.tsv
bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    rgr md stdin --num

```

| item    | count |
|---------|------:|
| strain  |   172 |
| species |    44 |
| genus   |     7 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus            | #species | #strains |
|------------------|---------:|---------:|
| Cladobotryum     |        2 |        3 |
| Escovopsis       |        2 |        7 |
| Hypomyces        |        2 |        2 |
| Mycogone         |        1 |        1 |
| Saccharomyces    |        1 |        1 |
| Sphaerostilbella |        1 |        1 |
| Trichoderma      |       35 |      157 |

### Download and check

* When `rsync.sh` is interrupted, run `check.sh` before restarting
* For projects that have finished downloading, but have renamed strains, you can run `reorder.sh` to
  avoid re-downloading
    * `misplaced.tsv`
    * `remove.list`
* The parameters of `n50.sh` should be determined by the distribution of the description statistics
* `collect.sh` generates a file of type `.tsv`, which is intended to be opened by spreadsheet
  software.
    * Information of assemblies are collected from *_assembly_report.txt *after* downloading
    * **Note**: `*_assembly_report.txt` have `CRLF` at the end of the line.
* `finish.sh` generates the following files
    * `omit.lst` - no annotations
    * `collect.pass.tsv` - passes the n50 check
    * `pass.lst` - passes the n50 check
    * `rep.lst` - representative or reference strains
    * `counts.tsv`

```shell
cd ~/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --ass

# Run
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
# LEN_N50   N_CONTIG    LEN_SUM
bash ASSEMBLY/n50.sh 100000 1000 1000000

# Adjust parameters passed to `n50.sh`
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50" --max "C" --min "S"
#N50_min C_max   S_min
#579860  533     33215161

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
#N50_pct10       103504.6
#N50_pct50       1412965.5
#C_pct50 157
#C_pct90 1044.8
#S_pct10 32206756.6
#S_pct50 37298829.5

# After the above steps are completed, run the following commands.

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```

| #item            | fields | lines |
|------------------|-------:|------:|
| url.tsv          |      3 |   172 |
| check.lst        |      1 |   172 |
| collect.tsv      |     20 |   173 |
| n50.tsv          |      4 |   173 |
| n50.pass.tsv     |      4 |   151 |
| collect.pass.tsv |     23 |   151 |
| pass.lst         |      1 |   150 |
| omit.lst         |      1 |   128 |
| rep.lst          |      1 |    48 |
| sp.lst           |      1 |    15 |

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

```shell
cd ~/data/Trichoderma

ulimit -n `ulimit -Hn`

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --bs

bash BioSample/download.sh

# Ignore rare attributes
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#170 lines, 39 fields

cp BioSample/attributes.lst summary/
cp BioSample/biosample.tsv summary/

```

## MinHash

Estimate nucleotide divergences among strains.

* Abnormal strains
    * This [paper](https://doi.org/10.1038/s41467-018-07641-9) showed that >95% intra-species and
      and <83% inter-species ANI values.
    * If the maximum value of ANI between strains within a species is greater than *0.05*, the
      median and maximum value will be reported. Strains that cannot be linked by the median
      ANI, e.g., have no similar strains in the species, will be considered as abnormal strains.
    * It may consist of two scenarios:
        1. Wrong species identification
        2. Poor assembly quality

* Non-redundant strains
    * If the ANI value between two strains within a species is less than *0.005*, the two strains
      are considered to be redundant.
    * Need these files:  representative.lst and omit.lst

* MinHash tree
    * A rough tree is generated by k-mean clustering.

* These abnormal strains should be manually checked to determine whether to include them in the
  subsequent steps.

```shell
cd ~/data/Trichoderma

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005

# Compute assembly sketches
bash MinHash/compute.sh

# Non-redundant strains within species
bash MinHash/nr.sh

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
#  80 summary/NR.lst
#  51 summary/redundant.lst

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
#13

# Distances between all selected sketches, then hierarchical clustering
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/redundant.lst \
    --height 0.4

bash MinHash/dist.sh

```

### Condense branches in the minhash tree

* This phylo-tree is not really formal/correct, and shouldn't be used to interpret phylogenetic
  relationships
* It is just used to find more abnormal strains

```shell
mkdir -p ~/data/Trichoderma/tree
cd ~/data/Trichoderma/tree

nw_reroot ../MinHash/tree.nwk Sa_cer_S288C |
    nwr order stdin --nd --an \
    > minhash.reroot.newick

nwr pl-condense --map -r species \
    minhash.reroot.newick ../MinHash/species.tsv |
    nwr comment stdin -r "(S|member)=" |
    nwr comment stdin -r "^\d+$" |
    nwr order stdin --nd --an \
    > minhash.condensed.newick

mv condensed.tsv minhash.condensed.tsv

nwr tex minhash.condensed.newick --bl -o Trichoderma.minhash.tex

tectonic Trichoderma.minhash.tex

```

## Count valid species and strains

### For *genomic alignments*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh

cat Count/taxa.tsv |
    rgr md stdin --num

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    rgr md stdin --num

# Can accept N_COUNT
bash Count/lineage.sh 1

cat Count/lineage.count.tsv |
    rgr md stdin --num

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv

```

| item    | count |
|---------|------:|
| strain  |   137 |
| species |    38 |
| genus   |     5 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
|---------------|---------:|---------:|
| Cladobotryum  |        2 |        3 |
| Escovopsis    |        2 |        7 |
| Hypomyces     |        2 |        2 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       31 |      124 |

| #family            | genus         | species                     | count |
|--------------------|---------------|-----------------------------|------:|
| Hypocreaceae       | Cladobotryum  | Cladobotryum mycophilum     |     2 |
|                    |               | Cladobotryum protrusum      |     1 |
|                    | Escovopsis    | Escovopsis sp.              |     5 |
|                    |               | Escovopsis weberi           |     2 |
|                    | Hypomyces     | Hypomyces perniciosus       |     1 |
|                    |               | Hypomyces rosellus          |     1 |
|                    | Trichoderma   | Trichoderma afroharzianum   |     5 |
|                    |               | Trichoderma aggressivum     |     1 |
|                    |               | Trichoderma arundinaceum    |     4 |
|                    |               | Trichoderma asperelloides   |     3 |
|                    |               | Trichoderma asperellum      |    18 |
|                    |               | Trichoderma atrobrunneum    |     1 |
|                    |               | Trichoderma atroviride      |    10 |
|                    |               | Trichoderma breve           |     1 |
|                    |               | Trichoderma brevicrassum    |     1 |
|                    |               | Trichoderma citrinoviride   |     4 |
|                    |               | Trichoderma cornu-damae     |     1 |
|                    |               | Trichoderma endophyticum    |     4 |
|                    |               | Trichoderma erinaceum       |     2 |
|                    |               | Trichoderma gamsii          |     2 |
|                    |               | Trichoderma gracile         |     1 |
|                    |               | Trichoderma guizhouense     |     1 |
|                    |               | Trichoderma hamatum         |     1 |
|                    |               | Trichoderma harzianum       |    10 |
|                    |               | Trichoderma koningii        |     1 |
|                    |               | Trichoderma koningiopsis    |     5 |
|                    |               | Trichoderma lentiforme      |     1 |
|                    |               | Trichoderma lixii           |     1 |
|                    |               | Trichoderma longibrachiatum |     7 |
|                    |               | Trichoderma orchidacearum   |     1 |
|                    |               | Trichoderma polysporum      |     1 |
|                    |               | Trichoderma reesei          |    19 |
|                    |               | Trichoderma semiorbis       |     1 |
|                    |               | Trichoderma simmonsii       |     1 |
|                    |               | Trichoderma sp.             |     4 |
|                    |               | Trichoderma virens          |     9 |
|                    |               | Trichoderma viride          |     3 |
| Saccharomycetaceae | Saccharomyces | Saccharomyces cerevisiae    |     1 |

### For *protein families*

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv and taxa.tsv
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
|---------|------:|
| strain  |    35 |
| species |    20 |
| genus   |     4 |
| family  |     2 |
| order   |     2 |
| class   |     2 |

| genus         | #species | #strains |
|---------------|---------:|---------:|
| Cladobotryum  |        1 |        1 |
| Escovopsis    |        1 |        1 |
| Saccharomyces |        1 |        1 |
| Trichoderma   |       17 |       32 |

## Collect proteins

```shell
cd ~/data/Trichoderma/

nwr template ~/Scripts/genomes/assembly/Trichoderma.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# collect proteins
bash Protein/collect.sh

# clustering
# It may need to be run several times
bash Protein/cluster.sh

rm -fr Protein/tmp/

# info.tsv
bash Protein/info.sh

# counts
bash Protein/count.sh

cat Protein/counts.tsv |
    tsv-summarize -H --count --sum 2-7 |
    sed 's/^count/species/' |
    datamash transpose |
    (echo -e "#item\tcount" && cat) |
    rgr md stdin --fmt

```

| #item      |   count |
|------------|--------:|
| species    |      21 |
| strain_sum |      39 |
| total_sum  | 363,402 |
| dedup_sum  | 363,402 |
| rep_sum    | 284,914 |
| fam88_sum  | 258,055 |
| fam38_sum  | 219,984 |

## Phylogenetics with fungi61

```shell
cd ~/data/Trichoderma/

mkdir -p HMM

# The Fungi HMM set
tar xvfz ~/data/HMM/fungi61/fungi61.tar.gz --directory=HMM
cp HMM/fungi61.lst HMM/marker.lst

```

## Phylogenetics with BUSCO

```shell
cd ~/data/Trichoderma/

rm -fr BUSCO

curl -L https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2024-01-08.tar.gz |
    tar xvz
mv fungi_odb10/ BUSCO

#curl -L https://busco-data.ezlab.org/v5/data/lineages/ascomycota_odb10.2024-01-08.tar.gz |
#    tar xvz
#mv ascomycota_odb10/ BUSCO

```

### Find corresponding representative proteins by `hmmsearch`

```shell
cd ~/data/Trichoderma

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

#fd --full-path "Protein/.+/busco.tsv" -X rm

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

fd --full-path "Protein/.+/busco.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-summarize --quantile 2:0.25,0.5,0.75
#23      24      26

# There are 21 species and 39 strains
fd --full-path "Protein/.+/fungi61.tsv" -X cat |
    tsv-summarize --group-by 1 --count |
    tsv-filter --invert --ge 2:20 --le 2:30 |
    cut -f 1 \
    > Protein/marker.omit.lst

cat BUSCO/scores_cutoff |
    parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 1 "
        echo {1}
    " \
    > Protein/marker.lst

wc -l Protein/marker.lst Protein/marker.omit.lst
# 758 Protein/marker.lst
#   0 Protein/marker.omit.lst

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

### Domain related protein sequences

```shell
cd ~/data/Trichoderma

mkdir -p Domain

# each assembly
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

cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

### Align and concat marker genes to create species tree

```shell
cd ~/data/Trichoderma

# Extract proteins
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

# Align each marker
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

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/busco.aln.fas stdin -o Domain/busco.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Domain/busco.aln.fa -out Domain/busco.trim.fa -automated1

hnsm size Domain/busco.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
#762750
#399438

# To make it faster
FastTree -fastest -noml Domain/busco.trim.fa > Domain/busco.trim.newick

```

### The protein tree

```shell
cd ~/data/Trichoderma/tree

nwr reroot  ../Domain/busco.trim.newick -n Sa_cer_S288C |
    nwr order stdin --nd --an \
    > busco.reroot.newick

nwr pl-condense --map -r species \
    busco.reroot.newick ../Count/species.tsv |
    nwr comment stdin -r "(S|member)=" |
    nwr comment stdin -r "^\d+$" |
    nwr order stdin --nd --an \
    > busco.condensed.newick

mv condensed.tsv busco.condense.tsv

nwr tex minhash.condensed.newick --bl -o Trichoderma.busco.tex

tectonic Trichoderma.busco.tex

```
