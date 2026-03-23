## 文件处理

```shell
cd ~/yeast/genomes
for file in *.Final.fasta; do
    mv "$file" "${file%.Final.fasta}.fna"
done
```

```shell
cd ~/yeast/PEP_fasta
for file in *.nuclear_genome.Final.pep.fa; do
    mv "$file" "${file%.nuclear_genome.Final.pep.fa}.faa"
done
sed -i 's/|.*//' *.faa
```

```shell
cd ~/yeast/CDS_fasta
for file in *.nuclear_genome.Final.cds.fa; do
    mv "$file" "${file%.nuclear_genome.Final.cds.fa}_cds_from_genomic.fna"
done
```

```shell
cd ~/yeast
mkdir summary
ls ./genomes/*.fna | sed 's|.*/||' | sed 's/\.fna$//' > name1.list
for name in $(cat ./genomes/name1.list); do
    mkdir -p ASSEMBLY/ASSEMBLY/$name
    cp ./PEP_fasta/${name}.faa ASSEMBLY/ASSEMBLY/$name/
    cp ./genomes/${name}.fna ASSEMBLY/ASSEMBLY/$name/
    cp ./CDS_fasta/${name}_cds_from_genomic.fna ./ASSEMBLY/ASSEMBLY/$name/
done
```

```shell
# 下载外群 Saccharomyces_paradoxus 基因组和蛋白文件
cd ~/yeast
NAME="Saccharomyces_paradoxus"
FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_paradoxus/latest_assembly_versions/GCF_002079055.1_ASM207905v1/GCF_002079055.1_ASM207905v1"

mkdir -p "ASSEMBLY/ASSEMBLY/$NAME"

# 下载基因组 (fna)
wget -O "ASSEMBLY/ASSEMBLY/$NAME/${NAME}.fna.gz" "${FTP_BASE}_genomic.fna.gz"

# 下载蛋白序列 (faa)
wget -O "ASSEMBLY/ASSEMBLY/$NAME/${NAME}.faa.gz" "${FTP_BASE}_protein.faa.gz"

# 下载 cds 文件
wget -O "ASSEMBLY/ASSEMBLY/$NAME/${NAME}_cds_from_genomic.fna.gz" "${FTP_BASE}_cds_from_genomic.fna.gz"

gzip -d "ASSEMBLY/ASSEMBLY/$NAME/"*.gz

cp ./ASSEMBLY/ASSEMBLY/Saccharomyces_paradoxus/Saccharomyces_paradoxus.fna ./genomes/Saccharomyces_paradoxus.fna

cd ~/yeast/ASSEMBLY/ASSEMBLY/Saccharomyces_paradoxus
# 在 > 后面插入 Spar_ 前缀
sed -i 's/^>/>Spar_/' Saccharomyces_paradoxus.faa
# 删掉第一个空格及其后面的所有内容
sed -i '/^>/ s/ .*//' Saccharomyces_paradoxus.faa
```

```shell
# 压缩文件
cd ~/yeast
find ASSEMBLY/ASSEMBLY/ -name "*.fna" | while read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file" .fna)
    mv "$file" "${dir}/${base}_genomic.fna"
    gzip "${dir}/${base}_genomic.fna"
done

find ASSEMBLY/ASSEMBLY/ -name "*.faa" | while read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file" .faa)
    mv "$file" "${dir}/${base}_protein.faa"
    gzip "${dir}/${base}_protein.faa"
done

find ASSEMBLY/ASSEMBLY/ -name "*_cds_from_genomic.fna" | while read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file" _cds_from_genomic.fna)
    gzip "${dir}/${base}_cds_from_genomic.fna"
done
```

```shell
cd ~/yeast/summary
ls ../genomes/*.fna | sed 's|.*/||' | sed 's/\.fna$//' > name.list

# 生成 yeast.assembly.tsv
echo -e "#name\tftp_path\tbiosample\tspecies\tassembly_level" > yeast.assembly.tsv
echo "==> 开始处理酵母文件..."

for file in ../genomes/*.fna; do
    strain=$(basename "$file" .fna)
    echo -e "${strain}\tlocal\t-\tSaccharomyces_cerevisiae\tChromosome" >> yeast.assembly.tsv
done

echo "yeast.assembly.tsv 已生成，共有 $(wc -l yeast.assembly.tsv) 行"
```

## 检查

```shell
cd ~/yeast
nwr template ./summary/yeast.assembly.tsv \
    --ass

rm ./ASSEMBLY/url.tsv
cp ./summary/name.list ./ASSEMBLY
awk '{print $1"\t.\tASSEMBLY"}' ./ASSEMBLY/name.list > ./ASSEMBLY/url.tsv

bash ASSEMBLY/n50.sh 100000 1000 1000000

cat ASSEMBLY/n50.tsv |
    tsv-filter -H |
    tsv-summarize -H --min "N50" --max "C" --min "S"
# N50_min C_max   S_min
# 151786  136     11359891

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "N50:0.1,0.5" --quantile "C:0.5,0.9" --quantile "S:0.1,0.5" |
    datamash transpose
# N50_pct10       823913.9
# N50_pct50       910124.5
# C_pct50 16
# C_pct90 18
# S_pct10 11700068.2
# S_pct50 11846490

bash ASSEMBLY/collect.sh
cat ./ASSEMBLY/name.list >> ./ASSEMBLY/collect.tsv

bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    rgr md stdin --fmt

```
| #item        | fields | lines |
| ------------ | -----: | ----: |
| url.tsv      |      3 |   860 |
| check.lst    |      0 |     0 |
| n50.tsv      |      4 |   861 |
| n50.pass.tsv |      4 |   861 |
| pass.lst     |      1 |   860 |
| omit.lst     |      0 |     0 |
| rep.lst      |      1 |   860 |

## MinHash

```shell
cd ~/yeast
nwr template ./summary/yeast.assembly.tsv \
    --mh \
    --parallel 8 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.001

sed -i 's/Saccharomyces_cerevisiae/ASSEMBLY/g' ./MinHash/species.tsv
bash MinHash/compute.sh

# 将 cluster 一行改为 hnsm clust cc stdin
bash MinHash/nr.sh

cp MinHash/ASSEMBLY/NR.lst summary/NR.lst
cp MinHash/ASSEMBLY/redundant.lst summary/redundant.lst
wc -l summary/NR.lst summary/redundant.lst
#  441 summary/NR.lst
#  419 summary/redundant.lst

# 将 cluster 一行改为 hnsm clust cc stdin
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst | wc -l
# 11（含 Saccharomyces_paradoxus ）

cd ~/yeast
nwr template ./summary/yeast.assembly.tsv \
    --mh \
    --parallel 8 \
    --not-in summary/redundant.lst \
    --height 0.01

# 将用到 species.tsv {2} 改为 ASSEMBLY
bash MinHash/dist.sh
```

## minhash tree

```shell
mkdir -p ~/yeast/tree
cd ~/yeast/tree

# 将树的根设置在 Saccharomyces_paradoxus 上
nw_reroot ../MinHash/tree.nwk Saccharomyces_paradoxus |
    nwr ops order stdin --nd --an \
    > minhash_S_cerevisiae_strain_tree.nwk

```

## 收集蛋白

```shell
cd ~/yeast
sed -i '/Saccharomyces_paradoxus/d' MinHash/abnormal.lst

nwr template ./summary/yeast.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in ASSEMBLY/omit.lst

# 修改 collect.sh
bash Protein/collect.sh

# 修改 cluster.sh
bash Protein/cluster.sh
rm -fr Protein/ASSEMBLY/tmp/

bash Protein/info.sh

```

## 用 BUSCO 进行系统发育分析

```shell
cd ~/yeast
rm -fr BUSCO

curl -L https://busco-data.ezlab.org/v5/data/lineages/saccharomycetes_odb10.2024-01-08.tar.gz | tar xvz
mv saccharomycetes_odb10/ BUSCO

```

## 通过`hmmsearch`筛选对应的代表性蛋白

```shell
cd ~/yeast

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 \
    > Protein/species-f.tsv

if [[ ! -s "Protein/ASSEMBLY/rep_seq.fa.gz" ]]; then
    exit 1
fi

if [[ -s "Protein/ASSEMBLY/busco.tsv" ]]; then
    log_info "busco.tsv already exists, skipping."
else
    cat BUSCO/scores_cutoff |
        parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 4 "
            gzip -dcf Protein/ASSEMBLY/rep_seq.fa.gz |
                hmmsearch -T {2} --domT {2} --noali --notextw BUSCO/hmms/{1}.hmm - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({1}), \$1; '
        " \
        > Protein/ASSEMBLY/busco.tsv
fi

sqlite3 Protein/ASSEMBLY/seq.sqlite "ALTER TABLE rep ADD COLUMN busco TEXT;"
nwr seqdb -d Protein/ASSEMBLY --rep busco=Protein/ASSEMBLY/busco.tsv

sqlite3 -tabs "Protein/ASSEMBLY/seq.sqlite" <<EOF > "Protein/ASSEMBLY/busco_stats.tsv"
    SELECT 
        rep.busco AS busco_id, 
        COUNT(DISTINCT asm_seq.asm_id) AS strain_count
    FROM asm_seq
    JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
    JOIN rep ON rep_seq.rep_id = rep.id
    WHERE rep.busco IS NOT NULL
    GROUP BY rep.busco;
EOF

# 在至少 95% 的菌株中存在且为单拷贝
sqlite3 -tabs Protein/ASSEMBLY/seq.sqlite "
    SELECT rep.busco
    FROM asm_seq
    JOIN rep_seq ON asm_seq.seq_id = rep_seq.seq_id
    JOIN rep ON rep_seq.rep_id = rep.id
    WHERE rep.busco IS NOT NULL
    GROUP BY rep.busco
    HAVING COUNT(asm_seq.seq_id) != COUNT(DISTINCT asm_seq.asm_id)
        OR COUNT(DISTINCT asm_seq.asm_id) < 817;
" > Protein/marker.omit.lst

cat BUSCO/scores_cutoff |
    parallel --colsep '\s+' --no-run-if-empty --linebuffer -k -j 1 "
        echo {1}
    " \
    > Protein/marker.lst

wc -l Protein/marker.lst Protein/marker.omit.lst
# 2137 Protein/marker.lst
# 1567 Protein/marker.omit.lst

# 去掉需要去除的基因，生成单拷贝基因列表
if [[ ! -s "Protein/ASSEMBLY/busco.tsv" ]]; then
    exit 1
fi

if [[ -s "Protein/ASSEMBLY/seq.sqlite" ]]; then
cat Protein/ASSEMBLY/busco.tsv |
    grep -v -Fw -f Protein/marker.omit.lst \
    > Protein/ASSEMBLY/busco.sc.tsv
    nwr seqdb -d Protein/ASSEMBLY --rep f3=Protein/ASSEMBLY/busco.sc.tsv
fi

```

## 结构域相关的蛋白质序列

```shell
cd ~/yeast
mkdir -p Domain

if [[ -f "Protein/ASSEMBLY/seq.sqlite" ]]; then
    # 1. 从数据库中提取满足 f3 条件的序列与菌株映射关系
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
    sqlite3 -tabs "Protein/ASSEMBLY/seq.sqlite" \
    > "Protein/ASSEMBLY/seq_asm_f3.tsv"

    # 2. 提取序列并压缩
    hnsm some "Protein/ASSEMBLY/pro.fa.gz" <(
        tsv-select -f 1 "Protein/ASSEMBLY/seq_asm_f3.tsv" |
            rgr dedup stdin
    ) |
    hnsm dedup stdin |
    hnsm gz stdin -o "Domain/busco.fa"
    
else
    echo "Error: Protein/ASSEMBLY/seq.sqlite not found!" >&2
    exit 1
fi

cp Protein/ASSEMBLY/seq_asm_f3.tsv Domain/seq_asm_f3.tsv

# 去冗余
cat Domain/seq_asm_f3.tsv |
    tsv-join -e -d 2 -f summary/redundant.lst -k 1 \
    > Domain/seq_asm_f3.NR.tsv

```

## 比对并串联标记基因以构建物种树

```shell
cd ~/yeast

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

    if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    cat Domain/seq_asm_f3.NR.tsv |
        tsv-filter --str-eq "3:${marker}" |
        tsv-select -f 1-2 |
        hnsm replace -s Domain/${marker}/${marker}.aln.fa stdin \
        > Domain/${marker}/${marker}.replace.fa
done

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

    echo
done \
    > Domain/busco.aln.fas

cat Domain/seq_asm_f3.NR.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/busco.aln.fas stdin -o Domain/busco.aln.fa

trimal -in Domain/busco.aln.fa -out Domain/busco.trim.fa -automated1

# 统计长度（上：原始串联长度，下：修剪后长度）
hnsm size Domain/busco.*.fa |
    rgr dedup stdin -f 2 |
    cut -f 2
# 343703
# 269218

```

## 建树

```shell
# 测试
# FastTree -fastest -noml Domain/busco.trim.fa > Domain/busco.trim.test.newick
# 正式建树，使用超算
# FastTree Domain/busco.trim.fa > Domain/busco.trim.ML.newick
FastTree /share/home/wangq/chenli/busco.trim.fa > /share/home/wangq/chenli/busco.trim.ML.newick

cd ~/yeast/tree
nw_reroot  ../Domain/busco.trim.ML.newick Saccharomyces_paradoxus |
    nwr ops order stdin --nd --an \
    > busco_S_cerevisiae_strain_tree.newick

```