# 文件处理

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

# 检查

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

# MinHash
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
# ANI_VALUE=0.001
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

# 合并 minhash tree 的分支
```shell
mkdir -p ~/yeast/tree
cd ~/yeast/tree

# 将树的根设置在 Saccharomyces_paradoxus 上
nw_reroot ../MinHash/tree.nwk Saccharomyces_paradoxus |
    nwr ops order stdin --nd --an \
    > minhash.reroot.newick

nwr ops order --nd --an minhash.reroot.newick > minhash_S_cerevisiae_strain_tree.nwk
```

# 统计有效的物种和菌株数量