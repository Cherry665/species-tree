#!/usr/bin/env bash

BASH_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "${BASH_DIR}"

#----------------------------#
# Colors in term
#----------------------------#
GREEN=
RED=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    NC='\033[0m' # No Color
fi

log_warn () {
    echo >&2 -e "${RED}==> $@ <==${NC}"
}

log_info () {
    echo >&2 -e "${GREEN}==> $@${NC}"
}

log_debug () {
    echo >&2 -e "==> $@"
}

export -f log_warn
export -f log_info
export -f log_debug

#----------------------------#
# helper functions
#----------------------------#
set +e

# set stacksize to unlimited
if [[ "$OSTYPE" != "darwin"* ]]; then
    ulimit -s unlimited
fi

signaled () {
    log_warn Interrupted
    exit 1
}
trap signaled TERM QUIT INT

readlinkf () {
    perl -MCwd -l -e 'print Cwd::abs_path shift' "$1";
}

export -f readlinkf

#----------------------------#
# Usage
#----------------------------#
USAGE="
Usage: $0 [STR_IN_FLD] ...

Default values:
    STR_IN_FLD  ''

$ bash collect.sh Klebsiella Stutzerimonas

"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo $USAGE
    exit 0
fi

#----------------------------#
# Run
#----------------------------#
log_warn Protein/collect.sh

#----------------------------#
# filtered species.tsv
#----------------------------#
log_info "Protein/species-f.tsv"
cat species.tsv |
    sort |
    if [ "$#" -gt 0 ]; then
        # Initialize an string to store the cmd
        result="tsv-filter --or"

        # Iterate over each argument and prepend the fixed string
        for arg in "$@"; do
            result+=" --str-in-fld '2:$arg'"
        done

        # Remove the trailing space from the result string
        result=${result% }

        # Execute the result string as a Bash command
        eval "$result"
    else
        rgr dedup stdin
    fi |
tsv-join -f ../ASSEMBLY/pass.lst -k 1 |
tsv-join -e -f ../ASSEMBLY/omit.lst -k 1 |
cat \
    > species-f.tsv

#----------------------------#
# Unique proteins
#----------------------------#
log_info "Unique proteins"

# 定义一个统一的输出目录
OUT_DIR="ASSEMBLY"
mkdir -p "${OUT_DIR}"

# 直接将所有过滤后的样本作为输入
cp species-f.tsv "${OUT_DIR}/strains.tsv"
rm -f "${OUT_DIR}/detail.tsv"
rm -f "${OUT_DIR}/detail.tsv.gz"

cat "${OUT_DIR}/strains.tsv" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [[ ! -d "../ASSEMBLY/ASSEMBLY/{1}" ]]; then
            echo "Warning: ../ASSEMBLY/ASSEMBLY/{1} not found" >&2
            exit
        fi

        gzip -dcf ../ASSEMBLY/ASSEMBLY/{1}/*_protein.faa.gz |
            grep "^>" |
            sed "s/^>//" |
            awk -v strain={1} "{print \$1 \"\t\" strain \"\t\" \"no_desc\"}" \
            >> "'${OUT_DIR}'"/detail.tsv

        gzip -dcf ../ASSEMBLY/ASSEMBLY/{1}/*_protein.faa.gz
    ' |
    hnsm filter stdin -u |
    hnsm gz stdin -p 4 -o "${OUT_DIR}/pro.fa"

    # 生成下游所需的映射文件
    tsv-select -f 1,3 "${OUT_DIR}/detail.tsv" | rgr dedup stdin | gzip > "${OUT_DIR}/anno.tsv.gz"
    tsv-select -f 1,2 "${OUT_DIR}/detail.tsv" | rgr dedup stdin | gzip > "${OUT_DIR}/asmseq.tsv.gz"
    rm -f "${OUT_DIR}/detail.tsv"

log_info "Protein extraction completed in ${OUT_DIR}"

log_info Done.

exit 0
