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

$ bash info.sh Klebsiella Stutzerimonas

"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo $USAGE
    exit 0
fi

#----------------------------#
# Run
#----------------------------#
log_warn Protein/info.sh

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
# seq.sqlite
#----------------------------#
OUT_DIR="ASSEMBLY"
log_info "seq.sqlite"
if [[ -f "${OUT_DIR}/seq.sqlite" ]]; then
    log_info "seq.sqlite already exists, skipping."
else
    log_debug "Initializing seqdb in ${OUT_DIR}"

    # 初始化数据库
    nwr seqdb -d "${OUT_DIR}" --init --strain

    # 导入序列长度信息
    nwr seqdb -d "${OUT_DIR}" \
        --size <(
            hnsm size "${OUT_DIR}"/pro.fa.gz
        ) \
        --clust

    # 导入注释与菌株映射关系 
    nwr seqdb -d "${OUT_DIR}" \
        --anno <(
            gzip -dcf "${OUT_DIR}"/anno.tsv.gz
        ) \
        --asmseq <(
            gzip -dcf "${OUT_DIR}"/asmseq.tsv.gz
        )

    # 导入不同层级的聚类结果
    nwr seqdb -d "${OUT_DIR}" --rep f1="${OUT_DIR}"/fam88_cluster.tsv
    nwr seqdb -d "${OUT_DIR}" --rep f2="${OUT_DIR}"/fam38_cluster.tsv
fi

log_info Done.

exit 0
