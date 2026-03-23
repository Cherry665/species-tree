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

$ bash cluster.sh Klebsiella Stutzerimonas

"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo $USAGE
    exit 0
fi

#----------------------------#
# Run
#----------------------------#
log_warn Protein/cluster.sh

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

OUT_DIR="ASSEMBLY"
#----------------------------#
# Clustering .95 .95
#----------------------------#
# The min sequence identity for clustering
# The min coverage of query and target for clustering
log_info "Clustering .95 .95"

if [[ -s "${OUT_DIR}/pro.fa.gz" ]] && [[ ! -s "${OUT_DIR}/rep_seq.fa.gz" ]]; then
    log_debug "Running mmseqs easy-cluster 0.95"
    
    mmseqs easy-cluster "${OUT_DIR}"/pro.fa.gz "${OUT_DIR}"/rep "${OUT_DIR}"/tmp \
        --threads 16 --remove-tmp-files -v 0 \
        --min-seq-id 0.95 -c 0.95

    rm -f "${OUT_DIR}"/rep_all_seqs.fasta
    hnsm gz "${OUT_DIR}"/rep_rep_seq.fasta -o "${OUT_DIR}"/rep_seq.fa
    rm -f "${OUT_DIR}"/rep_rep_seq.fasta
fi

#----------------------------#
# Family .8 .8 
#----------------------------#
log_info "Family .8 .8"

if [[ -s "${OUT_DIR}/rep_seq.fa.gz" ]] && [[ ! -s "${OUT_DIR}/fam88_cluster.tsv" ]]; then
    log_debug "Running mmseqs easy-cluster 0.8"

    mmseqs easy-cluster "${OUT_DIR}"/rep_seq.fa.gz "${OUT_DIR}"/fam88 "${OUT_DIR}"/tmp \
        --threads 16 --remove-tmp-files -v 0 \
        --min-seq-id 0.8 -c 0.8

    rm -f "${OUT_DIR}"/fam88_all_seqs.fasta
    rm -f "${OUT_DIR}"/fam88_rep_seq.fasta
fi

#----------------------------#
# Family .3 .8 
#----------------------------#
log_info "Family .3 .8"

if [[ -s "${OUT_DIR}/rep_seq.fa.gz" ]] && [[ ! -s "${OUT_DIR}/fam38_cluster.tsv" ]]; then
    log_debug "Running mmseqs easy-cluster 0.3"

    mmseqs easy-cluster "${OUT_DIR}"/rep_seq.fa.gz "${OUT_DIR}"/fam38 "${OUT_DIR}"/tmp \
        --threads 16 --remove-tmp-files -v 0 \
        --min-seq-id 0.3 -c 0.8

    rm -f "${OUT_DIR}"/fam38_all_seqs.fasta
    rm -f "${OUT_DIR}"/fam38_rep_seq.fasta
fi

log_info Done.

exit 0
