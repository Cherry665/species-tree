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

$ bash rsync.sh Klebsiella Stutzerimonas

"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo $USAGE
    exit 0
fi

#----------------------------#
# Run
#----------------------------#
log_warn rsync.sh

touch check.lst

cat url.tsv |
    tsv-join -f check.lst -k 1 -e |
    if [ "$#" -gt 0 ]; then
        # Initialize an string to store the cmd
        result="tsv-filter --or"

        # Iterate over each argument and prepend the fixed string
        for arg in "$@"; do
            result+=" --str-in-fld '3:$arg'"
        done

        # Remove the trailing space from the result string
        result=${result% }

        # Execute the result string as a Bash command
        eval "$result"
    else
        rgr dedup stdin
    fi |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        RAW_URL="{2}"
        CLEAN_URL=$(echo "$RAW_URL" | tr -d "\r\047\042 " | sed "s|^ftp.ncbi.nlm.nih.gov::|https://ftp.ncbi.nlm.nih.gov/|")
        echo >&2 "Debug URL: [${CLEAN_URL}]"

        mkdir -p "{3}/{1}"
    
        wget -r -l1 -np -nd -e robots=off \
            --user-agent="Mozilla/5.0" \
            -P "{3}/{1}/" \
            --reject "*_genomic.gbff.gz,*_genomic.gtf.gz,*_protein.gpff.gz,annotation_hashes.txt,assembly_status.txt" \
            "${CLEAN_URL}/"
    '

log_info Done.

exit 0
