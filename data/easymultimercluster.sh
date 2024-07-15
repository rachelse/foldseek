#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

exists() {
	[ -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

mapCmplName2ChainKeys() {
    awk -F"\t" 'FNR==1 {++fIndex}
        fIndex==1 {
            repName[$1]=1
            if (match($1, /MODEL/)){
                tmpName[$1]=1
            }else{
                tmpName[$1"_MODEL_1"]=1 
            }
            next
        }
        fIndex==2{
            if (match($2, /MODEL/)){
                if ($2 in tmpName){
                repId[$1]=1
                }else{
                    ho[1]=1
                }
            }else{
                if ($2 in repName){
                repId[$1]=1
                }
            }
            next
        }
        {
            if ($3 in repId){
                print $1
            }
        }
    ' "${1}" "${2}.source" "${2}.lookup" > "${3}"
}

postprocessFasta() {
    awk ' BEGIN {FS=">"}
    $0 ~/^>/ {
        # match($2, /(.*).pdb*/)
        split($2,parts,"_")
        complex=""
        for (j = 1; j < length(parts); j++) {
            complex = complex parts[j]
            if (j < length(parts)-1){
                complex=complex"_" 
            }
        }
        if (!(complex in repComplex)) {
            print "#"complex
            repComplex[complex] = ""
        }
    }
    {print $0}
    ' "${1}" > "${1}.tmp" && mv "${1}.tmp" "${1}"
}

if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "input createdb died"
fi

if notExists "${TMP_PATH}/complex_clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimercluster "${TMP_PATH}/input" "${TMP_PATH}/complex_clu" "${TMP_PATH}" ${MULTIMERCLUSTER_PAR} \
        || fail "Multimercluster died"
fi

SOURCE="${TMP_PATH}/input"
INPUT="${TMP_PATH}/latest/complex_db"
if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clu" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "Convert Alignments died"
fi

if notExists "${TMP_PATH}/complex_rep_seqs.dbtype"; then
    mapCmplName2ChainKeys "${TMP_PATH}/cluster.tsv" "${SOURCE}" "${TMP_PATH}/rep_seqs.list" 
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/rep_seqs.list" "${SOURCE}" "${TMP_PATH}/complex_rep_seqs" ${CREATESUBDB_PAR} \
        || fail "createsubdb died"
fi

if notExists "${TMP_PATH}/complex_rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${SOURCE}" "${SOURCE}"  "${TMP_PATH}/complex_rep_seqs" "${TMP_PATH}/complex_rep_seq.fasta" ${VERBOSITY_PAR} \
            || fail "result2flat died"
    postprocessFasta "${TMP_PATH}/complex_rep_seq.fasta"
fi

#TODO: generate fasta file for all sequences
# if notExists "${TMP_PATH}/complex_all_seqs.fasta"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" createseqfiledb "${INPUT}" "${TMP_PATH}/complex_clust" "${TMP_PATH}/complex_clust_seqs" ${THREADS_PAR} \
#             || fail "Result2repseq  died"

#     # shellcheck disable=SC2086
#     "$MMSEQS" result2flat "${INPUT}" "${INPUT}" "${TMP_PATH}/complex_clust_seqs" "${TMP_PATH}/complex_all_seqs.fasta" ${VERBOSITY_PAR} \
#             || fail "result2flat died"
# fi

# mv "${TMP_PATH}/complex_all_seqs.fasta"  "${RESULT}_all_seqs.fasta"
mv "${TMP_PATH}/complex_rep_seq.fasta"  "${RESULT}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"  "${RESULT}_cluster.tsv"

if [ -n "${REMOVE_TMP}" ]; then
    rm "${INPUT}.0"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_db" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    # "$MMSEQS" rmdb "${TMP_PATH}/complex_clu_seqs" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_rep_seqs" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_rep_seqs_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/complex_clu" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${INPUT}" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${INPUT}_h" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_ca" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_ss" ${VERBOSITY_PAR}
    rm "${TMP_PATH}/rep_seqs.list"
    rm -rf "${TMP_PATH}/latest"
    rm -f "${TMP_PATH}/easymultimercluster.sh"
fi