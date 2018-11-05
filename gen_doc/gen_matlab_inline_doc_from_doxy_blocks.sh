#!/bin/bash

# this script takes a matlab script documented as defined by doxygen matlab plugin
# parses the doxy doc block and generates a new file .m with inline doc corresponding to the doxygen doc (filtering doxygen tags)
# if a inner code comment block is already present in a function (in the start, just after its name) it will be replaced


function usage {
        echo "USAGE: $0 <input_path.m> <output_path.m>"
}

function usage_exit {
        usage
        exit $1
}

function parse_doxy_block {
        local FUNC=$(echo $1 | tr '[a-z]' '[A-Z]')
        local FUNC_LINE=$(echo "$3" | sed -e 's/^[[:blank:]]\{1,\}//')
        echo -e "$2" | sed -e "{
                                                        s/^[[:blank:]]*%=\{4,\}//;
                                                        s/%>[[:blank:]]*@brief/%% $FUNC/;
                                                        s/^\([[:blank:]]\{1,\}\)%> @param F/\1% SPECS: $FUNC_LINE:\n\1\t% - Parameter F: /;
                                                        s/%>/%/;
                                                         #s/^/\\t/;
                                                        s/Usage/Usage:/
                                                        s/@param[[:blank:]]\([^[:blank:]]\{1,\}\)/- Parameter \1: /;
                                                        s/@code/=========== code ===========/;
                                                        s/@endcode/========= end code =========/;
                                                        s/\(Examples\{0,1\}\)/\1:/
                                                        s/<[^>]*>//g;
                                                        s/&nbsp;/ /g;
                                                        s/@retval/Returns/g;
                                                        s/@[^[:blank:]]* \{0,1\}//g;
                                                        s/\\\./\./g
                                                }"
}

# parse input args

[[ $# -ge 2 && -r "$1" ]] && INPUT_FILE="$1" || usage_exit 1
[[ -z $(dirname "$2") || -d "$(dirname $2)" ]] && OUTPUT_FILE="$2" || usage_exit 2

# get function names
FUNCTIONS=$(sed -ne 's/^[[:blank:]]*function[[:blank:]]\{1,\}\([^=]\{1,\}=\)\{0,1\}[[:blank:]]*\([^([:blank:]]\{1,\}\).*$/\2/p' "$INPUT_FILE")

[[ -z "${FUNCTIONS}" ]] && echo "ERROR in getting the list of functions in $INPUT_FILE" && exit 3


cp "$INPUT_FILE" "$OUTPUT_FILE"

[ ! "$?" = 0 ] && echo "ERROR in copying file." >&2 && exit 4

DOXY_BLOCKS_LIMITS=($(sed -ne '/[[:blank:]]*\%=\{4,\}[[:blank:]]*/=' "$INPUT_FILE"))
typeset -i DOXY_BLOCK_ID=0

for FUNC in $FUNCTIONS
do
        # get the line numbers of the current inline code block
        # INLINE_DOC_LINE_NUMS=($(sed -ne '/^[[:blank:]]*function[[:blank:]]\{1,\}\([^=]\{1,\}=\)\{0,1\}[[:blank:]]*\('$FUNC'\)[[:blank:]]*(.*$/,/^[[:blank:]]*[^%]\{1,\}$/=' $OUTPUT_FILE))
        INLINE_DOC_START_NUM=$(sed -ne '/^[[:blank:]]*function[[:blank:]]\{1,\}\([^=]\{1,\}=\)\{0,1\}[[:blank:]]*\('$FUNC'\)[[:blank:]]*(.*$/='  $OUTPUT_FILE)
        NLINES=$(($(sed -ne '/^[[:blank:]]*function[[:blank:]]\{1,\}\([^=]\{1,\}=\)\{0,1\}[[:blank:]]*\('$FUNC'\)[[:blank:]]*(.*$/,/^[[:blank:]]*[^%]\{1,\}$/p' $OUTPUT_FILE | wc -l)-1))
        INLINE_DOC_END_NUM=$(($INLINE_DOC_START_NUM+$NLINES))
        FILE_NLINES=$(wc -l ${OUTPUT_FILE} | awk '{print $1}')
        [[ "$DEBUG" = "INLINE_DOC_BLOCK" ]] && echo "INLINE_DOC_START_NUM=$INLINE_DOC_START_NUM NLINES=$NLINES INLINE_DOC_END_NUM=$INLINE_DOC_END_NUM"
        # delete not empty inline code doc block
        # if [ ${#INLINE_DOC_LINE_NUMS[*]} -gt  2 ]
        # then
                # INLINE_DOC_START_NUM=${INLINE_DOC_LINE_NUMS[1]}
                # INLINE_DOC_END_NUM=${INLINE_DOC_LINE_NUMS[-2]}
                #[[ $INLINE_DOC_START_NUM =~ [[:digit:]]+ && $INLINE_DOC_END_NUM =~ [[:digit:]]+ ]] && sed -i "${INLINE_DOC_START_NUM},${INLINE_DOC_END_NUM}d" ${OUTPUT_FILE}
                [[ $INLINE_DOC_START_NUM =~ [[:digit:]]+ && $INLINE_DOC_END_NUM =~ [[:digit:]]+ ]] && sed -ne "1,${INLINE_DOC_START_NUM}p;${INLINE_DOC_END_NUM},${FILE_NLINES}p" "${OUTPUT_FILE}" > ${OUTPUT_FILE}_tmp && mv ${OUTPUT_FILE}_tmp ${OUTPUT_FILE}
                [  ! "$?" = 0 ] && echo "ERROR sed deleting old inline doc block." >&2 && exit 5
        # fi

        # get the doxy block for the function
        DOXY_BLOCK=$(sed -ne "$((${DOXY_BLOCKS_LIMITS[$DOXY_BLOCK_ID]}+1)),$((${DOXY_BLOCKS_LIMITS[$(($DOXY_BLOCK_ID+1))]}-1))p" "$INPUT_FILE")
        [[ "$DEBUG" = "DOXY_BLOCK" ]] && echo "DOXY_BLOCK=$DOXY_BLOCK"


        # get the new line position of the function
        FUNC_LINE=$(sed -ne '/^[[:blank:]]*function[[:blank:]]\{1,\}\([^=]\{1,\}=\)\{0,1\}[[:blank:]]*\('$FUNC'\)[[:blank:]]*(.*$/=' "$OUTPUT_FILE")
        [[ "$DEBUG" = "FUNC_LINE" ]] && echo FUNC_LINE=$FUNC_LINE FUNC=$FUNC >&2
        [[ ! $FUNC_LINE =~ [[:digit:]]+ ]] && echo "ERROR in getting function $FUNC line number in $OUTPUT_FILE." >&2 && DOXY_BLOCK_ID+=2 && continue

        # insert parsed block
        FUNC_DEF_LINE=$(sed -ne $FUNC_LINE'p' $OUTPUT_FILE)
    sed -ne '1,'$FUNC_LINE'p' $OUTPUT_FILE > ${OUTPUT_FILE}_tmp
            # parse doxy block for inline block
        INLINE_BLOCK=$(parse_doxy_block $FUNC "$DOXY_BLOCK" "$FUNC_DEF_LINE")
        echo -e "$INLINE_BLOCK" >> ${OUTPUT_FILE}_tmp
        sed -ne $(($FUNC_LINE+1))','$(wc -l ${OUTPUT_FILE} | awk '{print $1}')'p' $OUTPUT_FILE >> ${OUTPUT_FILE}_tmp
        [ "$?" = 0 ] && mv ${OUTPUT_FILE}_tmp ${OUTPUT_FILE}
        DOXY_BLOCK_ID+=2
done



# filter header for class doc
HEADER_START_END=($(sed -ne '/^[[:blank:]]*% ==*/=' $OUTPUT_FILE))
OUT_NUM_LINES=$(wc -l ${OUTPUT_FILE}| awk '{print $1}')
DOXY_HEADER=$(sed -ne "$((${HEADER_START_END[0]}+1)),$((${HEADER_START_END[1]}-1))p"  ${OUTPUT_FILE})
[[ "$DEBUG" = "DOXY_HEADER" ]] && echo "DOXY_HEADER=$DOXY_HEADER"
INLINE_HEADER=$(parse_doxy_block "FAUST" "$DOXY_HEADER")
[[ "$DEBUG" = "INLINE_HEADER" ]] && echo "INLINE_HEADER=$INLINE_HEADER"
echo -e "$INLINE_HEADER" | sed 's/^[[:blank:]]\{1,\}//' > ${OUTPUT_FILE}_tmp
echo -e '\n\n' >> ${OUTPUT_FILE}_tmp
sed -ne "${HEADER_START_END[0]},${OUT_NUM_LINES}p" $OUTPUT_FILE >> ${OUTPUT_FILE}_tmp
[ "$?" = 0 ] && mv ${OUTPUT_FILE}_tmp ${OUTPUT_FILE}
[[ -z "$DEBUG" ]] && rm  ${OUTPUT_FILE}_tmp
