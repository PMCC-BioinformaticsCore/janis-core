cwlVersion: v1.0
class: CommandLineTool
id: codex_filter
label: codex_filter

requirements:
  ShellCommandRequirement: {}

baseCommand: [bash, -c]

inputs:
  script:
    type: string?
    inputBinding:
      position: 1
    default: |
      INPUT_FILE=$0
      SAMPLE_FILE=$1
      MIN_LEN=$2
      MAX_LEN=$3
      MIN_LRATIO=$4
      TEMP_NAME=${INPUT_FILE/_Chrom.txt/.filtered.txt}
      TEMP_FILE=$(echo "$TEMP_NAME" | sed "s/.*\///")
      OUTPUT_NAME=${TEMP_FILE/txt/bed}
      OUTPUT_FILE=$(echo "$OUTPUT_NAME" | sed "s/.*\///")

      tail -n +2 $INPUT_FILE | awk -v maxLen=$MAX_LEN -v minLen=$MIN_LEN -v minLRatio=$MIN_LRATIO '{ \
      chr=$2; \
      start=$4; \
      end=$5; \
      sample=$1; \
      len=$6; \
      lratio=$12; \
      st_exon=$7; \
      ed_exon=$8; \
      type=$3=="del"?"DEL":"DUP"; \
      tool="codex"; \
      if(len>minLen && \
         len<maxLen && \
         len/(ed_exon-st_exon)<50 && \
         lratio>=minLRatio){
           print chr"\t"start"\t"end"\t"sample"\t"type"\t"lratio"\t"tool}}' > ${TEMP_FILE}

      # mapping case_id with sample_id
      join -1 4 -2 2 -o 1.1,1.2,1.3,2.1,1.5,1.6,1.7 <(sort -k 4 $TEMP_FILE) <(sort -k 2 $SAMPLE_FILE) > ${OUTPUT_FILE}

  input:
    type: File
    inputBinding:
      position: 2
  samples:
    type: File
    inputBinding:
      position: 3
  min_len:
    type: string
    inputBinding:
      position: 4
  max_len:
    type: string
    inputBinding:
      position: 5
  min_lratio:
    type: string
    inputBinding:
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bed"
