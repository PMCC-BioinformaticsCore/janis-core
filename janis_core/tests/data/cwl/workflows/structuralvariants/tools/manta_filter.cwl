cwlVersion: v1.0
class: CommandLineTool
id: manta_filter
label: manta_filter

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
      MIN_Q=$4
      OUTPUT_NAME=${INPUT_FILE/raw/filtered}
      OUTPUT_FILE=$(echo "$OUTPUT_NAME" | sed "s/.*\///")

      # mapping case_id with sample_id
      SAMPLE_ID=(${OUTPUT_FILE//./ })
      CASE_ID=$(awk -v search="$SAMPLE_ID" '$0 ~ search{print $1}' "$SAMPLE_FILE")

      tail -n +2 $INPUT_FILE | awk -v maxLen=$MAX_LEN -v minLen=$MIN_LEN -v minQ=$MIN_Q '{ \
      chr=$1; \
      start=$2; \
      end=$6; \
      sample="'$CASE_ID'"; \
      q=$8; \
      len=(end-start)/1000; \
      type=$11; \
      tool="manta"; \
      if(len>minLen && \
        len<maxLen && \
        q>=minQ && \
        (type=="DEL" || type=="DUP")){
          print chr"\t"start"\t"end"\t"sample"\t"type"\t"q"\t"tool}}' > ${OUTPUT_FILE}

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
  min_q:
    type: string
    inputBinding:
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: "*"
