cwlVersion: v1.0
class: CommandLineTool
id: exomedepth_filter
label: exomedepth_filter

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
      MIN_BF=$4
      TEMP_NAME=${INPUT_FILE/csv/filtered.csv}
      TEMP_FILE=$(echo "$TEMP_NAME" | sed "s/.*\///")
      OUTPUT_NAME=${TEMP_FILE/csv/bed}
      OUTPUT_FILE=$(echo "$OUTPUT_NAME" | sed "s/.*\///")

      tail -n +2 $INPUT_FILE | sed 's/"//g' | awk -F"," -v maxLen=$MAX_LEN -v minLen=$MIN_LEN -v minBf=$MIN_BF '{ \
      chr=$7; \
      start=$5; \
      end=$6; \
      sample=$13; \
      bf=$9; \
      len=(end-start)/1000; \
      type=$3=="deletion"?"DEL":"DUP"; \
      tool="exomeDepth"; \
      if(len>minLen && \
         len<maxLen && \
         bf>=minBf){
           print chr"\t"start"\t"end"\t"sample"\t"type"\t"bf"\t"tool}}' > ${TEMP_FILE}

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
  min_bf:
    type: string
    inputBinding:
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bed"
