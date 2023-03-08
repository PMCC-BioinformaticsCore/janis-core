cwlVersion: v1.0
class: CommandLineTool
id: collapse
label: collapse

doc: |
  Merge CNVs predictions for each tool

requirements:
  ShellCommandRequirement: {}

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bedtools:2.26.0gx--he513fc3_4'


baseCommand: [bash, -c]

inputs:
  script:
    type: string?
    inputBinding:
      position: 1
    default: |
      INPUT_FILE=$0
      OUTPUT_NAME=${INPUT_FILE/bed/sorted.merged.bed}
      OUTPUT_FILE=$(echo "$OUTPUT_NAME" | sed "s/.*\///")

      for i in $(cat $INPUT_FILE | sed -e "s/[[:space:]]\+/\t/g" | cut -f4 | sort -u)  # FIXME: remove sed -e "s/[[:space:]]\+/\t/g" try better solution
      do
        for j in "DEL" "DUP"
        do
          grep $i $INPUT_FILE | grep $j | sed -e "s/[[:space:]]\+/\t/g" | bedtools sort | bedtools merge -c 6,7 -o max,distinct | awk -v sample=$i -v type=$j '{ \
            print $1,$2,$3,sample,type,$4,$5}' OFS="\t"
        done >> ${OUTPUT_FILE}
      done
  input:
    type: File
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bed"
