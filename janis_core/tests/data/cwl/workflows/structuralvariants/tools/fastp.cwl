cwlVersion: v1.0
class: CommandLineTool
id: fastp
label: fastp

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.fastq1)
      - $(inputs.fastq2)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0'

baseCommand: [fastp]

arguments:
  - prefix: '-o'
    valueFrom: $(inputs.fastq1.nameroot.replace(/\b.fastq\b/g, '')).trimmed.fastq
  - |
    ${
      if (inputs.fastq2){
        return '-O';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.fastq2){
        return inputs.fastq2.nameroot.replace(/\b.fastq\b/g, '') + '.trimmed.fastq';
      } else {
        return '';
      }
    }
  - prefix: '--unpaired1'
    valueFrom: $(inputs.fastq1.nameroot.replace(/\b.fastq\b/g, '')).unpaired.trimmed.fastq
  - |
    ${
      if (inputs.fastq2){
        return '--unpaired2';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.fastq2){
        return inputs.fastq2.nameroot.replace(/\b.fastq\b/g, '') + ".unpaired.trimmed.fastq";
      } else {
        return '';
      }
    }

inputs:
  fastq1:
    type: File
    inputBinding:
      prefix: -i
  fastq2:
    type: File
    inputBinding:
      prefix: -I
  threads:
    type: int?
    default: 12
    inputBinding:
      position: 1
      prefix: '--thread'
  cut_right:
    type: boolean?
    default: true
    inputBinding:
      prefix: '--cut_right'
  cut_right_window_size:
    type: int?
    default: 5
    inputBinding:
      prefix: '--cut_right_window_size'
  cut_right_mean_quality:
    type: int?
    default: 24
    inputBinding:
      prefix: '--cut_right_mean_quality'
  trim_tail1:
    type: int?
    default: 1
    inputBinding:
      prefix: '--trim_tail1'
  length_required:
    type: int?
    default: 70
    inputBinding:
      prefix: '--length_required'

outputs:
  paired_fastq1:
    type: File
    outputBinding:
      glob: $(inputs.fastq1.nameroot.replace(/\b.fastq\b/g, '')).trimmed.fastq
  paired_fastq2:
    type: File
    outputBinding:
      glob: $(inputs.fastq2.nameroot.replace(/\b.fastq\b/g, '')).trimmed.fastq
  unpaired_fastq1:
    type: File
    outputBinding:
      glob: $(inputs.fastq1.nameroot.replace(/\b.fastq\b/g, '')).unpaired.trimmed.fastq
  unpaired_fastq2:
    type: File
    outputBinding:
      glob: $(inputs.fastq2.nameroot.replace(/\b.fastq\b/g, '')).unpaired.trimmed.fastq
  html_report:
    type: File
    outputBinding:
      glob: 'fastp.html'
      outputEval: |
        ${
          self[0].basename = inputs.fastq1.nameroot.split('_', 2)[0] + ".fastp.html";
          return self[0]
        }
  json_report:
    type: File
    outputBinding:
      glob: 'fastp.json'
      outputEval: |
        ${
          self[0].basename = inputs.fastq1.nameroot.split('_', 2)[0] + ".fastp.json";
          return self[0]
        }
