#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: GenerateGenomeFileForBedtoolsCoverage
doc: ''

requirements:
- class: InitialWorkDirRequirement
  listing:
  - entryname: GenerateGenomeFileForBedtoolsCoverage-script.py
    entry: |2

      import argparse, json, sys
      from typing import Optional, List, Dict, Any
      cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")
      cli.add_argument("--reference", type=str, required=True)
      cli.add_argument("--output_filename", type=str, help='Filename to output to')

      FastaDict = str
      String = str
      Filename = str
      Int = int
      Float = float
      Double = float
      Boolean = str
      File = str
      Directory = str
      Stdout = str
      Stderr = str
      Array = List
      class PythonTool:
          File = str
          Directory = str



      def code_block(
          reference: FastaDict, output_filename: str = "genome_file.txt"
      ) -> Dict[str, Any]:
          """
          :param reference: Reference file to generate genome for (must have ^.dict) pattern
          :param output_filename: Filename to output to
          """
          from re import sub

          ref_dict = sub("\.fa(sta)?$", ".dict", reference)

          with open(ref_dict) as inp, open(output_filename, "w+") as out:
              for line in inp:
                  if not line.startswith("@SQ"):
                      continue
                  pieces = line.split("\t")
                  chrom = pieces[1].replace("SN:", "")
                  length = pieces[2].replace("LN:", "")

                  out.write(f"{chrom}\t{length}\n")

              return {"out": output_filename}


      try:
          args = cli.parse_args()
          result = code_block(reference=args.reference, output_filename=args.output_filename)

          from os import getcwd, path
          cwd = getcwd()
          def prepare_file_or_directory_type(file_or_directory, value):
              if value is None:
                  return None
              if isinstance(value, list):
                  return [prepare_file_or_directory_type(file_or_directory, v) for v in value]
              return {
                  "class": file_or_directory,
                  "path": path.join(cwd, value)
              }
          result["out"] = prepare_file_or_directory_type("File", result["out"])
          print(json.dumps(result))
      except Exception as e:
          print(str(e), file=sys.stderr)
          raise
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: python:3.8.1

inputs:
- id: reference
  label: reference
  type: File
  secondaryFiles:
  - ^.dict
  inputBinding:
    prefix: --reference
- id: output_filename
  label: output_filename
  doc: Filename to output to
  type: string
  default: genome_file.txt
  inputBinding:
    prefix: --output_filename

outputs:
- id: out
  label: out
  doc: Genome file for BedToolsCoverage
  type: File
stdout: cwl.output.json
stderr: python-capture.stderr

baseCommand:
- python
- GenerateGenomeFileForBedtoolsCoverage-script.py
id: GenerateGenomeFileForBedtoolsCoverage
