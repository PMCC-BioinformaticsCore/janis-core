version development

task GenerateVardictHeaderLines {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File reference
    File reference_dict
    String? output_filename
  }
  command <<<
    
cat <<EOT >> GenerateVardictHeaderLines-script.py
    
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
    reference: FastaDict, output_filename: str = "output.txt"
) -> Dict[str, Any]:
    """
    :param reference: Reference file to generate vardict header lines for (must have ^.dict) pattern
    :param output_filename: Filename to output to
    """
    from re import sub

    ref_dict = sub("\.fa(sta)?$", ".dict", reference)

    with open(output_filename, "w+") as out, open(ref_dict) as inp:
        out.write("##source=vardict\n")
        for line in inp:
            if not line.startswith("@SQ"):
                continue
            pieces = line.split("\t")
            chrom = pieces[1].replace("SN:", "")
            length = pieces[2].replace("LN:", "")

            out.write(f"##contig=<ID={chrom},length={length}>\n")

        return {"out": output_filename}


try:
    args = cli.parse_args()
    result = code_block(reference=args.reference, output_filename=args.output_filename)

    print(json.dumps(result))
except Exception as e:
    print(str(e), file=sys.stderr)
    raise

EOT
    python GenerateVardictHeaderLines-script.py \
      --reference '~{reference}' \
      ~{if defined(select_first([output_filename, "output.txt"])) then ("--output_filename '" + select_first([output_filename, "output.txt"]) + "'") else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "python:3.8.1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = read_json(stdout())["out"]
  }
}