version development

task ParseFastqcAdaptors {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] fastqc_datafiles
    File? cutadapt_adaptors_lookup
  }
  command <<<
    
cat <<EOT >> ParseFastqcAdaptors-script.py
    
import argparse, json, sys
from typing import Optional, List, Dict, Any
cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")
cli.add_argument("--fastqc_datafiles", nargs='+', default=[], type=str)
cli.add_argument("--cutadapt_adaptors_lookup", type=str, help='Specifies a file which contains the list of adapter sequences which will\nbe explicity searched against the library. The file must contain sets of named adapters in\nthe form name[tab]sequence. Lines prefixed with a hash will be ignored.')

Array = List
File = str
String = str
Filename = str
Int = int
Float = float
Double = float
Boolean = str
Directory = str
Stdout = str
Stderr = str
class PythonTool:
    File = str
    Directory = str



def code_block(
    fastqc_datafiles: List[File], cutadapt_adaptors_lookup: Optional[File]
):
    """

    :param fastqc_datafiles:

    :param cutadapt_adaptors_lookup: Specifies a file which contains the list of adapter sequences which will
        be explicity searched against the library. The file must contain sets of named adapters in
        the form name[tab]sequence. Lines prefixed with a hash will be ignored.
    :return:
    """
    if not cutadapt_adaptors_lookup:
        return {"adaptor_sequences": []}

    import mmap, re, csv
    from io import StringIO
    from sys import stderr

    def get_overrepresented_text(f):
        """
        Get the table "Overrepresented sequences" within the fastqc_data.txt
        """
        adapt_section_query = (
            br"(?s)>>Overrepresented sequences\t\S+\n(.*?)>>END_MODULE"
        )
        # fastqc_datafile could be fairly large, so we'll use mmap, and then
        with open(f) as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as fp:
            overrepresented_sequences_match = re.search(adapt_section_query, fp)
            if overrepresented_sequences_match is None:
                raise Exception(
                    f"Couldn't find query ('{adapt_section_query.decode('utf8')}') in {fastqc_datafiles}"
                )

            return overrepresented_sequences_match.groups()[0].decode("utf8")

    def parse_tsv_table(tbl: str, skip_headers):
        """
        Parse a TSV table from a string using csvreader
        """

        rd = csv.reader(StringIO(tbl), delimiter="\t", quotechar='"')
        ret = list(rd)
        if len(ret) == 0:
            return ret
        if skip_headers:
            ret.pop(0)  # discard headers
        return ret

    def get_cutadapt_map():
        """
        Helper method to parse the file 'cutadapt_adaptors_lookup' with
        format: 'name[tab]sequence' into the dictionary: '{ sequence: name }'
        """
        cutadapt_map = {}
        with open(cutadapt_adaptors_lookup) as fp:
            for row in fp:
                st = row.strip()
                if not st or st.startswith("#"):
                    continue

                # In reality, the format is $name[\t+]$seqence (more than one tab)
                # so we'll just split on a tab, and remove all the empty elements.
                split = [f for f in st.split("\t") if bool(f) and len(f) > 0]

                # Invalid format for line, so skip it.
                if len(split) != 2:
                    print(
                        f"Skipping cutadapt line '{st}' as irregular elements ({len(split)})",
                        file=stderr,
                    )
                    continue

                # reverse the order from name[tab]sequence to { sequence: tab }
                cutadapt_map[split[1]] = split[0]
        return cutadapt_map

    # Start doing the work
    adaptor_ids = set()
    for fastqcfile in fastqc_datafiles:
        text = get_overrepresented_text(fastqcfile)
        adaptor_ids = adaptor_ids.union(
            set(a[0] for a in parse_tsv_table(text, skip_headers=True))
        )

    adaptor_sequences = []

    if adaptor_ids:
        cutadapt_map = get_cutadapt_map()
        for aid in adaptor_ids:
            if aid in cutadapt_map:
                print(
                    f"Identified sequence '{aid}' as '{cutadapt_map.get(aid)}' in lookup",
                    file=stderr,
                )
                adaptor_sequences.append(aid)

            else:
                print(
                    f"Couldn't find a corresponding sequence for '{aid}' in lookup map",
                    file=stderr,
                )

    return {"adaptor_sequences": adaptor_sequences}


try:
    args = cli.parse_args()
    result = code_block(fastqc_datafiles=args.fastqc_datafiles, cutadapt_adaptors_lookup=args.cutadapt_adaptors_lookup)

    print(json.dumps(result))
except Exception as e:
    print(str(e), file=sys.stderr)
    raise

EOT
    python ParseFastqcAdaptors-script.py \
      ~{if length(fastqc_datafiles) > 0 then "--fastqc_datafiles '" + sep("' '", fastqc_datafiles) + "'" else ""} \
      ~{if defined(cutadapt_adaptors_lookup) then ("--cutadapt_adaptors_lookup '" + cutadapt_adaptors_lookup + "'") else ""}
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
    Array[String] adaptor_sequences = read_json(stdout())["adaptor_sequences"]
  }
}