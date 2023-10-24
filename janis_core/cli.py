

import argparse
import sys 

from janis_core.ingestion import SupportedIngestion 
from janis_core.translation_deps.supportedtranslations import SupportedTranslation
from janis_core.ingestion import ingest
from janis_core.translations import translate


def main() -> None:
    sysargs = sys.argv[1:]
    args_namespace = parse_args(sysargs)
    args_dict = interpret_args(args_namespace)
    do_translate(args_dict)

def do_translate(args: dict[str, str]) -> None:
    internal = ingest(args['infile'], args['from']) 
    return translate(internal, args['to'], mode=args['mode'], export_path=args['outdir'])

def interpret_args(args: argparse.Namespace) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, val in args._get_kwargs():  # workaround for '--from' name: usually a python error.
        if key == 'from':
            out['from'] = val
        elif key == 'to':
            out['to'] = val
        elif key == 'mode':
            out['mode'] = val
        elif key == 'infile':
            out['infile'] = val
        elif key == 'output_dir':
            out['outdir'] = val
    return out

def parse_args(sysargs: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Translate a janis workflow to CWL, WDL or Nextflow')

    parser.add_argument(
        "infile", 
        help="Path to input file",
    )
    parser.add_argument(
        "--from",
        help="Language of infile. Will be autodetected if not supplied",
        choices=SupportedIngestion.all(),
        type=str
    )
    parser.add_argument(
        "--to",
        help="Language to translate to.",
        choices=SupportedTranslation.all(),
        type=str
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory to write output to (default: translated).",
        type=str,
        default="translated"
    )
    parser.add_argument(
        "--mode",
        help="Translate mode (default: regular). Controls extent of tool translation\n\
        - skeleton: ignores inputs which aren't used in workflow. no CLI command generation.\n\
        - regular: ignores inputs which aren't used in workflow. \n\
        - extended: full translation of all inputs & CLI command",
        type=str,
        choices=["skeleton", "regular", "extended"],
        default="regular"
    )

    return parser.parse_args(sysargs)

