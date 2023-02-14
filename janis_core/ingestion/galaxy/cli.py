

import argparse
from typing import Optional


useage_str = '''
gxtool2janis command [OPTIONS]

Commands:
   tool         Parse galaxy tool (.xml format)
   workflow     Parse galaxy workflow (.ga format)

'''

class CLIparser:
    def __init__(self, argv: list[str]):
        parser = argparse.ArgumentParser(
            description='gxtool2janis.py',
            usage=useage_str
        )
        parser.add_argument('command', help='Command')
        args = parser.parse_args(argv[1:2])
        if not hasattr(self, args.command):
            print(f'Unrecognized command {argv[1]}')
            parser.print_help()
            exit(1)
        # dispatch pattern (calls the function with same name as args.command)
        self.command = args.command
        getattr(self, args.command)(argv)

    def tool(self, argv: list[str]):
        parser = argparse.ArgumentParser(
            description='Parse single galaxy tool'
        )
        parser.add_argument("infile",
                            help="path to tool.xml file to parse.", 
                            type=str,
                            )
        parser.add_argument("-r",
                            "--remote", 
                            help="comma separated list of [owner],[repo],[toolid],[revision] if downloading wrapper from toolshed. example: devteam,picard,picard_MarkDuplicates,29:1aac2a13842a",
                            type=str
                            )
        parser.add_argument("-o",
                            "--outdir",
                            help="output folder to place translation", 
                            type=str,
                            )
        args = parser.parse_args(argv[2:])
        out: dict[str, Optional[str]] = args.__dict__
        out['command'] = self.command
        self.args = out

    def workflow(self, argv: list[str]):
        parser = argparse.ArgumentParser(
            description='') 
        parser.add_argument("infile", 
                            help="path to workflow.ga file to parse.", 
                            type=str)
        parser.add_argument("-o",
                            "--outdir",
                            help="output folder to place translation", 
                            type=str,
                            )
        parser.add_argument("--dev-partial-eval", 
                            help="turn off partial cheetah evaluation when identifying tool values",
                            default=False,
                            action='store_true')
        args = parser.parse_args(argv[2:])
        out: dict[str, Optional[str]] = args.__dict__
        out['command'] = self.command
        self.args = out

