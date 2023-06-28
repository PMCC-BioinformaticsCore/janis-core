#!/usr/bin/env python3
"""
This script generates/returns the parameter file for next iteration.
"""

import argparse
import os
import json
import yaml

class WriteParam:
    """
    Class to write parameter file for next iterations
    """
    def __init__(self, pfam_xmap, cath_xmap):
        self.pfam_id = pfam_xmap
        self.cath_id = cath_xmap

    @staticmethod
    def create_attributes(out_loc, parm_solo):
        """
        This assigns the attribute path/location of outdir instead of tmp_dir (default)
        """
        attribute_list = ['class', 'basename', 'format', 'location']
        tmp_new = {k:v for k, v in parm_solo.items() if k in attribute_list}
        tmp_new['location'] = os.path.join(out_loc, parm_solo['basename'])
        return tmp_new


    def write_parmfile(self, in_parm, iter_parm, out_parm):
        """
        Reads the iteration file for current iteration and writes the next iteration parameter file
        """
        with open(in_parm) as wrty:
            inp_parm = yaml.load(wrty, Loader=yaml.FullLoader)
            # globals().update(inp_parm.get("arguments"))

        inp_parm['pfam'] = self.pfam_id
        inp_parm['cath'] = self.cath_id
        inp_parm['iteration'] += 1
        loc = inp_parm['outputdir']['location'] # getting outdir location
        
        # Current iteration parameters
        iter_parm = json.load(open(iter_parm, 'r'))

        inp_parm['filename'] = self.create_attributes(loc, iter_parm['filename'])
        inp_parm['true_domain_file'] = self.create_attributes(loc, iter_parm['true_domain_file'])
        inp_parm['prev_crossMapped_pfam'] = self.create_attributes(loc, iter_parm['crossmap_pfam'])
        inp_parm['prev_crossMapped_cath'] = self.create_attributes(loc, iter_parm['crossmap_cath'])
        inp_parm['core_domain_struct'] = self.create_attributes(loc, iter_parm['core_domain_struct'])
        inp_parm['domain_like'] = self.create_attributes(loc, iter_parm['domain_like'])
        inp_parm['failed_domain'] = self.create_attributes(loc, iter_parm['failed_domain'])
        inp_parm['next_parmfile'] = inp_parm['next_parmfile'].replace('.yml', '1.yml')
        if isinstance(iter_parm['pfam_lost_structs'], dict):
            in_parm['pfam_lost'] = self.create_attributes(loc, iter_parm['pfam_lost_structs'])
        if isinstance(iter_parm['cath_lost_structs'], dict):
            in_parm['cath_lost'] = self.create_attributes(loc, iter_parm['cath_lost_structs'])

        if (len(self.pfam_id) == 0) and (len(self.cath_id) == 0):
            inp_parm['last_iteration'] = -1

        # Writing new param File for next iteration
        with open(out_parm, 'w') as wrty:
            yaml.dump(inp_parm, wrty, sort_keys=True, indent=2)

        return out_parm


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-i", "--input", type=str, required=True, \
        help="Parameter file (yml) for current iteration")
    P.add_argument("-px", "--pfam_xmap", type=str, \
        help="Cross-mapped Pfam families passing the threshold")
    P.add_argument("-cx", "--cath_xmap", type=str, \
        help="Cross-mapped CATH superfamilies passing the threshold")
    P.add_argument("-o", "--outparam", type=str, default="new_param.yml", \
        help="Output filename for next iteration parameters")
    ARGS = P.parse_args()

    curr_iter_parm = 'runtime_parm.json'
    cath = json.load(open(ARGS.pfam_xmap, "r")).keys() if ARGS.pfam_xmap is not None else []
    pfam = json.load(open(ARGS.cath_xmap, "r")).keys() if ARGS.cath_xmap is not None else []

    WP = WriteParam(list(pfam), list(cath))
    OUT_FILE = WP.write_parmfile(ARGS.input, curr_iter_parm, ARGS.outparam)
    print(OUT_FILE)
