#!/usr/bin/env python3
"""
The script adds crossMapped structures from previous iteration to resMapped structures
from current iteration.
"""
import argparse
import json

class CombineFiles:
    """
    Combines the content of two files
    """
    def __init__(self):
        self.colnames = 'PDB_id,Chain_id,Domain,Family_id,PDB_start,PDB_end,Domain_pos,' \
                        + 'UNP_id,UNP_start,UNP_end'

    @staticmethod
    def csv2list(filecsv):
        """Function to convert a csv file into a list"""
        data = open(filecsv).read().split('\n')
        data = data[1:] if data[0].startswith('PDB') else data
        while '' in data:
            data.remove('')
        return data

    @staticmethod
    def json2list(filejson, data):
        """Function to convert a json file into a list"""
        cross_mapped = json.load(open(filejson, 'r'))
        for fam in cross_mapped:
            if fam == 'unmapped':
                continue
            data.extend(cross_mapped[fam])
        return data


    # def combine_file_content(self):
    #     """
    #     Adds the crossMapped structures from last iteration to the resMapped structures from
    #     current iteration
    #     """
    #     data = open(self.filecsv).read().split('\n')
    #     while '' in data:
    #         data.remove('')

    #     cross_mapped = json.load(open(self.filejson, 'r'))

    #     for fam in cross_mapped:
    #         if fam == 'unmapped':
    #             continue
    #         data.extend(cross_mapped[fam])
    #     return data

    def __str__(self):
        return "Adding crossMapped structures from last iteration to resMapped from current "\
            "iteration, if any."



##########################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-p", "--pfam", type=str, \
        help="File containing res-mapped domain StIs from Pfam family(ies)")
    P.add_argument("-c", "--cath", type=str, \
        help="File containing res-mapped domain StIs from CATH family(ies)")
    P.add_argument("-px", "--pfamX", type=str, \
        help="File containing cross-mapped domain StIs from Pfam family(ies)")
    P.add_argument("-cx", "--cathX", type=str, \
        help="File containing cross-mapped domain StIs from CATH family(ies)")
    P.add_argument("-pr", "--pfamresult", type=str, default="pfam_res_crossMapped.csv", \
        help="Output file with cross-mapped and res-mapped domain StIs from Pfam family(ies)")
    P.add_argument("-cr", "--cathresult", type=str, default="cath_res_crossMapped.csv", \
        help="Output file with cross-mapped and res-mapped domain StIs from CATH family(ies)")
    ARGS = P.parse_args()

    COLS = 'PDB_id,Chain_id,Domain,Family_id,PDB_start,PDB_end,Domain_pos,UNP_id,UNP_start,UNP_end'

    PF = CombineFiles()
    CROSSED = PF.csv2list(ARGS.pfam) if ARGS.pfam is not None else []
    CROSSED = PF.json2list(ARGS.pfamX, CROSSED) if ARGS.pfamX is not None else CROSSED
    CROSSED.insert(0, COLS)

    # CROSSED = PF.combine_file_content()
    with open(ARGS.pfamresult, 'w') as wrt:
        wrt.write('\n'.join(CROSSED))

    CA = CombineFiles()
    CROSSED = PF.csv2list(ARGS.cath) if ARGS.cath is not None else []
    CROSSED = PF.json2list(ARGS.cathX, CROSSED) if ARGS.cathX is not None else CROSSED
    CROSSED.insert(0, COLS)

    # CROSSED = CA.combine_file_content()
    with open(ARGS.cathresult, 'w') as wrt:
        wrt.write('\n'.join(CROSSED))
