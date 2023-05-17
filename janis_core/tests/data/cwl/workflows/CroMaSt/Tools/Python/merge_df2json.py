#!/usr/bin/env python3
"""
Script to collect all passed instances together.
"""
import argparse
import json
import pandas as pd

def df2list(tmp_df, lst):
    """Returns a list with elements as a row from dataframe"""
    for col in tmp_df.columns:
        tmp_df[col] = tmp_df[col].astype(str)
    return lst + [','.join(row[1].to_list()) for row in tmp_df.iterrows()]

def dict2list(tmp_dict, lst):
    """Returns a list with all values from given dictionary"""
    for key in tmp_dict:
        lst.extend(tmp_dict[key])
    return lst

###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-p", "--pfam_unmap", type=str, \
        help=" Un-mapped domain StIs from Pfam (passing or failing the threshold)")
    P.add_argument("-c", "--cath_unmap", type=str, \
        help=" Un-mapped domain StIs from Pfam (passing or failing the threshold)")
    P.add_argument("-px", "--pfam_xmap", type=str, \
        help="Cross-mapped domain StIs from Pfam only for failed StIs")
    P.add_argument("-cx", "--cath_xmap", type=str, \
        help="Cross-mapped domain StIs from CATH only for failed StIs")
    P.add_argument("-o", "--outfile", type=str, default='absd.json',\
        help="Output filename for either Domain-like or failed domain StIs")

    ARGS = P.parse_args()

    all_unmap = json.load(open(ARGS.outfile, 'r'))
    MERGE = []
    MERGE = df2list(pd.read_csv(ARGS.pfam_unmap), MERGE) if ARGS.pfam_unmap is not None else MERGE
    MERGE = df2list(pd.read_csv(ARGS.cath_unmap), MERGE) if ARGS.cath_unmap is not None else MERGE

    MERGE = dict2list(json.load(open(ARGS.pfam_xmap, "r")), MERGE) if ARGS.pfam_xmap is not None \
        else MERGE
    MERGE = dict2list(json.load(open(ARGS.cath_xmap, "r")), MERGE) if ARGS.cath_xmap is not None \
        else MERGE

    all_unmap[str(len(all_unmap))] = MERGE
    with open(ARGS.outfile, 'w') as f:
        json.dump(all_unmap, f, indent=2, sort_keys=True)
