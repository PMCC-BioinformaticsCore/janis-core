#!/usr/bin/env python3
"""
This script filters the structural instances passing and failing the given threshold from chosen
Kpax score.
"""

import argparse
import json
import pandas as pd

class FilterStructInsta:
    """
    Class with functions to retrieve the structural instances for their respective averaged
    UniProt domains"""
    def __init__(self, aln_df, aln_score, aln_thresh):
        self.aln_score = aln_score
        self.aln_thresh = aln_thresh
        self.aln_columns = aln_df.columns.to_list()
        self.passed = aln_df[aln_df[self.aln_score] >= self.aln_thresh].reset_index(drop=True)
        self.failed = aln_df[aln_df[self.aln_score] < self.aln_thresh].reset_index(drop=True)
        self.res_cols = ['PDB_id', 'Chain_id', 'Domain', 'Family_id', 'PDB_start', 'PDB_end', \
            'Domain_pos', 'UNP_id', 'UNP_start', 'UNP_end']
        self.origin_struct = pd.DataFrame(columns=self.res_cols)


    def single_keydict_2dataframe(self, indict, query_key):
        """
        Funtion to convert a dictionary into a datagrame
        """
        for insta in indict[query_key]:
            self.origin_struct.loc[len(self.origin_struct)] = insta.split(',')
        return self.origin_struct

    def unp2struct_instances(self, aln_df, flt_df_orig):
        """
        Function to retrieve structural instances based on Avgeraged UNP domains from Kpax result.
        """
        final_df = pd.DataFrame(columns=self.res_cols)
        for row in aln_df.iterrows():
            query_unp = row[1][self.aln_columns[1]].split('_', 1)[1].rsplit('_', 2)[:-1]
            tmp_group = flt_df_orig[(flt_df_orig[self.res_cols[-4]] == query_unp[0]) &\
                (flt_df_orig[self.res_cols[-3]] == query_unp[1])]

            final_df = final_df.append(tmp_group, ignore_index=True)
        return final_df

    def wrapper_unp2struct(self):
        """
        Wrapper around unp2struct_instances; runs this function separately to retrieve the
        passed and failed structural instances for given alignment threshold
        """
        passed = self.unp2struct_instances(self.passed, self.origin_struct)
        failed = self.unp2struct_instances(self.failed, self.origin_struct)
        return passed, failed

    def __str__(self):
        return "Retrieves the structural instances belonging to the averaged domains that passes"\
            " and fails the given threshold."


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-i", "--alnres", type=str, help="Result from Kpax pairwise alignment")
    P.add_argument("-p", "--passed", type=str, default='passed_structs.csv', \
        help="Output filename for the domain StIs passing the given threshold")
    P.add_argument("-f", "--failed", type=str, default='failed_structs.csv', \
        help="Output filename for the domain StIs failed to pass the given threshold")
    P.add_argument("-x", "--inifile", type=str, \
        help="File having the whole list domain StIs before alignment")
    P.add_argument("-s", "--score", type=str, choices=['Kscore', 'Gscore', 'Jscore', 'Mscore', \
        'Tscore', 'Nanchor', 'RMSD-aligned', 'RMSD-matched', 'RMSD-anchor', 'Naligned'], \
            default="Mscore", help="Alignment score to chose from for the evaluation")
    P.add_argument("-t", "--threshold", type=float, default=0.6, \
        help="The threshold value for the chosen score")

    ARGS = P.parse_args()
    score_type = '#' + ARGS.score

    # Reading alignment results as dataframe and all structural instances as dict
    all_structs = pd.read_csv(ARGS.alnres)
    whl_dict = json.load(open(ARGS.inifile))

    FSI = FilterStructInsta(all_structs, score_type, ARGS.threshold)
    all_unmapped = FSI.single_keydict_2dataframe(whl_dict, list(whl_dict.keys())[0])
    mapped_struct, failed_struct = FSI.wrapper_unp2struct()
    mapped_struct.to_csv(ARGS.passed, index=False)
    failed_struct.to_csv(ARGS.failed, index=False)
