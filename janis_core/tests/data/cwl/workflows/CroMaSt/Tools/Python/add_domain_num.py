#!/usr/bin/env python3
"""This script  collects multiple files (.csv) from subworkflow after residue mapping and
adds domain position labels to each structural instance."""

import argparse
import pandas as pd

class AddDomainNum:
    """
    This class contains functions for collects multiple files (.csv), merging them and
    adding domain number to each structural instance.
    """
    def __init__(self, inputs):
        self.files = inputs

    def combine_files(self):
        """Combines all the input (.csv) files into a dataframe"""
        whl_df = ''
        for afile in self.files:
            if isinstance(whl_df, str):
                whl_df = pd.read_csv(afile, index_col=False)
                continue
            tmp_df = pd.read_csv(afile, index_col=False)
            whl_df = pd.concat([whl_df, tmp_df], ignore_index=True)
        return whl_df

    @staticmethod
    def domain_posi_per_unp(mult_query, one_target, compare_col):
        """
        This function takes all the domain structural entries corresponding to one UniProt ID;
        and assignes the domain position/number to each of the entry
        """
        prev_up = None
        tmp_row = []
        for index, row in mult_query.iterrows():
            if index == 0:
                prev_up = row
                dom_position = 0

            if (prev_up[compare_col[2]] < row[compare_col[1]]) or \
                    (abs(prev_up[compare_col[2]] - row[compare_col[1]]) < 20): #ThresholdForNxt_dom
                dom_position += 1

            tmp_row = row.to_list() + ['domain_{0}'.format(dom_position + 1)]
            one_target.loc[len(one_target)] = tmp_row
            prev_up = row
        return one_target

    @staticmethod
    def wrapper_domain_pos(final_df):
        """wrapper around the function that adds domain position for each structural instance"""
        sort_col = ['UNP_id', 'UNP_start', 'UNP_end']
        for num_db in sort_col[1:] + ['PDB_start', 'PDB_end']:
            final_df[num_db] = final_df[num_db].astype(int)
        res_map = final_df.sort_values(sort_col, ascending=True).reset_index(drop=True)

        target = pd.DataFrame(columns=list(res_map.columns) + ['Domain_pos'])
        done_unp = []
        for row in res_map.iterrows():
            if row[1][sort_col[0]] in done_unp:
                continue
            done_unp.append(row[1][sort_col[0]])
            query = res_map[res_map[sort_col[0]] == row[1][sort_col[0]]].reset_index(drop=True)
            target = D.domain_posi_per_unp(query, target, sort_col)

        target = target[[c for c in target if c not in sort_col] + sort_col]
        return target


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-i", "--infiles", type=str, nargs="+", required=True, \
        help="List of files with residue mapped domain structural instances")
    P.add_argument("-o", "--outfile", type=str, default="resmapped_dom_posi.csv", \
        help="Output filename for domain StIs with the domain position labels")
    ARGS = P.parse_args()

    D = AddDomainNum(ARGS.infiles)
    COMBINED = D.combine_files()
    DOM_POS = D.wrapper_domain_pos(COMBINED)
    DOM_POS.to_csv(ARGS.outfile, index=False)
