#!/usr/bin/env python3
"""
The script compares the two set (CATH & Pfam) of structures with the given length threshold.
"""

import argparse
import json
import os
import pandas as pd

class ComparingStructures:
    """
    Class to compare the structures and return common structures along with unique to each set """
    def __init__(self, cath_df, pfam_df, true_doms, dom_length):
        self.cath = cath_df
        self.pfam = pfam_df
        self.true_domain = true_doms
        self.dom_size = dom_length

    def compare_cath_pfam(self):
        """
        This function takes two dataframes with the Uniprot numbering for all the structures from
        CATH and Pfam respectively. And one dictionary of true domain StIs to find out common and
        unique structures in both dataframes(CATH and Pfam)
        """
        # rearranging the columns of pfam dataframe
        helping_cols = ['UNP_id', 'UNP_start', 'UNP_end']
        self.pfam = self.pfam[[c for c in self.pfam if c not in helping_cols] + helping_cols]
        self.cath = self.cath[[c for c in self.cath if c not in helping_cols] + helping_cols]
        pfam_x = self.pfam
        index = self.pfam.index
        matched_pfam, unmatched_cath = [], []
        num_true = len(self.true_domain)

        for i in range(len(self.cath)):
            common = False
            one_cath = list(self.cath.loc[i])

            condition = (self.pfam['PDB_id'] == one_cath[0]) & \
                        (self.pfam['Chain_id'] == one_cath[1])
            positions = index[condition]
            for j in positions:
                one_pfam = list(self.pfam.loc[j])
                if one_pfam[0] == 'PDB_id':
                    continue

                if (abs(int(one_pfam[-2]) - int(one_cath[-2])) < self.dom_size) and \
                    (abs(int(one_pfam[-1]) - int(one_cath[-1])) < self.dom_size):
                    tmp_dict = {'Pfam': ','.join(one_pfam), 'CATH': ','.join(one_cath)}
                    self.true_domain[num_true] = tmp_dict
                    num_true += 1
                    matched_pfam.append(j)
                    common = True
                    break
                elif (int(one_pfam[-2]) > int(one_cath[-2])) and (
                        int(one_pfam[-1]) < int(one_cath[-1])):
                    tmp_dict = {'Pfam': ','.join(one_pfam), 'CATH': ','.join(one_cath)}
                    self.true_domain[num_true] = tmp_dict
                    num_true += 1
                    matched_pfam.append(j)
                    common = True
                    break
            if not common:
                unmatched_cath.append(','.join(one_cath))

        unmatched_pfam = []
        for pos in range(len(pfam_x)):
            if pos not in matched_pfam:
                unmatched_pfam.append(','.join(list(pfam_x.loc[pos])))

        return self.true_domain, unmatched_pfam, unmatched_cath

    def __str__(self):
        return "Compares two sets of domain structural instances from Pfam and CATH based on the" \
               " UniProt numbering."


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-p", "--Pfam", type=str, required=True, \
        help="File with residue-mapped Pfam domain StIs from given family(ies)")
    P.add_argument("-c", "--CATH", type=str, required=True, \
        help="File with residue-mapped CATH domain StIs from given family(ies)")
    P.add_argument("-f", "--truedom", default="true_domains.json", type=str, \
        help="Output filename for TRUE domain StIs, keep it same for all iterations")
    P.add_argument("-uq_pf", "--unique_pfam", default="unique_pfam.csv", type=str, \
        help="Output filename for domain StIs that are only in Pfam & not in CATH list")
    P.add_argument("-uq_ca", "--unique_cath", default="unique_cath.csv", type=str, \
        help="Output filename for domain StIs that are only in CATH & not in Pfam list")
    P.add_argument("-l", "--len_dom", type=int, default=31, \
        help="Minimum domain length criteria to filter/map the structural instances")
    ARGS = P.parse_args()

    true_domains = {}
    cath = pd.read_csv(ARGS.CATH, dtype=object)
    pfam = pd.read_csv(ARGS.Pfam, dtype=object)

    CMPR = ComparingStructures(cath, pfam, true_domains, ARGS.len_dom)
    true_domains, unique_pfam, unique_cath = CMPR.compare_cath_pfam()

    if not os.path.exists(ARGS.truedom.split('/')[-1]):
        common_domains = {}
    else:
        common_domains = json.load(open(ARGS.truedom.split('/')[-1], 'r'))

    iter_n = str(len(common_domains))
    common_domains[iter_n] = true_domains

    # Writing all the results to separate files
    with open(ARGS.truedom, "w") as outfile:
        json.dump(common_domains, outfile, indent=2)

    with open(ARGS.unique_pfam, "w") as outfile:
        outfile.write('\n'.join(unique_pfam))

    with open(ARGS.unique_cath, "w") as outfile:
        outfile.write('\n'.join(unique_cath))
