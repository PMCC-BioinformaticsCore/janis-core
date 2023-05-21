#!/usr/bin/env python3
"""
This script  maps the residues from PDB to UniPort numbering for all the CATH structural instances
using SIFTS mapping.
"""

import argparse
import pandas as pd
from residue_mapping import ResidueMapper

class ResMapping:
    """
    Wrapper around the ResidueMapper from residue_mapping to map all instances from a file
    """
    def __init__(self, in_file, sifts_path, lost):
        self.filename = in_file
        self.sifts = sifts_path
        self.lost_file = lost

    def mapping_cath2up(self):
        """
        This function maps the PDB numbering to UniProt numbering for all the structures from CATH
        """
        cath_rrm = pd.read_csv(self.filename, dtype=object)
        lost_cath = []

        cath_colname = ['PDB_id', 'Chain_id', 'Domain', 'Family_id', 'PDB_start', 'PDB_end', \
            'UNP_id', 'UNP_start', 'UNP_end']
        result_df = pd.DataFrame(columns=cath_colname)

        for i in range(len(cath_rrm)):
            entry = list(cath_rrm.loc[i])
            entry = [str(x) for x in entry]

            M = ResidueMapper(entry[0].lower(), [int(entry[-2]), int(entry[-1])], self.sifts)
            if (not hasattr(M, 'mapped')) or (not hasattr(M, 'mapped')):
                lost_cath.append(','.join(entry))
                continue

            tot = M.resmapper_pdb2unp(chain=entry[1])
            all_good = True
            x_factor = [None, 'null']

            while len(tot) < 1 or tot[0][2] in x_factor or tot[1][2] in x_factor:
                tmp_tot = tot
                # loop for calling the same function within
                if M.res_pos[0] >= M.res_pos[1]:
                    all_good = False
                    break
                if (len(tot) < 1) or (tot[0][2] in x_factor and tot[1][2] in x_factor):
                    M.res_pos = [M.res_pos[0]+1, M.res_pos[1]-1]
                    try:
                        tot = M.resmapper_pdb2unp(chain=entry[1])
                    except:
                        continue
                elif (len(tot[0]) < 1) or (tot[0][2] in x_factor):
                    M.res_pos[0] += 1
                    try:
                        tot = M.resmapper_pdb2unp(chain=entry[1])
                    except:
                        continue
                elif (len(tot[1]) < 1) or (tot[1][2] in x_factor):
                    M.res_pos[1] -= 1
                    try:
                        tot = M.resmapper_pdb2unp(chain=entry[1])
                    except:
                        continue
                tot = tot if len(tot) == 2 else tmp_tot

            if (not all_good) or (len(tot) < 2) or (tot[0][1] != tot[1][1]):
                lost_cath.append(','.join(entry))
                continue

            entry.extend([tot[0][1], tot[0][2], tot[1][2]])
            result_df.loc[len(result_df)] = entry

        with open(self.lost_file, 'w') as wrt:
            wrt.write('\n'.join(lost_cath))
        return result_df

    def __str__(self) -> str:
        return "Mapping of UniProt Residues to corresponding PDBs from {0}".format(self.filename)



###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--cathfile", type=str, required=True, \
        help="File with all the filtered Pfam domain StIs")
    P.add_argument("-s", "--sifts", type=str, required=True, \
        help="Directory where you have/want to store SIFTS files")
    P.add_argument("-m", "--resmapped", type=str, default="cath_resMapped-X.csv", \
        help="Output filename for residue-mapped Pfam domain StIs")
    P.add_argument("-l", "--reslost", type=str, default="lost_cath-X.txt", \
        help="Output filename for obsolete or inconsistent Pfam domain StIs")
    ARGS = P.parse_args()

    CATH = ResMapping(ARGS.cathfile, ARGS.sifts, ARGS.reslost)
    OUT_DF = CATH.mapping_cath2up()
    OUT_DF.to_csv(ARGS.resmapped, index=False)
