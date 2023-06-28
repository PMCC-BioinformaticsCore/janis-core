#!/usr/bin/env python3
"""This script  maps the residue numbering from UniPort to PDB for all the Pfam structural
instances using SIFTS mapping."""

import argparse
import pandas as pd
from residue_mapping import ResidueMapper

class Pfam2pdb:
    """
    Wrapper around the ResidueMapper from residue_mapping to map all instances from a file"""
    def __init__(self, in_file, sifts_path, lost):
        self.filename = in_file
        self.sifts = sifts_path
        self.lost_file = lost

    def mapping_pfam2pdb(self):
        """
        Main function mapping the residue numbering from a dataframe.
        """
        pfam_rrm = pd.read_csv(self.filename, dtype=object)
        lost_pfam = []

        pfam_colnames = ['PDB_id', 'Chain_id', 'Domain', 'Family_id', 'UNP_id', 'UNP_start', \
            'UNP_end', 'PDB_start', 'PDB_end']
        result_df = pd.DataFrame(columns=pfam_colnames)

        for i in range(len(pfam_rrm)):
            entry = list(pfam_rrm.loc[i])
            tmp = entry[-1]
            entry = entry[:-1] + tmp.split('-')

            M = ResidueMapper(entry[4], [int(entry[-2]), int(entry[-1])], path=self.sifts)
            if (not hasattr(M, 'mapped')) or (not hasattr(M, 'mapped')):
                lost_pfam.append(','.join(entry))
                continue

            tot = M.resmapper_unp2pdb(pdb=entry[0], chain=entry[1])
            all_good = True
            x_factor = [None, 'null']
            while(len(tot) < 1) or (tot[0][-2] in x_factor or tot[1][-2] in x_factor):
                # loop for calling the same function within
                if M.res_pos[0] >= M.res_pos[1]:
                    all_good = False
                    break
                if (len(tot) < 1) or (tot[0][-2] in x_factor and tot[1][-2] in x_factor):
                    M.res_pos = [M.res_pos[0]+1, M.res_pos[1]-1]
                    try:
                        tot = M.resmapper_unp2pdb(pdb=entry[0], chain=entry[1])
                    except:
                        continue
                elif (len(tot[0]) < 1) or (tot[0][-2] in x_factor):
                    M.res_pos[0] += 1
                    try:
                        tot = M.resmapper_unp2pdb(pdb=entry[0], chain=entry[1])
                    except:
                        continue
                elif (len(tot[1]) < 1) or (tot[1][-2] in x_factor):
                    M.res_pos[1] -= 1
                    try:
                        tot = M.resmapper_unp2pdb(pdb=entry[0], chain=entry[1])
                    except:
                        continue

            if not all_good:
                lost_pfam.append(','.join(entry))
                continue

            entry.extend([tot[0][-2], tot[1][-2]])
            result_df.loc[len(result_df)] = entry

        with open(self.lost_file, 'w') as file:
            file.write('\n'.join(lost_pfam))
        return result_df

    def __str__(self):
        return "Mapping of UniProt Residues to corresponding PDBs from {0}".format(self.filename)


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--filename", type=str, required=True, \
        help="File with all the filtered Pfam domain StIs")
    P.add_argument("-s", "--sifts", type=str, required=True, \
        help="Directory where you have/want to store SIFTS files")
    P.add_argument("-m", "--resmapped", type=str, default="pfam_resMapped.csv", \
        help="Output filename for residue-mapped Pfam domain StIs")
    P.add_argument("-l", "--reslost", type=str, default="lost_pfam.txt", \
        help="Output filename for obsolete or inconsistent Pfam domain StIs")
    ARGS = P.parse_args()

    PF = Pfam2pdb(ARGS.filename, ARGS.sifts, ARGS.reslost)
    OUT_DF = PF.mapping_pfam2pdb()
    OUT_DF.to_csv(ARGS.resmapped, index=False)
