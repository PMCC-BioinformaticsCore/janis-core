#!/usr/bin/env python3
"""
This script filters all the structural instances available for given Pfam families.
"""

import argparse
import pandas as pd
from get_family_ids import FamilyIds

class SeparateStructs:
    """
    Class used to separate PDB structures provided family ids (CATH & Pfam)
    """
    def __init__(self, fam_ids, fam_raw, obs_file, dom_length):
        self.obs_tot = open(obs_file).read()[1:-1].replace('"', '').replace(' ', '').split(',')
        self.src_file = fam_raw
        self.fam_ids = fam_ids
        self.obs_query = []
        self.dom_size = dom_length

    @staticmethod
    def write_structs(structs, struct_f, obs_fam, obs_f):
        """
        Writes separated and obsolete structures in separate files
        """
        structs.to_csv(struct_f, index=False)
        with open(obs_f, 'w') as wrt:
            wrt.write('\n'.join(obs_fam))

    def split_dataframe(self, df_whole, outname):
        """
        This function will split and writes the dataframe into multiple files for parallel
        processing in CWL.
        """
        splited_files = []
        for n in range(0, len(df_whole), 100):
            x = n + 100 if (n + 100) < len(df_whole) else len(df_whole)
            tmp_sep = df_whole.iloc[n:x]
            splited_files.append("{0}_{1}".format(len(splited_files), outname))
            tmp_sep.to_csv(splited_files[len(splited_files)-1], index=False)

        return splited_files

    def __str__(self):
        return "Separating structures for following families {0}".format(', '.join(self.fam_ids))


class FilterPfamStructs(SeparateStructs):
    """
    Inherited class used to filter PDB structures from Pfam
    """
    def filter_pfam(self):
        """
        Separate all PDB instances/structures corresponding to given pfam_ids from the pdbmap file
        """
        all_struct = open(self.src_file).read().split('\n')
        pfam_df = pd.DataFrame(columns=['PDB_id', 'Chain_id', 'Domain', \
            'Family_id', 'UNP_id', 'UNP_range'])

        for line in all_struct:
            if line == '':
                continue
            line = ''.join(line.split()).split(';')
            if line[3] not in self.fam_ids:
                continue

            if line[0].upper() in self.obs_tot:
                self.obs_query.append(','.join(line))
            else:
                res_range = line[5].split('-')
                if abs(int(res_range[1]) - int(res_range[0])) < self.dom_size:
                    continue
                pfam_df = pfam_df.append({'PDB_id': line[0], 'Chain_id': line[1], \
                    'Domain': line[2], 'Family_id': line[3], 'UNP_id': line[4], \
                        'UNP_range': line[5]}, ignore_index=True)

        return pfam_df, self.obs_query


############################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--file", type=str, required=True, \
        help="File with the family IDs per iteration")
    P.add_argument("-d", "--obs_pdb", type=str, default="Data/obsolete_PDB_entry_ids.txt", \
        help="File with all the obsolete PDB IDs")
    P.add_argument("-p", "--pfam_raw", type=str, default="Data/pdbmap", \
        help="Raw file from Pfam containing all the domain information")
    P.add_argument("-n", "--filtered_pfam", type=str, default="Filtered_Pfam.csv", \
        help="User-defined filename for filtered structures from Pfam")
    P.add_argument("-o", "--obsolete_pfam", type=str, default="obsolete_pfam.txt", \
        help="User-defined filename for obsolete Pfam structs from the given list of Pfam IDs")
    P.add_argument("-s", "--split_files", type=str, default="part.csv", \
        help="User-defined suffix for splitted files (.csv)")
    P.add_argument("-l", "--len_dom", type=int, default=31, \
        help="Minimum domain length criteria to filter/map the structural instances")
    ARGS = P.parse_args()

    # Pre-processing
    FAM = FamilyIds(ARGS.file)
    PF_X, PF = FAM.read_ids('Pfam')

    # Filtering PDB structures/instances from list of given IDs
    FLT = FilterPfamStructs(PF, ARGS.pfam_raw, ARGS.obs_pdb, ARGS.len_dom)
    FT_PF, OBS_PF = FLT.filter_pfam()
    FLT.write_structs(FT_PF, ARGS.filtered_pfam, OBS_PF, ARGS.obsolete_pfam)
    splitted_filnames = FLT.split_dataframe(FT_PF, ARGS.split_files)
