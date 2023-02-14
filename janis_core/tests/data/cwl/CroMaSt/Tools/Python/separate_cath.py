#!/usr/bin/env python3
"""
This script filters all the structural instances available for given CATH superfamilies.
"""

import argparse
import pandas as pd
from get_family_ids import FamilyIds
from separate_pfam import SeparateStructs

class FilterCATHStructs(SeparateStructs):
    """
    Inherited class used to filter PDB structures from CATH
    """
    def filter_cath_with_num(self):
        """
        Separate all PDB instances/structures corresponding to given superfamil(y)ies
        from file cath-domain-list.txt
        """
        all_dom = open(self.src_file).read().split('//')
        ca_df = pd.DataFrame(columns=['PDB_id', 'Chain_id', 'Domain', 'Family_id', \
            'PDB_start', 'PDB_end'])

        for domain in all_dom:
            for sup in self.fam_ids:
                if any(sup in fam for fam in domain.split('\n')):
                    one = self.filter_cath_helper(domain)
                    if len(sup) != len(one[1]):
                        break
                    elif one[0][:4].upper() in self.obs_tot:
                        one = [str(z) for z in one]
                        self.obs_query.append(','.join(one))
                        break

                    if abs(int(one[3]) - int(one[2])) < self.dom_size:  # domain length threshold
                        continue
                    ca_df = ca_df.append({'PDB_id': one[0][:4].upper(), 'Chain_id': one[0][4:5],\
                        'Domain': one[0][5:], 'Family_id': one[1], 'PDB_start': one[2], \
                            'PDB_end': one[3]}, ignore_index=True)
                    break
        return ca_df, self.obs_query


    def filter_cath_helper(self, descr):
        """
        Returns necessary information from description of each domain
        """
        info = []
        for line in descr.split('\n'):
            if line.startswith('DOMAIN'):
                info.append(line.split()[1])
            elif line.startswith('CATHCODE'):
                info.append(line.split()[1])
            elif line.startswith('SRANGE'):
                line = line.split()
                info.append(int(line[1][6:]))
                info.append(int(line[2][5:]))

        if len(info) > 4:    # for multiple segments
            if abs(int(info[4]) - int(info[3])) < self.dom_size:
                del info[3:5]
            elif abs(int(info[3]) - int(info[2])) < self.dom_size:
                del info[2:4]

        return info



############################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--file", type=str, required=True, \
        help="File with the family IDs per iteration")
    P.add_argument("-d", "--obs_pdb", type=str, default="Data/obsolete_PDB_entry_ids.txt", \
        help="File with all the obsolete PDB IDs")
    P.add_argument("-c", "--cath_raw", default="Data/cath-domain-description-file.txt",\
        type=str, help="Raw file from CATH containing all the domain information")
    P.add_argument("-n", "--filtered_cath", type=str, default="Filtered_CATH.csv", \
        help="User-defined filename for filtered structures from CATH")
    P.add_argument("-o", "--obsolete_cath", type=str, default="obsolete_cath.txt", \
        help="User-defined filename for obsolete CATH structures from the given list of CATH families")
    P.add_argument("-s", "--split_files", type=str, default="part.csv", \
        help="User-defined suffix for splitted files (.csv)")
    P.add_argument("-l", "--len_dom", type=int, default=31, \
        help="Minimum domain length criteria to filter/map the structural instances")
    ARGS = P.parse_args()

    # Pre-processing
    FAMC = FamilyIds(ARGS.file)
    CATH_X, CATH = FAMC.read_ids('CATH')

    # Filtering PDB structures/instances from list of given IDs
    FLT = FilterCATHStructs(CATH, ARGS.cath_raw, ARGS.obs_pdb, ARGS.len_dom)
    FT_CA, OBS_CA = FLT.filter_cath_with_num()
    FLT.write_structs(FT_CA, ARGS.filtered_cath, OBS_CA, ARGS.obsolete_cath)
    splitted_filnames = FLT.split_dataframe(FT_CA, ARGS.split_files)
