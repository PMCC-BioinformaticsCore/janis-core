#!/usr/bin/env python3
"""
Script to cross-map structural instances from Pfam to CATH.
"""

import argparse
import json

class StructsPfam2CATH:
    """
    Class to cross-map Pfam structural instances to CATH
    """
    def __init__(self, pfam, inp, dom_length):
        self.pfam = pfam
        self.inp = inp
        self.dom_size = dom_length

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
        # info = cath_map_num(info)
        return info

    def map_pfam_struct2cath(self, pfam_struct):
        """
        Map the PDB structure from Pfam to associated CATH superfamily
        """
        all_dom = open(self.inp).read().split('//')
        mapping = {}
        pfam_struct = pfam_struct.split(',')

        for domain in all_dom:
            if any(pfam_struct[0].lower()+pfam_struct[1] in fam for fam in domain.split('\n')):
                one = self.filter_cath_helper(domain)
                if (abs(int(pfam_struct[4]) - one[2]) < self.dom_size) and (abs(int(pfam_struct[5]) - \
                    one[3]) < self.dom_size):
                    mapping[one[1]] = ','.join(pfam_struct)
                    return mapping
                # elif (int(pfam_struct[-4]) > one[2]) and (int(pfam_struct[-3]) < one[3]):
                #     mapping[one[1]] = ','.join(pfam_struct)
                #     return mapping

        if len(mapping) < 1:
            mapping['unmapped'] = ','.join(pfam_struct)
        return mapping


    def wrapper_struct_pfam2cath(self):
        """
        wrapper for mapping pfam structures to CATH superfamilies
        """
        pfam_map_cath = {}
        for one in self.pfam:
            if one == '':
                continue
            mapped_inst = self.map_pfam_struct2cath(one) #filename
            tmp_ls = list(mapped_inst.values())[0]

            if list(mapped_inst.keys())[0] in pfam_map_cath:
                pfam_map_cath[list(mapped_inst.keys())[0]].append(tmp_ls)
            else:
                pfam_map_cath[list(mapped_inst.keys())[0]] = [tmp_ls]

        return pfam_map_cath


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-p", "--Pfam", type=str, \
        help="File with CATH domain StIs to cross-map against whole Pfam db")
    P.add_argument("-c", "--cathraw", default="Data/cath-domain-description-file.txt", \
        type=str, help="Raw file from Pfam with all the domain information")
    P.add_argument("-x", "--name_out", type=str, default="pfam_crossMapped_cath.jsonx", \
        help="Output filename to store all cross-mapped domain StIs")
    P.add_argument("-u", "--unmapped", type=str, default="pfam_unq_unmapped.jsonx", \
        help="Output filename to store all un-mapped domain StIs")
    P.add_argument("-l", "--len_dom", type=int, default=31, \
        help="Minimum domain length or size criteria to filter/map the domain StIs")
    ARGS = P.parse_args()

    unique_pfam = open(ARGS.Pfam).read().split('\n')
    PF2CA = StructsPfam2CATH(unique_pfam, ARGS.cathraw, ARGS.len_dom)
    mapped_pfam = PF2CA.wrapper_struct_pfam2cath()

    if 'unmapped' in mapped_pfam:
        with open(ARGS.unmapped, "w") as f:
            json.dump({'unmapped': mapped_pfam['unmapped']}, f, indent=2, sort_keys=True)
        del mapped_pfam['unmapped']

    for fam in mapped_pfam:
        new_dict = {fam: mapped_pfam[fam]}
        with open(fam + '.json', "w") as f:
            json.dump(new_dict, f, indent=2)

    all_mapped = [value for key, value in mapped_pfam.items() if key not in ['unmapped']]
    only_mapped = {'mapped': all_mapped}
    with open(ARGS.name_out, "w") as outfile:
        json.dump(mapped_pfam, outfile, indent=2, sort_keys=True)
