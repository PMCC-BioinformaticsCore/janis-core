#!/usr/bin/env python3
"""
Script to cross-map structural instances from CATH to Pfam.
"""

import argparse
import json

class StructsCATH2Pfam:
    """
    Class to cross-map CATH structural instances to Pfam
    """
    def __init__(self, cath, inp, dom_length):
        self.cath = cath
        self.inp = inp
        self.dom_size = dom_length

    def map_cath_struct2pfam(self, cath_struct):
        """
        Map the PDB structure from CATH to associated Pfam family
        """
        all_struct = open(self.inp).read().split('\n')
        mapping = {}

        for line in all_struct:
            if line == '':
                continue
            line = line.split()
            line = ''.join(line).split(';')

            if (line[:2] == cath_struct[:2]) and (line[4] == cath_struct[-3]):
                if (abs(int(line[5].split('-')[0]) - int(cath_struct[-2])) < self.dom_size) and \
                    (abs(int(line[5].split('-')[1]) - int(cath_struct[-1])) < self.dom_size):
                    mapping[line[3]] = ','.join(cath_struct)
                    return mapping
                elif (int(line[5].split('-')[0]) > int(cath_struct[-2])) and \
                    (int(line[5].split('-')[1]) < int(cath_struct[-1])) and \
                    (abs(int(line[5].split('-')[1]) - int(line[5].split('-')[0])) > self.dom_size):
                    mapping[line[3]] = ','.join(cath_struct)
                    return mapping

        if len(mapping) < 1:
            mapping['unmapped'] = ','.join(cath_struct)
        return mapping


    def wrapper_struct_cath2pfam(self):
        """
        wrapper for mapping CATH structures to Pfam family
        """
        cath_map_pfam = {}
        for one in self.cath:
            mapped_inst = self.map_cath_struct2pfam(one.split(','))   #filename
            tmp_ls = list(mapped_inst.values())[0]

            if list(mapped_inst.keys())[0] in cath_map_pfam:
                cath_map_pfam[list(mapped_inst.keys())[0]].append(tmp_ls)
            else:
                cath_map_pfam[list(mapped_inst.keys())[0]] = [tmp_ls]

        return cath_map_pfam


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-c", "--CATH", type=str, \
        help="File with CATH domain StIs to cross-map against whole Pfam db")
    P.add_argument("-p", "--pfamraw", type=str, default="Data/pdbmap", \
        help="Raw file from Pfam with all the domain information")
    P.add_argument("-x", "--name_out", type=str, default="cath_crossMapped_pfam.jsonx", \
        help="Output filename to store all cross-mapped domain StIs")
    P.add_argument("-u", "--unmapped", type=str, default="cath_unq_unmapped.jsonx", \
        help="Output filename to store all un-mapped domain StIs")
    P.add_argument("-l", "--len_dom", type=int, default=31, \
        help="Minimum domain length criteria to filter/map the domain StIs")
    ARGS = P.parse_args()

    unique_cath = open(ARGS.CATH).read().split('\n')
    CA2PF = StructsCATH2Pfam(unique_cath, ARGS.pfamraw, ARGS.len_dom)
    mapped_cath = CA2PF.wrapper_struct_cath2pfam()

    if 'unmapped' in mapped_cath:
        with open(ARGS.unmapped, "w") as f:
            json.dump({'unmapped':mapped_cath['unmapped']}, f, indent=2, sort_keys=True)
        del mapped_cath['unmapped']

    for fam in mapped_cath:
        new_dict = {fam: mapped_cath[fam]}
        with open(fam + '.json', "w") as f:
            json.dump(new_dict, f, indent=2)

    all_mapped = [value for key, value in mapped_cath.items() if key not in ['unmapped']]
    only_mapped = {'mapped': all_mapped}
    with open(ARGS.name_out, "w") as outfile:
        json.dump(mapped_cath, outfile, indent=2, sort_keys=True)
