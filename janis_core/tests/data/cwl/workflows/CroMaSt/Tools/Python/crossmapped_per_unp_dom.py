#!/usr/bin/env python3
"""
This script takes a json file per crossmapped family and divide a json per domain instance.
"""
import argparse
import json

class UnpDomain:
    """
    Class divides a json with structs per family to structs per domain instance per family.
    """
    def __init__(self, in_file):
        self.fam_file = in_file

    def structs_unp_dom(self):
        """
        Divides structures per doamin instance
        """
        all_data = json.load(open(self.fam_file, 'r'))
        done_dom = {}
        fam_in = []
        for fam in all_data:
            for struct in all_data[fam]:
                struct = struct.split(',')
                if fam + '_' + '_'.join(struct[6:8]) in done_dom:
                    done_dom[fam + '_' + '_'.join(struct[6:8])].append(','.join(struct))
                else:
                    done_dom[fam + '_' + '_'.join(struct[6:8])] = [','.join(struct)]
                fam_in.append(fam)
        return done_dom, fam_in

    def __str__(self):
        return "Dividing all structures from {0} according to their domain instances"\
            "".format(self.fam_file)


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--file", type=str, required=True, \
        help="Cross-mapped domain StIs per family")
    ARGS = P.parse_args()

    DOM = UnpDomain(ARGS.file)
    domain_structs, fam_names = DOM.structs_unp_dom()
    for domain_unp in domain_structs:
        new_dict = {domain_unp: domain_structs[domain_unp]}
        with open(domain_unp + '.json', "w") as f:
            json.dump(new_dict, f, indent=2)
    print(fam_names[0])
