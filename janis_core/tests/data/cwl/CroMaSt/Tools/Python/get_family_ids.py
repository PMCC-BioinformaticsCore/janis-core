#!/usr/bin/env python3
"""
This script takes the input arguments for family ids from Pfam and CATH with the no of iteration
and returns a file nby adding this information."""
import argparse
import json

class FamilyIds:
    """
    Reads family ids provided by user and writes it to out file per iteration"""
    def __init__(self, file):
        self.file = file

    def write_ids(self, pfam, cath, num):
        """
        Write family IDs provided by user in a json file per iteration.
        """
        filename = self.file.rsplit('/', 1)[-1]
        data = json.load(open(filename, 'r'))
        assert str(num) not in data, """You are trying to re-run the same iteration. Are you sure?
            If you want to proceed please remove the key for the iteration from file {0} 
            and try again.""".format(self.file)

        data[str(num)] = {'Pfam': pfam, 'CATH': cath}
        with open(self.file, 'w') as out:
            json.dump(data, out, indent=2, sort_keys=True)

        return filename


    def read_ids(self, src_db):
        """
        Reads json file and returns already processed and to be processed family ids
        """
        data = json.load(open(self.file, 'r'))
        fam_done = []
        for iter_n in range(len(data)-2, -1, -1):
            fam_done.extend(data[str(iter_n)][src_db])

        fam = data[str(len(data) - 1)][src_db]

        target = [*set(fam)]
        for one in [*set(fam)]:
            if one in fam_done:
                target.remove(one)
            else:
                fam_done.append(one)
        return fam_done, target


############################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Writes the user-given family Ids to a file')
    P.add_argument("-p", "--pfam", type=str, nargs="+", default=[], \
        help="Pfam ID(s) to start the current iteration with")
    P.add_argument("-c", "--cath", type=str, nargs="+", default=[], \
        help="CATH ID(s) to start the current iteration with")
    P.add_argument("-f", "--file", type=str, required=True, \
        help="Output filename to store all family IDs per iteration")
    P.add_argument("-n", "--number", type=int, required=True, help="Current iteration number")
    ARGS = P.parse_args()

    WRT_ID = FamilyIds(ARGS.file)
    WRT_ID.write_ids(ARGS.pfam, ARGS.cath, ARGS.number)
