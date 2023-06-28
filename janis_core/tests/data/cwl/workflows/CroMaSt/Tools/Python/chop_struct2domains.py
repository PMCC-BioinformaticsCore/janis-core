#!/usr/bin/env python3
"""
This script chops the PDB files into given domains according to the mapping provided in a different
directory and aligns all the chooped structures using Kpax. Then the aligned structures are used to
compute the average structure.
This script requires a json file where key is family id with and value is list of domain structures
with residue numbering.
"""

import argparse
import os
import json
import warnings
from sh import gunzip
import wget
import Bio.PDB as bpdb

warnings.filterwarnings("ignore")


class PDButil:
    """Class to provide some utilities for PDB files"""
    def __init__(self, pdb_code, dir_pdb):
        self.pdb_code = pdb_code.lower()
        self.dir_pdb = dir_pdb

    def get_pdb_file(self):
        """
        Download and unzips the PDB file from rcsb in given directory
        """
        url = 'https://files.rcsb.org/download/' + self.pdb_code + '.pdb.gz'
        wget.download(url, out=self.dir_pdb)
        gunzip(os.path.join(self.dir_pdb, self.pdb_code + '.pdb.gz'))

    def get_cif_file(self):
        """
        Download and unzips the PDB file from rcsb in given directory
        """
        url = 'https://files.rcsb.org/download/' + self.pdb_code + '.cif'
        if self.pdb_code.lower() + '.cif' not in os.listdir(self.dir_pdb):
            wget.download(url, out=self.dir_pdb)

    def pdb_chopper(self, chain_id, start_res, end_res, chopped, stype):
        """
        Chops the PDB structure/file
        """
        class ResSelect(bpdb.Select):
            """Re-Using the class from Bio.PDB"""
            def accept_residue(self, res):
                """Actual chopping of the structure"""
                return bool((res.id[1] >= int(start_res)) and (res.id[1] <= int(end_res)) and \
                    (res.parent.id == chain_id))

        if stype == 'pdb':
            wholepdb = os.path.join(self.dir_pdb, self.pdb_code + '.pdb')
            stparse = bpdb.PDBParser().get_structure('temp', wholepdb)
        elif stype in ('cif', 'mmcif'):
            wholepdb = os.path.join(self.dir_pdb, self.pdb_code + '.cif')
            stparse = bpdb.MMCIFParser().get_structure('temp', wholepdb)

        if len(chain_id) > 1:
            stparse, chain_id = self.multi_letter_chains(stparse, chain_id)

        inop = bpdb.PDBIO()
        inop.set_structure(stparse)
        inop.save(chopped, ResSelect())
        return chopped

    @staticmethod
    def multi_letter_chains(structure, chain_id):
        """
        Assigns a single letter to the target chain, as only first letter of original id.
        """
        for chain in structure.get_chains():
            if chain.id != chain_id[0]:
                continue
            chain.id = chain_id[-1]

        for chain in structure.get_chains():
            if chain.id != chain_id:
                continue
            chain.id = chain_id[0]

            return structure, chain_id[0]


class UtilityDir:
    """Class to provide some utilities for PDB files"""
    def __init__(self, adir):
        self.dir = adir

    def recursively_remove_files(self, dirpath=None):
        """USE WITH CAUTION:
        Removes the directory even if it's not empty"""
        dirpath = self.dir if dirpath is None else dirpath
        if os.path.isfile(dirpath):
            os.remove(dirpath)
        elif os.path.isdir(dirpath):
            for sub in os.listdir(dirpath):
                self.recursively_remove_files(dirpath=os.path.join(self.dir, sub))
            os.rmdir(dirpath)

    def empty_dir(self):
        """Remove all files from a directory"""
        for afile in os.listdir(self.dir):
            try:
                os.remove(os.path.join(self.dir, afile))
            except:
                print('Unable to delete ', afile)


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,\
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--file", type=str, required=True, \
        help="File with list of StIs to chop into domains")
    P.add_argument("-p", "--pdbdir", type=str, required=True, \
        help="Directory where all the PDB files are/will be stored")
    P.add_argument("-s", "--splitdir", type=str, default="split_PDB/", \
        help="Directory where all the chopped structures will be stored")
    P.add_argument("-k", "--kpax_res", type=str, default="KPAX_RESULTS/", \
        help="Directory path to store results from  kpax alignments")

    ARGS = P.parse_args()

    file = ARGS.file
    util_d = UtilityDir(ARGS.splitdir)
    if os.path.exists(ARGS.splitdir):
        util_d.recursively_remove_files()

    for a in [ARGS.pdbdir, ARGS.splitdir]:
        if not os.path.exists(a):
            os.makedirs(a)

    no_pdb, avg_structures = [], []
    data = json.load(open(file, 'r'))

    for a in data:
        util_d.empty_dir() # removing files from previous iterations #split_dir

        # For all mapped structures
        for b in data[a]:
            # reading each structure from given json file
            struct = b.replace('"', '').split(',')
            pdb_id = struct[0].lower()
            pdb = PDButil(pdb_id, ARGS.pdbdir)

            out_file = os.path.join(ARGS.splitdir, '_'.join(struct[6:8]) + '_' + \
                '_'.join(struct[:2]) + '_' + '_'.join(struct[4:6]) + '.pdb')
            if pdb_id + '.pdb' not in os.listdir(ARGS.pdbdir):
                try:
                    pdb.get_pdb_file()
                except:
                    try:
                        pdb.get_cif_file()
                        filename = pdb.pdb_chopper(struct[1], struct[4], struct[5], out_file, 'cif')
                    except:
                        no_pdb.append(pdb_id)
                    continue

            filename = pdb.pdb_chopper(struct[1], struct[4], struct[5], out_file, 'pdb')
        print(a)
