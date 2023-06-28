#!/usr/bin/env python3
"""
This script aligns all the chooped structures using Kpax. Then the aligned structures are used
to compute the average structure. This script requires a json file where key is family id with
and value is list of domain structures with residue numbering.
"""

import argparse
import warnings
import os
import sys
import Bio.PDB as bpdb
import pandas as pd
warnings.filterwarnings("ignore")

class AvgStructComputation:
    """Class to compute average structure from given n number of 3-D structures,
    without the need to have same number of atoms or residues"""

    def __init__(self, fasta, a_dir, out):
        self.fasta_file = fasta
        self.struct_dir = a_dir
        self.out_file = out
        self.backbone = ['N', 'CA', 'C', 'O']


    def compute_avg_struct(self):
        """
        Computes average structure from a given set of structures
        """
        dict3_1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', \
                   'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', \
                   'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', \
                   'TYR': 'Y', 'MET': 'M', 'MSE': 'M'}
        dict1_3 = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE', \
                   'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 'G': 'GLY', 'H': 'HIS', \
                   'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V': 'VAL', 'E': 'GLU', \
                   'Y': 'TYR', 'M': 'MET'}

        seq = self.fasta_reader(self.fasta_file)  # Read fasta file

        # setting of the position to be considered based on percentage of gaps;
        # value ranges from 0 to 1
        fraction = 0.5
        for col in range(1, len(seq.loc[0])):
            pdb_names = {}
            list_x = list(seq.loc[0])
            for pos in range(len(seq)):
                # getting residue position in PDB
                row_m = list(seq.loc[pos])[:col + 1]

                if row_m[-1] == '-':
                    continue
                if '-' in row_m:
                    row_m = [val for val in row_m if val != '-']

                pdb_names[row_m[0]] = {len(row_m) - 1: row_m[-1]}
            coor_res = self.get_coord_mean(pdb_names, self.struct_dir)

            if self.most_frequent(list(seq[col])) == '-':
                if list(seq[col]).count('-') > fraction * len(list(seq[col])):
                    continue
                else:
                    res_x = list(seq[col])
                    res_x = [val for val in res_x if val != '-']
                    res_x = dict1_3[self.most_frequent(res_x)]
            else:
                res_x = dict1_3[self.most_frequent(list(seq[col]))]

            with open(self.out_file, 'a') as wrtav:  # Writing the average structure
                line = self.get_pdb_line(col, res_x, coor_res)
                wrtav.write('\n' + '\n'.join(line))
        return self.out_file

    @staticmethod
    def fasta_reader(filename):
        """
        Reads the fasta file and returns all sequneces in a dataframe with first column as headers
        """
        all_seq = open(filename).read().split('>')
        all_seq = [nn.replace('\n', '') for nn in all_seq[1:]]
        all_seq = [nn.split() for nn in all_seq]
        with open(filename.replace('.fasta', '.tmp'), 'w') as wrtfs:
            for pos in all_seq:
                wrtfs.write('\n')
                wrtfs.write(pos[0] + ',' + ','.join(list(pos[1])[1:]))

        align_df = pd.read_csv(filename.replace('.fasta', '.tmp'), header=None)
        return align_df

    def get_coord_mean(self, pdb_dict, path):
        """
        Computes average position for all backbone atoms (N, CA, C, O) & returns average positions
        """
        bb_coor = {}
        parsepdb = bpdb.PDBParser()

        for pdb in pdb_dict:
            # kpax retunrs .pdb extension with name only if only 2 structures are passed to align
            pdbname = pdb if pdb.endswith('.pdb') else pdb + '.pdb'
            # reading each structure
            structure = parsepdb.get_structure("struct", path + pdbname)
            for chain in structure[0]:
                num_n = 1
                for residue in chain:
                    if list(pdb_dict[pdb].keys())[0] == num_n:
                        x_res = residue
                        break
                    num_n += 1

                for atom in x_res:
                    if atom.get_name() in self.backbone:
                        if atom.get_name() in bb_coor.keys():
                            bb_coor[atom.get_name()].append(list(atom.get_vector()))
                        else:
                            bb_coor[atom.get_name()] = [list(atom.get_vector())]

        avg_coor = {}
        for atom in bb_coor:
            atom_coor = bb_coor[atom]
            avg_coor[atom] = [sum(nn)/len(nn) for nn in zip(*atom_coor)]
        return avg_coor

    @staticmethod
    def most_frequent(query_list):
        """Returns most frequent value from query list"""
        return max(set(query_list), key=query_list.count)


    def get_pdb_line(self, n, res, avg_coor):
        """Create the lines to be written for Average PDB structure """
        lines = []
        for x in self.backbone:
            if x not in avg_coor.keys():
                continue
            num = (n-1)*4
            record = 'ATOM' + ' '*(7 - len(str(num + 1 + self.backbone.index(x))))
            record += str(num + 1 + self.backbone.index(x)) + ' '*2 + x
            record += ' '*(17 - len(record)) + res + ' ' + 'A'
            record += ' '*(4 - len(str(n))) + str(n) + ' '*(12 - len("{:.3f}".format(avg_coor[x][0])))
            record += "{:.3f}".format(avg_coor[x][0]) + ' '*(8 - len("{:.3f}".format(avg_coor[x][1])))
            record += "{:.3f}".format(avg_coor[x][1]) + ' '*(8 - len("{:.3f}".format(avg_coor[x][2])))
            record += "{:.3f}".format(avg_coor[x][2])

            lines.append(record)
        return lines


class UtilityDir:
    """Class to provide some utilities for PDB files"""
    def __init__(self, a_dir):
        self.dir = a_dir

    def recursively_remove_files(self):
        """USE WITH CAUTION:
        Removes the directory even if it's not empty"""
        if os.path.isfile(self.dir):
            os.remove(self.dir)
        elif os.path.isdir(self.dir):
            for sub in os.listdir(self.dir):
                self.recursively_remove_files(os.path.join(self.dir, sub))
            os.rmdir(self.dir)

    def empty_dir(self):
        """Remove all files from a directory"""
        for afile in os.listdir(self.dir):
            try:
                os.remove(self.dir + afile)
            except:
                print('Unable to delete ', afile)


class Kpax:
    """Class for aligning structurs with Kpax; Currently only supports multiple alignment"""
    def __init__(self, a_dir):
        self.dir = a_dir

    def kpax_multi(self):
        """
        Multiple alignment using Kpax
        """
        log = 'kpax_log1.tmp'
        if (len(os.listdir(self.dir)) < 3) and (len(os.listdir(self.dir)) > 1):
            # For only 2 PDB files
            cmd = 'kpax5.1.3.x64 -broken -nowrite -pdb -fasta ' + \
                  os.path.join(self.dir, '*.pdb') + ' > ' + log
        else:
            # For more than two 2 PDB files
            cmd = 'kpax5.1.3.x64 -multi -broken -nowrite -pdb -fasta ' + \
                  os.path.join(self.dir, '*.pdb') + ' > ' + log
        os.system(cmd)
        resdir, fasta = '', ''

        for line in open(log, 'r'):
            if line.startswith('Creating RESULTS directory'):
                resdir = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line.startswith('Writing FASTA multiple alignment file'):
                fasta = line.split(':')[1].replace(' ', '').replace('\n', '')
            if (len(os.listdir(self.dir)) < 3) and line.startswith('Writing TARGET pdb file'):
                fasta = line.split(':')[1].replace(' ', '').replace('\n', '')
                fasta = fasta.replace('.pdb', '.fasta')

            if (resdir != '') and (fasta != ''):
                break
        os.remove(log)
        return resdir + '/', fasta

    def __str__(self):
        return "Aligning all the structures from directory: {0}".format(self.dir)


###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-f", "--file", type=str, required=True, \
        help="The file with instance name to be assigned to average structure")
    P.add_argument("-s", "--splitdir", type=str, default="split_PDB/", \
        help="The directory path for the PDB structures to be averaged")
    P.add_argument("-k", "--kpax_res", type=str, default="KPAX_RESULTS/", \
        help="The directory path where all the kpax alignment results will be stored")

    ARGS = P.parse_args()

    a = open(ARGS.file).read().replace('\n', '')

    # Default directory for KPAX_RESULTS
    if 'KPAX_RESULTS' not in os.environ:
        os.environ['KPAX_RESULTS'] = ARGS.kpax_res
    os.environ['KPAX_RESULTS'] = os.environ['KPAX_RESULTS'] if os.environ['KPAX_RESULTS'].endswith('/') else os.environ['KPAX_RESULTS'] + '/'

    if os.path.exists(ARGS.kpax_res):
        util_d = UtilityDir(ARGS.kpax_res)
        util_d.recursively_remove_files()

    if len(os.listdir(ARGS.splitdir)) > 1:
        kp = Kpax(ARGS.splitdir)
        kpax_dr, fastafile = kp.kpax_multi()
    else:
        if '.' in str(a):       # The change is bcz of kpax; Sometimes it replaces . with ~
            a = str(a).replace('.', '-')            # To avoid the change in name
        pdbout = str(a) + '_avgStruct.pdb'              # Name of the average PDB structure file

        if len(os.listdir(ARGS.splitdir)) < 1:
            print("\n\nNo PDB files in following Dir", ARGS.splitdir, ARGS.file, a)
            sys.exit()
        os.system('mv '+ os.path.join(ARGS.splitdir, os.listdir(ARGS.splitdir)[0]) + ' '+ pdbout)
        sys.exit()

    if '.' in str(a):                # The change is bcz of kpax; Sometimes it replaces . with ~
        a = str(a).replace('.', '-')            # To avoid the change in name

    pdbout = str(a) + '_avgStruct.pdb'              # Name of the average PDB structure file

    print(a, len(os.listdir(ARGS.splitdir)), pdbout, "IF error happens")

    avg = AvgStructComputation(fastafile, kpax_dr, pdbout)
    file = avg.compute_avg_struct()
