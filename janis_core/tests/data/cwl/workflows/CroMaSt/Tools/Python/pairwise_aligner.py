#!/usr/bin/env python3
"""
This script takes several structures as input and align them against core_avgStruct.pdb.
Returns a table containing RMSD, number of aligned residues and M-score.
"""

import argparse
import os
import pandas as pd

class KpaxPair:
    """
    Pairwise alignment of structures from a directory
    """
    def __init__(self, core, listS, resdir):
        self.core = core
        self.list_structs = listS
        self.kpaxres = resdir

    def kpax_pairwise(self):
        """
        Pairwise alignment"""
        os.environ['KPAX_RESULTS'] = self.kpaxres
        if os.path.exists(self.kpaxres):
            self.recursively_remove_files(self.kpaxres)

        log = 'kpax_log1.tmp'
        cmd = 'kpax5.1.3.x64 -show=1500 -broken -nowrite -pdb -kalign -fasta -knames ' + \
              self.core + '  ' + self.list_structs + ' > ' + log
        os.system(cmd)
        resdir, fasta_file = '', ''

        for line in open(log, 'r'):
            if line.startswith('Creating RESULTS directory'):
                resdir = line.split(':')[1].replace(' ', '').replace('\n', '')
            elif line.startswith('Writing FASTA multiple alignment file'):
                fasta_file = line.split(':')[1].replace(' ', '').replace('\n', '')

            if (resdir != '') and (fasta_file != ''):
                break
        os.remove(log)
        return resdir + '/', fasta_file

    def recursively_remove_files(self, adir):
        """
        Recurs"""
        if os.path.isfile(adir):
            os.remove(adir)
        elif os.path.isdir(adir):
            for sub in os.listdir(adir):
                self.recursively_remove_files(os.path.join(adir, sub))
            os.rmdir(adir)


    #kalign file reader in loop for all the structures
    @staticmethod
    def processing_kalign(query_name, path, outfile):
        """
        query is name of core_avg_struct and path is where all result files are
        """
        kalign = open(path + query_name.replace('.pdb', '.knames')).read().split('\n')

        characteristics = ['#Query In', '#Target In', '#Query chain', '#Target chain', \
        '#Query Out', '#Target Out', '#Kscore', '#Gscore', '#Jscore', '#Mscore', \
        '#Tscore', '#Nquery', '#Ntarget', '#Nsegments', '#Ncolumns', '#Ncover', \
        '#Naligned', '#Nmatched', '#Nanchor', '#Nidentity', '#RMSD-aligned', \
        '#RMSD-matched', '#RMSD-anchor', '#MaxDist', '#Contact']

        kaln_df = pd.DataFrame(columns=characteristics)
        for struct in kalign:
            if struct == '':
                continue
            one = {}
            kaln_file = path + struct + '_' + query_name.replace('.pdb', '.kalign')

            data = open(kaln_file).read().split('\n')                # reading every kalign file
            for line in data:                                          # Scanning every .kalign file
                if not line.startswith('#'):
                    continue

                aln_scores = [nn for nn in characteristics if line.startswith(nn)]
                if len(aln_scores) > 0:
                    line = line.split(':')
                    one[line[0]] = line[1].replace(' ', '')

            # create series from dictionary and append it to dataframe as row
            one = pd.Series(one, index=characteristics)
            kaln_df = kaln_df.append(one, ignore_index=True)

        kaln_df.to_csv(outfile, index=False)                 # Write dataframe to csv file
        return outfile



###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-c", "--queryCath", type=str, nargs="+", \
        help="PDB structures to align as query from CATH")
    P.add_argument("-p", "--queryPfam", type=str, nargs="+", \
        help="PDB structures to align as query from Pfam")
    P.add_argument("-t", "--target", type=str, required=True, \
        help="The Core average structure")
    P.add_argument("-d", "--directory", type=str, \
        help="Directory with PDB files to align with target structure")
    P.add_argument("-r", "--result_file", type=str, default="align_Struct_analysis.csv", \
        help="Output filename for the alignment results")
    ARGS = P.parse_args()

    queryCath = ARGS.queryCath if ARGS.queryCath else []
    queryPfam = ARGS.queryPfam if ARGS.queryPfam else []

    if isinstance(ARGS.target, list):
        ARGS.target = ARGS.target[0]

    query = ''
    if (len(queryCath) == len(queryPfam) == 0) and (ARGS.directory is None):
        query = ARGS.target + ' ' + ARGS.target
    elif len(queryCath) == 0:
        query = ' '.join(queryPfam)
    elif len(queryPfam) == 0:
        query = ' '.join(queryCath)
    else:
        query = ' '.join(queryCath) + ' ' + ' '.join(queryPfam)

    if ARGS.directory is not None:
        queryDir = ARGS.directory if ARGS.directory.endswith('/') else ARGS.directory + '/'
        query += ' ' + queryDir + '*.pdb'

    # Default directory for KPAX_RESULTS
    kpaxRes = 'KPAX_RESULTS/'
    kp = KpaxPair(ARGS.target, query, kpaxRes)
    directory, fasta = kp.kpax_pairwise()

    file = kp.processing_kalign(ARGS.target.split('/')[-1], directory, ARGS.result_file)
