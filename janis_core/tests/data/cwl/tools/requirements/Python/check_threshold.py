#!/usr/bin/env python3
"""
This script reads a csv file containing information from Kpax alignment
and creates different set for structures passing and failing the threshold.
"""

import argparse
import json
import pandas as pd

class ThresholdChecker:
    """
    Class to check the thresholds scores for the aligned structures
    """
    def __init__(self, analysis_file, score_x, threshold_used):
        self.query_df = pd.read_csv(analysis_file, index_col=None)
        self.param = score_x
        self.thresh = threshold_used
        self.passed, self.failed = [], []

    def check_threshold(self):
        """
        Checks the threshold for given score for each alignment and returns a list of targets
        passing threshold value for given score
        """
        for row in self.query_df.iterrows():
            target_id = row[1]['#Target In'].rsplit('_', 1)[0]
            target_id = target_id.replace('-', '.')       # Because of Kpax file naming problem
            if row[1][self.param] >= self.thresh:
                self.passed.append(target_id)
            else:
                self.failed.append(target_id)

        return self.passed, self.failed

    @staticmethod
    def prev_fams(fam_file, src_db):
        """Returns all family ids from the fam_file for given src_db"""
        db = ['Pfam', 'CATH']
        assert src_db in db, 'src_db parameter must be from {0}'.format(', '.join(db))
        data = json.load(open(fam_file, 'r'))
        all_fam = []
        for iterate in data:
            all_fam.extend(data[iterate][src_db])
        return [*set(all_fam)]

    @staticmethod
    def compare_fams(prev_fam, next_fam):
        """Returns the list next_fam without any element from prev_fam"""
        return [fam for fam in next_fam if fam not in prev_fam]



###################################################################################################
if __name__ == '__main__':
    P = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, \
        epilog='Enjoy the program! :)')
    P.add_argument("-a", "--analysis", type=str, required=True, \
        help="Result from Kpax alignements.")
    P.add_argument("-f", "--fam_ids", type=str, required=True, \
        help="File with family IDs per iteration.")
    P.add_argument("-s", "--score", type=str, choices=['Kscore', 'Gscore', 'Jscore', 'Mscore', \
        'Tscore', 'Nanchor', 'RMSD-aligned', 'RMSD-matched', 'RMSD-anchor', 'Naligned'], \
            default="Mscore", help="Alignment score to chose from for the evaluation")
    P.add_argument("-t", "--threshold", type=float, default=0.6, \
        help="The threshold value for the chosen score")
    P.add_argument("-px", "--crosspfam", type=str, \
        help="All Pfam domain StIs corresponding to cross-mapped family")
    P.add_argument("-cx", "--crosscath", type=str, \
        help="All CATH domain StIs corresponding to cross-mapped family")

    ARGS = P.parse_args()

    score = '#' + ARGS.score if ARGS.score else '#Mscore'
    TC = ThresholdChecker(ARGS.analysis, score, ARGS.threshold)
    passed_ids, failed_ids = TC.check_threshold()

    pf_passed = [fam for fam in passed_ids if fam.startswith('PF')]
    ca_passed = [fam for fam in passed_ids if not fam.startswith('PF')]

    pf_passed = TC.compare_fams(TC.prev_fams(ARGS.fam_ids, 'Pfam'), pf_passed)
    with open('pfam.json', 'w') as f:
        json.dump({'Pfam': pf_passed}, f, indent=2)

    ca_passed = TC.compare_fams(TC.prev_fams(ARGS.fam_ids, 'CATH'), ca_passed)
    with open('cath.json', 'w') as f:
        json.dump({'CATH': ca_passed}, f, indent=2)

    # if crossmapped instances are passed
    if ARGS.crosspfam is not None:
        pf_cx = json.load(open(ARGS.crosspfam, 'r'))
        pf_pass = {key: val for key, val in pf_cx.items() if key in ca_passed}
        with open('crossmapped_pfam_passed.json', 'w') as f:
            json.dump(pf_pass, f, indent=2)

        pf_fail = {key: val for key, val in pf_cx.items() if key not in ca_passed}
        with open('crossmapped_pfam_failed.json', 'w') as f:
            json.dump(pf_fail, f, indent=2)

    if ARGS.crosscath is not None:
        ca_cx = json.load(open(ARGS.crosscath, 'r'))
        ca_pass = {key: val for key, val in ca_cx.items() if key in pf_passed}
        with open('crossmapped_cath_passed.json', 'w') as f:
            json.dump(ca_pass, f, indent=2)

        ca_fail = {key: val for key, val in ca_cx.items() if key not in pf_passed}
        with open('crossmapped_cath_failed.json', 'w') as f:
            json.dump(ca_fail, f, indent=2)
