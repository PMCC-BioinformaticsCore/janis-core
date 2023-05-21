#!/usr/bin/env python3
"""
The class for residue mapping from PDB to UniProt and Vice-versa.
"""

from urllib.request import urlopen
import os
import json
import urllib.response
import wget
from bs4 import BeautifulSoup
from sh import gunzip

class ResidueMapper():
    """
    Class for mapping the residue numbering between PDB <--> UniProt residues
    """

    def __init__(self, src_id, res_pos, path=None):
        self.src_id = src_id
        self.res_pos = [res_pos] if isinstance(res_pos, int) else res_pos
        assert isinstance(self.res_pos[0], int), "The residue positions passed are not integers"
        # self.data = pd.read_csv('pdb_chain_uniprot.csv') # File from SIFTS for UniProt <-> PDB
        self.path = os.getcwd() if path is None else path

        if len(self.src_id) > 4:  # to get numbering for PDB residues from UniProt
            try:
                self.mapped = self.unp2pdb_api()
            except urllib.error.HTTPError as err: #response_data
                print("The mappings for given UniProt ID are {0}: {1}".format(err.reason, self.src_id))
        else:  # to get numbering for UniProt residues from PDB
            self.mapped = self.pdb2unp_api()

    def unp2pdb_api(self):
        """
        Function to map uniprot ids to corresponding pdb ids with chains
        :return: list of tuples :[('1B7F', 'A', 'P19339'), ('1B7F', 'B', 'P19339')]
        """
        url = "https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/" + self.src_id
        data_dict = json.loads(urlopen(url).read())[self.src_id.upper()]['PDB']
        new_dict = []
        for pdb in list(data_dict.keys()):
            if '-' in pdb:
                continue
            new_dict.extend([(self.src_id.upper(), pdb, x['chain_id']) for x in data_dict[pdb]])
        assert len(new_dict) > 0, "The given UniProt ID can not be mapped to any PDB Ids"
        return new_dict

    def pdb2unp_api(self):
        """
        Function to get all chain ids and map the pdb ids to corresponding uniprot ids
        :return: list of tuples :[('1B7F', 'A', 'P19339'), ('1B7F', 'B', 'P19339')]
        """
        url = "https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/" + self.src_id
        data_dict = json.loads(urlopen(url).read())[self.src_id.lower()]['UniProt']
        new_dict = []
        for unp in list(data_dict.keys()):
            if '-' in unp:
                continue
            new_dict.extend([(self.src_id, x['chain_id'], unp) for x in data_dict[unp]['mappings']])
        assert len(new_dict) > 0, "The given PDB ID can not be mapped to any UniProt Ids"
        return new_dict

    def get_sifts_file(self, pdb_code):
        """
        Download residue mapping files from SIFTS (.xml)
        :param pdb_code: Protein Data Bank (PDB) identifier
        :return: Nothing
        """
        # prev_url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'+pdb_code.upper()+'.xml.gz'
        url = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/' + pdb_code.lower() + '.xml.gz'
        wget.download(url, out=self.path)
        gunzip(os.path.join(self.path, pdb_code.lower() + '.xml.gz'))

    def gather_residues_xml(self, pdb=None):
        """
        Download SIFTS file if not already & collect all the residue entities together
        :param pdb: Protein Data Bank (PDB) identifier
        :return: all residue elements gathered from SIFTS file
        """
        pdb = pdb if pdb is not None else self.src_id
        if pdb + '.xml' not in os.listdir(self.path):
            self.get_sifts_file(pdb)

        data = open(os.path.join(self.path, pdb + '.xml'), 'r').read()  # Reading the xml file
        data = BeautifulSoup(data, "xml")
        entity = data.find_all('entity')

        all_res = []
        for one in entity:  # get all residues at one place to iterate over
            part = one.find_all('residue')
            if len(all_res) == 0:
                all_res = part
                continue
            all_res.extend(part)
        return all_res

    def resmapper_unp2pdb(self, pdb=None, chain=None):
        """
        Residue mapping is done using this function (UniProt --> PDB)
        :param pdb: Protein Data Bank (PDB) identifier
        :param chain: chain identifier from the given PDB structure
        :return: list of tuples
        """
        if pdb is not None:
            self.mapped = [x for x in self.mapped if x[1].lower() == pdb.lower()]
        uni_res_pos_list = [str(a) for a in self.res_pos]
        query_pdbs = list({x[1].upper() for x in self.mapped})

        final = []
        for pdb_id in query_pdbs:
            residues = self.gather_residues_xml(pdb=pdb_id.lower())

            # looking for list of residue numbers (int) from example
            for residue in residues:
                crossref = residue.find_all('crossRefDb')
                pdb = [aa for aa in crossref if aa.get('dbSource') == 'PDB']
                pdb = None if len(pdb) < 1 else pdb[0]
                uniprot = [aa for aa in crossref if aa.get('dbSource') == 'UniProt']
                uniprot = None if len(uniprot) < 1 else uniprot[0]

                if (uniprot is not None) and (uniprot.get('dbAccessionId') == self.src_id) and \
                        (uniprot.get('dbResNum') in uni_res_pos_list):
                    tmp_tup = (uniprot.get('dbSource'), uniprot.get('dbAccessionId'), \
                               uniprot.get('dbResNum'), uniprot.get('dbResName'))
                    tmp_tup += ('PDB', pdb_id, None, None, None) if pdb is None else \
                        (pdb.get('dbSource'), pdb.get('dbAccessionId'), pdb.get('dbChainId'), \
                         pdb.get('dbResNum'), pdb.get('dbResName'))

                    final.append(tmp_tup)

        assert len(final) > 0, "The given residue position(s) is not present in any of the " \
                               "PDB structures {0}".format(', '.join(query_pdbs))
        if chain is not None:
            tmp_tup = [x for x in final if x[-3] == chain]
            return tmp_tup

        return final

    def resmapper_pdb2unp(self, chain=None):
        """
        Residue mapping is done using this function (PDB --> UniProt)
        :param chain: chain identifier from the given PDB structure
        :return: list of tuples[('UniProt', 'P19339', '253', 'G', 'PDB', '1b7f', 'A', '253', 'GLY')]
        """
        self.src_id = self.src_id.lower()
        pdb_res_pos_list = [str(x) for x in self.res_pos]

        residues = self.gather_residues_xml()
        final = []
        # looking for list of residue numbers (int) from example
        for residue in residues:
            crossref = residue.find_all('crossRefDb')

            pdb = [aa for aa in crossref if aa.get('dbSource') == 'PDB']
            pdb = None if len(pdb) < 1 else pdb[0]
            uniprot = [aa for aa in crossref if aa.get('dbSource') == 'UniProt']
            uniprot = None if len(uniprot) < 1 else uniprot[0]
            if (pdb is not None) and (pdb.get('dbResNum') != 'null') and (pdb.get('dbResNum') \
                                                                          in pdb_res_pos_list):
                tmp_tup = ('UniProt', None, None, None) if uniprot is None else \
                    (uniprot.get('dbSource'), uniprot.get('dbAccessionId'), \
                     uniprot.get('dbResNum'), uniprot.get('dbResName'))

                tmp_tup += (pdb.get('dbSource'), pdb.get('dbAccessionId'), pdb.get('dbChainId'), \
                            pdb.get('dbResNum'), pdb.get('dbResName'))
                final.append(tmp_tup)

        assert len(final) > 0, "The given residue position(s) is not present in any chain: {0}" \
                               "".format(', '.join(list({x[1].upper() for x in self.mapped})))
        if chain is not None:
            tmp_tup = [x for x in final if x[-3] == chain]
            return tmp_tup

        return final

    def __str__(self):
        return '{0} can be mapped to as follows: {1}'.format(self.src_id, self.mapped)


if __name__ == "__main__":
    M = ResidueMapper('P19339', 122)
    if not hasattr(M, 'mapped'):
        print('Yes')
    MAP = M.resmapper_unp2pdb(pdb='1b7f')
    if 'null' in MAP[0]:
        print('Map', MAP, M.mapped, M.res_pos)
