'''
Functions for running queries to find homologous proteins, for excluding
from benchmark
'''

import pymol
import traceback
import os
import requests
import yaml
import prody
import json
import helix.workspace as ws
from helix.utils import utils

query = {
    "query":{
        "type": "terminal",
        "service": "sequence",
        "parameters":{
            "identity_cutoff": 0.84,
            "target": "pdb_protein_sequence",
            "value": "NLKASLSLTLKHYVPLSGNLLMPIKSGEMPKFQYSRDEGDSITLIVAESDQDFDYLKGHQLVDSNDLHGLFYVMPRVIRTMQDYKVIPLVAVQVTVFPNRGIAVALTAHHSIADAKSFVMFINAWAYINKFGKDADLLSANLLPSFDRSIIKDLYGLEETFWNEMQDVLEMFSRFGSKPPRFNKVRATYVLSLAEIQKLKNKVLNLRGSEPTIRVTTFTMTCGYVWTCMVKSKDDVVSEESSNDENELEYFSFTADCRGLLTPPCPPNYFGNCLA"
        }
    },
    "return_type": "entry"
}


def find_homologs(pdbid, id_cutoff, chain=None):
    polymers = prody.parsePDBHeader(pdbid, 'polymers')
    seq = ''
    for polymer in polymers:
        if chain:
            if polymer.chid == chain:
                seq += polymer.sequence
        else:
            seq += polymer.sequence

    query['query']['parameters']['value'] = seq
    query['query']['parameters']['identity_cutoff'] = id_cutoff / 100
    query['query']['parameters']['evalue_cutoff'] = 1
    print('Sequence:')
    print(query['query']['parameters']['value'])
    print('ID cutoff:')
    print(query['query']['parameters']['identity_cutoff'])

    response = requests.get("https://search.rcsb.org/rcsbsearch/v1/query",
            {"json": json.dumps(query)})
    if str(response) == '<Response [204]>':
        print("No results found for target {}.".format(pdbid))
        results = []
    else:
        results = [result["identifier"] for result in response.json()["result_set"]]

    print('Results found:')
    print(results)
    return results


def test():
    pdbid = '4hkr'
    find_homologs(pdbid, 95, chain='A')


if __name__=='__main__':
    test()