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
    try:
        atoms = utils.pose_from_wynton(pdbid, use_prody=True)
    except:
        atoms = utils.pose_from_rcsb(pdbid, use_prody=True)
    if chain:
        atoms = atoms.select('chain {}'.format(chain))
    atoms = atoms.select('name CA')
    seq = atoms.getSequence()

    query['query']['parameters']['value'] = seq
    query['query']['parameters']['identity_cutoff'] = id_cutoff / 100
    print('Sequence:')
    print(query['query']['parameters']['value'])

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
