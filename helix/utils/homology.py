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
import helix.workspace as ws
from helix.utils import utils

query_template = {
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
    atoms = utils.pose_from_wynton(pdbid, prody=True)
    if chain:
        atoms = atoms.select('chain {}'.format(chain))
    atoms = atoms.select('name CA')
    seq = atoms.getSequence()

    query['query']['parameters']['value'] = seq
    query['query']['parameters']['identity_cutoff'] = id_cutoff

    response = response.get("https://search.rcsb.org/rcsbsearch/v1/query",
            {"json": json.dumps(query_dict)})
    print('Response:')
    print(response)
    if "result_set" in response.json():
        results = [result["identifier"] for result in response.json()["result_set"]]
    else:
        print("No results found for target {}.".format(tarpath))
        results = []

    return results
