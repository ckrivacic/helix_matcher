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


def find_homologs(workspace, target, id_cutoff):
    tarpath = workspace.target_path_clean
    pymol.cmd.load(tarpath, 'complex')
    fastaseq = pymol.cmd.get_fastastr()
    seqlist = fastaseq.split('\n')

    print('Finding homologs for {}'.format(tarpath))
    print(fastaseq)
    seq = ''

    for subseq in seqlist:
        if not subseq.startswith('>'):
            seq += subseq
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
