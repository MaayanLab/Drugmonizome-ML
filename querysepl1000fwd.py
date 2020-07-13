""" Class for downloading data from the SEP-L1000 and L1000FWD projects
    and querying drugs.
"""

import json, requests
import pandas as pd
import numpy as np

L1000FWD_URL = 'http://amp.pharm.mssm.edu/L1000FWD/'
L1000FWD_METADATA = 'data/L1000FWD/Drugs_metadata.csv'

def _convert_pert_id_to_InChI(pert_ids):
    """ Given a list of pert_id drug identifiers from the L1000FWD project,
        converts them to a list of InChI keys by reading the L1000FWD metadata.
        All ids in pert_ids must be present in the metadata.
    """
    l1000meta_df = pd.read_csv(L1000FWD_METADATA, index_col=0)
    return list(map(lambda s: s.replace('InChIKey=', '') if isinstance(s, str) else s, l1000meta_df['inchi_key'].loc[pert_ids]))

def _get_drugs_in_metadata(pert_ids):
    l1000meta_df = pd.read_csv(L1000FWD_METADATA, index_col=0)
    return set(pert_ids).intersection(l1000meta_df.index)

def query_drug_names(names, verbose=0):
    """ Given a list of drug names, queries them through the L1000FWD API and converts
        them to a corresponding set of InChI keys. Drug names without matches are
        excluded, and all matches to a drug name are included if there are multiple
        matches. (Adapted from http://amp.pharm.mssm.edu/l1000fwd/api_page).
    """
    all_pert_ids = set()
    for query_string in names:
        query_string = query_string.replace(' ', '-').upper()
        response = requests.get(L1000FWD_URL + 'synonyms/' + query_string)
        found_match = False
        if response.status_code == 200:
            for result in response.json():
                if query_string == result['Name'].upper():
                    all_pert_ids.add(result['pert_id'])
                    found_match = True
        if verbose and not found_match:
            print(query_string + ' not found')

    all_pert_ids = _get_drugs_in_metadata(all_pert_ids)
    return _convert_pert_id_to_InChI(all_pert_ids)

def get_drug_names(keys):
    """ Given a list of drug InChI keys, converts them to a corresponding list of drug names.
    """
    l1000meta_df = pd.read_csv('data/L1000FWD/Drugs_metadata.csv', index_col=5)
    l1000meta_df.index = l1000meta_df.index.map(lambda s: s.replace('InChIKey=', '') if isinstance(s, str) else s)
    l1000meta_df = l1000meta_df.iloc[np.logical_not(l1000meta_df.index.duplicated())]

    return list(l1000meta_df['pert_iname'].reindex(keys))


### TODO: convert probe id to gene name