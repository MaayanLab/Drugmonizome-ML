import requests
import time
from ExponentialBackoff import ExponentialBackoff

PUBCHEM_BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'

class PubChemQuery():
    """
    Class for making queries through the PubChem REST API.
    Uses exponential backoff to limit server request rates when 503 error is encountered.
    """
    
    backoff = ExponentialBackoff()
    
    @classmethod
    def make_query(cls, url):
        while True:
            time.sleep(cls.backoff.value())
            r = requests.get(url)
            if r.status_code != 503:
                cls.backoff.halve()
                break
            cls.backoff.double()
        return r

    @classmethod
    def query_by_name(cls, name, property):
        return cls.make_query(f'{ PUBCHEM_BASE_URL }/name/{ name }/property/{ property }/TXT')

    @classmethod
    def query_by_smiles(cls, smiles, property):
        return cls.make_query(f'{ PUBCHEM_BASE_URL }/smiles/{ smiles }/property/{ property }/TXT')

    @classmethod
    def name_to_inchi_keys(cls, name):
        r = cls.query_by_name(name, 'InChIKey')
        return set(r.content.decode('utf-8').split('\n')[:-1])

    @classmethod
    def name_to_smiles(cls, name, isomeric=True):
        if isomeric:
            r = cls.query_by_name(name, 'isomericSMILES')
        else:
            r = cls.query_by_name(name, 'canonicalSMILES')
        return set(r.content.decode('utf-8').split('\n')[:-1])