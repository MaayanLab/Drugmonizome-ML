from PubChemQuery import PubChemQuery
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial

class DrugNameConverter():
    """
    Class for converting drug names to InChI keys using PubChem API to query drug names and RDKit for generating InChI keys.
    Includes options for using isomeric forms and for removing salts from drugs.
    """

    query = PubChemQuery()
    remover = SaltRemover()

    @classmethod
    def to_inchi_keys(cls, name, isomeric=True, strip_salts=True):
        inchi_keys = set()
        for smiles in cls.query.name_to_smiles(name, isomeric=isomeric):
            mol = Chem.MolFromSmiles(smiles)
            inchi_keys.add(Chem.MolToInchiKey(mol))
            if strip_salts:
                stripped_mol = cls.remover.StripMol(mol, dontRemoveEverything=True)
                inchi_keys.add(Chem.MolToInchiKey(stripped_mol))
        return inchi_keys

    @classmethod
    def batch_to_inchi_keys_single_thread(cls, names, verbose=0, **kwargs):
        all_inchi_keys = {}
        names = set(names)
        for name in names:
            inchi_keys = cls.to_inchi_keys(name, **kwargs)
            all_inchi_keys[name] = inchi_keys

            if verbose:
                print(f'Completed { len(all_inchi_keys) }/{ len(names) } drugs...', end='\r')
        return all_inchi_keys

    @classmethod
    def batch_to_inchi_keys(cls, names, num_cores=None, **kwargs):
        if num_cores is None:
            num_cores = min(max(mp.cpu_count(), 1), 12)
        names = list(set(names))
        with Pool(num_cores) as p:
            res = p.map(partial(cls.to_inchi_keys, **kwargs), names)
        return dict(zip(names, res))
