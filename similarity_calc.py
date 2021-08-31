import enum
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit import DataStructs

class DataLoader(object):
    def __init__(self, file_url):
        self.file_url = file_url

    def csv_loader(file_url):
        _df = pd.read_csv(file_url, sep=',')
        return _df
    

class Similarity(object):
    def __init__(self, dataframe: pd.DataFrame, reference: str):
        self.dataframe = dataframe
        self.reference = reference

    def rejection_similarity(dataframe, reference=None):
        '''
        Returns the highest Tanimoto similarity of a molecule based on a set of molecules.
        '''
        _ref = rdkit.Chem.RDKFingerprint(rdkit.Chem.MolFromSmiles(reference))
        _mols = [rdkit.Chem.MolFromSmiles(smile) for smile in dataframe['smiles']]
        _fps = [rdkit.Chem.RDKFingerprint(x) for x in _mols]

        def _sim(_ref, _mol):
            return DataStructs.FingerprintSimilarity(_ref, _mol, metric=DataStructs.TanimotoSimilarity)
        similarity_list = [(enum, _sim(_fp, _ref)) for enum, _fp in enumerate(_fps)]
        similarity_list.sort(key=lambda x: x[1])
        return similarity_list[-1], _mols[similarity_list[-1][0]], dataframe.iloc[similarity_list[-1][0]]