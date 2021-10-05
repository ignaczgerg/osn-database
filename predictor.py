import joblib
import pandas as pd
import numpy as np
from mordred import Calculator, descriptors
from rdkit import Chem
from multiprocessing import freeze_support


descr_used = ['NssssN', 'SssssN', 'ETA_shape_y', 'RNCG', 'AATSC1are',
    'GhoseFilter', 'AATSC1se', 'AATSC1pe', 'JGI4', 'AATSC1i', 'NssssC',
    'SssssC', 'AETA_eta_B', 'piPC8', 'NddsN', 'SIC3', 'JGI3', 'NssNH',
    'MATS1i', 'AATS4v', 'BIC3', 'ATSC6s', 'GATS4d', 'n10FHRing',
    'VE2_Dt', 'VE2_Dzm', 'VE2_DzZ', 'VE2_A', 'VE2_Dzi', 'VE2_Dzare',
    'VE2_Dzpe', 'VE2_Dzse', 'VE2_D', 'AXp-5d', 'MATS1are', 'BCUTs-1l',
    'VE2_Dzv', 'AETA_eta_BR', 'VE2_Dzp', 'NsNH2', 'AETA_eta_FL',
    'SlogP_VSA5', 'AATSC2v', 'ATSC3s', 'SsNH2', 'MATS1pe', 'nBase',
    'piPC7', 'NdCH2', 'Xch-7dv', 'ATSC8Z', 'ATSC4s', 'NssO',
    'BCUTc-1l', 'ATSC7m', 'fMF', 'AATSC3s', 'AETA_eta_RL', 'SssO',
    'GATS2c', 'ATSC8m', 'MPC10', 'ATS8s', 'HybRatio', 'ATSC7Z',
    'FCSP3', 'GATS2se', 'ATSC5are', 'AETA_eta_L', 'Xch-6dv', 'AXp-3d',
    'MATS4pe', 'MPC8', 'C4SP3', 'MPC9', 'AETA_beta', 'n10FAHRing',
    'GATS4pe', 'AATSC0d', 'SsssN', 'AETA_beta_ns_d', 'AETA_beta_ns',
    'SIC0', 'SsBr', 'GGI6', 'RotRatio', 'GATS2are', 'BIC0', 'MWC08',
    'MWC06', 'MWC10', 'NdsN', 'MWC07', 'MWC09', 'GATS2s', 'NdssS',
    'MWC05', 'MWC04', 'VE3_A', 'MPC7', 'SRW08', 'SRW06', 'SRW10',
    'NtCH', 'Xpc-6d', 'AETA_dBeta', 'GATS2pe', 'Xc-4dv', 'VE3_Dzp',
    'MWC03', 'VR3_Dzi', 'ATS8m']


def descripter(smiles, descriptors_used=descr_used):
    freeze_support()
    calc = Calculator(descriptors, ignore_3D=False)
    mol = Chem.MolFromSmiles(smiles)
    descr_lst = calc(mol)
    descr_lst = { descr_to_keep: descr_lst[descr_to_keep] for descr_to_keep in descr_used }
    if descr_lst['GhoseFilter'] == False:
        descr_lst['GhoseFilter'] = 0.0
    else:
        descr_lst['GhoseFilter'] = 1.0
    descr_lst = {k:(v if isinstance(v, np.floating)==True 
                    or isinstance(v, int)== True 
                    or isinstance(v, float)== True else 0) for (k,v) in descr_lst.items()}
    return descr_lst

def predictor(descriptors_lst):
    loaded_model = joblib.load("pls_model.sav")
    prediction = loaded_model.predict(np.fromiter(descriptors_lst.values(), dtype=float).reshape(1,-1).astype('float64'))
    if prediction > 1:
        return [[1.0]]
    elif prediction < 0:
        return [[0.0]]
    else: 
        return prediction




