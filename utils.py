
from mordred import descriptors
from numpy import zeros

descriptors_used =  ['NssssN', 'SssssN', 'ETA_shape_y', 'RNCG', 'AATSC1are',
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

descriptors_dict = dict(zip(descriptors_used, zeros(len(descriptors_used))))
# print(descriptors_dict)
'''
Iterate over all the key value pairs in dictionary and call the given
callback function() on each pair. Items for which callback() returns True,
add them to the new dictionary. In the end return the new dictionary.

def filter(dictObj, callback):
    newDict = dict()
    # Iterate over all the items in dictionary
    for (key, value) in dictObj.items():
        # Check if item satisfies the given condition then add to new dict
        if callback((key, value)):
            newDict[key] = value
    return newDict


# Filter a dictionary to keep elements only whose keys are even
newDict = filter(descriptors_dict, lambda elem : elem[0] % 2 == 0)
print('Filtered Dictionary : ')
print(newDict)
''' 