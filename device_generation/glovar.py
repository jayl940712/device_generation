# GLOBAL VARIABLES
# All units are in um

class tsmc40_glovar:
    GRID = 0.005
    min_w = dict()
    min_w['CO'] = 0.17
    min_w['LI'] = 0.17
    min_w['M1'] = 0.33
    min_w['SP'] = 0.33
# Layer Map
    layer = dict()
    layer['LVSDMY'] = 208
    layer['SUB'] = 236
    layer['MOMDMY'] = 155 
    layer['DMEXCL'] = 150 
    layer['RH'] = 117
    layer['RPDMY'] = 115
    layer['VTH'] = 68
    layer['VTL'] = 67
    layer['VIA4'] = 54
    layer['VIA3'] = 53
    layer['VIA2'] = 52
    layer['VIA1'] = 51
    layer['M5'] = 35
    layer['M4'] = 34
    layer['M3'] = 33
    layer['M2'] = 32
    layer['M1'] = 31
    layer['CO'] = 30
    layer['RPO'] = 29
    layer['LCO'] = 28
    layer['LI'] = 27
    layer['NP'] = 26
    layer['PP'] = 25
    layer['NPC'] = 24
    layer['OD_25'] = 18
    layer['PO'] = 17
    layer['NT_N'] = 11
    layer['OD'] = 6
    layer['TAP'] = 5
    layer['NW'] = 3
    layer['CAP'] = 89
    datatype = dict()
    datatype['M1'] = 0
    datatype['M2'] = 0
    datatype['M3'] = 0
    datatype['M4'] = 0
    datatype['M5'] = 0
    datatype['VIA1'] = 0
    datatype['VIA2'] = 0
    datatype['VIA3'] = 0
    datatype['VIA4'] = 0
    datatype['PO'] = 0
    datatype['CO'] = 0
    datatype['LCO'] = 0
    datatype['LI'] = 0
    datatype['NPC'] = 0
    datatype['NP'] = 0
    datatype['PP'] = 0
    datatype['OD'] = 0
    datatype['VTL'] = 0
    datatype['VTH'] = 0
    datatype['TAP'] = 0
    datatype['NW'] = 0
    datatype['SUB'] = 0
    datatype['CAP'] = 0
    sp = dict()
    sp['CO'] = dict()
    sp['CO']['CO'] = 0.19
    sp['CO']['PO'] = 0.06
    sp['CO']['RPO'] = 0.2
    sp['M1'] = dict()
    sp['M1']['M1'] = 0.15 
    en = dict()
    en['M1'] = dict()
    en['M1']['CO'] = 0.03
    en['OD'] = dict()
    en['OD']['PO'] = 0.10
    en['OD']['CO'] = 0.06
    en['PO'] = dict()
    en['PO']['CO'] = 0.08
    en['NP'] = dict()
    en['NP']['PO'] = 0.10
    en['NW'] = dict()
    en['NW']['OD'] = [0.05, 0.15] 
    en['PP'] = en['NP']  
    en['NT_N'] = dict()
    en['NT_N']['OD'] = 0.25
    en['RH'] = dict()
    en['RH']['PO'] = 0.2
    en['NPC'] = dict()
    en['NPC']['CO'] = 0.1
    en['LI'] = dict()
    en['LI']['LCO'] = 0.08
    ex = dict()
    ex['PO'] = dict()
    ex['PO']['OD'] = 0.13
    ex['NP'] = dict()
    ex['NP']['OD'] = 0.125
    ex['PP'] = ex['NP']
    ex['RPO'] = dict()
    ex['RPO']['PO'] = 0.5
    NWELL_GR = True
    NW_OD = 0.10
    NP_OD = 0.05
    OD_W = 0.15
    SUB_GR = False
