# GLOBAL VARIABLES
# All units are in um
# MOCK PDK

# Manufactoring Grid
GRID = 0.005
# Minimum Width Rule
min_w = dict()
min_w['CO'] = 0.05
min_w['M1'] = 0.1
min_w['SP'] = 0.1
# Layer Map
layer = dict()
layer['XXX0'] = 1
layer['XXX1'] = 2
layer['XXX2'] = 3
layer['XXX3'] = 4
layer['XXX6'] = 5
layer['VTH_P'] = 6
layer['VTH_N'] = 7
layer['VIA8'] = 8
layer['VIA7'] = 9
layer['VIA6'] = 10
layer['VIA5'] = 11
layer['VIA4'] = 12 
layer['VIA3'] = 13
layer['VIA2'] = 14
layer['VIA1'] = 15
layer['M8'] = 16
layer['M7'] = 17
layer['M6'] = 18
layer['M5'] = 19
layer['M4'] = 20 
layer['M3'] = 21
layer['M2'] = 22
layer['M1'] = 23
layer['CO'] = 24
layer['XXX10'] = 25
layer['NP'] = 26
layer['PP'] = 27
layer['OD_25'] = 28
layer['PO'] = 29
layer['VTL_P'] = 30
layer['VTL_N'] = 31
layer['NT_N'] = 32
layer['OD'] = 33
layer['NW'] = 34
# Datatype map for metal
datatype = dict()
datatype['M1'] = 0
datatype['M2'] = 0
datatype['M3'] = 0
datatype['M4'] = 0
datatype['M5'] = 0
datatype['M6'] = 0
datatype['M7'] = 0
datatype['M8'] = 0
datatype['VIA1'] = 0
datatype['VIA2'] = 0
datatype['VIA3'] = 0
datatype['VIA4'] = 0
datatype['VIA5'] = 0
datatype['VIA6'] = 0
datatype['VIA7'] = 0
# Minimum Spacing Rule
sp = dict()
sp['CO'] = dict()
sp['CO']['CO'] = 0.1
sp['CO']['PO'] = 0.05
sp['CO']['XXX10'] = 0.2
sp['M1'] = dict()
# Normally M1 spacing is min_w, but this is for spacing of gate horizontal contact
sp['M1']['M1'] = 0.1
# Enclosure Rule
en = dict()
en['M1'] = dict()
en['M1']['CO'] = 0.05
en['OD'] = dict()
en['OD']['PO'] = 0.1
en['OD']['CO'] = 0.05
en['PO'] = dict()
en['PO']['CO'] = 0.05
en['NP'] = dict()
en['NP']['PO'] = 0.1
en['NW'] = dict()
en['NW']['OD'] = [0.1, 0.2] 
en['PP'] = en['NP'] 
en['NT_N'] = dict()
en['NT_N']['OD'] = 0.25
en['XXX3'] = dict()
en['XXX3']['PO'] = 0.2
# Extension Rule
ex = dict()
ex['PO'] = dict()
ex['PO']['OD'] = 0.1
ex['NP'] = dict()
ex['NP']['OD'] = 0.1
ex['PP'] = ex['NP']
ex['XXX10'] = dict()
ex['XXX10']['PO'] = 0.5
# NWELL GuardRing All Hard Encoded
NWELL_GR = True
NW_OD = 0.1
NP_OD = 0.05
OD_W = 0.2
# Sub GR
SUB_GR = False
# This is the M1 gate_hori spacing to S/D fingers
# Currently this feature has been removed
# Should be >= 0
KR_SP = 0
