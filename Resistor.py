import gdspy
import basic
import glovar
from Pin import Pin

# Standard Rules from glovar.py
min_w = glovar.min_w
layer = glovar.layer
sp = glovar.sp
en = glovar.en
ex = glovar.ex
NWELL_GR = glovar.NWELL_GR
NW_OD = glovar.NW_OD
NP_OD = glovar.NP_OD
OD_W = glovar.OD_W
GRID = glovar.GRID

# Special Rule for XXX6 layers
en_pp_po_wo = 0.15

class Resistor:
    def __init__(self, series, name, w, l, seg_num, seg_space=0.18, attr=[]):
    # seg_space will be hornoured, users should check if spacing satisfy grid 
        self.series = series
        self.name = name
        self.w = w
        self.l = l
        self.seg_num = seg_num
        self.seg_space = seg_space
        self.plus = Pin('PLUS')
        self.minus = Pin('MINUS')
        self.m = False
        self.wo = False
        if 'm' in attr:
            self.m = True
        if 'wo' in attr:
            self.wo = True
            en_pp_po = en['PP']['PO']
            en['PP']['PO'] = en_pp_po_wo
        self.cell = gdspy.Cell(name, True) 
        self.res_core()
        if 'wo' in attr:
            self.xxx6_layer()
            en['PP']['PO'] = en_pp_po
        self.flatten()
        self.print_pins()

    def pin(self):
        return [self.plus, self.minus]

    def flatten(self):
        if self.origin:
            self.origin = [self.origin[0] + 0.5*min_w['M1'], self.origin[1] + 0.5 * min_w['M1']]
            temp = gdspy.CellReference(self.cell, (-self.origin[0],-self.origin[1]))
            self.cell = gdspy.Cell(self.name, True)
            self.cell.add(temp)
            self.plus.adjust(self.origin)
            self.minus.adjust(self.origin)
            self.origin = None
        self.cell.flatten()

    def to_gds(self, multiplier):
        #self.flatten()
        return self.cell.to_gds(multiplier)

    def xxx6_layer(self):
        xxx6_shape = gdspy.Rectangle((self.xxx6_x1, -ex['XXX6']['PO']), (self.xxx6_x2, self.cell_poly_w+ex['XXX6']['PO']), layer['XXX6'])
        self.cell.add(xxx6_shape)
        #self.flatten()

    def res_core(self):
# Poly Resistor Core
    # Poly Shape Core
        x_pos1 = en['PO']['CO'] - 0.5 * (min_w['M1'] - min_w['CO'])
        poly_l = 2 * (en['PO']['CO'] + sp['CO']['XXX6'] + min_w['CO']) + self.l
        m1_vert_space = poly_l - 2 * (min_w['M1'] + x_pos1)
        poly_cell = gdspy.Cell('POLY', True)
        self.origin = [x_pos1, 0]
    # Contact Shape
        m1_vert = basic.metal_vert(min_w['M1'], self.w)
        m1_shape_h = self.w
        m1_shape = gdspy.Rectangle((0,0),(min_w['M1'],m1_shape_h),layer['M1'])
        m1_vert.add(m1_shape)
        m1_vert_ref1 = gdspy.CellReference(m1_vert, (x_pos1, 0))
        x_pos2 = basic.legal_coord([x_pos1 + m1_vert_space + min_w['M1'],0],self.origin,3)[0]
        m1_vert_ref2 = gdspy.CellReference(m1_vert, (x_pos2, 0))
        # Moved add poly shape here
        poly_l = poly_l + x_pos2 - (x_pos1 + m1_vert_space + min_w['M1'])
        poly_shape = gdspy.Rectangle((0, 0), (poly_l, self.w), layer['PO'])
        poly_cell.add(poly_shape)
        poly_cell.add(m1_vert_ref1)
        poly_cell.add(m1_vert_ref2)
    # XXX6 Layer
        self.xxx6_x1 = en['PO']['CO'] + min_w['CO'] + sp['CO']['XXX6']
        self.xxx6_x2 = self.xxx6_x1 + self.l
        xxx6_datatype = 0
        if self.m:
            xxx6_datatype = 1
        xxx6_shape = gdspy.Rectangle((self.xxx6_x1, 0), (self.xxx6_x2, self.w), layer['XXX6'], datatype=xxx6_datatype)
        poly_cell.add(xxx6_shape)
    # Poly Array for seg_num
        poly_space = self.w + self.seg_space
        poly_array = gdspy.CellArray(poly_cell, 1, self.seg_num, [0, poly_space])
        self.cell.add(poly_array)
    # M1 Connection for series/parallel
        if self.series and self.seg_num > 1:
            m1_num_1 = int(self.seg_num/2) + 1
            m1_num_2 = int((self.seg_num+1)/2)
            for i in range(m1_num_1):
                if i == m1_num_1 - 1 and self.seg_num % 2 == 0:
                    y1 = 2*(self.w + self.seg_space)*(i-1) + self.w + self.seg_space
                    m1_ll_y = basic.legal_coord((x_pos1,y1),self.origin,1)[1]
                    m1_ur_y = basic.legal_len(self.w+y1-m1_ll_y) + m1_ll_y
                    m1_connect = gdspy.Rectangle((x_pos1,m1_ll_y),(x_pos1+min_w['M1'],m1_ur_y),layer['M1'])
                    m1_ll_y_1 = m1_ll_y
                    m1_ur_y_1 = m1_ur_y
                elif i == 0:
                    m1_ur_y = basic.legal_len(self.w)
                    m1_connect = gdspy.Rectangle((x_pos1,0),(x_pos1+min_w['M1'],m1_ur_y),layer['M1'])
                else:
                    y1 = 2*(self.w + self.seg_space)*(i-1) + self.w + self.seg_space
                    m1_ll_y = basic.legal_coord((x_pos1,y1),self.origin,1)[1]
                    m1_ur_y = basic.legal_len(2*self.w+self.seg_space+y1-m1_ll_y) + m1_ll_y
                    m1_connect = gdspy.Rectangle((x_pos1,m1_ll_y),(x_pos1+min_w['M1'],m1_ur_y),layer['M1'])
                    if i == m1_num_1 - 1:
                        m1_ll_y_1 = m1_ll_y
                        m1_ur_y_1 = m1_ur_y
                self.cell.add(m1_connect)
            for i in range(m1_num_2):
                if i == m1_num_2 - 1 and self.seg_num % 2 != 0:
                    y1 = 2*(self.w + self.seg_space)*i
                    m1_ll_y = basic.legal_coord((x_pos2,y1),self.origin,1)[1]
                    m1_ur_y = basic.legal_len(self.w+y1-m1_ll_y) + m1_ll_y
                    m1_connect = gdspy.Rectangle((x_pos2,m1_ll_y),(x_pos2+min_w['M1'],m1_ur_y),layer['M1'])
                    m1_ll_y_2 = m1_ll_y
                    m1_ur_y_2 = m1_ur_y
                else:
                    y1 = 2*(self.w + self.seg_space)*i
                    m1_ll_y = basic.legal_coord((x_pos1,y1),self.origin,1)[1]
                    m1_ur_y = basic.legal_len(2*self.w+self.seg_space+y1-m1_ll_y) + m1_ll_y
                    m1_connect = gdspy.Rectangle((x_pos2,m1_ll_y),(x_pos2+min_w['M1'],m1_ur_y),layer['M1'])
                    if i == m1_num_2 - 1:
                        m1_ll_y_2 = m1_ll_y
                        m1_ur_y_2 = m1_ur_y
                self.cell.add(m1_connect)
        #elif not self.series and self.seg_num > 1:
        else:
            y1 = (self.seg_num-1)*(self.w+self.seg_space)+self.w
            m1_ur_y = basic.legal_len(y1)
            m1_connect1 = gdspy.Rectangle((x_pos1,0),(x_pos1+min_w['M1'],m1_ur_y),layer['M1'])
            self.cell.add(m1_connect1)
            m1_connect2 = gdspy.Rectangle((x_pos2,0),(x_pos2+min_w['M1'],m1_ur_y),layer['M1'])
            self.cell.add(m1_connect2)
               # PP Layer
        self.cell_poly_w = self.w*self.seg_num+self.seg_space*(self.seg_num-1)
        pp_p1 = [-en['PP']['PO'], -en['PP']['PO']]
        pp_p2 = [poly_l+en['PP']['PO'], self.cell_poly_w+en['PP']['PO']]
        pp_shape = gdspy.Rectangle(pp_p1, pp_p2, layer['PP'])
        self.cell.add(pp_shape)
    # XXX3 Layer
        xxx3_p1 = [-en['XXX3']['PO'], -en['XXX3']['PO']]
        xxx3_p2 = [poly_l+en['XXX3']['PO'], self.cell_poly_w+en['XXX3']['PO']]
        xxx3_shape = gdspy.Rectangle(xxx3_p1, xxx3_p2, layer['XXX3'])
        self.cell.add(xxx3_shape)
        #self.flatten()
    # Adding Pins
        if not self.series or self.seg_num == 1:
            self.minus.add_shape('M1', [[x_pos1, 0], [x_pos1+min_w['M1'], m1_ur_y]])
            self.plus.add_shape('M1', [[x_pos2, 0], [x_pos2+min_w['M1'], m1_ur_y]])
        else:
            self.minus.add_shape('M1', [[x_pos1, 0], [x_pos1+min_w['M1'], basic.legal_len(self.w)]])
            if self.seg_num % 2 == 0:
                self.plus.add_shape('M1', [[x_pos1, m1_ll_y_1], [x_pos1+min_w['M1'], m1_ur_y_1]])
            else:
                self.plus.add_shape('M1', [[x_pos2, m1_ll_y_2], [x_pos2+min_w['M1'], m1_ur_y_2]])

    def print_pins(self):
        if not (self.plus.check() and self.minus.check()):
            print "Pin location not legal"
        #print self.plus, self.minus

    def flip_vert(self):
        flip_cell = gdspy.Cell(self.cell.name, True)
        bounding_box = self.cell.get_bounding_box()
        x_sym_axis = bounding_box[0][0] + bounding_box[1][0]
        polydict = self.cell.get_polygons(by_spec=True)
        for key in polydict:
            layer, datatype = key
            for shape in polydict[key]:
                x_min = shape[0][0]
                y_min = shape[0][1]
                x_max = shape[2][0]
                y_max = shape[2][1]
                x_min_s = x_sym_axis - x_max
                x_max_s = x_sym_axis - x_min
                new_shape = gdspy.Rectangle([x_min_s,y_min], [x_max_s,y_max], layer, datatype=datatype)
                flip_cell.add(new_shape)
        self.cell = flip_cell
        self.plus.flip_vert(x_sym_axis)
        self.minus.flip_vert(x_sym_axis)

    def bounding_box(self):
        if self.origin:
            assert self.origin == [0, 0], "Cell origin not reset."
        bounding_box = self.cell.get_bounding_box()
        ll = list(basic.legal_coord(bounding_box[0],[0,0],1))
        ur = list(basic.legal_coord(bounding_box[1],[0,0],3))
        ll[0] = ll[0] - min_w['M1'] - min_w['SP']
        ll[1] = ll[1] - min_w['M1'] - min_w['SP']
        ur[0] = ur[0] + min_w['M1'] + min_w['SP']
        ur[1] = ur[1] + min_w['M1'] + min_w['SP']
        return [ll,ur]
