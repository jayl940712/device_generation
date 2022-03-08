import gdspy
from .basic import basic
from .glovar import tsmc40_glovar as glovar
from .Pin import Pin

# Standard Rules from glovar.py
min_w = glovar.min_w
layer = glovar.layer
datatype = glovar.datatype
sp = glovar.sp
en = glovar.en
ex = glovar.ex
NWELL_GR = glovar.NWELL_GR
NW_OD = glovar.NW_OD
NP_OD = glovar.NP_OD
OD_W = glovar.OD_W
GRID = glovar.GRID

# Special Rules for Capacitor
SP_CAP = 0.2
SP_BOT = 0.8
EN_VIA = 0.2
SP_VIA = 0.2
W_VIA = 0.2

class Capacitor:
    def __init__(self, name, w, l, m_bot=3, attr=[], flatten=True):
        self.name = name
        self.w = w
        self.l = l
        self.m_bot = m_bot
        self.plus = Pin('PLUS')
        self.minus = Pin('MINUS')
        self.cell = gdspy.Cell(name, True)
        self.finger_core()
        if flatten:
            self.flatten()
        #self.print_pins()

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

# Compatible with gdspy
    def to_gds(self, *args):
        """
        @param first: outfile
        @param second: multiplier
        The reason for have varidic length of variables:
        different versions of gdspy have different interfaces for this callback function.
        What is silly is that the previous version is ( multiplier)
        but the new version becomes (outfile, multiplier).
        So it is not very suitable to just give the second argument a default None.
        """
        if len(args) == 1:
            return self.cell.to_gds(args[0])
        elif len(args) == 2:
            return self.cell.to_gds(args[0], args[1])

    def finger_core(self):
        self.origin = [0, 0]
        # Add top plate
        w = self.w - 2 * SP_CAP
        h = self.l - 2 * SP_CAP
        contact_num_w = int((w-2*EN_VIA+SP_VIA)/(W_VIA+SP_VIA))
        contact_num_h = int((h-2*EN_VIA+SP_VIA)/(W_VIA+SP_VIA))
        contact_space_w = (w-2*EN_VIA-contact_num_w*W_VIA)/(contact_num_w-1)+W_VIA
        contact_space_h = (h-2*EN_VIA-contact_num_h*W_VIA)/(contact_num_h-1)+W_VIA
        if contact_space_w < SP_VIA + W_VIA:
            contact_num_w = contact_num_w - 1 
            contact_space_w = (h-2*EN_VIA-contact_num_w*W_VIA)/(contact_num_w-1)+W_VIA
        if contact_space_h < SP_VIA + W_VIA:
            contact_num_h = contact_num_h - 1 
            contact_space_h = (h-2*EN_VIA-contact_num_h*W_VIA)/(contact_num_h-1)+W_VIA
        contact_space_h = round(contact_space_h/basic.GRID)*basic.GRID  
        contact_space_w = round(contact_space_w/basic.GRID)*basic.GRID  
        x_offset = (w-W_VIA-contact_space_w*(contact_num_w-1))*0.5 
        y_offset = (h-W_VIA-contact_space_h*(contact_num_h-1))*0.5 
        x_offset = round(x_offset/basic.GRID)*basic.GRID
        y_offset = round(y_offset/basic.GRID)*basic.GRID
        met_layer = 'M' + str(self.m_bot+1)
        self.top_shape = gdspy.Rectangle((0, 0),(w,h),basic.layer[met_layer], basic.datatype[met_layer])
        self.cell.add(self.top_shape)
        contact_cell = basic.contact(self.m_bot)
        contact_array = gdspy.CellArray(contact_cell, contact_num_w, contact_num_h, [contact_space_w, contact_space_h], [x_offset, y_offset])
        self.cell.add(contact_array)
        cap_layer = gdspy.Rectangle((-SP_CAP, -SP_CAP), (self.w-SP_CAP, self.l-SP_CAP), layer['CAP'], datatype['CAP'])
        self.cell.add(cap_layer)
        # Add top plate pin shape
        legal_x, legal_y = basic.legal_coord((w, h), self.origin)
        self.plus.add_shape(met_layer, [[0,0], [legal_x, legal_y]])
        # Add bottom plate 
        offset_x = legal_x+min_w['M1']+min_w['SP']
        offset_y = -min_w['M1']-min_w['SP']
        bot_contact_y = max(legal_y+2*(min_w['M1']+min_w['SP']), basic.legal_len(w+SP_BOT-offset_y-min_w['M1']))
        bot_contact = basic.metal_vert(min_w['M1'], bot_contact_y, lay=self.m_bot+1)
        while offset_x - w < min_w['M1']: # 0.4
            offset_x += min_w['M1'] + min_w['SP']
        bot_contact_ref = gdspy.CellReference(bot_contact, (offset_x, offset_y))
        self.cell.add(bot_contact_ref)
        self.minus.add_shape(met_layer, [[offset_x, offset_y], [offset_x+min_w['M1'], offset_y+legal_y+2*(min_w['M1']+min_w['SP'])]])
        met_layer = 'M' + str(self.m_bot)
        bot_y = max(bot_contact_y + offset_y, w + SP_BOT)
        bot_met = gdspy.Rectangle((-SP_BOT, -SP_BOT), (offset_x+min_w['M1'], bot_y), layer[met_layer], datatype[met_layer])
        self.cell.add(bot_met)
        pr_layer = gdspy.Rectangle((-SP_BOT, -SP_BOT), (offset_x, bot_y), layer['PR'], datatype['PR'])
        self.cell.add(pr_layer)

 
    def print_pins(self):
        if not (self.plus.check() and self.minus.check()):
            print("Pin location not legal")

    def flip_vert(self):
        flip_cell = gdspy.Cell(self.cell.name, True)
        bounding_box = self.cell.get_bounding_box()
        x_sym_axis = bounding_box[0][0] + bounding_box[1][0]
        # Floating point error 
        # Since gdsii precision is 5nm, here we only round to 1nm precision
        #x_sym_axis = round(x_sym_axis * 10000) / 10000.0
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
        #self.flatten()
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
