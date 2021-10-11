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
SUB_GR = glovar.SUB_GR

class Mosfet:
# Currently supported Mosfet types: nch, pch
# Supported attribute keywords: od25, ud, od, na, lvt, hvt, mac
# Note: na devices should only be compatible with nch
# Unsupported attributes: x, dnw, hv25, snw, na25
# GATE Connection metal check legalization if min_w['M1'] not 0.07
    def __init__(self, nch, name, w, l, nf, attr=[], spectre=True, pinConType=None, bulkCon=[]):
        self.nmos = nch
        self.name = name
        if spectre:
            self.w = w/nf
        else:
            self.w = w 
        self.l = l
        self.nf = nf
        self.drain = Pin('D')
        self.gate = Pin('G')
        self.source = Pin('S')
        self.bulk = Pin('B')
        self.cell = gdspy.Cell(name, True)
        self.mos_core()
        self.doping_layer()
        if 'lvt' in attr:
            self.vth_layer(True)
        if 'hvt' in attr:
            assert(not self.nmos), "HVT only valid for PMOS"
            self.vth_layer(False)
        self.bulkCon = bulkCon
        self.nwell_gr()
        self.connect_pin(pinConType)
        self.connect_pin_bulk()
        self.flatten()
        self.print_pins()

    def pin(self):
        # TODO Future version should include all shapes
        if not self.nmos and NWELL_GR:
            return [self.drain, self.gate, self.source, self.bulk]
        else:
            return [self.drain, self.gate, self.source]
       
    def flatten(self):
        if self.origin:
            self.origin = [self.origin[0] + 0.5*min_w['M1'], self.origin[1] + 0.5 * min_w['M1']]
            temp = gdspy.CellReference(self.cell, (-self.origin[0],-self.origin[1]))
            self.cell = gdspy.Cell(self.name, True)
            self.cell.add(temp)
            self.drain.adjust(self.origin)
            self.gate.adjust(self.origin)
            self.source.adjust(self.origin)
            self.bulk.adjust(self.origin)
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

    def vth_layer(self, lvt):
        if lvt:
            vth_layer = layer['VTL']
            vth_datatype = datatype['VTL']
        else:
            vth_layer = layer['VTH']
            vth_datatype = datatype['VTH']
        vth_shape = gdspy.Rectangle(self.nw_ll, self.nw_ur, vth_layer, vth_datatype)
        self.cell.add(vth_shape)
        #self.flatten()    

    def mos_core(self):
        self.nf_change = False
        self.li_m1 = (min_w['M1'] - min_w['LI'])*0.5
    ### Source Drain Metal Contact
        if self.nf == 1:
            m1_cell = basic.metal_vert(min_w['M1'], self.w)
        else:
            m1_cell = basic.metal_vert(min_w['M1'], self.w, 0)
        m1_y_legal = -basic.legal_len(self.w) + self.w
        m1_cell_space = basic.legal(self.l + min_w['CO'] + 2 * sp['CO']['PO'])
        m1_array_offset = self.l + 0.5*(m1_cell_space - self.l - min_w['M1']) - m1_cell_space
        if self.nf == 1:
            m1_legal_shape = gdspy.Rectangle((0,m1_y_legal), (min_w['M1'], self.w), layer['M1'], datatype['M1'])
            m1_cell.add(m1_legal_shape)
        else:
            pass
        m1_array = gdspy.CellArray(m1_cell, (self.nf+1), 1, [m1_cell_space, 0], [m1_array_offset, 0])
        self.cell.add(m1_array)
        self.origin = [m1_array_offset, m1_y_legal]
    ### Gate Poly
        gate_cell = gdspy.Cell('GATE', True)
        gate_shape = gdspy.Rectangle((0, -ex['PO']['OD']), (self.l, self.w+ex['PO']['OD']), layer['PO'], datatype['PO'])
        gate_cell.add(gate_shape)
        self.gate_space = m1_cell_space
        gate_array = gdspy.CellArray(gate_cell, self.nf, 1, [self.gate_space, 0])
        self.cell.add(gate_array)
    ### OD Layer
        self.x_od1 = -(self.gate_space - min_w['M1'] - self.l)/2 - min_w['CO'] - 0.5*(min_w['M1'] - min_w['CO']) - en['OD']['CO']
        if self.nf % 2 == 0 or self.nf == 1:
            self.x_od2 = (self.nf-1)*self.gate_space + self.l - self.x_od1
            od_shape = gdspy.Rectangle((self.x_od1, 0), (self.x_od2, self.w), layer['OD'], datatype['OD'])
        else:
            self.x_od2 = self.nf*self.gate_space + self.l - self.x_od1
            od_shape = gdspy.Rectangle((self.x_od1, 0), (self.x_od2 - self.gate_space, self.w), layer['OD'], datatype['OD'])
        self.cell.add(od_shape)
    ### GATE Connection
        # Single Finger, no need to connect source/drain
        if self.nf == 1: 
            # Gate Contact
            width = min_w['CO'] + 2 * en['M1']['CO']
            self.gate_ext_len = min_w['SP']
            #self.gate_ext_len = basic.legal(ex['PO']['OD'] + 0.5 * (min_w['CO'] + 2*en['PO']['CO'] - min_w['M1'])) - min_w['M1']
            if (width+en['PO']['CO']-en['M1']['CO'])*2 >= self.l:
                x_pos = 0.5*(width-self.l) 
                y_pos = self.w + self.gate_ext_len - 0.5 * (min_w['CO'] + 2*en['PO']['CO'] - min_w['M1'])
                # Legalization added con_shape
                gate_con_shape = gdspy.Rectangle((0, self.w), (self.l, y_pos+en['PO']['CO']), layer['PO'], datatype['PO'])
                self.cell.add(gate_con_shape)
                # Legalization for top metal
                y_pos = y_pos + en['PO']['CO']-0.5*(min_w['M1']-min_w['CO'])
                #x_pos_legal1 = basic.legal_coord((-x_pos, y_pos),self.origin,1)[0]
                #x_pos_legal2 = basic.legal_coord((-x_pos+width, y_pos+min_w['M1']),self.origin,3)[0]
                x_pos_legal1 = m1_array_offset
                x_pos_legal2 = m1_array_offset + self.gate_space 
                gate_m1_shape = gdspy.Rectangle((x_pos_legal1, y_pos),(x_pos_legal2+min_w['M1'], y_pos+min_w['M1']), layer['M1'], datatype['M1'])
                self.cell.add(gate_m1_shape)
                pm_cell = basic.poly_metal_hori(width*2, width)
                pm_cell_ref = gdspy.CellReference(pm_cell, (-width+0.5*self.l, y_pos))
                self.cell.add(pm_cell_ref)
                # Adding Pin:G
                self.gate.add_shape('M1', gate_m1_shape.get_bounding_box())
            else:
                y_pos = self.w + self.gate_ext_len + 0.5*(min_w['M1']+min_w['CO']) + en['PO']['CO']
                gate_con_shape = gdspy.Rectangle((0, self.w), (self.l, y_pos), layer['PO'], datatype['PO'])
                self.cell.add(gate_con_shape)
                y_pos = self.w + self.gate_ext_len
                m1_gate = basic.poly_metal_hori(self.l-2*(en['PO']['CO']-en['M1']['CO']), min_w['M1'])
                m1_gate_ref = gdspy.CellReference(m1_gate, (en['PO']['CO']-en['M1']['CO'], y_pos))
                self.cell.add(m1_gate_ref)
                # Legalization M1
                x_pos_legal1 = m1_array_offset
                x_pos_legal2 = m1_array_offset + self.gate_space 
                gate_m1_shape = gdspy.Rectangle((x_pos_legal1, y_pos),(x_pos_legal2+min_w['M1'], y_pos+min_w['M1']), layer['M1'], datatype['M1'])
                self.cell.add(gate_m1_shape)
                # Adding Pin:G
                self.gate.add_shape('M1', gate_m1_shape.get_bounding_box())
        # Multiple Finger, need to connect source/drain
        else:
            # Source/Drain M1 Conection
            # Gate Connection
            gate_extension = gdspy.Cell('GATE_EXT', True)
            # Modified for legal
            self.gate_ext_len = 2 * min_w['SP'] + min_w['M1']
            gate_ext_shape = gdspy.Rectangle((0, self.w), (self.l, self.w+self.gate_ext_len+en['PO']['CO']), layer['PO'], datatype['PO'])
            gate_extension.add(gate_ext_shape)
            gate_extension_array = gdspy.CellArray(gate_extension, self.nf, 1, [self.gate_space, self.gate_space])
            self.cell.add(gate_extension_array)
            #end else
            y_pos = self.w + self.gate_ext_len - ( en['PO']['CO'] - 0.5 * ( min_w['M1'] - min_w['CO'] ) )
            gate_contact = basic.poly_metal_hori(self.gate_space*(self.nf-1)+self.l, min_w['M1'])
            gate_contact_ref = gdspy.CellReference(gate_contact, (0, self.w+self.gate_ext_len))
            self.cell.add(gate_contact_ref)
            # Generate symmetry source connections
            # self.nf changes here!!!
            if self.nf % 2 == 1:
                self.nf_change = True
                self.nf += 1
                m1_dummy_x = m1_array_offset + m1_cell_space * self.nf + self.li_m1
                m1_dummy_source = gdspy.Rectangle((m1_dummy_x, 0), (m1_dummy_x+min_w['LI'],self.w), layer['LI'], datatype['LI'])
                self.cell.add(m1_dummy_source)
            m1_square_s = gdspy.Cell('M1_SQUARE', True)
            m1_square_d = gdspy.Cell('M1_SQUARE', True)
            m1_sq_shape_s = gdspy.Rectangle((self.li_m1, 0), (min_w['LI']+self.li_m1, min_w['SP']+self.li_m1), layer['LI'], datatype['LI'])
            m1_sq_shape_d = gdspy.Rectangle((self.li_m1, -self.li_m1), (min_w['LI']+self.li_m1, min_w['SP']+basic.legal_len(self.w)-self.w), layer['LI'], datatype['LI'])
            m1_square_s.add(m1_sq_shape_s)
            m1_square_d.add(m1_sq_shape_d)
            source_count = int(self.nf/2+1)
            drain_count = int((self.nf+1)/2)
            m1_source = gdspy.CellArray(m1_square_s, source_count, 1, [2*self.gate_space, self.gate_space], [m1_array_offset, self.w])
            m1_drain = gdspy.CellArray(m1_square_d, drain_count, 1, [2*self.gate_space, self.gate_space], [m1_array_offset+self.gate_space, -min_w['M1']+self.origin[1]])
            m1_source_hori_bb = [ [m1_array_offset, self.w+min_w['SP']], [m1_array_offset+(2*source_count-2)*self.gate_space+min_w['M1'], self.w+min_w['SP']+min_w['M1']] ]
            m1_drain_hori_bb = [ [m1_array_offset+self.gate_space, self.origin[1]-min_w['SP']-min_w['M1']], [m1_array_offset+(2*drain_count-1)*self.gate_space+min_w['M1'], self.origin[1]-min_w['SP']] ]
            m1_source_hori = basic.metal_hori(m1_source_hori_bb[1][0]-m1_source_hori_bb[0][0], m1_source_hori_bb[1][1]-m1_source_hori_bb[0][1], licon=False)
            m1_source_hori_ref = gdspy.CellReference(m1_source_hori, m1_source_hori_bb[0])
            m1_drain_hori= basic.metal_hori(m1_drain_hori_bb[1][0]-m1_drain_hori_bb[0][0], m1_drain_hori_bb[1][1]-m1_drain_hori_bb[0][1], licon=False)
            m1_source_hori_ref = gdspy.CellReference(m1_source_hori, m1_source_hori_bb[0])
            m1_drain_hori_ref = gdspy.CellReference(m1_drain_hori, m1_drain_hori_bb[0])
            self.cell.add(m1_source)
            self.cell.add(m1_source_hori_ref)
            if drain_count > 1:
                self.cell.add(m1_drain)
                self.cell.add(m1_drain_hori_ref)
            else:
                m1_drain_vert_bb = [[self.origin[0]+self.gate_space, 0], [self.origin[0]+self.gate_space+min_w['M1'], self.w]]
                m1_drain_vert = basic.metal_vert(m1_drain_vert_bb[1][0]-m1_drain_vert_bb[0][0], m1_drain_vert_bb[1][1]-m1_drain_vert_bb[0][1])
                m1_drain_vert_ref = gdspy.CellReference(m1_drain_vert, m1_drain_vert_bb[0])
                m1_drain_vert_m1_shape = gdspy.Rectangle((m1_drain_vert_bb[0][0], self.origin[1]),m1_drain_vert_bb[1],layer['M1'],datatype['M1'])
                self.cell.add(m1_drain_vert_ref)
                self.cell.add(m1_drain_vert_m1_shape)
            # Legalization M1 Gate
            y_pos = self.w + self.gate_ext_len
            x_pos_legal1 = m1_array_offset
            x_pos_legal2 = m1_array_offset + self.nf*self.gate_space 
            gate_m1_shape = gdspy.Rectangle((x_pos_legal1, y_pos),(x_pos_legal2+min_w['M1'], y_pos+min_w['M1']), layer['M1'], datatype['M1'])
            self.cell.add(gate_m1_shape)
            # Adding Pin:G
            self.gate.add_shape('M1', gate_m1_shape.get_bounding_box())
        # Adding Pin:S/D
        if self.nf == 1:
            self.source.add_shape('M1', [[m1_array_offset, self.origin[1]], [m1_array_offset+min_w['M1'], self.w]])
            self.drain.add_shape('M1', [[m1_array_offset+self.gate_space, self.origin[1]], [m1_array_offset+min_w['M1']+self.gate_space, self.w]])
        else:
            self.source.add_shape('M1', m1_source_hori_ref.get_bounding_box())
            if drain_count > 1:
                self.drain.add_shape('M1', m1_drain_hori_ref.get_bounding_box())
            else:
                self.drain.add_shape('M1', m1_drain_vert_m1_shape.get_bounding_box())
            for i in range(self.nf+1):
                if i % 2 == 0:
                    self.source.add_shape('LI', [[self.origin[0]+i*self.gate_space, self.origin[1]], [self.origin[0]+i*self.gate_space+min_w['M1'], self.w]])
                else:
                    self.drain.add_shape('LI', [[self.origin[0]+i*self.gate_space, self.origin[1]], [self.origin[0]+i*self.gate_space+min_w['M1'], self.w]])
        #self.flatten()
    
    def nwell_gr(self):
        if not self.nmos and NWELL_GR:
            nwell_gr, self.bulk = basic.nwell_GR(self.cell.get_bounding_box()[0], self.cell.get_bounding_box()[1],self.origin, True)
            #self.add_bulk_shape(nwell_gr.get_polygons(True)[(layer['M1'],0)])
            nwell_gr_ref = gdspy.CellReference(nwell_gr)
            self.cell.add(nwell_gr_ref)
            #self.flatten()
        if self.nmos and SUB_GR:
            sub_gr, self.bulk = basic.sub_GR(self.cell.get_bounding_box()[0], self.cell.get_bounding_box()[1],self.origin, True)
            #self.add_bulk_shape(sub_gr.get_polygons(True)[(layer['M1'],0)])
            sub_gr_ref = gdspy.CellReference(sub_gr)
            self.cell.add(sub_gr_ref)
            #self.flatten()

    def add_bulk_shape(self, shapes):
        # TODO future should include all shapes
        shape = shapes[3] # This is the vertical top strip of M1
        self.bulk.add_shape('M1', [shape[0], shape[2]])

    def doping_layer(self):
        if self.nmos:
            doping_layer = layer['NP']
            doping_datatype = datatype['NP']
        else:
            doping_layer = layer['PP']
            doping_datatype = datatype['PP']
    # Define NP/PP and NW Shapes
        self.dope_ll = [self.x_od1-ex['NP']['OD'], -ex['NP']['OD']]
        self.dope_ur = [self.x_od2+ex['NP']['OD'], self.w+ex['NP']['OD']]
        self.nw_ll = [self.x_od1-en['NW']['OD'][0], -en['NW']['OD'][1]]
        self.nw_ur = [self.x_od2+en['NW']['OD'][0], self.w+en['NW']['OD'][1]]
        # Shrink dope due to self.nf change
        # Draw NP/PP
        doping_shape = gdspy.Rectangle(self.dope_ll, self.dope_ur, doping_layer, doping_datatype)
        self.cell.add(doping_shape)
        # For PMOS
        if not self.nmos:
            # Draw NW
            nw_shape = gdspy.Rectangle(self.nw_ll, self.nw_ur, layer['NW'], datatype['NW'])
            self.cell.add(nw_shape)
        #self.flatten()

    def print_pins(self):
        if not (self.drain.check() and self.gate.check() and self.source.check() and self.bulk.check()):
            print(self.name, "Pin location not legal")
        #print(self.drain, self.gate, self.source, self.bulk)

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
        self.drain.flip_vert(x_sym_axis)
        self.gate.flip_vert(x_sym_axis)
        self.source.flip_vert(x_sym_axis)
        self.bulk.flip_vert(x_sym_axis)

    def connect_pin_bulk(self):
        # This is only valid for device with guard ring
        if 1 in self.bulkCon:
            x = self.origin[0] + (self.nf / 2) * self.gate_space + self.li_m1
            _, ll, ur = self.gate.shape[0]
            _, _, ur_bulk = self.bulk.shape[1]
            con = gdspy.Rectangle((x, ur_bulk[1]-min_w['M1']+self.li_m1),(x+min_w['LI'], ll[1]+self.li_m1), layer['LI'], datatype['LI'])
            self.cell.add([con])
        if 2 in self.bulkCon:
            _, _, ur = self.source.shape[0]
            _, _, ur_bulk = self.bulk.shape[0]
            con = gdspy.Rectangle((ur[0]-min_w['M1']+self.li_m1, ur_bulk[1]-min_w['M1']+self.li_m1),(ur[0]-self.li_m1, ur[1]), layer['LI'], datatype['LI'])
            self.cell.add(con)
        if 0 in self.bulkCon:
            _, _, ur = self.drain.shape[0]
            _, _, ur_bulk = self.bulk.shape[0]
            con = gdspy.Rectangle((ur[0]-min_w['M1']+self.li_m1, ur_bulk[1]-min_w['M1']+self.li_m1),(ur[0]-self.li_m1, ur[1]), layer['LI'], datatype['LI'])
            self.cell.add(con)

    def connect_pin(self, conType):
        # conType, 0: G,S; 1: G,D; 2: S,D;
        # special case G,D; S,B; call connect_pin then connect_pin_bulk
        if conType in ['GS','SG']:
            _, ll, _ = self.source.shape[0]
            _, _, ur = self.gate.shape[0]
            if self.nf == 1:
                con = gdspy.Rectangle((ll[0]+self.li_m1, 0),(ll[0]+self.li_m1+min_w['LI'], ur[1]-self.li_m1), layer['LI'], datatype['LI'])
                patch = gdspy.Rectangle((ll[0]+self.li_m1, ur[1]-min_w['M1']+en['M1']['CO']), (ll[0]+min_w['M1'], ur[1]-min_w['M1']+en['M1']['CO']+min_w['LI']), layer['LI'], datatype['LI'])
                self.cell.add(patch)
            else:
                x = self.origin[0] + (self.nf / 2) * self.gate_space
                con = gdspy.Rectangle((x, ll[1]+self.li_m1),(x+min_w['M1'], ur[1]-self.li_m1), layer['LI'], datatype['LI'])
            self.cell.add(con)
        elif conType in ['GD','DG']: 
            _, ll, _ = self.source.shape[0]
            _, _, ur = self.gate.shape[0]
            if self.nf == 1:
                con = gdspy.Rectangle((ll[0]+self.li_m1, 0),(ll[0]+self.li_m1+min_w['LI'], ur[1]-self.li_m1), layer['LI'], datatype['LI'])
                patch = gdspy.Rectangle((ll[0]+self.li_m1, ur[1]-min_w['M1']+en['M1']['CO']), (ll[0]+min_w['M1'], ur[1]-min_w['M1']+en['M1']['CO']+min_w['LI']), layer['LI'], datatype['LI'])
                self.cell.add(patch)
            else:
                x = self.origin[0] + (self.nf / 2) * self.gate_space
                con = gdspy.Rectangle((x, ll[1]+self.li_m1),(x+min_w['M1'], ur[1]-self.li_m1), layer['LI'], datatype['LI'])
            self.cell.add(con)
            # swap source/drain
            tempPin = self.source
            self.source = self.drain
            self.drain = tempPin
        elif conType in ['SD','DS']:
            if self.nf == 1:
                _, ll, _ = self.source.shape[0]
                _, _, ur = self.drain.shape[0]
                con1 = gdspy.Rectangle((ll[0], ll[1]),(ur[0], ll[1]+min_w['M1']), layer['M1'], datatype['M1'])
                con2 = gdspy.Rectangle((ll[0], ur[1]-min_w['M1']),(ur[0], ur[1]), layer['M1'], datatype['M1'])
                self.cell.add([con2])#,con2])
            else:
                _, ll, _ = self.source.shape[1]
                _, _, ur = self.source.shape[-1]
                con = gdspy.Rectangle((ll[0]+self.li_m1, 0),(ur[0]-self.li_m1, min_w['LI']), layer['LI'], datatype['LI'])
                self.cell.add(con)

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
