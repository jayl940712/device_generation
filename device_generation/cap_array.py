# Slightly modified from Capacitor for uniform connections
import gdspy
from .Capacitor import Capacitor
from .glovar import tsmc40_glovar as glovar
from .basic import basic
from .Pin import Pin

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

# Special Rules for Capacitor
W_CON = 0.14 # Cap contact width
W_VIA = 0.07 # VIA width
SP_VIA = 0.09 # VIA to VIA spacing
EN_VIA = 0.05 # Metal Enclose VIA Rule
VIA_OFF = (W_CON - W_VIA)/2
VIA_DIST = (SP_VIA - W_VIA)/2
SP_SUB = 0.5 # Substrate contact to finger space
W_SUB = min_w['M1'] #0.12 # Sub width

class custom_cap:
    def __init__(self, name, w, sp, nf, l, m_bot=3, m_top=5, attr=[], f_tip=0.14, flatten=True):
        self.name = name
        self.w = w
        self.l = l
        self.sp = sp
        self.nf = nf
        self.metal = []
        self.via = []
        self.m_bot = m_bot
        self.m_top = m_top
        self.plus = Pin('PLUS')
        self.minus = Pin('MINUS')
        self.bulk = Pin('BULK')
        self.t2 = False
        if '2t' in attr:
            self.t2 = True
        # Reversing pin order to top metal first, fix for routing
        for i in range(m_bot, m_top+1):
        #for i in range(m_top, m_bot-1, -1):
            self.metal.append('M'+str(i))
            if i != m_top:
                self.via.append(i)
        self.f_tip = f_tip
        self.cell = gdspy.Cell(name, True)
        self.finger_core()
        if self.t2:
            self.t2_layer()
        else:
            self.lvs_layer()
        self.via_layer()
        if flatten:
            self.flatten()
        #self.print_pins()

    def pin(self):
        if self.t2:
            return [self.plus, self.minus]
        else:
            return [self.plus, self.minus, self.bulk]

    def flatten(self):
        if self.origin:
            self.origin = [self.origin[0] + 0.5*min_w['M1'], self.origin[1] + 0.5 * min_w['M1']]
            temp = gdspy.CellReference(self.cell, (-self.origin[0],-self.origin[1]))
            self.cell = gdspy.Cell(self.name, True)
            self.cell.add(temp)
            self.plus.adjust(self.origin)
            self.minus.adjust(self.origin)
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

    def finger_core(self):
# MOMCAP Finger CORE
    # Define CAP finger shapes
    # This should be legal by default if f_tip=170n, sp=70n, w=70n
        self.W_CON = (basic.legal_len(2*W_CON + 2*self.f_tip + self.l) - 2*self.f_tip - self.l) * 0.5
        self.VIA_OFF = int((self.W_CON - W_VIA)*100)/200.0
        self.finger_y1 = self.W_CON + self.f_tip
        self.finger_y2 = self.finger_y1 + self.l
        self.finger_sp = self.sp + self.w
        con_len = basic.legal_len(self.finger_sp*(self.nf-1)+self.w)
        self.con1 = [[0, 0], [con_len, self.W_CON]]
        self.con2 = [[0, self.finger_y2+self.f_tip], [self.con1[1][0], self.con1[1][1]+self.finger_y2+self.f_tip]]
        self.origin = [0, 0]
        for metal in self.metal:
        # Cap core finger cell
            finger_cell = gdspy.Cell("FINGER", True)
            finger_shape = gdspy.Rectangle((0, self.finger_y1), (self.w, self.finger_y2), layer[metal])
            finger_cell.add(finger_shape)
            finger_array = gdspy.CellArray(finger_cell, self.nf, 1, [self.finger_sp, 0])
            self.cell.add(finger_array)
        # Finger extension
            ext_cell = gdspy.Cell("EXT", True)
            ext_shape = gdspy.Rectangle((0, self.W_CON), (self.w, self.finger_y1), layer[metal])
            ext_cell.add(ext_shape)
            bot_ext = gdspy.CellArray(ext_cell, self.nf/2, 1, [2*self.finger_sp, 0])
            bot_cell = gdspy.Cell("BOT_EXT", True)
            bot_cell.add(bot_ext)
            top_ext = gdspy.CellReference(bot_cell, (self.finger_sp, self.finger_y2-self.W_CON))
            self.cell.add([bot_ext, top_ext])
        # Finger connection metal
            con_bot_shape = gdspy.Rectangle(self.con1[0], self.con1[1], layer[metal])
            con_top_shape = gdspy.Rectangle(self.con2[0], self.con2[1], layer[metal])
            self.cell.add([con_bot_shape, con_top_shape])
        # Adding Pins
            self.plus.add_shape(metal, [self.con1[0],self.con1[1]])
            self.minus.add_shape(metal, [self.con2[0], self.con2[1]])
            #self.plus.add_shape(metal, [self.con1[0],[self.con1[1][0],self.con1[0][1]+min_w['M1']]])
            #self.minus.add_shape(metal, [[self.con2[0][0],self.con2[1][1]-min_w['M1']],self.con2[1]])
        # MOMDMY Layer
            momdmy_shape = gdspy.Rectangle((-GRID, self.W_CON), (self.con2[1][0]+GRID, self.con2[0][1]), layer['MOMDMY'], datatype=int(metal[1]))
            self.cell.add(momdmy_shape)
        #self.flatten()

    def t2_layer(self):
        for metal in self.metal:
        # DMEXCL Layer
            dmexcl_shape = gdspy.Rectangle((0, 0), self.con2[1], layer['DMEXCL'], datatype=int(metal[1]))
            self.cell.add(dmexcl_shape)
    # MOMDMY test0 and dummy8 layer for LVS
    # Layer datatype hard encoded
        shape = gdspy.Rectangle((-GRID, self.W_CON), (self.con2[1][0]+GRID, self.con2[0][1]), layer['MOMDMY'], datatype=100)
        self.cell.add(shape)
        shape = gdspy.Rectangle((-GRID, self.W_CON), (self.con2[1][0]+GRID, self.con2[0][1]), layer['MOMDMY'], datatype=27)
        self.cell.add(shape)
        #self.flatten()

    def lvs_layer(self):
        # Substrate contact
        temp = sp['CO']['CO']
        sp['CO']['CO'] = 0.11
        self.sub_x1 = basic.legal_coord((-SP_SUB-W_SUB,0),(0,0),1)[0]
        self.sub_x2 = basic.legal_coord((self.con2[1][0]+SP_SUB,0),(0,0),2)[0]
        #self.sub_x1 = -SP_SUB-W_SUB
        #self.sub_x2 = self.con2[1][0]+SP_SUB
        for i in [1]:
            sub_cell = basic.metal_vert(W_SUB, self.con2[1][1],lay=i)
            sub_cell1 = gdspy.CellReference(sub_cell, (self.sub_x1, 0))
            sub_cell2 = gdspy.CellReference(sub_cell, (self.sub_x2, 0))
            self.cell.add(sub_cell1)
            self.cell.add(sub_cell2)
        sp['CO']['CO'] = temp
        self.bulk.add_shape('M6', sub_cell1.get_bounding_box())
        self.bulk.add_shape('M6', sub_cell2.get_bounding_box())
        # PO Dummy Layer
        podmy_shape = gdspy.Rectangle((self.sub_x1, 0), (self.sub_x2+W_SUB, self.con2[1][1]), layer['PO'], datatype=7)
        self.cell.add(podmy_shape)
        # OD25 Layer
        self.od25_p1 = [self.sub_x1-en['OD']['PO'], -en['OD']['PO']]
        self.od25_p2 = [self.sub_x2+en['OD']['PO']+W_SUB, self.con2[1][1]+en['OD']['PO']]
        od25_shape = gdspy.Rectangle(self.od25_p1, self.od25_p2, layer['OD_25'])
        self.cell.add(od25_shape)
        # MOMDMY test0 
        shape = gdspy.Rectangle((-GRID, self.W_CON), (self.con2[1][0]+GRID, self.con2[0][1]), layer['MOMDMY'], datatype=100)
        self.cell.add(shape)
        # MOMDMY dummy2
        shape = gdspy.Rectangle(self.od25_p1, self.od25_p2, layer['MOMDMY'], datatype=21)
        self.cell.add(shape)
        # DMEXCL Layer
        for metal in self.metal:
            dmexcl_shape = gdspy.Rectangle((-GRID, self.W_CON), (self.con2[1][0]+GRID, self.con2[0][1]), layer['DMEXCL'], datatype=int(metal[1]))
            self.cell.add(dmexcl_shape)
        # OD/PO BLK
        shape = gdspy.Rectangle(self.od25_p1, self.od25_p2, layer['DMEXCL'], datatype=20)
        self.cell.add(shape)
        shape = gdspy.Rectangle(self.od25_p1, self.od25_p2, layer['DMEXCL'], datatype=21)
        self.cell.add(shape)
        #self.flatten()


    def via_layer(self):
        width = self.con1[1][0]
        via_count = int((width - 2 * EN_VIA + SP_VIA)/(W_VIA + SP_VIA))
        via_offset = (width - W_VIA - (SP_VIA + W_VIA) * (via_count - 1))/2
        delta = 0.5 * (SP_VIA + W_VIA)
        for via in self.via:
            via_cell = gdspy.Cell("VIA", True)
            via_layer = layer['VIA'+str(via)]
            via_shape = gdspy.Rectangle((0, self.VIA_OFF), (W_VIA, self.VIA_OFF+W_VIA), via_layer)
            via_cell.add(via_shape)
            if (via - self.m_bot) % 2 == 0:
                via_array1 = gdspy.CellArray(via_cell, via_count, 1, [SP_VIA+W_VIA, 0], (via_offset, 0))
                via_array2 = gdspy.CellArray(via_cell, via_count, 1, [SP_VIA+W_VIA, 0], (via_offset, self.con2[0][1]))
            else:
                via_array1 = gdspy.CellArray(via_cell, via_count-1, 1, [SP_VIA+W_VIA, 0], (via_offset+delta, 0))
                via_array2 = gdspy.CellArray(via_cell, via_count-1, 1, [SP_VIA+W_VIA, 0], (via_offset+delta, self.con2[0][1]))
            self.cell.add(via_array1)
            self.cell.add(via_array2)

class unit_cap(custom_cap):
    def __init__(self, name, w=0.07, sp=0.07, nf=8, l=1, m_bot=3, m_top=5, attr=[], f_tip=0.14, column=True):
        super(unit_cap, self).__init__(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, False)
        self.column = column
        self.extend_bulk()

    def extend_bulk(self):
        # Extended bulk connection
        sub_cell1 = gdspy.Rectangle((self.sub_x1, self.od25_p1[1]),(self.sub_x1+W_SUB, self.od25_p2[1]), layer['M1'])
        sub_cell2 = gdspy.Rectangle((self.sub_x2, self.od25_p1[1]),(self.sub_x2+W_SUB, self.od25_p2[1]), layer['M1'])
        self.cell.add(sub_cell1)
        self.cell.add(sub_cell2)
    
    def bottom_plate(self, width=0.15, lay=4, dummy=False):
        assert lay >= self.m_bot and lay <= self.m_top
        if self.column:
            self.botx = self.sub_x1 + width + W_SUB
            vert_con = gdspy.Rectangle((self.sub_x1+width+W_SUB, self.od25_p1[1]), (self.sub_x1+2*width+W_SUB, self.od25_p2[1]), layer['M'+str(lay)])
            self.cell.add(vert_con)
            if not dummy:
                hori_con = gdspy.Rectangle((self.sub_x1+width+W_SUB, self.con1[0][1]), (self.con1[0][0],self.con1[1][1]), layer['M'+str(lay)])
                self.cell.add(hori_con)
        else:
            # dummy not implemented for column=False
            hori_con = gdspy.Rectangle((self.od25_p1[0],self.con1[0][1]), (self.od25_p2[0],self.con1[1][1]), layer['M'+str(lay)])
            self.cell.add(hori_con)

    def top_ext(self):
        sub_cell = basic.metal_hori(self.con2[1][0]-self.con2[0][0], self.con2[1][1]-self.con2[0][1],lay=self.m_top+1)
        sub_cell2 = gdspy.CellReference(sub_cell, (self.con2[0][0], self.con2[0][1]))
        self.cell.add(sub_cell2)
        ext_dist = 2 * self.od25_p2[0] - self.od25_p1[0]
        hori_ext = gdspy.Rectangle((self.con2[0][0], self.con2[0][1]), (ext_dist, self.con2[1][1]), layer['M'+str(self.m_top+1)])
        self.cell.add(hori_ext)
        
    def top_plate(self, width=0.15, lay=4, dummy=False):
        assert lay >= self.m_bot and lay <= self.m_top
        if self.column:
            self.topx = self.sub_x2 - 2 * width
            vert_con = gdspy.Rectangle((self.sub_x2-2*width, self.od25_p1[1]), (self.sub_x2-width, self.od25_p2[1]), layer['M'+str(lay)])
            self.cell.add(vert_con)
            if not dummy:
                hori_con = gdspy.Rectangle((self.con2[1][0], self.con2[0][1]), (self.sub_x2-2*width,self.con2[1][1]), layer['M'+str(lay)])
                self.cell.add(hori_con)
        else:
            # dummy not implemented for column=False
            hori_con = gdspy.Rectangle((self.od25_p1[0],self.con2[0][1]), (self.od25_p2[0],self.con2[1][1]), layer['M'+str(lay)])
            self.cell.add(hori_con)

class cap_array:
    def __init__(self, name, w=0.07, sp=0.07, nf=8, l=1, m_bot=3, m_top=5, attr=[], f_tip=0.14, column=True, route_width=0.15):
        # Params
        self.route_width = route_width
        self.cell = gdspy.Cell(name, True)
        self.origin = [0,0]
        self.name = name
        # Array
        self.array_cell = unit_cap(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, column)
        self.array_cell.bottom_plate(route_width)
        self.array_cell.top_plate(route_width)
        self.array_cell.cell.flatten()
        # Column
        self.column_cell = unit_cap(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, column)
        self.column_cell.bottom_plate(route_width)
        self.column_cell.top_ext()
        self.column_cell.cell.flatten()
        # Peripheral dummies
        self.dummy_cell_t = unit_cap(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, column)
        self.dummy_cell_t.top_plate(route_width, dummy=True)
        self.dummy_cell_t.cell.flatten()
        self.dummy_cell_b = unit_cap(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, column)
        self.dummy_cell_b.bottom_plate(route_width, dummy=True)
        self.dummy_cell_b.cell.flatten()
        self.dummy_cell = unit_cap(name, w, sp, nf, l, m_bot, m_top, attr, f_tip, column)
        self.dummy_cell.cell.flatten()

    def init_array(self, bit=8, col_bit=4, ref=True, x_overlap=True):
        # Generate array of self.row * self.column
        # Currently ref column only True works
        # Implement Cell Array
        self.col_bit = col_bit
        self.row_bit = bit - col_bit
        self.column = 2**col_bit
        self.row = 2**(bit - col_bit) - 1  # Minus 1 because of column_cell
        if ref:
            self.row = self.row + 1
            self.row_bit = self.row_bit + 1
        bound = self.array_cell.cell.get_bounding_box()
        self.x_dist = bound[1][0] - bound[0][0]
        self.y_dist = bound[1][1] - bound[0][1]
        if x_overlap:
            self.overlap_offset = 2 * (self.array_cell.sub_x1 - self.array_cell.od25_p1[0]) + min_w['M1']
            self.x_dist = self.x_dist - self.overlap_offset
        self.dist = [self.x_dist, self.y_dist]
        self.array = gdspy.CellArray(self.array_cell.cell, self.row, self.column, self.dist)
        self.col = gdspy.CellArray(self.column_cell.cell, 1, self.column, self.dist, [self.x_dist*self.row, 0])
        self.dummy_hori_1 = gdspy.CellArray(self.dummy_cell_t.cell, self.row, 1, self.dist, [0, -self.y_dist])
        self.dummy_hori_2 = gdspy.CellArray(self.dummy_cell_b.cell, self.row+1, 1, self.dist, [0, self.y_dist*self.column])
        self.dummy_hori_3 = gdspy.CellReference(self.dummy_cell.cell, (self.x_dist*self.row ,-self.y_dist))
        self.dummy_vert_1 = gdspy.CellArray(self.dummy_cell.cell, 1, self.column+2, self.dist, [-self.x_dist, -self.y_dist])
        self.dummy_vert_2 = gdspy.CellArray(self.dummy_cell.cell, 1, self.column+2, self.dist, [self.x_dist*(self.row+1), -self.y_dist])
        self.cell.add(self.array)
        self.cell.add(self.col)
        self.cell.add(self.dummy_vert_1)
        self.cell.add(self.dummy_vert_2)
        self.cell.add(self.dummy_hori_1)
        self.cell.add(self.dummy_hori_2)
        self.cell.add(self.dummy_hori_3)
        # Calculate crutial coordinates
        self.sub_x = []
        self.top_x = []
        self.top_y = []
        self.bot_x = []
        for i in range(self.column):
            self.top_y.append(self.array_cell.con2[0][1] + i * self.y_dist)
        for i in range(-1, self.row+3):
            self.sub_x.append(self.array_cell.sub_x1 + i * self.x_dist)
        for i in range(self.row):
            self.top_x.append(self.array_cell.topx + i * self.x_dist)
        for i in range(self.row+1):
            self.bot_x.append(self.array_cell.botx + i * self.x_dist)
        self.cell.flatten()
        self.bound = self.cell.get_bounding_box()

    def connect_bot(self, width=0.2, lay=4):
        self.ref_p1 = [self.bound[0][0], self.bound[1][1]+width]
        self.ref_p1 = basic.legal_coord(self.ref_p1, self.origin, 3)
        self.ref_p2 = [self.bot_x[-1]+self.route_width, self.ref_p1[1]+width]
        ref_shape = gdspy.Rectangle(self.ref_p1, self.ref_p2, layer['M'+str(lay)])
        self.cell.add(ref_shape)
        for x in self.bot_x:
            con_shape = gdspy.Rectangle([x, self.bound[1][1]], [x+self.route_width, self.ref_p1[1]], layer['M'+str(lay)])
            self.cell.add(con_shape)

    def connect_sub(self, width=0.3):
        # Connect substrate on top of array
        self.sub_p1 = [self.sub_x[0], self.ref_p2[1]+2*min_w['M1']]
        self.sub_p2 = [self.sub_x[-1]+W_SUB, self.sub_p1[1]+width]
        od_shape = gdspy.Rectangle(self.sub_p1, self.sub_p2, layer['OD'])
        self.cell.add(od_shape)
        m_hori = basic.metal_hori(self.sub_p2[0]-self.sub_p1[0], width)
        m_hori = gdspy.CellReference(m_hori, (self.sub_p1[0],self.sub_p1[1]))
        self.cell.add(m_hori)
        for x in self.sub_x:
            con_shape = gdspy.Rectangle((x, self.bound[0][1]), (x+W_SUB,self.sub_p1[1]), layer['M1'])
            self.cell.add(con_shape)
        nw_shape = gdspy.Rectangle(self.bound[0], (self.bound[1][0]-self.overlap_offset, self.sub_p2[1]), layer['NW'])
        np_shape = gdspy.Rectangle(self.bound[0], (self.bound[1][0]-self.overlap_offset, self.sub_p2[1]), layer['NP'])
        self.cell.add(nw_shape)
        self.cell.add(np_shape)

    def interdigitation_code(self, bit):
        seed = [[0]]
        for i in range(1,bit-1):
            new_seed = list(range(1,2**i,2))
            for j in range(len(seed)):
                for k in range(len(seed[j])):
                    seed[j][k] = seed[j][k] * 2 
            seed.append(new_seed)
        ref = False
        tot = 2 ** (bit - 2) - 1
        set_coord = []
        for s in seed:
            if not ref:
                set_coord.append([tot-s[0]])
                set_coord.append([tot+s[0]+1])
                ref = True
            else:
                new_set = []
                for n in s:
                    new_set.append(tot-n)
                for n in s:
                    new_set.append(tot+n+1)
                set_coord.append(new_set)
        return set_coord

    def connect_top_1(self, width=0.2, con_l=0.6, lay=4):
        # MSB bit below the array
        self.row_p1 = [self.bound[0][0], self.bound[0][1]-2*width]
        self.row_p1 = basic.legal_coord(self.row_p1, self.origin, 3)
        self.row_p2 = [self.top_x[-1]+self.route_width+2*width, self.row_p1[1]+width]
        row_shape = gdspy.Rectangle(self.row_p1, self.row_p2, layer['M'+str(lay-1)])
        self.cell.add(row_shape)
        for i in range(1,self.row_bit):
            row_shape = gdspy.Rectangle([self.row_p1[0],self.row_p1[1]-2*i*width], [self.row_p2[0],self.row_p2[1]-2*i*width], layer['M'+str(lay-1)])
            self.cell.add(row_shape)
        for x in self.top_x:
            con_shape = gdspy.Rectangle([x, self.row_p1[1]-2*(self.row_bit-1)*width-width], [x+self.route_width, self.bound[0][1]], layer['M'+str(lay)])
            self.cell.add(con_shape)
        # Connect according to interdigitation
        row_coord = self.interdigitation_code(self.row_bit)
        hori_cell = basic.metal_hori(con_l,width,lay)
        bit = 0
        diff = 0.5 * con_l - 0.5 * self.route_width
        for coord in row_coord:
            for x in coord:
                temp_con = gdspy.CellReference(hori_cell, (self.top_x[x]-diff, self.row_p1[1]-bit*2*width))
                self.cell.add(temp_con)
            bit = bit + 1

    def connect_top_2(self, width=0.2, con_l=0.6):
        # LSB bit on the right of array
        self.col_p1 = [self.bound[1][0]+width, self.top_y[0]-2*width]
        self.col_p1 = basic.legal_coord(self.col_p1, self.origin, 2)
        self.col_p2 = [self.col_p1[0]+width, self.bound[1][1]]
        col_shape = gdspy.Rectangle(self.col_p1, self.col_p2, layer['M'+str(self.array_cell.m_top)])
        self.cell.add(col_shape)
        for i in range(1, self.col_bit+1):
            col_shape = gdspy.Rectangle([self.col_p1[0]+2*i*width, self.col_p1[1]], [self.col_p2[0]+2*i*width, self.col_p2[1]], layer['M'+str(self.array_cell.m_top)]) 
            self.cell.add(col_shape)
        for y in self.top_y:
            con_shape = gdspy.Rectangle([self.bound[1][0], y], [self.col_p2[0]+2*self.col_bit*width+2*width, y+self.array_cell.con2[1][1]-self.array_cell.con2[0][1]], layer['M'+str(self.array_cell.m_top+1)])
            self.cell.add(con_shape)
        # Connect according to interdigitation
        col_coord = self.interdigitation_code(self.col_bit+1)
        vert_cell = basic.metal_vert(width,con_l,self.array_cell.m_top+1)
        bit = 0
        diff = 0.5 * con_l - 0.5 * (self.array_cell.con2[1][1] - self.array_cell.con2[0][1])
        for coord in col_coord:
            for y in coord:
                temp_con = gdspy.CellReference(vert_cell, (self.col_p1[0]+2*bit*width, self.top_y[y]-diff))
                self.cell.add(temp_con)
            bit = bit + 1

    def flatten(self):
        if self.origin:
            self.origin = [self.origin[0] + 0.5*min_w['M1'], self.origin[1] + 0.5 * min_w['M1']]
            temp = gdspy.CellReference(self.cell, (-self.origin[0],-self.origin[1]))
            self.cell = gdspy.Cell(self.name, True)
            self.cell.add(temp)
            #self.plus.adjust(self.origin)
            #self.minus.adjust(self.origin)
            #self.bulk.adjust(self.origin)
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

