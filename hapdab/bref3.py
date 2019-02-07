import struct

class bref3:
    def __init__(self, filename):
        self.stream = open(filename, 'rb')
        self.snvPerms = ['A','C','G','T']

    def readRecords(self):
        # read the magin number
        if self.read_int() != 2055763188:
            raise ValueError('file is not in bref3 format')
        program = self.read_utf()
        samples = self.read_string_array()
        nHaps = 2*len(samples)

        recList = []
        nRecs = self.read_int()
        while(nRecs != 0):
            self.readDataBlock(samples, recList, nRecs)
            nRecs = self.read_int()

        return recList

    def readDataBlock(self,samples, recList, nRecs):
        chrom = self.read_utf()
        nSeqs = self.read_unsigned_short()

        hap2Seq = []
        for j in range(0,2*len(samples)):
            hap2Seq.append(self.read_unsigned_short()) 
        
        for j in range(0,nRecs):
            rec = self.readRecord(chrom,samples,nSeqs,hap2Seq)
        breakpoint()


    def read_string_array(self):
        length = self.read_int()
        entries = [self.read_utf() for _ in range(length)]
        return entries

    def read_boolean(self):
        return struct.unpack('?', self.stream.read(1))[0]

    def read_byte(self):
        return struct.unpack('b', self.stream.read(1))[0]

    def read_unsigned_byte(self):
        return struct.unpack('B', self.stream.read(1))[0]

    def read_char(self):
        return chr(struct.unpack('>H', self.stream.read(2))[0])

    def read_double(self):
        return struct.unpack('>d', self.stream.read(8))[0]

    def read_float(self):
        return struct.unpack('>f', self.stream.read(4))[0]

    def read_short(self):
        return struct.unpack('>h', self.stream.read(2))[0]

    def read_unsigned_short(self):
        return struct.unpack('>H', self.stream.read(2))[0]

    def read_long(self):
        return struct.unpack('>q', self.stream.read(8))[0]

    def read_utf(self):
        utf_length = struct.unpack('>H', self.stream.read(2))[0]
        return self.stream.read(utf_length).decode('utf-8')

    def read_int(self):
        return struct.unpack('>i', self.stream.read(4))[0]
