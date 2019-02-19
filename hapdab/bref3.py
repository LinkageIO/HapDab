import struct
from itertools import permutations

class bref3:
    def __init__(self, filename):
        self.stream = open(filename, 'rb')
        self.snvPerms = list(permutations(['A','C','G','T']))

    def readRecords(self):
        # read the magic number
        if self.read_int() != 2055763188:
            raise ValueError('file is not in bref3 format')
        program = self.read_utf()
        samples = self.read_string_array()
        nHaps = 2*len(samples)

        recList = []
        nRecs = self.read_int()
        print(f'Reading {nRecs} records!')
        while(nRecs != 0):
            self.readDataBlock(samples, recList, nRecs)
            nRecs = self.read_int()

        return recList

    def read_string_array(self):
        length = self.read_int()
        entries = [self.read_utf() for _ in range(length)]
        return entries

    def readByteLengthStringArray(self):
        length = self.read_unsigned_byte()
        array = []
        for j in range(0,length):
            array.append(self.read_utf())
        return array

    def readDataBlock(self,samples, recList, nRecs):
        # Chrom for all records in data block
        chrom = self.read_utf()
        # Number of distinct allele sequences in sequence coded records
        nSeqs = self.read_unsigned_short()
        print(f'On chrom {chrom}')

        # index of sequence carried by each haplotype at sequence-coded records
        hap2Seq = []
        for j in range(0,2*len(samples)):
            hap2Seq.append(self.read_unsigned_short()) 

        for j in range(0,nRecs):
            rec = self.readRecord(chrom,samples,nSeqs,hap2Seq)
            recList.append(rec)

    def readRecord(self, chrom, samples, nSeqs, hap2Seq):
        marker = self.readMarker(chrom)
        coding = self.read_byte()
        if coding == 0:
            print(f"{marker['id']}:seq coded")
            return self.readSeqCodedRecord(samples,marker,nSeqs,hap2Seq)
        elif coding == 1:
            print(f"{marker['id']}:allele coded")
            return self.readAlleleCodedRecord(samples, marker)

    def readMarker(self, chrom):
        marker = dict()
        marker['pos'] = self.read_int()
        marker['id'] = self.readByteLengthStringArray()

        alleleCode = self.read_byte()
        if alleleCode == -1:
            marker['alleles'] = self.read_string_array()
            marker['end'] = self.read_int()
        else:
            marker['nAlleles'] = 1 + (alleleCode & 0b11)
            permIndex = (alleleCode >> 2)
            marker['alleles'] = self.snvPerms[permIndex][0:marker['nAlleles']]
            marker['end'] = -1
        return marker

    def readSeqCodedRecord(self,samples,marker,nSeqs,hap2Seq):
        seq2Allele = []
        for _ in range(nSeqs):
            seq2Allele.append(self.read_unsigned_byte())
        breakpoint()
        hap2Allele = []
        for x in hap2Seq:
            hap2Allele.append(seq2Allele[x])
        record = dict()
        record['marker'] = marker
        record['samples'] = samples
        record['hap2Allele'] = hap2Allele
        return record

    def readAlleleCodedRecord(self,samples,marker):
        nHaps = 2*len(samples)
        nAlleles = len(marker['alleles'])
        hapIndices = []
        majorAllele = -1
        for j in range(0,nAlleles):
            hapIndices.append(self.readIntArray())
            if hapIndices[j] is None:
                majorAllele = j

        hap2Allele = []
        for j in range(0,nHaps):
            hap2Allele.append(majorAllele)
        for j in range(0,len(hapIndices)):
            if hapIndices[j] != None:
                for hap in hapIndices[j]:
                    hap2Allele[hap] = j
        record = dict()
        record['marker'] = marker
        record['samples'] = samples
        record['hapToAllele'] = hap2Allele
        return record

    def readIntArray(self):
        length = self.read_int()
        if length == -1:
            return None
        else:
            array = []
            for j in range(0,length):
                array.append(self.read_int())
            return array

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
