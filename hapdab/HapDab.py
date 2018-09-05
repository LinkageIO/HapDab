#!/bin/env/python3

import os
import sys

from minus80 import Freezable, Cohort, Accession
from minus80.RawFile import RawFile
from locuspocus import Locus, RefLoci, Fasta


class HapDab(Freezable):

    def __init__(self,name,parent=None):
        super().__init__(name)
        self._initialize_tables()
        if 'Fasta' in self._dict:
            self._fasta = Fasta(self._dict['Fasta'],parent=self)

    def _initialize_tables(self):
        pass


    @classmethod
    def from_files(cls,name,vcffile,fastafile):
        self = cls(name)
        self._add_fasta(fastafile)
        self._add_vcf(vcffile)
        return self

    def _add_fasta(self,fastafile):
        if 'Fasta' in self._dict:
            raise ValueError('A Fasta has already been assigned to this database!')
        f = Fasta.from_file(self._m80_name,fastafile,parent=self)
        self._fasta = f
        # and remember for next time
        self._dict['Fasta'] = self._m80_name

    def _add_vcf(self,vcf_file):
        sorted_vcf = self_sort(vcf_file) 
        # conform
        # phase
        # java -jar ../../hapdab/include/beagle/beagle.08Jun17.d8b.jar gt=2M_CHIP_1000.vcf.gz out=test.phased nthreads=30
        pass

    def _sort_vcf(self,vcf_file): 
        '''
            Sorts an input vcf based on the fasta

            Parameters
            ----------
            vcf_file : path
                The path the VCF file

            Returns
            -------
            A named temp file containing the sorted vcf
        '''
        def log(message,*formatting):
            print(message.format(*formatting),file=sys.stderr)
        headers = list()
        variants = list()
        cur_byte = 0
        chroms = list() 
        temps  = dict()
    
        log("Sorting {}",vcf_file)
        # Get the chromosome order
        # Iterate through the chromosome keys and open temp files
        try:
            for chrom in self._fasta.chrom_names():
                temps[chrom] = self._tmpfile(suffix=f'sorted_vcf_{chrom}') 
                chroms.append(chrom)
        except Exception as e:
            log('{}: try increasing the open file limit on your system',e)
        # Get headers and extract positions with file byte offsets
        log("Reading in VCF: {}",vcf_file)
        with RawFile(vcf_file) as VCF:
            for i,line in enumerate(VCF):
                if line.startswith("#"):
                    headers.append(line.strip())
                else:
                    chrom,pos,*junk = line.split()
                    temps[chrom].write(line)
        # close all temp files
        for key,val in temps.items():
            log("flushing tmp file: {}",key)
            val.flush()
        log("Outputting sorted chroms")
        out = self._tmpfile(suffix='sorted_vcf')
        # print headers
        print("\n".join(headers),file=out)
        for chrom in chroms:
            # read in that chroms bullshit
            with open(temps[chrom].name,'r') as CHROM:
                variants =  CHROM.readlines()
                # sort by position
                variants.sort(key=lambda x: int(x.split()[1]))
                log("printing chrom {}",chrom)
                print("".join(variants),file=out,end="")
        return out
