#!/bin/env/python3

from minus80 import Freezable, Cohort, Accession
from minus80.RawFile import RawFile
from locuspocus import Locus, RefLoci, Fasta


class HapDab(Freezable):

    def __init__(self,name):
        super().__init__(name)
        self._initialize_tables()

    def add_fasta(self,fasta):
        if not isinstance(fasta,Fasta):
            raise ValueError(f'The fasta must be a locuspocus Fasta object')
        self.fasta = fasta

    def add_vcf(self,vcffile):
        pass

    def _initialize_tables():
        pass


    @classmethod
    def from_files(cls,vcffile,fastafile)
