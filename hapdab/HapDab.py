#!/bin/env/python3

import os
import sys
import pkg_resources
import subprocess
import shutil

from minus80 import Freezable, Cohort, Accession
from minus80.RawFile import RawFile
from locuspocus import Locus, RefLoci, Fasta


from .Variant import Variant


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
        self._add_ref_vcf(vcffile)
        return self

    def _add_fasta(self,fastafile):
        if 'Fasta' in self._dict:
            raise ValueError('A Fasta has already been assigned to this database!')
        f = Fasta.from_file(self._m80_name,fastafile,parent=self)
        self._fasta = f
        # and remember for next time
        self._dict['Fasta'] = self._m80_name

    def _add_ref_vcf(self,vcf_file,skip_conform=False,skip_phase=False):
        original_vcf = vcf_file
        tmp_files = set()
        if not skip_conform:
            conformed = self._conform_vcf(vcf_file)
            tmp_files.add(conformed)
            vcf_file = conformed.name
        if not skip_phase:
            phased = self._phase_vcf(vcf_file)
            tmp_files.add(phased)
            vcf_file = phased.name+'.vcf.gz'
        # Copy the file 
        dest_path = os.path.join(self._basedir,'reference.vcf.gz')
        shutil.copy(vcf_file,dest_path)
        

    def _phase_vcf(self,vcf_file):
        # Get a temp file for the output
        imputed_vcf = self._tmpfile(delete=False)
        # Get the path of BEAGLE
        beagle_path = pkg_resources.resource_filename(
            'hapdab',
            'include/beagle/beagle.08Jun17.d8b.jar'
        )
        # Create a command
        cmd = f"java -jar {beagle_path} gt={vcf_file} out={imputed_vcf.name}".split()
        phaser = subprocess.run(cmd)
        if phaser.returncode != 0:
            raise ValueError("Phasing failed!")
        return imputed_vcf

    def _conform_vcf(self,vcf_file): 
        '''
            Conforms and sorts an input vcf based on the fasta

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
                    var = Variant.from_str(line)
                    #var.conform(self._fasta[var.chrom][var.pos].upper())
                    temps[var.chrom].write(str(var)+'\n')
        # close all temp files
        for key,val in temps.items():
            log("flushing tmp file: {}",key)
            val.flush()
        log("Outputting sorted chroms")
        out = self._tmpfile(suffix='_sorted_vcf')
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
