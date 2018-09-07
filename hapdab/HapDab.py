#!/bin/env/python3

import os
import io
import sys
import pkg_resources
import subprocess
import shutil
import logging
import gzip

import minus80 as m80

from minus80 import Freezable, Cohort, Accession
from minus80.RawFile import RawFile
from locuspocus import Locus, RefLoci, Fasta


from .Variant import Variant


class HapDab(Freezable):
    
    # Init method
    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
                )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    def __init__(self,name,parent=None):
        super().__init__(name)
        if 'Fasta' in self._dict:
            self._fasta = Fasta(self._dict['Fasta'],parent=self)

    # --------------------------------------------
    # Class Methods
    # --------------------------------------------

    @classmethod
    def from_files(cls,name,vcffile,fastafile,skip_phase=False):
        self = cls(name)
        self._add_fasta(fastafile)
        self._add_ref_vcf(vcffile,skip_phase=skip_phase)
        return self

    # --------------------------------------------
    # Public Methods
    # --------------------------------------------
    
    def impute(self,vcf_file):
        '''
            Impute a VCF file to the density of the 
            Reference VCF
        '''
        # Get a temp file for the output
        imputed_vcf = self._tmpfile()
        # Beagle needs a file with one common samples to exclude
        #self.log.info('Finding common samples to exclude')
        #common_sample = self._tmpfile()
        #with open(vcf_file,'r') as IN:
        #    for line in IN:
        #        if line.startswith('#CHROM'):
        #            fields = line.split()
        #            self.log.info(f'Found {len(fields[9:])} samples in common')
        #            common_sample.write("\n".join(fields[9:]))
        #            common_sample.flush()
        #            break
        #        elif line.startswith('##'):
        #            continue
        #        else:
        #            raise Exception('Could not find sample IDS in VCF header')
        # Get the path of BEAGLE
        beagle_path = pkg_resources.resource_filename(
            'hapdab',
            #'include/beagle/beagle.08Jun17.d8b.jar'
            'include/beagle/beagle.03Jul18.40b.jar'
        )
        ref_vcf = os.path.join(self._basedir,'reference.vcf.gz')
        # Create a command
        cmd = (f"java -jar {beagle_path} impute=true ref={ref_vcf} "
               f"gt={vcf_file} out={imputed_vcf.name}").split()
        phaser = subprocess.run(cmd)
        if phaser.returncode != 0:
            raise ValueError("Imputation failed!")
        return imputed_vcf.name + '.vcf.gz'
        

    # --------------------------------------------
    # Private Methods
    # --------------------------------------------

    def _add_fasta(self,fasta):
        '''
            Add a reference geneome sequence to the database.
            This reference sequence will be used to sort and 
            conform the genotype (VCF) files used for imputation.
            A fasta can only be attached once to a hapdab.

            Parameters
            ----------
            fasta : str, path-like or a locuspocus.Fasta object
                Path to the fasta file

            Returns
            -------
            None if successful. See the API for accessing the
            fasta object.

            Raises
            ------
            ValueError if a fasta has already been assigned to the
            database.

        '''
        if 'Fasta' in self._dict:
            raise ValueError('A Fasta has already been assigned to this database!')
        if os.path.exists(fasta):
            f = Fasta.from_file(self._m80_name,fasta,parent=self)
        elif m80.Tools.available('Fasta',fasta):
            f = Fasta(fasta)
            os.symlink(
                f._basedir, 
                os.path.join(self._basedir,f'Fasta.{f._m80_name}')
            )
        else:
            raise ValueError(f'Unable to determine the type of fasta')
        self._fasta = f
        # and remember for next time
        self._dict['Fasta'] = f._m80_name

    def _add_ref_vcf(self, vcf_file, skip_conform=False, skip_phase=False):
        '''
            Add a reference VCF to the database.
        '''
        original_vcf = vcf_file
        if not skip_conform:
            conformed = self._conform_vcf(vcf_file)
            vcf_file = conformed.name
        if not skip_phase:
            phased = self._phase_vcf(vcf_file)
            vcf_file = phased.name+'.vcf.gz'
        # Copy the file 
        dest_path = os.path.join(self._basedir,'reference.vcf.gz')
        try:
            # peek will return if gzipped
            gzip.open(vcf_file).peek(1)
            shutil.copy(vcf_file,dest_path)
        except OSError as e:
            # gzip will throw an AttributeError excpetion if not gzipped
            self.log.info('compressing the output file')
            with open(vcf_file,'rb') as IN:
                with gzip.open(dest_path,'wb',compresslevel=6) as OUT:
                    shutil.copyfileobj(IN,OUT,length=512*(2**20))
            self.log.info('done compressing')
         

    def _phase_vcf(self,vcf_file):
        '''
            Phase a VCF file using beagle

            Parameters
            ----------
            vcf_file : str, path-like
                A path to a VCF file

            Returns
            -------
            A handle to a temporary file containing phased data
        '''
        # Get a temp file for the output
        phased_vcf = self._tmpfile()
        # Get the path of BEAGLE
        beagle_path = pkg_resources.resource_filename(
            'hapdab',
            #'include/beagle/beagle.08Jun17.d8b.jar'
            'include/beagle/beagle.03Jul18.40b.jar'
        )
        # Create a command
        cmd = f"java -jar {beagle_path} gt={vcf_file} out={phased_vcf.name}".split()
        phaser = subprocess.run(cmd)
        if phaser.returncode != 0:
            raise ValueError("Phasing failed!")
        return phased_vcf

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
        headers = list()
        variants = list()
        cur_byte = 0
        chroms = list() 
        temps  = dict()
    
        self.log.info(f"Sorting {vcf_file}")
        # Get the chromosome order
        # Iterate through the chromosome keys and open temp files
        try:
            for chrom in self._fasta.chrom_names():
                temps[chrom] = self._tmpfile(suffix=f'sorted_vcf_{chrom}') 
                chroms.append(chrom)
        except Exception as e:
            self.log.info(f'{e}: try increasing the open file limit on your system')
        # Get headers and extract positions with file byte offsets
        self.log.info(f"Reading in VCF: {vcf_file}")
        with RawFile(vcf_file) as VCF:
            for i,line in enumerate(VCF):
                if line.startswith('#CHROM'):
                    fields = line.strip().split()
                    for i in range(9,len(fields)):
                        fields[i] = fields[i] + '_ref'
                    headers.append('\t'.join(fields))
                elif line.startswith("#"):
                    headers.append(line.strip())
                else:
                    var = Variant.from_str(line)
                    #var.conform(self._fasta[var.chrom][var.pos].upper())
                    temps[var.chrom].write(str(var)+'\n')
        # close all temp files
        for key,val in temps.items():
            self.log.info(f"flushing tmp file: {key}")
            val.flush()
        self.log.info("Outputting sorted chroms")
        out = self._tmpfile(suffix='_sorted_vcf')
        # print headers
        print("\n".join(headers),file=out)
        for chrom in chroms:
            # read in that chroms bullshit
            with open(temps[chrom].name,'r') as CHROM:
                variants =  CHROM.readlines()
                # sort by position
                variants.sort(key=lambda x: int(x.split()[1]))
                self.log.info(f"printing chrom {chrom}")
                print("".join(variants),file=out,end="")
        return out
