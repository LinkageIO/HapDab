#!/bin/env python3

import numpy as np
import scipy.stats
import pandas as pd
import time as time
import logging
import pkg_resources
import sys

from collections import namedtuple
from subprocess import Popen, PIPE

from minus80 import Freezable, Cohort, Accession
from itertools import chain,repeat
from hapdab.RawFile import RawFile 

from locuspocus import Locus,Loci,Fasta

from apsw import ConstraintError

log = logging.getLogger('VarDab')
# Setup some logging information
handler = logging.StreamHandler()
formatter = logging.Formatter(
            '%(asctime)s %(name)-8s %(levelname)-5s %(message)s'
        )
handler.setFormatter(formatter)
if not len(log.handlers):
    log.addHandler(handler)
    log.setLevel(logging.INFO)



def HW_chi(n_AA,n_AB,n_BB):
    '''
        Calculates a chi-square for allele frequencies
    '''
    n = sum([n_AA,n_AB,n_BB])
    p = (2*n_AA+n_AB)/(2*n)
    q = 1-p
    exp_AA = (p**2)*n
    exp_AB = 2*p*q*n
    exp_BB = (q**2)*n
    return scipy.stats.chisquare([n_AA,n_AB,n_BB],[exp_AA,exp_AB,exp_BB],ddof=1).pvalue

class VarDab(Freezable):

    '''
        A VarDab is not your parents variant database.
    '''
    genoRecord = namedtuple(
        'genoRecord',
        ['varid','chrom','pos','ref','alt','sample','dosage','flag']
    )

    genoVCF = namedtuple(
        'genoVCF', # A genotype record from a VCF file 
        ['fileid','varid','sampleid','flag','dosage']
    )


    def __init__(self,name,fasta=None):
        '''
            Initialize a database.

            Parameters
            ----------
            name : str
                The name of the database you are initializing

            fasta : locuspocus.Fasta object
                The reference sequence for the database you are
                initializing. If you've build this database before
                this option can be left blank (None).
        '''
        super().__init__(name)
        # Attach the Loci database
        self.loci = Loci(name)     
        self._db.cursor().execute(
            'ATTACH DATABASE ? AS loci;',(self.loci._dbfilename(),)
        )
        # Attach the Cohort database
        self.cohort = Cohort(name)
        self._db.cursor().execute(
            'ATTACH DATABASE ? AS cohort;',(self.cohort._dbfilename(),)
        )
        # Attach the fasta
        if fasta is None:
            # pull the stored fasta
            try:
                self.fasta = Fasta.from_minus80(self._dict('fasta'))
            except ValueError as e:
                raise ValueError(f'Provide a valied Fasta for {name}')
        else:
            fasta.to_minus80(name)
            self._dict('fasta',name)
            self.fasta = fasta
        # Initialize the tables
        self._initialize_tables()

    def hw_disequilibrium(self, loci=None):
        '''
            Calculate Hardy Weinberg disequilibrium

            Parameters
            ----------
            loci : an iterable of locuspocus.Locus
                If not specified, will iterate over all
                loci in the Vardab.

            Returns
            -------
            A Dataframe with counts and chi square pvalues

        '''
        if loci == None:
            loci = self.loci
        cur = self._db.cursor()
        allele_counts = cur.execute('''
                SELECT id,dosage,COUNT(dosage) 
                FROM genotypes 
                JOIN loci ON VARIANTID = loci.rowid
                GROUP BY id,dosage
        ''').fetchall()
        allele_counts = pd.DataFrame(
                allele_counts,
                columns=['id','allele','count']
        )  
        allele_counts = pd.pivot_table(
            allele_counts,
            index='id',
            columns='allele',
            values='count'
        ).fillna(0)  
        allele_counts.columns = ['AA','AB','BB']
        allele_counts['pval'] = [HW_chi(x['AA'],x['AB'],x['BB']) for i,x in allele_counts.iterrows()] 
        return allele_counts

    def conform(self):
        '''
            Conform the Ref genotype dosages to a reference
            Fasta Object. This only works when the ALT allele 
            matches the FASTA. This swaps the ALT and REF allele.

            Parameters
            ----------

            Returns
            -------
            loci that were conformed

        '''
        offenders = []
        for snp in self.loci:
            # If the REF doesn't match the FASTA 
            if snp['ref'].upper() != fasta[snp.chrom][snp.start].upper(): 
                if snp['alt'].upper() == fasta[snp.chrom][snp.start].upper():
                    offenders += [snp]
        # Returns the offenders
        self._swap_locus_alt_dosage(offenders)
        log.info(f'Conforming {len(offenders)} variants')
        return offenders


    def _call_beagle_phase(self,filename=None):
        if filename == None:
            tmp = self._tmpfile(suffix='.vcf',delete=False)
            filename = tmp.name
            # write the VCF to tmp file
            log.info(f'Ouputting genotypes to VCF: {filename}...')
            self.to_VCF(filename)
            log.info("done")
        # Get a temp file for the output
        imputed_vcf = self._tmpfile(delete=False)
        # Get the path of BEAGLE
        beagle_path = pkg_resources.resource_filename(
            'hapdab',
            'include/beagle/beagle.08Jun17.d8b.jar'
        )
        # Create a command
        cmd = f"java -jar {beagle_path} gt={filename} out=out"
        try:
            log.info(f'Phasing {filename} into {imputed_vcf.name}')
            p = Popen(cmd, stdout=PIPE, stderr=sys.stderr, shell=True)
            sout = p.communicate()[0]
            p.wait()
            if p.returncode == 0:
                # Read in the results
                import ipdb; ipdb.set_trace()
                pass
            else:
                raise ValueError('BEAGLE failed to phase')
        except FileNotFoundError as e:
            raise e
        
    def genotypes(self,cohort=None,loci=None,as_dataframe=True):
        '''
            Returns genotypes from a list of cohort
            and a list of loci
        
            Parameters
            ----------
            cohort : iterable of cohort

            loci : an iterable of loci 

            Returns
            -------
            a DataFrame if as_dataframe is true
            else an iterable of genoRecords (a named tuple)

        '''
        # Build the basic query
        query = 'SELECT * FROM sample_genotypes '

        # add on options for cohort ...
        if cohort != None:
            # NOTE: This is susceptible to SQL Injection ...
            names = [x.name for x in cohort] 
            query += " WHERE name IN ('{}')".format("','".join(names))
        # ... and loci
        if loci != None:
            if 'WHERE' in query:
                query += ' AND'
            else:
                query += ' WHERE'
            ids = [x.id for x in loci]
            query += " id IN ('{}')".format("','".join(ids))
        # Execute the query
        data = self._db.cursor().execute(
            query 
        )
        # Package up the results
        if as_dataframe == True:
            data = pd.DataFrame(
                list(data),
                columns = ['ID','chrom','pos','ref','alt','sample','dosage','flag']
            )
            #data['flag'] = [bin(x) for x in data['flag']]
        else:
            data = (self.genoRecord(*x) for x in data)
        return data

    @property
    def accessions(self):
        return [ x[0] for x in self._db.cursor().execute(
            "SELECT name FROM accessions ORDER BY name"        
        ).fetchall()]
    
    @property
    def num_files(self):
        return self._db.cursor().execute(
            'SELECT COUNT(DISTINCT(filename)) FROM files'  
        ).fetchone()[0]

    def files(self):
        return [ x[0] for x in self._db.cursor().execute(
            'SELECT filename FROM files'
        ).fetchall()]

    @property
    def num_genotypes(self):
        return self._db.cursor().execute(
          'SELECT COUNT(*) FROM genotypes'       
        ).fetchone()[0]

    @property
    def num_accessions(self):
        return len(self.cohort)

    @property
    def num_snps(self):
        return len(self.loci)

    @property
    def shape(self):
        return (self.num_snps,self.num_accessions)

    def _add_header(self, file_id, line_id, items, cur=None):
        if cur == None:
            cur = self._db.cursor()
        for item in items:
            cur.execute(
                "INSERT INTO headers VALUES ({},{},?,?,?)".format(file_id,line_id),
                item
            )

    def _dump_VCF_records_to_db(self,cur,variants,genotypes,start_time):
        '''
            A convenience method to add_VCF
        '''
        # conform the variants to the reference fasta
        variants_to_conform = set()
        for var in variants:
            true_ref = self.fasta[var.chrom][var.start].upper()
            if var['ref'].upper() != true_ref:
                if var['alt'].upper() ==  true_ref:
                    log.info(f'{var.id} needed to be confmed')
                    variants_to_conform.add(var.id)
                    # Swap!
                    var['ref'],var['alt'] = var['alt'],var['ref']

        self.loci.add_loci(variants)
        GID_map = {
            x.id : self.loci.rowid(x.id) for x in variants 
        }
        # swap out the IDS for the GIDS
        for i,geno in enumerate(genotypes):
            if geno.varid in variants_to_conform:
                # first update the dosage of the local variable geno
                geno = geno._replace(dosage=2-geno.dosage)
            # replace the original with the updated version
            genotypes[i] = geno._replace(varid=GID_map[geno.varid])
        cur.executemany('''
            INSERT OR REPLACE INTO genotypes (FILEID,VARIANTID,SAMPLEID,flag,dosage) VALUES (?,?,?,?,?) 
        ''',genotypes)
        
        # Track the elapsed time for processing all this data
        elapsed = time.time() - start_time 
        var_rate = int(len(variants) / elapsed)
        geno_rate = int(len(genotypes) / elapsed)
        log.info(f"Processed {len(variants)} variants ({var_rate}/second) containing {len(genotypes)} genotypes ({geno_rate}/second)")
        # Delete the processed variants and genotypes
        del genotypes[:]
        del variants[:]


    async def async_add_VCF(self, filename):
        '''
            Asyncronous version of add_VCF
        '''
        pass

    def add_VCF(self, filename, force=False):
        '''
            Adds variants from a VCF file to the database.

            Parameters
            ----------
            filename : str
                The VCF filename
            fasta : a locuspocus Fasta object
                Alleles will be conformed according to the reference
                genotype in the fasta. This save TONS of headaches down
                the line

            Returns
            -------
            None if successful
        '''
        log.info(f'Importing genotypes from {filename}')
        with self.bulk_transaction() as cur:
            #cur.execute('PRAGMA synchronous = off')
            cur.execute('PRAGMA journal_mode = memory')
            start_time = time.time()
            # Add the filename to the Database
            try:
                cur.execute('''
                    INSERT INTO files (filename) VALUES (?)
                ''',(filename,))
            except ConstraintError as e:
                log.warn(f'{filename} was already added.')
                if force == False:
                    return
            file_id, = cur.execute('SELECT FILEID FROM files WHERE filename = ?;',(filename,)).fetchone()
            # Iterate over the file and build the pieces of the database
            with RawFile(filename) as IN:
                variants = []
                genotypes = []
                for line_id,line in enumerate(IN):
                    if len(variants) >= 100_000:
                        self._dump_VCF_records_to_db(cur,variants,genotypes,start_time)
                        start_time = time.time()
                    line = line.strip()
                    # Case 1: Header line
                    if line.startswith('##'):
                        items = self.parse_header(line)
                        self._add_header(file_id,line_id,items,cur=cur)
                    # Case 2: Sample Line
                    elif line.startswith('#CHROM'):
                        cur_samples = [Accession(name,files=[filename]) for name in line.split()[9:]]
                        self.cohort.add_accessions(cur_samples)
                        cur_samples = self.cohort.AID_mapping
                    # Case 3: Genotypes
                    else:
                        # Split the line into its parts
                        chrom,pos,id,ref,alt,qual,fltr,info,fmt,*genos = line.split()
                        info = [field.split('=') for field in info.split(';')]
                        # Make a variant
                        variant = Locus(
                           chrom, pos,
                           id=id, ref=ref, alt=alt,
                        )
                        variants.append(variant)
                        var_id = variant.id
                        # Insert the observed QUAL score
                        cur.execute('''
                            INSERT OR REPLACE INTO variant_qual (FILEID,VARIANTID,qual,filter) VALUES (?,?,?,?)
                        ''',(file_id,var_id,qual,fltr))
                        # Find the Genotype Field Index
                        GT_ind = fmt.split(':').index('GT')
                        # Insert the genotypes 
                        for g,sample in zip(genos,cur_samples.values()):
                            GT = g.split(':')[GT_ind]
                            dosage = self._GT_to_dosage(GT)
                            if not np.isnan(dosage):
                                genotypes.append(self.genoVCF(
                                    file_id,
                                    var_id,
                                    sample,
                                    self._GT_to_flag(GT),
                                    dosage
                                ))
            self._dump_VCF_records_to_db(cur,variants,genotypes,start_time)
            log.info('Import Successful.')

    def to_VCF(self,filename, cohort=None, loci=None):
        '''
            Outputs genotypes to a VCF file

            Parameters
            ----------
            filename : string
                outputs the VCF to the filename provided
        '''
        
        genos = []
        tab = '\t'          # Line delimiter
        cvar = None         # Current var
        if cohort is None:
            cohort = self.cohort
        num_accessions = len(cohort)
        with open(filename,'w') as OUT:
            print('##fileformat=VCFv4.1',file=OUT) 
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format('\t'.join([x.name for x in cohort])),file=OUT)
            # Iterate over the records 
            genotypes = self.genotypes(
                cohort=cohort,
                loci=loci,
                as_dataframe=False
            )
            for record in genotypes: 
                if len(genos) > 0 and cvar.varid != record.varid:
                    # Emit a line when you find a new variant
                    print(
                        "{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT\t{}".format(
                            cvar.chrom,
                            cvar.pos,
                            cvar.varid,
                            cvar.ref,
                            cvar.alt,
                            '\t'.join([self._dosage_to_string(cvar.dosage) for x in genos])
                        ),
                        file=OUT
                    )
                    genos = []
                genos.append(record)
                cvar = record 
            else:
                #Print out the last one
                print(  
                    "{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT\t{}".format(
                       cvar.chrom,
                       cvar.pos,
                       cvar.varid,
                       cvar.ref,
                       cvar.alt,
                       '\t'.join([self._dosage_to_string(cvar.dosage) for x in genos])
                    ),
                    file=OUT
                )

    def _drop_tables(self):
        tables = [x[0] for x in self._db.cursor().execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall() ]
        for table in tables:
            if table != 'sqlite_sequence':
                self._db.cursor().execute('''
                    DROP TABLE IF EXISTS {};
                '''.format(table))
        views = [x[0] for x in self._db.cursor().execute(
            "SELECT name FROM sqlite_master WHERE type='view'"
        ).fetchall() ]
        for table in tables:
            if table != 'sqlite_sequence':
                self._db.cursor().execute('''
                    DROP VIEW IF EXISTS {};
                '''.format(table))

    def _reset_tables(self):
        self._drop_tables()
        self._initialize_tables()

    # -----------------------------------------------------------------
    # Static Methods
    # -----------------------------------------------------------------

    @staticmethod
    def _dosage_to_string(dosage):
        if dosage == 0:
            return '0/0'
        elif dosage == 1:
            return '0/1'
        elif dosage == 2:
            return '1/1'
        elif dosage == None:
            return './.'
        else:
            raise ValueError(f'{dosage} could not be converted to a string')

    @staticmethod
    def _GT_to_flag(alleles):
        flag = 0b0000000
        # Phased
        if '|' in alleles:
            flag = flag | 0b00000001
        # Heterozygous and 1|0
        if alleles == '1|0':
            flag = flag | 0b00000010
        return int(flag)
        
    @staticmethod
    def _GT_to_dosage(alleles):
        alleles = alleles.replace('|','').replace('/','')
        if '.' in alleles:
            return np.nan
        elif alleles == '00':
            return 0.0
        elif alleles == '01' or alleles == '10':
            return 1.0
        elif alleles == '11':
            return 2.0
        else:
            return np.nan

    @staticmethod
    def _flag_is_phased(flag):
        if flag & 0b000001 == 1:
            return True
        else:
            return False

    @staticmethod
    def parse_header(header_string):
        '''
        Parses a header line into a list of tuples: (tag, key, val)
        Example
        -------
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        gets parsed into:
            [('FORMAT', 'ID', 'GQ'),
             ('FORMAT', 'Number', '1'),
             ('FORMAT', 'Type', 'Integer'),
             ('FORMAT', 'Description', '"Genotype Quality"')]

        '''
        header_string = header_string.strip('#')
        tag = None
        keys = []
        vals = []
        current = ''
        in_str = False
        str_delim = None
        for i,char in enumerate(header_string):
            if in_str:
                if (char == '"' or char == '"') and header_string[i-1] != '\\': 
                    current += char
                    vals.append(current)
                    current = ''
                    in_str = False
                else:
                    current += char
            elif char == '"' or char == "'":
                current += char
                in_str = True
            elif char == '=':
                if tag == None:
                    tag = current
                else:
                    if len(keys) == len(vals):
                        keys.append(current)
                    else:
                        raise ValueError('fuck')
                current = ''
            elif char == ',' or char == '>' or char == '<':
                if in_str:
                    current += char
                elif len(vals) == len(keys)-1:
                    vals.append(current)
                    current = ''
            else:
                current += char 
        # Clean up the loose ends
        if current != '':
            vals.append(current)
        if len(vals) == 1:
            keys.append(None)
        return [ (tag,x,y) for x,y in zip(keys,vals)]

    def _initialize_tables(self):
        cur = self._db.cursor()
        # Files -- Meta Data
        cur.execute('''
        CREATE TABLE IF NOT EXISTS files (
            FILEID INTEGER PRIMARY KEY AUTOINCREMENT,
            filename TEXT NOT NULL UNIQUE,
            added datetime DEFAULT CURRENT_TIMESTAMP
        );
        ''')
        # Headers
        cur.execute('''
        CREATE TABLE IF NOT EXISTS headers (
            FILEID INT,  -- Maps to files
            LINEID INT,  -- Records the line number
            tag TEXT,
            key TEXT DEFAULT NULL,
            val TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID)
        );
        ''')
        # Variant Quality
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variant_qual (
            FILEID INTEGER,
            VARIANTID INTEGER,
            qual REAL,
            filter TEXT,
            PRIMARY KEY(FILEID, VARIANTID),
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID)
        );
        ''')
        # Variant Info
        cur.execute('''
        CREATE TABLE IF NOT EXISTS info (
            INFOID INTEGER PRIMARY KEY AUTOINCREMENT,
            info TEXT
        );
        ''')
        # Genotypes
        cur.execute('''
        CREATE TABLE IF NOT EXISTS genotypes (
            FILEID INT,
            VARIANTID INT,
            SAMPLEID INT,  -- Accession Name
            flag INT,      -- bit flag containing information for the variant 
            dosage FLOAT,  -- ALT Dosage
            PRIMARY KEY(FILEID,VARIANTID,SAMPLEID),
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(VARIANTID) REFERENCES loci(ROWID),
            FOREIGN KEY(SAMPLEID) REFERENCES accessions(AID)
        ); 
        CREATE INDEX IF NOT EXISTS genotype_VARIANTID ON genotypes (VARIANTID);
        CREATE INDEX IF NOT EXISTS genotype_FILEID ON genotypes (FILEID);
        CREATE INDEX IF NOT EXISTS genotype_SAMPLEID ON genotypes (SAMPLEID);
        ''')

        # Create a View where variants also have ref and alt alleles
        cur.execute('''
            CREATE TEMP VIEW loci_alleles AS
            SELECT 
                loci.rowid as LID,
                loci.id, 
                chromosome, 
                start, 
                end, 
                ref.val AS rAllele, 
                alt.val AS aAllele
            FROM loci
            LEFT JOIN loci_attrs ref 
                ON loci.id = ref.id 
                AND ref.key = "ref" 
            LEFT JOIN loci_attrs alt 
                ON loci.id = alt.id 
                AND alt.key = "alt"
        ''')

        # Create some views
        cur.execute(''' 
            CREATE TEMP VIEW sample_genotypes AS
            SELECT 
                LA.LID,
                LA.chromosome,
                LA.start,
                LA.rAllele,
                LA.aAllele,
                A.name,
                G.dosage,
                G.flag
            FROM loci_alleles LA 
                CROSS JOIN accessions A 
            LEFT OUTER JOIN genotypes G 
                ON  LA.LID=VARIANTID 
                AND A.AID=SAMPLEID;
        ''')
    

    def _sqliteshell_debug(self):
        '''
            Print out the commands to recreate the database environment 
            seen within this class from the command line. 
        '''
        print(f'''
            rlwrap sqlite3 {self._dbfilename()}
            ATTACH DATABASE "{self.loci._dbfilename()}" as loci;
            ATTACH DATABASE "{self.cohort._dbfilename()}" as cohort;

            CREATE TEMP VIEW loci_alleles AS
            SELECT 
                loci.rowid as LID,
                loci.id, 
                chromosome, 
                start, 
                end, 
                ref.val AS rAllele, 
                alt.val AS aAllele
            FROM loci
            LEFT JOIN loci_attrs ref 
                ON loci.id = ref.id 
                AND ref.key = "ref" 
            LEFT JOIN loci_attrs alt 
                ON loci.id = alt.id 
                AND alt.key = "alt";

            CREATE TEMP VIEW sample_genotypes AS
            SELECT 
             loci_alleles.id,
             chromosome,
             start,
             loci_alleles.rAllele,
             loci_alleles.aAllele,
             cohort.accessions.name, 
             dosage, 
             flag
            FROM genotypes 
            CROSS JOIN loci_alleles on genotypes.VARIANTID = loci_alleles.AID
            CROSS JOIN cohort.accessions on genotypes.SAMPLEID = cohort.accessions.AID;

        '''
        )

