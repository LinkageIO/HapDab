#!/bin/env python3

import numpy as np
import scipy.stats
import pandas as pd
import time as time
import logging
import pkg_resources
import sys
import asyncio
import concurrent

from pysam import VariantFile
from collections import namedtuple
from subprocess import Popen, PIPE
from apsw import ConstraintError
from contextlib import contextmanager

# Linkage Imports
from minus80 import Freezable, Cohort, Accession
from itertools import chain,repeat,product
from locuspocus import Locus,RefLoci,Fasta

# Internal Imports
from .RawFile import RawFile 

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
        Calculates a chi-square for allele frequencies.

        This is the canonoical calculation for HardyWeinberg 
        disequilibrium tell you in allele frequencies are
        likely constant between generations.
    '''
    n = sum([n_AA,n_AB,n_BB])
    p = (2*n_AA+n_AB)/(2*n)
    q = 1-p
    exp_AA = (p**2)*n
    exp_AB = 2*p*q*n
    exp_BB = (q**2)*n
    return scipy.stats.chisquare(
            [n_AA,n_AB,n_BB],
            [exp_AA,exp_AB,exp_BB],
            ddof=1
           ).pvalue


def timed(method):
    def timed_dec(*args, **kwargs):
        start = time.time()
        result = method(*args,**kwargs)
        end = time.time()
        elap = end - start

        units = ['days','hours','minutes','seconds']
        days = elap // 86400
        hours = elap // 3600 % 24
        minutes = elap // 60 % 60
        seconds = elap % 60
        elap = [f"{e}:{u}" for e,u in filter(lambda x: x[0], zip([days,hours,minutes,seconds],units))]
        if len(elap) == 0:
            print('Took less than a second')
        else:
            print("Took: " + ','.join(elap))
        return result
    return timed_dec

class VarDab(Freezable):

    '''
        A VarDab is not your parents variant database.
    '''

    # Class variables holding records for various tuples generated
    # by VarDab
    genoRecord = namedtuple(
        'genoRecord',
        ['varid','chrom','pos','ref','alt','sample','dosage','flag']
    )

    VCFRecord = namedtuple(
        'VCFRecord',
        ['chrom','start','id','ref','alt']
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
        self.loci = RefLoci(name)     
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
                fasta_name = self._dict('fasta')
                self.fasta = Fasta(fasta_name)
            except ValueError as e:
                raise ValueError(f'Provide a valied Fasta for {name}')
        else:
            self._dict('fasta',fasta._m80_name)
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
        LIDS = [self.loci.LID(x) for x in loci]
        cur = self._db.cursor()
        allele_counts = cur.executemany('''
                SELECT id,dosage,COUNT(dosage) 
                FROM genotypes 
                JOIN loci ON VARIANTID = loci.rowid
                WHERE loci.rowid  = ?
                GROUP BY id,dosage
        ''',((x,) for x in LIDS)).fetchall()
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


    def _call_beagle_phase(self,loci=None,cohort=None,filename=None):
        # If not filename, them create a temp file
        if filename == None:
            tmp = self._tmpfile(suffix='.vcf')
            filename = tmp.name
            # write the VCF to tmp file
            log.info(f'Ouputting genotypes to VCF: {filename}...')
            self.to_VCF(filename,loci=loci,cohort=cohort)
            log.info("done")
        # Get a temp file for the output
        imputed_vcf = self._tmpfile(delete=False)
        # Get the path of BEAGLE
        beagle_path = pkg_resources.resource_filename(
            'hapdab',
            'include/beagle/beagle.08Jun17.d8b.jar'
        )
        # Create a command
        cmd = f"java -jar {beagle_path} gt={filename} out={imputed_vcf.name}"
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
            cohort : minus80 Cohort object

            loci : an iterable of loci objects

            Returns
            -------
            a DataFrame if as_dataframe is true
            else an iterable of genoRecords (a named tuple)

        '''
        # set default for cohort and loci
        if cohort is None:
            cohort = self.cohort
        if loci is None:
            loci = self.loci
        # Build the basic query
        cohort_ids = [x['AID'] for x in cohort]
        data = []
        for locus in loci:
            lid = self.loci.LID(locus)
            locus = self.VCFRecord(*self._db.cursor().execute('''
                SELECT chromosome,start,id,rAllele,aAllele
                FROM loci_alleles WHERE LID = ?
            ''',(lid,)).fetchone())
            genotypes = self._db.cursor().executemany('''
                SELECT 
                    CASE WHEN COUNT(dosage) > 0 
                    THEN dosage 
                    ELSE NULL 
                    END 
                AS dosage FROM genotypes 
                WHERE VARIANTID = ? 
                AND SAMPLEID = ?;
            ''',product([lid],cohort_ids))
            genos = [x[0] for x in genotypes.fetchall()]
            data.append(list(locus)+genos)
        return pd.DataFrame(data)   

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

    def _add_VCF_header(self, file_id, line_id, line, cur=None,force=False):
        items = self.parse_header(line)
        if cur == None:
            cur = self._db.cursor()
        for item in items:
            cur.execute(
                "INSERT OR REPLACE INTO headers VALUES ({},{},?,?,?)".format(file_id,line_id),
                item
            )

   
    def add_VCF(self, filename, force=False):
        '''
            Adds variants from a VCF file to the database.

            Parameters
            ----------
            filename : str
                The VCF filename

            Returns
            -------
            None if successful
        '''
        log.info(f'Importing genotypes from {filename}')
        with self.bulk_transaction() as cur, self.suspend_sqlite_indices():
            cur.execute('PRAGMA journal_mode = memory')
            # Add the filename to the Database
            FID = self._add_VCF_filename(filename,cur=cur)
            # Iterate over the file and build the pieces of the database
            variants = []
            genotypes = []
            vfile = VariantFile(filename)
            # Handle Accessions
            sample_names = list(vfile.header.samples)
            cur_samples = [Accession(name,files=[filename]) for name in sample_names]
            self.cohort.add_accessions(cur_samples)
            # Genotype
            genotypes = []
            variants = []
            for line_id,rec in enumerate(vfile):
                if len(rec.alts) > 1:
                    log.warn(f'{rec.chrom}:{rec.pos} is not bi-allelic')
                else:
                    if self.fasta[rec.chrom][rec.pos].upper() != rec.ref.upper():
                        ref = rec.alts[0]
                        alt = rec.ref
                        conform = True
                    else:
                        ref = rec.ref
                        alt = rec.alts[0]
                        conform = False
                    variant = Locus(
                       rec.chrom, rec.pos,
                       id=rec.id, ref=rec.ref, alt=rec.alts[0],
                    )
                    variants.append(variant)
                    #LID = line_id # self.loci.add_locus(variant)
                    #cur.execute('''
                    #    INSERT OR REPLACE INTO variant_qual (FILEID,VARIANTID,qual,filter) VALUES (?,?,?,?)
                    #''',(FID,LID,rec.qual,'PASS'))
                    dosages = []
                    for i,alleles in enumerate(rec.samples.values()):
                        if None in alleles.allele_indices:
                            continue
                        dose = int(sum(alleles.allele_indices))
                        if conform:
                            dose = 2 - dose
                        AID = self.cohort._get_AID(sample_names[i])
                        dosages.append((FID,AID,0,dose))
                    genotypes.append(dosages)

            # Take care of the database
            self.loci.add_loci(variants)
            LIDs = [self.loci.LID(x) for x in variants ]
            cur.executemany('''
                INSERT OR REPLACE INTO genotypes (FILEID,VARIANTID,SAMPLEID,flag,dosage) VALUES (?,?,?,?,?) 
            ''',((FID,LID,AID,flag,dose) for LID,row in zip(LIDs,genotypes) for FID,AID,flag,dose in row))

            log.info('Import Successful.')

    def _add_VCF_filename(self,filename,cur=None,force=False):
        '''
            Helper method to add_VCF. Adds a filename to the
            provided cursor.
        '''
        if cur is None:
            cur = self._db.cursor()
        try:
           cur.execute('''
                INSERT INTO files (filename) VALUES (?)
            ''',(filename,))
        except ConstraintError as e:
            log.warn(f'{filename} was already added.')
            if force == False:
                pass
        file_id, = cur.execute(
            'SELECT FILEID FROM files WHERE filename = ?;',
            (filename,)
        ).fetchone()
        return file_id

    def to_VCF(self,filename,cohort=None,loci=None):
        '''
            Outputs genotypes to a VCF file

            Parameters
            ----------
            filename : string
                outputs the VCF to the filename provided
            cohort : iterable of m80.Accessions
                An iterable of the cohort that will
                be sent to file
            loci : iterable of locuspocus.Loci
                An iterable of the loci that will be 
                sent to file
        '''
        if loci is None:
            loci = self.loci
        if cohort is None:
            cohort = self.cohort
        loci.sort()
        # Get loci ids
        async def produce_records(queue,loci,cohort):
            cohort_ids = [x['AID'] for x in cohort]
            # Iterate and produce records
            for locus in loci:
                lid = self.loci.LID(locus)
                locus = self.VCFRecord(*self._db.cursor().execute('''
                    SELECT chromosome,start,id,rAllele,aAllele
                    FROM loci_alleles WHERE LID = ?
                ''',(lid,)).fetchone())
                genotypes = self._db.cursor().executemany('''
                    SELECT 
                        CASE WHEN COUNT(dosage) > 0 
                        THEN dosage 
                        ELSE NULL 
                        END 
                    AS dosage FROM genotypes 
                    WHERE VARIANTID = ? 
                    AND SAMPLEID = ?;
                ''',product([lid],cohort_ids))
                genos = genotypes.fetchall()
                genos = [self._dosage_to_string(x[0]) for x in genos]
                await queue.put((locus,genos))
            await queue.put((None,None))

        async def consume_records(queue):
            with open(filename,mode='w') as OUT:
                print('##fileformat=VCFv4.1',file=OUT) 
                print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
                    '\t'.join([x.name for x in cohort]
                    )
                ),file=OUT)
                while True:            
                    (locus,genos) = await queue.get()
                    if locus is None:
                        break
                    print("{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT\t{}".format(
                            locus.chrom,
                            locus.start,
                            locus.id,
                            locus.ref,
                            locus.alt,
                            '\t'.join(genos)
                        ),file=OUT)

        loop = asyncio.get_event_loop()
        queue = asyncio.Queue(loop=loop)
        # Create some tasks
        prod = produce_records(queue,loci,cohort)
        cons = consume_records(queue)
        loop.run_until_complete(asyncio.gather(prod,cons))
        loop.close()

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
        if dosage == 0.0:
            return '0/0'
        elif dosage == 1.0:
            return '0/1'
        elif dosage == 2.0:
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
    def _process_GT(alleles):
        flag = 0b0000000
        # Phased
        if '|' in alleles:
            flag = flag | 0b00000001
        # Heterozygous and 1|0
        if alleles == '1|0':
            flag = flag | 0b00000010
        alleles = alleles.replace('|','').replace('/','')
        if '.' in alleles:
            dosage = np.nan
        elif alleles == '00':
            dosage = 0.0
        elif alleles == '01' or alleles == '10':
            dosage = 1.0
        elif alleles == '11':
            dosage = 2.0
        else:
            dosage = np.nan
        return (flag,dosage)

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


    @contextmanager
    def suspend_sqlite_indices(self):
        cur = self._db.cursor()
        log.info('Dropping INDEXES')
        cur.execute('''
            DROP INDEX IF EXISTS genotype_VARIANTID;
            DROP INDEX IF EXISTS genotype_SAMPLE_VARIANT;
            DROP INDEX IF EXISTS genotype_FILEID;
            DROP INDEX IF EXISTS genotype_SAMPLEID;
        ''')
        yield
        log.info('Rebuilding INDEXES')
        cur.execute('''
            CREATE INDEX IF NOT EXISTS genotype_VARIANTID ON genotypes (VARIANTID);
            CREATE INDEX IF NOT EXISTS genotype_SAMPLE_VARIANT on genotypes (VARIANTID,SAMPLEID);
            CREATE INDEX IF NOT EXISTS genotype_FILEID ON genotypes (FILEID);
            CREATE INDEX IF NOT EXISTS genotype_SAMPLEID ON genotypes (SAMPLEID);
        ''')


    def _initialize_tables(self):
        # Initialize the cassandra namespace
        cas =  self._cassandra()
        default_keyspace = f'{self._m80_name}_{self._m80_type}'
        cas.execute(f'''
            CREATE KEYSPACE {default_keyspace}
            WITH replication = {{'class':'SimpleStrategy','replication_factor':2}};
        ''')
        # Create the table for variants

        cur = self._db.cursor()
        self._schema = []
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
        CREATE INDEX IF NOT EXISTS genotype_SAMPLE_VARIANT on genotypes (VARIANTID,SAMPLEID);
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
    

