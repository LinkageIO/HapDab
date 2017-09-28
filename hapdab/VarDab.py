#!/bin/env python3

import numpy as np
import pandas as pd
import time as time
import logging
from collections import namedtuple

from minus80 import Freezable, Cohort, Accession
from itertools import chain,repeat
from hapdab.RawFile import RawFile 

from locuspocus import Locus,Loci

class VarDab(Freezable):

    genoRecord = namedtuple(
        'genoRecord',
        ['varid','chrom','pos','ref','alt','sample','dosage','flag']
    )

    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
                )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)


    def __init__(self,name):
        super().__init__(name)
        self.loci = Loci(name+'.VarDab')     
        self._db.cursor().execute(
            'ATTACH DATABASE ? AS loci;',(self.loci._dbfilename(),)
        )
        self.cohort = Cohort(name+'.VarDab')
        self._db.cursor().execute(
            'ATTACH DATABASE ? AS cohort;',(self.cohort._dbfilename(),)
        )
        self._initialize_tables()

    def conform(self, Fasta):
        '''
            Conform the Ref and Alt genotypes to a reference
            Fasta Object
        '''
        raise NotImplementedError()
        
    def genotypes(self,accessions=None,variants=None,as_dataframe=True):
        '''
            Returns genotypes from a list of accessions
            and a list of loci
        
            Parameters
            ----------
            accessions : iterable of accessions

            variants : an iterable of variants 

            Returns
            -------
            a DataFrame if as_dataframe is true
            else an iterable of genoRecords (a named tuple)

        '''
        # Build the basic query
        query = 'SELECT * FROM sample_genotypes '

        # add on options for accessions ...
        if accessions != None:
            # NOTE: This is susceptible to SQL Injection ...
            names = [x.name for x in accessions] 
            query += " WHERE name IN ('{}')".format("','".join(names))
        # ... and variants
        if variants != None:
            if 'WHERE' in query:
                query += ' AND'
            else:
                query += ' WHERE'
            ids = [x.id for x in variants]
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
            data['flag'] = [bin(x) for x in data['flag']]
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
            #if item[0] == 'FORMAT' and item[1] == 'ID': 
            #    cur.execute('INSERT OR IGNORE INTO format (format) VALUES (?)',(item[2],))
            #elif item[0] == 'INFO' and item[1] == 'ID':
            #    cur.execute('INSERT OR IGNORE INTO info (info) VALUES (?)',(item[2],))


    def _dump_VCF_records_to_db(self,cur,variants,genotypes,start_time):
        '''
            A convenience method to add_VCF
        '''
        self.loci.add_loci(variants)
        GID_map = {
            x.id : self.loci.rowid(x.id) for x in variants 
        }
        # swap out the IDS for the GIDS
        for geno in genotypes:
            geno[1] = GID_map[geno[1]]
        cur.executemany('''
            INSERT OR REPLACE INTO genotypes (FILEID,VARIANTID,SAMPLEID,flag,dosage) VALUES (?,?,?,?,?) 
        ''',genotypes)
        elapsed = time.time() - start_time 
        var_rate = int(len(variants) / elapsed)
        geno_rate = int(len(genotypes) / elapsed)
        self.log.info(f"Processed {len(variants)} variants ({var_rate}/second) containing {len(genotypes)} genotypes ({geno_rate}/second)")
        del genotypes[:]
        del variants[:]

    def add_VCF(self, filename):
        self.log.info(f'Importing genotypes from {filename}')
        cur = self._db.cursor()
        cur.execute('PRAGMA synchronous = off')
        cur.execute('PRAGMA journal_mode = memory')
        cur.execute('SAVEPOINT add_vcf')
        start_time = time.time()
        try:
            # Add the filename to the Database
            cur.execute('''
                INSERT OR IGNORE INTO files (filename) VALUES (?)
            ''',(filename,))
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
                        cur_samples = [Accession(name,files=[filename]) for name in line.split('\t')[9:]]
                        self.cohort.add_accessions(cur_samples)
                        cur_samples = self.cohort.AID_mapping
                    # Case 3: Genotypes
                    else:
                        # Split the line into its parts
                        chrom,pos,id,ref,alt,qual,fltr,info,fmt,*genos = line.split('\t')
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
                                genotypes.append([
                                    file_id,
                                    var_id,
                                    sample,
                                    self._GT_to_flag(GT),
                                    dosage
                                ])
            self._dump_VCF_records_to_db(cur,variants,genotypes,start_time)
            cur.execute('RELEASE SAVEPOINT add_vcf')
            self.log.info('Import Successful.')

        except Exception as e:
            cur.execute('ROLLBACK TO SAVEPOINT add_vcf')
            cur.execute('RELEASE SAVEPOINT add_vcf')
            raise e

    def to_VCF(self,filename, accessions=None):
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
        with open(filename,'w') as OUT:
            print('##fileformat=VCFv4.1',file=OUT) 
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format('\t'.join(self.accessions)),file=OUT)
            # Iterate over the records 
            for record in self.genotypes(as_dataframe=False): 
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
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID),
            FOREIGN KEY(SAMPLEID) REFERENCES accessions(AID)
        ); 
        ''')

        # Create a View where variants also have ref and alt alleles
        cur.execute('''
            CREATE TEMP VIEW loci_alleles AS
            SELECT 
                loci.rowid AID, 
                loci.id, 
                chromosome, 
                start, 
                end, 
                ref.val AS rAllele, 
                alt.val AS aAllele
            FROM loci
            LEFT OUTER JOIN loci_attrs ref 
                ON loci.id = ref.id 
                AND ref.key = "ref" 
            LEFT OUTER JOIN loci_attrs alt 
                ON loci.id = alt.id 
                AND alt.key = "alt"

        ''')

        # Create some views
        cur.execute(''' 
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
        ''')
    

