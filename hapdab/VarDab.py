from minus80 import Freezable
from itertools import chain,repeat

from hapdab.RawFile import RawFile 

import numpy as np
import pandas as pd
import time as time

from locuspocus import Locus,Loci

class VarDab(Freezable):

    def __init__(self,name):
        super().__init__(name)
        self.loci = Loci(name)     
        self._initialize_tables()


    def conform(self, Fasta):
        '''
            Conform the Ref and Alt genotypes to a reference
            Fasta Object
        '''

    def genotypes(self,samples=None,variants=None,as_dataframe=True):
        '''
            Returns genotypes from a list of samples
            and a list of loci
        
            Parameters
            ----------
            samples : iterable of sample ids

            variants : an iterable of variants 

        '''
        query = 'SELECT * FROM sample_genotypes '
        if samples != None:
            query += 'WHERE sample in ("{}")'.format('","'.join(samples))

        data = self._db.cursor().execute(
            query
        )
        if as_dataframe == True:
            data = pd.DataFrame(
                data.fetchall(),
                columns = ['chrom','pos','alt','ref','sample','dosage','flag']
            )
            data['flag'] = [bin(x) for x in data['flag']]
        return data

    @property
    def samples(self):
        return [ x[0] for x in self._db.cursor().execute(
            "SELECT sample FROM samples ORDER BY sample"        
        ).fetchall()]

    @property
    def num_accessions(self):
        return self._db.cursor().execute(
            "SELECT COUNT(*) FROM samples"        
        ).fetchone()[0]

    @property
    def num_snps(self):
        return self._db.cursor().execute(
            "SELECT COUNT(*) FROM variants" 
        ).fetchone()[0]

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
            if item[0] == 'FORMAT' and item[1] == 'ID': 
                cur.execute('INSERT OR IGNORE INTO format (format) VALUES (?)',(item[2],))
            elif item[0] == 'INFO' and item[1] == 'ID':
                cur.execute('INSERT OR IGNORE INTO info (info) VALUES (?)',(item[2],))

    def _add_sample(self,samples,cur=None):
        if cur == None:
            cur = self._db.cursor()
        cur.executemany(
                'INSERT INTO samples (FILEID,sample) VALUES ({},?)'.format(file_id),
                [(x,) for x in samples ]
            )
        # Parse out the INFO and FMT fields
        cur_info = {
            key:ID for key,ID in cur.execute('SELECT info, INFOID FROM info')        
        } 
        cur_format = {
            key:ID for key,ID in cur.execute('SELECT format, FMTID FROM format')        
        } 
        cur_samples = {
            key:ID for key,ID in cur.execute(
                '''SELECT sample,SAMPLEID from samples 
                   WHERE FILEID = "{}" 
                   ORDER BY SAMPLEID'''.format(file_id)
            )        
        } 
        return cur_info,cur_format,cur_samples

    def add_VCF(self, filename):
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
            # Get a list of current variants and their ids
            cur_vars = { 
                (chrom,pos):(VID,id) \
                    for VID,chrom,pos,id \
                    in cur.execute(
                        'SELECT VARIANTID, chrom, pos, id FROM variants'
                    )
            }
            # Iterate over the file and build the pieces of the database
            with RawFile(filename) as IN:
                genotypes = []
                for line_id,line in enumerate(IN):
                    if line_id % 100_000 == 0 and line_id > 0:
                        cur.executemany('''
                            INSERT INTO genotypes VALUES (?,?,?,?,?) 
                        ''',genotypes)
                        elapsed = time.time() - start_time 
                        rate = int(100_000 / elapsed)
                        print(f"Processed {line_id} variants ({rate}/second)")
                        genotypes = []
                        start_time = time.time()
                    line = line.strip()
                    # Case 1: Header line
                    if line.startswith('##'):
                        items = self.parse_header(line)
                        self._add_header(file_id,line_id,items,cur=cur)
                    # Case 2: Sample Line
                    elif line.startswith('#CHROM'):
                        cur_samples = line.split('\t')[9:]
                    # Case 3: Genotypes
                    else:
                        chrom,pos,id,ref,alt,qual,fltr,info,fmt,*genos = line.split('\t')
                        info = [field.split('=') for field in info.split(';')]
                        variant = Locus(
                           chrom, pos,
                           id=id, ref=ref, alt=alt,
                        )
                        if variant not in self.loci:
                            self.loci.add_locus(variant)     
                        var_id = self.loci.rowid(variant.id) 
                        # Insert the observed QUAL score
                        cur.execute('''
                            INSERT INTO variant_qual (FILEID,VARIANTID,qual,filter) VALUES (?,?,?,?)
                        ''',(file_id,var_id,qual,fltr))
                        # Insert the info fields
                        # Find the Genotype Field Index
                        GT_ind = fmt.split(':').index('GT')
                        # Insert the genotypes 
                        for g,sample in zip(genos,cur_samples.values()):
                            GT = g.split(':')[GT_ind]
                            genotypes.append((
                                file_id,
                                var_id,
                                sample,
                                self._GT_to_flag(GT),
                                self._GT_to_dosage(GT)
                            ))
            cur.executemany('''
                INSERT INTO genotypes VALUES (?,?,?,?,?) 
            ''',genotypes)
            cur.execute('RELEASE SAVEPOINT add_vcf')
            elapsed = time.time() - start_time 
            rate = int(100_000 / elapsed)
            print(f"Processed {line_id} variants ({rate}/second)")

        except Exception as e:
            cur.execute('ROLLBACK TO SAVEPOINT add_vcf')
            cur.execute('RELEASE SAVEPOINT add_vcf')
            raise e

    def to_VCF(self,filename, samples=None):
        '''
            Outputs genotypes to a VCF file

            Parameters
            ----------
            filename : string
                outputs the VCF to the filename provided
        '''
        genos = []
        tab = '\t'          # Line delimiter
        cvar = None         # current var
        with open(filename,'w') as OUT:
            print('##fileformat=VCFv4.1',file=OUT) 
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format('\t'.join(self.samples)),file=OUT)
            # Iterate over the records 
            for record in self.genotypes(as_dataframe=False): 
                if len(genos) > 0 and (cvar[0],cvar[1]) != (record[0],record[1]):
                    # Emit a line when you find a new variant
                    print(
                        "{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT\t{}".format(
                            cvar[0],
                            cvar[1],
                            '.',
                            cvar[2],
                            cvar[3],
                            '\t'.join([self._dosage_to_string(x[5]) for x in genos])
                        ),
                        file=OUT
                    )
                    genos = []
                genos.append(record)
                cvar = record 
            else:
                #Prin   t out the last one
                print(  
                    "{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT\t{}".format(
                        cvar[0],
                        cvar[1],
                        cvar[2],
                        cvar[3],
                        cvar[4],
                        '\t'.join([self._dosage_to_string(x[5]) for x in genos])
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
        else:
            return '1/1'

    
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
        # Samples
        cur.execute('''
        CREATE TABLE IF NOT EXISTS samples (
            FILEID INTEGER,  -- Maps to files
            SAMPLEID INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID)
        );
        ''')
        # Variants
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            VARIANTID INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT,
            pos INT,
            id TEXT,
            ref TEXT,
            alt TEXT
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
        # Variant Info
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variant_info(
            FILEID INTEGER,
            INFOID INTEGER,
            VARIANTID INTEGER,
            value text,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(INFOID) REFERENCES info(INFOID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID)
        );
        ''')
        # FORMAT
        cur.execute('''
        CREATE TABLE IF NOT EXISTS format (
            FMTID INTEGER PRIMARY KEY AUTOINCREMENT,
            format TEXT UNIQUE
        );
        ''')
        # Genotypes
        cur.execute('''
        CREATE TABLE IF NOT EXISTS genotypes (
            FILEID INT,
            VARIANTID INT,
            SAMPLEID INT,
            flag INT,      -- bit flag containing information for the variant 
            dosage FLOAT,
            PRIMARY KEY(FILEID,VARIANTID,SAMPLEID),
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID),
            FOREIGN KEY(SAMPLEID) REFERENCES samples(SAMPLEID)
        ); 
        ''')

        # Create some views
        cur.execute(''' 
            CREATE VIEW IF NOT EXISTS sample_genotypes AS
            SELECT 
             chrom,pos,ref,alt,
             sample, dosage, flag
            FROM genotypes 
            CROSS JOIN variants on genotypes.VARIANTID = variants.VARIANTID
            CROSS JOIN samples on genotypes.SAMPLEID = samples.SAMPLEID
            ORDER BY chrom, pos, sample
        ''')
