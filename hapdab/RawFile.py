import gzip
import bz2

class RawFile(object):
    def __init__(self,filename):
        self.filename = filename
        if filename.endswith('.gz'):
            self.handle = gzip.open(filename,'rt')
        elif filename.endswith('bz2'):
            self.handle = bz2.open(filename,'rt')
        else:
            self.handle = open(filename,'r')
    def __enter__(self):
        return self.handle
    def __exit__(self,type,value,traceback):
        self.handle.close()
