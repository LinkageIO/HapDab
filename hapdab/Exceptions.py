class HapDabError(Exception):
    pass

class TriAllelicError(HapDabError):
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = message

class MissingChromosomeError(HapDabError):
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = message

class InfoKeyError(HapDabError):
    def __init__(self,expr,message=''):
        self.expr = expr
        self.message = message
