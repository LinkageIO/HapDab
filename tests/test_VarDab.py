
def test_VarDab_shape(IncreasingAF,IncreasingAFSwapRef):
    assert IncreasingAF.shape == (11,5)
    assert IncreasingAFSwapRef.shape == (11,5)


def test_IncreasingAF_Conforms(IncreasingAF,IncreasingAFSwapRef):
    x = IncreasingAF.genotypes(as_dataframe=True)
    y = IncreasingAFSwapRef.genotypes(as_dataframe=True)
    z = x == y
    z = z.all()
    z = z.all()
    assert z

def test_MissingGenotypes(MissingGenotypes):
    x = MissingGenotypes
    (nrows,ncol) = x.genotypes().shape
    ncohort = ncol - len(['chr','pos','id','alt','red'])
    nloci = nrows
    assert (ncohort == len(x.cohort)) and (nloci == len(x.loci))

def test_MissingGenotypesCountNotCrossJoin(MissingGenotypes):
    x = MissingGenotypes
    (real,) = x._db.cursor().execute('SELECT COUNT(*) from genotypes').fetchone()
    assert real == 15
