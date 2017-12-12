
def test_VarDab_shape(IncreasingAF,IncreasingAFSwapRef):
    assert IncreasingAF.shape == (11,5)
    assert IncreasingAFSwapRef.shape == (11,5)


def test_IncreasingAF_Conforms(IncreasingAF,IncreasingAFSwapRef):
    assert all(IncreasingAF.genotypes(as_dataframe=True) == IncreasingAFSwapRef.genotypes(as_dataframe=True))

def test_MissingGenotypes(MissingGenotypes):
    x = MissingGenotypes
    assert len(x.genotypes()) == len(x.cohort) * len(x.loci)

def test_MissingGenotypesCountNotCrossJoin(MissingGenotypes):
    x = MissingGenotypes
    (real,) = x._db.cursor().execute('SELECT COUNT(*) from genotypes').fetchone()
    assert real == 15
