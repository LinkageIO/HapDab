
def test_VarDab_shape(IncreasingAF,IncreasingAFSwapRef):
    assert IncreasingAF.shape == (11,5)
    assert IncreasingAFSwapRef.shape == (11,5)


def test_IncreasingAF_Conforms(IncreasingAF,IncreasingAFSwapRef):
    assert all(IncreasingAF.genotypes(as_dataframe=True) == IncreasingAFSwapRef.genotypes(as_dataframe=True))

def test_MissingGenotypes(MissingGenotypes):
    assert True
