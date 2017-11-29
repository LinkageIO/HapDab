import pytest
import os

import hapdab as dab
import locuspocus as lp

import minus80.Tools as m80Tools

@pytest.fixture(scope='module')
def ACGTFasta():
    x = lp.Fasta.from_file('data/ACGT.fasta') 
    return x

@pytest.fixture(scope='module')
def IncreasingAF(ACGTFasta):
    m80Tools.delete('PytestIncAF',force=True)
    x = dab.VarDab('PytestIncAF',ACGTFasta)
    x.add_VCF('data/IncreasingAC.vcf')
    return x

@pytest.fixture(scope='module')
def IncreasingAFSwapRef(ACGTFasta):
    m80Tools.delete('PytestIncAFSwapRef',force=True)
    x = dab.VarDab('PytestIncAFSwapRef',ACGTFasta)
    x.add_VCF('data/IncreasingACSwappedRef.vcf')
    return x
