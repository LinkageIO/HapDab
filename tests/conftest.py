import pytest
import os

import hapdab as dab
import locuspocus as lp

import minus80.Tools as m80Tools

@pytest.fixture(scope='module')
def ACGTFasta():
    x = lp.Fasta.from_file('data/ACGT.fasta') 
    return x


