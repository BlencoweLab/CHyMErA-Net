import logging
import pandas as pd

import RNA
from Bio.SeqUtils import MeltingTemp as mt


logger = logging.getLogger(__name__)


def melting_temperature(seq, mt_table=mt.DNA_NN2):
    """ Compute melting temperatures of various segments of the guide

    Melting temperatures are computed using nearest neighbor thermodynamics:
    http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html
    """
    Tm = []

    # Whole sequence
    Tm.append(mt.Tm_NN(seq, nn_table=mt_table))

    # Seed 1-6
    Tm.append(mt.Tm_NN(seq[0:6], nn_table=mt_table))

    # Trunk 7-18
    Tm.append(mt.Tm_NN(seq[6:18], nn_table=mt_table))

    # Promiscuous 19-23
    Tm.append(mt.Tm_NN(seq[18:], nn_table=mt_table))

    return pd.Series(Tm, index=['Tm_NN.1_23', 'Tm_NN.1_6', 'Tm_NN.7_18',
                                'Tm_NN.19_23'])


def rnafold(seq):
    """ Compute minimum free energy of guide sequence
    """
    _, mfe = RNA.fold(seq)
    return pd.Series([mfe], index=["MFE"])
