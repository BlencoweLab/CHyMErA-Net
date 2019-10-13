import logging
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import chymeranet as chy

def test_preprocess(mocker, caplog):
    # If sequence has N, return None
    seqobj = SeqRecord(Seq('AAATTTCCCN'), id='test', name='test')
    with caplog.at_level(logging.DEBUG, logger='chymeranet'):
        result = chy.preprocess(seqobj)
    assert result is None
    assert "Sequence test contains N. Skipping." in caplog.text

    seqobj = SeqRecord(Seq('AAATTTCCC'), id='test', name='test')
    mock_fold = mocker.patch('chymeranet.aux.rnafold',
            return_value=pd.Series([0]))
    mock_melt = mocker.patch('chymeranet.aux.melting_temperature',
        return_value=pd.Series([1]))
    result = chy.preprocess(seqobj)
    mock_fold.assert_called_once_with(str(seqobj.seq))
    mock_melt.assert_called_once_with(seqobj.seq)
    assert isinstance(result, pd.Series)

def test_main(mocker, caplog):
    pass

    # set args

    # Mock loading fasta

    # Mock preprocess

    # Join results

    # StandardScaler


    # Prepare for prediction


    # Predict

    # Output


    
