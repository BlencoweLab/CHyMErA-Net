import pytest
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from chymeranet import utils


@pytest.fixture
def seq():
    return 'GATTACA'

@pytest.fixture
def onehotseq():
    ohs = [[0, 0, 0, 1],
           [1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 1, 0, 0],
           [1, 0, 0, 0],
           [0, 0, 1, 0],
           [1, 0, 0, 0]]
    return ohs

@pytest.fixture
def features():
    df = pd.DataFrame({
        'ID': ['A', 'B', 'C'],
        'Seq': ['GAT', 'CCA', 'GGG'],
        'F1': np.random.rand(3),
        'F2': np.random.randint(3)
        }, columns=['ID', 'Seq', 'F1', 'F2'])
    return df

def test_default_alphabet():
    # Default alphabet
    assert utils._ALPHABET == 'ATCG'

def test_onehot_to_char(onehotseq, seq):
    # Should convert one hot encoding np.array to sequence
    result = utils.onehot_to_char(np.array(onehotseq))
    assert result == seq

def test_onehot_to_char_non_default_alphabet(onehotseq):
    # Changing the alphabet changes the sequence
    result = utils.onehot_to_char(np.array(onehotseq), 'GCTA')
    assert result == 'AGCCGTG'

def test_onehot_to_char_bad_encoding():
    # one hot encoding cannot have more than one "1" in each element
    ohs = [[1, 1, 0, 0],
           [0, 1, 0, 0]]
    with pytest.raises(ValueError) as cm:
        utils.onehot_to_char(np.array(ohs))
    assert  'too many values' in str(cm.value)

def test_char_to_onehot(onehotseq, seq):
    # Should convert sequence to one hot encoding
    result = utils.char_to_onehot(seq)
    assert result == onehotseq

def test_char_to_onehot_non_default_alphabet(seq):
    # Changing the alphabet changes the encoding
    result = utils.char_to_onehot(seq, 'GCTA')
    expected = [[1, 0, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 1, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                [0, 1, 0, 0],
                [0, 0, 0, 1]]
    assert result == expected

def test_has_N():
    assert utils.has_N('G') == False
    assert utils.has_N('N') == True
    assert utils.has_N('ATN') == True
    assert utils.has_N('NTC') == True
    assert utils.has_N('NTN') == True

def test_encode_data(mocker, features):
    mock_scaler = mocker.patch.object(StandardScaler, 'transform')
    mock_scaler.return_value = np.random.rand(features.shape[0],
            features.shape[1]-2)
    sc = StandardScaler()
    ohe, aux = utils.encode_data(features, sc)
    assert isinstance(ohe, np.ndarray) == True
    assert isinstance(aux, pd.DataFrame) == True
    assert aux.columns.tolist() == ['F1', 'F2']
    assert aux.index.tolist() == ['A', 'B', 'C']
    assert mock_scaler.called is True

