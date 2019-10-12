import numpy as np
import pandas as pd
import re
import logging

logger = logging.getLogger(__name__)

_ALPHABET = 'ATCG'


def encode_data(X, scaler):
    '''
    Convenience function that one hot encodes sequences and scales auxiliary
    features
    '''
    # One-hot encode for cnn
    cnn_X = np.array([char_to_onehot(x) for x in X.Seq])

    # Standardize aux features
    aux_std = pd.DataFrame(scaler.transform(X.iloc[:, 2:]),
                           columns=X.columns[2:],
                           index=X.ID)

    return cnn_X, aux_std



def onehot_to_char(one_hot_seq, alphabet = _ALPHABET):
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    charseq = []
    for i in range(one_hot_seq.shape[0]):
        if sum(one_hot_seq[i]) != 1:
            raise ValueError("Element %d:%s has too many values" %
                                (i, one_hot_seq[i]))
        letter = int_to_char[np.where(one_hot_seq[i] == 1)[0][0]]
        charseq.append(letter)
    return ''.join(charseq)


def char_to_onehot(seq, alphabet = _ALPHABET):
    # https://machinelearningmastery.com/how-to-one-hot-encode-sequence-data-in-python/

    # define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))

    # integer encode input data
    integer_encoded = [char_to_int[char] for char in seq]
    onehot_encoded = []
    for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1
        onehot_encoded.append(letter)
    return onehot_encoded


def has_N(seq):
    return bool(re.search(r'N', seq))
           
