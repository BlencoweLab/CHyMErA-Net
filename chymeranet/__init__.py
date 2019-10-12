import sys
from pkg_resources import resource_filename
import argparse
import fileinput
import logging
from tqdm import tqdm

import pandas as pd

from Bio import SeqIO

import joblib
import keras.models

from chymeranet.version import __version__
from chymeranet import utils
from chymeranet import aux

logging.basicConfig(format='%(asctime)s - %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger(__name__)

_SCALER = resource_filename(__name__, "data/aux_scaler.human_mouse.full.pkl")
_MODEL = resource_filename(__name__, "data/CNN.human_mouse.full.h5")


def getoptions(args=None):
    desc = "Deep learning framework for predicting efficient Cas12a guides"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fasta_file', metavar='FASTA', nargs=1,
                        help="Input FASTA file")
    parser.add_argument('-s', '--seq', action='store_true',
                        help="Print guide sequence in output [%(default)s]")
    parser.add_argument('-a', '--aux', action='store_true',
                        help="Print MFE and melting temperatures in output [%(default)s]")
    parser.add_argument('-b', '--batch_size', default=100, type=int,
                        help="Batch size for prediction [%(default)s]")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args(args=args)

    return args

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Source: http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def parse_sequences(fasta_file):
    """Load FASTA sequence and return SeqRecord iterator
    """
    fin = fileinput.input(fasta_file, openhook=fileinput.hook_compressed)
    return SeqIO.parse(fin, 'fasta')


def preprocess(seqobj):
    """ Perform pre-processing steps of a given guide sequence

    This involves:
        - Removing guide examples containing an N
        - Computing secondary structure
        - Computing melting temperatures
    """

    # Remove examples containing an N
    if utils.has_N(str(seqobj.seq)):
        logger.debug("Sequence %s contains N. Skipping." % seqobj.id)
        return None

    # RNAfold
    mfe = aux.rnafold(str(seqobj.seq))

    # Compute melting temperatures
    mt = aux.melting_temperature(seqobj.seq)

    # Add ID and sequence
    meta = pd.Series([seqobj.id, str(seqobj.seq)], index=["ID", "Seq"])

    return pd.concat([meta, mt, mfe])


def main(args=None):
    args = getoptions(args)

    logger.info("Loading FASTA")
    seq_iter = parse_sequences(args.fasta_file)

    # Pre-processing
    logger.info("Computing auxiliary features")
    aux_features = []
    for seqobj in tqdm(seq_iter):
        aux_features.append(preprocess(seqobj))

    # Join results
    logger.info("Joining preprocessed data")
    candidates = pd.concat([x for x in aux_features if x is not None], axis=1).transpose()

    assert all(candidates.columns == ['ID', 'Seq', 'Tm_NN.1_23', 'Tm_NN.1_6',
        'Tm_NN.7_18', 'Tm_NN.19_23', 'MFE'])

    # Load StandardScaler
    logger.info("Loading standard scaler")
    scaler = joblib.load(_SCALER)

    # Prepare data for prediction
    logger.info("Preparing data for prediction")
    (cnn_X, aux_std) = utils.encode_data(candidates, scaler)

    # Load CNN model
    logger.info("Loading the ChymeraNet CNN")
    classifier = keras.models.load_model(_MODEL)

    # Predict
    logger.info("Predicting with ChymeraNet CNN")
    if not args.seq:
        candidates = candidates.drop(['Seq'], axis=1)
    if not args.aux:
        candidates = candidates.drop(candidates.columns[1:], axis=1)
    candidates["CNN"] = classifier.predict([cnn_X, aux_std],
                                           batch_size=args.batch_size)

    # Output
    logger.info("Writing output")
    candidates.to_csv(sys.stdout, sep="\t", index=False)

    logger.info("Done!")
