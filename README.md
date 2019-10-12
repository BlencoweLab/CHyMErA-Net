# ChymeraNet

A deep learning framework for predicting efficient Cas12a guides. This software accompanies the paper: "Systematic multi-site genomic editing using a hybrid CRISPR-Cas platform reveals gene paralog interactions and essential alternative exons" by Thomas Gonatopoulos-Pournatzis, Michael Aregger, Kevin R. Brown, Ulrich Braunschweig, Shaghayegh Farhangmehr, Henry Ward, Kevin C. H. Ha, Alexander Weiss, Tanja Durbic, Chad Myers, Benjamin J. Blencowe and Jason Moffat.

ChymeraNet is a Python software package that scores the efficacy of 39-nt Cas12a guide sequences (6 nt flanking upstream + 4 nt PAM + 23 nt guide + 6 nt flanking downstream).  The model is a convolutional neural network trained using Keras (2.2.4) and TensorFlow (1.13.1), which takes as input the guide sequence as well as computed auxiliary features: secondary structure (minimum free energy) and melting temperatures.

## Installation

The recommended installation method is via the conda package manager.  ChymeraNet requires Python 3.5 or higher.

1. Install [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)

2. Create a virtual environment using the provided `environment.yml` file, which contains the pinned dependencies used.

        conda env create -n cas12a -f environment.yml

3. Activate the virtual environment
   
        conda activate cas12a

4. Clone or download the ChymeraNet repository:

        git clone https://github.com/BlencoweLab/ChymeraNet.git
        cd ChymeraNet

5. Install the package

        python setup.py install


## Manual installation

If the `environment.yml` is not used (step 2 above), then the
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) package must be installed
manually:

        conda install -c bioconda viennarna

Then, install ChymeraNet as described in steps 4-5 above.
    

## Data files

### Input

A FASTA file with 39 nt guide sequences. Each guide sequence must consist of:

 - 6 nt upstream flanking sequence
 - 4 nt PAM sequence
 - 23 nt **guide sequence**
 - 6 nt downstream flanking sequence
  
Sequences cannot contain N characters.

## Usage

    chymeranet guides.fasta > scores.txt

The default output is a two column table:

    1. ID of the guide sequence (taken from the FASTA file)
    2. Prediction score between 0 and 1, where 0 is not effective and 1 is highly effective.

For additional options, open the help message by calling `chymeranet -h`.

### Model files

The trained models are saved in the `chymeranet/data` directory: 

1. The CNN Keras model
2. Scikit-learn scaling model: contains the scaling factors used in the
       training data

## License

