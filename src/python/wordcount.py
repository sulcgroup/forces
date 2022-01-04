"""
Word counting/dimer bias computing module
"""

import sys
from math import log, exp
import numpy as np
from scipy import linalg
import _wordcount

# some default parameters
MAX_NUM_ITER = 100  # maximum number of iterations when computing the bias
TOLN = 0.01  # tolerance for the count of dimers
TOLF = 0.001  # tolerance for the force correction

class WordCountException(Exception):
    pass


def _alphabet_code(alphabet):
    """
    Translate alphabet name to one of the codes accepted by the C extension
    """
    la = alphabet.lower()
    if la in ('n', 'na', 'nt', 'nucleotide'):
        # nucleotide alphabet
        return b'n'
    elif la in ('p', 'pp', 'purine-pyrimidine', 'reduced nucleotide'):
        # reduced purine/pyrimidine alphabet
        return b'p'
    elif la in ('a', 'aa', 'amino acid'):
        # amino acid
        return b'a'
    elif la in ('r', 'reduced aa', 'reduced amino acid'):
        # reduced amino acid
        return b'r'
    raise WordCountException('Unknown alphabet: {}'.format(alphabet))


def _normalize_count(X):
    """
    Normalize vector or matrix X so that the each row adds up to one
    """
    S = np.sum(X, axis=-1)
    X /= S[..., np.newaxis]


def _calc_force(F1, N2, L, i, j, tolerance_n, tolerance_f, eps, max_iter):
    """
    Auxiliary function to compute the bias (force) 
    for a particular dinucleotide (i, j).
    F1 = character frequencies (normalized)
    N2 = dimer counts (raw)
    L = original sequence length
    For other parameters see DimerForce help
    """
    if eps is None:
        eps = tolerance_f / 10.
    N = len(F1)  # matrix size (equals alphabet size)alph_size
    n_obs = N2[i * N + j]
    x = 0
    try:
        for counter in range(max_iter):
            TM = np.ones((N, N))
            TM[i, j] = exp(x)
            TM = TM * F1[:, np.newaxis]
            TMp = np.ones((N, N))
            TMp[i, j] = exp(x + eps)
            TMp = TMp * F1[:, np.newaxis]
            TMm = np.ones((N, N))
            TMm[i, j] = exp(x - eps)
            TMm = TMm * F1[:, np.newaxis]
            # partition function
            Z = np.sum(np.dot(np.linalg.matrix_power(TM, L-1), F1))
            Zp = np.sum(np.dot(np.linalg.matrix_power(TMp, L-1), F1))
            Zm = np.sum(np.dot(np.linalg.matrix_power(TMm, L-1), F1))
            lZ, lZp, lZm = log(Z), log(Zp), log(Zm)
            # estimate of n
            n = (lZp - lZm) / (2 * eps)
            dn = (lZp - 2 * lZ + lZm) / eps**2
            dx = (n_obs - n) / dn
            x += dx
            if abs(n_obs - n) <= tolerance_n and abs(dx) <= tolerance_f:
                break
        else:
            # we have exceeded MAX_NUM_ITER without reaching the stable solution
            x = float('NaN')
    except:
        x = float('NaN')
    return x


def count_words(seq, word_length=2, alphabet='n', normalize=True,
    pseudocount=None):
    """
    Count words of word_length within the sequence without overlapping.
    seq = sequence
    word_length = word (n-mer) length 
    alphabet = alphabet type
    normalize = normalize raw counts so that the add up to one
    pseudocount = value to be artificially added to the each count
    
    Return a numpy array of length alphabet_size**word_length
    """
    _alphabet = _alphabet_code(alphabet)
    count_vector = _wordcount.count_words(seq, word_length, _alphabet)
    if pseudocount is not None:
        count_vector += pseudocount
    if normalize:
        _normalize_count(count_vector)
    return count_vector


def count_overlapping_words(seq, word_length=2, alphabet='n', normalize=True,
    pseudocount=None):
    """
    Count words of word_length within the sequence with overlapping.
    seq = sequence
    word_length = word (n-mer) length 
    alphabet = alphabet type
    normalize = normalize raw counts so that the add up to one
    pseudocount = value to be artificially added to the each count
    
    Return a numpy array of length alphabet_size**word_length
    """
    _alphabet = _alphabet_code(alphabet)
    count_vector = _wordcount.count_overlapping_words(seq, word_length, 
        _alphabet)
    if pseudocount is not None:
        count_vector += pseudocount
    if normalize:
        _normalize_count(count_vector)
    return count_vector


def count_sliding_overlapping_words(seq, word_length, window_size, step, 
    alphabet='n', normalize=True, pseudocount=None):
    """
    Count words of word_length with overlapping within a window moving
    with some step along the sequence.
    seq = sequence
    word_length = word (n-mer) length 
    window_size = size of the window to count n-mers in
    step = step to move window by
    alphabet = alphabet type
    normalize = normalize raw counts so that the add up to one
    pseudocount = value to be artificially added to the each count
    
    Return a numpy array of shape 
    (#(windows that fit within seq), alphabet_size**word_length);
     each row is the n-mer count within the corresponding window.
     
    Use it if you want to study the fluctuations of the character 
    frequencies aling the sequence.
    For computing mere frequencies for the whole sequence use
    count_words or count_overlapping_words. 
    """
    _alphabet = _alphabet_code(alphabet)
    count_matrix = _wordcount.count_sliding_overlapping_words(seq, word_length, 
        window_size, step, _alphabet)
    if pseudocount is not None:
        count_matrix += pseudocount
    if normalize:
        _normalize_count(count_matrix)
    return count_matrix


def diffusion_matrix(word_length=2, alphabet='n', rate=0):
    """
    Compute the diffusion kernel for the unit time interval,
    assuming the total substitution rate equal to rate
    and assuming all possible substitutions equiprobable.
    Our conventions imply alphabet_size (not alphabet_size - 1)
    different possible equiprobable substitutions, thus counting
    a trivial substitution (e. g., A->A) towards the rate.
    """
    # Build the diffusion matrix D
    _alphabet = _alphabet_code(alphabet)
    N = _wordcount.alphabet_size(_alphabet)
    subs_prob = float(rate) / N  # probability of one particular change
    D = np.zeros((N**word_length, N**word_length))
    X = [0] * (word_length + 1)  # N-ary number with (word_length+1) digits
    for k in range(N**word_length):
        D[k, k] -= word_length * rate  # substitutions at all possible sites
        for n in range(word_length):
            kk = k - X[n] * N**n  # zero out n-th N-ary digit
            for d in range(N):
                D[k, kk + d * N**n] += subs_prob
        X[0] += 1
        for n in range(word_length):
            if X[n] == N:
                X[n] = 0
                X[n+1] += 1
    # exponentiate D to get the kernel
    K = np.real(linalg.funm(D, np.exp))
    return K


def DimerForces(seq, alphabet='n', pseudocount=None, tolerance_n=TOLN, 
    tolerance_f=TOLF, eps=None, max_iter=MAX_NUM_ITER):
    """
    Compute dimer force (bias)
    seq = sequence
    alphabet = alphabet to use
    pseudocount = pseudocount to add when computing the raw counts
    tolerance_n = stop iterating when the difference between the
        observed and predicted counts N2[i, j] does not exceed 
        this threshold
    tolerance_f = stop iterating when the shift in the force value
        at the current step does not exceed this threshold
    --- note that tolerance_f and tolerance_n must both hold for
        the iteration to stop
    eps = step to use when doing numerical differentiation
    """
    L = len(seq)
    F1 = count_overlapping_words(seq, 1, alphabet, normalize=True, 
        pseudocount=pseudocount)
    alph_size = len(F1)
    N2 = count_overlapping_words(seq, 2, alphabet, normalize=False, 
        pseudocount=pseudocount)
    dnf = np.array([ _calc_force(F1, N2, L, i, j, tolerance_n, tolerance_f, eps,
        max_iter) for i in range(alph_size) for j in range(alph_size)])
    dnf = dnf.reshape((alph_size, alph_size))
    return dnf

def DimerForce(seq, i, j, alphabet='n', pseudocount=None, tolerance_n=TOLN, 
    tolerance_f=TOLF, eps=None, max_iter=MAX_NUM_ITER):
    """
    Compute dimer force (bias)
    seq = sequence
    alphabet = alphabet to use
    pseudocount = pseudocount to add when computing the raw counts
    tolerance_n = stop iterating when the difference between the
        observed and predicted counts N2[i, j] does not exceed 
        this threshold
    tolerance_f = stop iterating when the shift in the force value
        at the current step does not exceed this threshold
    --- note that tolerance_f and tolerance_n must both hold for
        the iteration to stop
    eps = step to use when doing numerical differentiation
    """
    L = len(seq)
    F1 = count_overlapping_words(seq, 1, alphabet, normalize=True, 
        pseudocount=pseudocount)
    N2 = count_overlapping_words(seq, 2, alphabet, normalize=False, 
        pseudocount=pseudocount)
    dnf =  _calc_force(F1, N2, L, i, j, tolerance_n, tolerance_f, eps,
        max_iter)
    return dnf
