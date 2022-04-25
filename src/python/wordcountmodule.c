#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_3kcompat.h"
/* #include <numpy/arrayobject.h> */


/* enumerate the 4 nucleotides */
static Py_ssize_t nt2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
			return 0;
		case 'C':
		case 'c': 
			return 1;
		case 'G':
		case 'g': 
			return 2;
		case 'T':
		case 't': 
			return 3;
	}
	// invalid character
	return -1;
}


/* enumerate purines and pyrimidines */
static Py_ssize_t simple_nt2int(char c)
{
	switch (c) {
		// purine
		case 'A':
		case 'a': 
		case 'G':
		case 'g': 
			return 0;
		// pyrimidine
		case 'C':
		case 'c': 
		case 'T':
		case 't': 
			return 1;
	}
	// invalid character
	return -1;
}


/* enumerate the 20 standard amino acids */
static Py_ssize_t aa2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
			return 0;
		case 'C':
		case 'c': 
			return 1;
		case 'D':
		case 'd': 
			return 2;
		case 'E':
		case 'e': 
			return 3;
		case 'F':
		case 'f': 
			return 4;
		case 'G':
		case 'g': 
			return 5;
		case 'H':
		case 'h': 
			return 6;
		case 'I':
		case 'i': 
			return 7;
		case 'K':
		case 'k': 
			return 8;
		case 'L':
		case 'l': 
			return 9;
		case 'M':
		case 'm': 
			return 10;
		case 'N':
		case 'n': 
			return 11;
		case 'P':
		case 'p': 
			return 12;
		case 'Q':
		case 'q': 
			return 13;
		case 'R':
		case 'r': 
			return 14;
		case 'S':
		case 's': 
			return 15;
		case 'T':
		case 't': 
			return 16;
		case 'V':
		case 'v': 
			return 17;
		case 'W':
		case 'w': 
			return 18;
		case 'Y':
		case 'y': 
			return 19;
	}
	// invalid character
	return -1;
}


/* enumerate the amino acid types 
 * (reduced amino acid alphabet) */
static Py_ssize_t simple_aa2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
		case 'I':
		case 'i': 
		case 'L':
		case 'l': 
		case 'M':
		case 'm': 
		case 'F':
		case 'f': 
		case 'P':
		case 'p': 
		case 'W':
		case 'w': 
		case 'V':
		case 'v': 
			// non-polar
			return 0;
		case 'N':
		case 'n': 
		case 'C':
		case 'c': 
		case 'Q':
		case 'q': 
		case 'G':
		case 'g': 
		case 'S':
		case 's': 
		case 'T':
		case 't': 
		case 'Y':
		case 'y': 
			// polar
			return 1;
		case 'D':
		case 'd': 
		case 'E':
		case 'e':
			// acidic
			return 2;
		case 'R':
		case 'r': 
		case 'H':
		case 'h': 
		case 'K':
		case 'k': 
			// basic
			return 3;
	}
	// invalid character
	return -1;
}

static Py_ssize_t alphabet_size(char c)
{
	switch(c) {
		case 'p':
			// purine/pyrimidine
			return 2;
		case 'a':
			// regular amino acids
			return 20;
		case 'r':
			// amino acids reduced to types
			return 4;
		case 'n': 
			// regualr nucleotides
			// fall through
		default:
			return 4;
		}
}


static Py_ssize_t(*select_char2int(char c)) (char)
{
	switch(c) {
	case 'a':
		// amino acid
		return aa2int;
	case 'r':
		// reduced amino acid
		return simple_aa2int;
	case 'p':
		// purine/pyrimidine
		return simple_nt2int;
	case 'n':
		// nucleotide
		// fall through
	default:
		return nt2int;
	}
}


static Py_ssize_t power(Py_ssize_t base, Py_ssize_t exponent)
{
	Py_ssize_t res = 1, i;
	for (i=0; i<exponent; i++) {
		res *= base;
	}
	return res;
}


static PyObject* count_words(PyObject* self, PyObject* args)
{
    const char* seq;
    Py_ssize_t len;
    Py_ssize_t K;
    char alphabet;
 
    if (!PyArg_ParseTuple(args, "s#nc", &seq, &len, &K, &alphabet)) {
		return NULL;
	}
	
	Py_ssize_t ALPHABET_SIZE = alphabet_size(alphabet);
	Py_ssize_t (*char2int)(char) = select_char2int(alphabet);
	
	Py_ssize_t N = power(ALPHABET_SIZE, K); // distinct words ALPHABET_SIZE**K

	// build the array
	npy_intp dims[] = {N};
	PyArrayObject *rslt =  (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	npy_double * fv = (npy_double*) rslt->data;
	memset(fv, 0, N * sizeof(npy_double));

	// fill the vector
	Py_ssize_t i, j;
	for(i=0; i<=len-K; i+=K){
		npy_intp ix=0;
		for(j=0; j<K; ++j){
			char c = toupper(seq[i+j]);
			Py_ssize_t c_index = char2int(c);
			if (c_index == -1){
				break;
			}
			ix = ix*ALPHABET_SIZE + c_index;
		}
		if (j==K){
			fv[ix] += 1;
		}
	}
	
    return PyArray_Return(rslt);
}


static PyObject* count_overlapping_words(PyObject* self, PyObject* args)
{
    const char* seq;
    Py_ssize_t len;
    Py_ssize_t K;
    char alphabet;
 
    if (!PyArg_ParseTuple(args, "s#nc", &seq, &len, &K, &alphabet)) {
		return NULL;
	}

	Py_ssize_t ALPHABET_SIZE = alphabet_size(alphabet);
	Py_ssize_t (*char2int)(char) = select_char2int(alphabet);
	

	Py_ssize_t N = power(ALPHABET_SIZE, K); // distinct words
	Py_ssize_t MODULUS = N/ALPHABET_SIZE;
	
	// build the array
	npy_intp dims[] = {N};
	PyArrayObject *rslt =  (PyArrayObject *) PyArray_SimpleNew(1, dims, 
	NPY_DOUBLE);
	npy_double * fv = (npy_double*) rslt->data;
	memset(fv, 0, N * sizeof(npy_double));

    // implement counting
	Py_ssize_t pos = 0;
	char in_progress = 0; // set to 1 when extending a valid word
	npy_intp ix = 0;
	Py_ssize_t i;
	while (pos<len){
		if (!in_progress){
			if (pos>len-K){
				break;
			}
			ix = 0;
			in_progress = 1;
			for(i=pos; (i<len) && (i<pos+K); ++i){
				char c = toupper(seq[i]);
				Py_ssize_t c_index = char2int(c);
				if (c_index == -1){
					pos = i+1;
					in_progress = 0;
					break;
				}	
				ix = ALPHABET_SIZE*ix + c_index;
			}
			if (in_progress){
				fv[ix] += 1;
				pos += K;
			}
		} else { // if (!in_progress)
			char c = toupper(seq[pos]);
			Py_ssize_t c_index = char2int(c);
			if (c_index == -1){
				in_progress = 0;
				pos += 1;
				continue;
			}
			ix = ALPHABET_SIZE * (ix % MODULUS) + c_index;
			*(fv+ix) += 1;
			pos += 1;
		}
	}
   
    return PyArray_Return(rslt);
}


static PyObject* count_sliding_overlapping_words(PyObject* self, 
PyObject* args)
{
    const char* seq;
    Py_ssize_t len;
    Py_ssize_t K;
    Py_ssize_t win_size;
    Py_ssize_t step;
    char alphabet;
 
    if (!PyArg_ParseTuple(args, "s#nnnc", &seq, &len, &K, &win_size, &step, 
    &alphabet)) {
		return NULL;
	}

	Py_ssize_t ALPHABET_SIZE = alphabet_size(alphabet);
	Py_ssize_t (*char2int)(char) = select_char2int(alphabet);

	Py_ssize_t N = power(ALPHABET_SIZE, K); // distinct words
	Py_ssize_t MODULUS = N/ALPHABET_SIZE;
	
	// number of windows to process
	// this expression can handle short len giving 0 or a negative M 
	// when necessary
	Py_ssize_t M = (len - win_size - K + 1 + step) / step; 
	if (M <= 0) {
		Py_RETURN_NONE;
	}

	// build the array
	npy_intp dims[] = {M, N};
	PyArrayObject *rslt =  (PyArrayObject *) PyArray_SimpleNew(2, dims, 
	NPY_DOUBLE);
	npy_double * fv = (npy_double*) rslt->data;
	memset(fv, 0, M * N * sizeof(npy_double));
	if (!PyArray_ISCARRAY(rslt)) {
		// check that it is C-contiguous, writeable and has a native byte order
		return NULL;
	}

	// word at each position
	Py_ssize_t *words = (Py_ssize_t*) malloc((len - K + 1) * sizeof(Py_ssize_t)); 
	memset(words, -1, (len - K + 1) * sizeof(Py_ssize_t)); // initialize with -1
	// the value of -1 indicates the word is not a valid k-mer

    // identify words and store them in a vector
	Py_ssize_t pos = 0;
	char in_progress = 0; // set to 1 when extending a valid word
	Py_ssize_t ix = 0;
	Py_ssize_t i;
	while (pos<len){
		if (!in_progress){
			if (pos>len-K){
				break;
			}
			ix = 0;
			in_progress = 1;
			for(i=pos; i<pos+K; ++i){
				char c = toupper(seq[i]);
				Py_ssize_t c_index = char2int(c);
				if (c_index == -1){
					pos = i+1;
					in_progress = 0;
					break;
				}
				ix = ALPHABET_SIZE*ix + nt2int(c);
			}
			if (in_progress) {
				words[pos] = ix;
				pos += K;
			}
		} else { // if (!in_progress)
			char c = toupper(seq[pos]);
			Py_ssize_t c_index = char2int(c);
			if (c_index == -1){
				in_progress = 0;
				pos++;
				continue;
			}
			ix = ALPHABET_SIZE * (ix % MODULUS) + nt2int(c);
			pos++;
			words[pos-K] = ix;
		}
	}

	//compute array elements
	Py_ssize_t *starting_position = words; // starting poPy_ssize_t in the words vector
	Py_ssize_t m;
	for (m=0; m<M; m++) {
		// pointer to the m-th row
		npy_double* array_row = (npy_double*) PyArray_GETPTR2(rslt, m, 0);

		for (i=0; i<win_size; i++) {
			Py_ssize_t cur_word = starting_position[i];
			if (cur_word >= 0) {
				array_row[cur_word]++;
			}
		}
		starting_position += step;
	}
	free(words);

    return PyArray_Return(rslt);
}


static PyObject* alphabet_size_ext(PyObject* self, PyObject* args)
{
    char alphabet;
 
    if (!PyArg_ParseTuple(args, "c", &alphabet)) {
		return NULL;
	}
	
	Py_ssize_t ALPHABET_SIZE = alphabet_size(alphabet);

    return Py_BuildValue("n", ALPHABET_SIZE);
}


static PyObject* res2int_ext(PyObject* self, PyObject* args)
{
    char ch;
    char alphabet;
    
    if (!PyArg_ParseTuple(args, "cc", &ch, &alphabet)) {
		return NULL;
	}
    
    Py_ssize_t (*char2int)(char) = select_char2int(alphabet);
    Py_ssize_t numeric_repr = char2int(ch);
    
    return Py_BuildValue("n", numeric_repr);
}


static PyMethodDef WordCountMethods[] =
{
	{"count_words", count_words, METH_VARARGS, "Count n-mers in a string.\n"
		"Arguments: sequence, n-mer size, alphabet type"},
	{"count_overlapping_words", count_overlapping_words, METH_VARARGS,
		"Count n-mers in a string.\n" 
		"Arguments: sequence, n-mer size, alphabet type"},
	{"count_sliding_overlapping_words", count_sliding_overlapping_words,
		METH_VARARGS, "Count n-mers in a liding window over a string.\n" 
		"Arguments: sequence, n-mer size, window length, step, alphabet type"},
	{"alphabet_size", alphabet_size_ext,
		METH_VARARGS, "Size of the alphabet.\n" 
		"Arguments: alphabet type"},
    {"res2int", res2int_ext,
		METH_VARARGS, "Numeric index of the character.\n" 
        "Arguments: character, alphabet code"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef wordcountmodule_definition = {
    PyModuleDef_HEAD_INIT,
    "_wordcount",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    WordCountMethods
};

PyMODINIT_FUNC PyInit__wordcount(void) {
    Py_Initialize();
    import_array();
    return PyModule_Create(&wordcountmodule_definition);
}

