# RNA motif forces
A set of tools to calculate the "force" acting on the maximum length of complementary segments in a given transcript, dinucleotide motif, and sequence complexity

## Installation

### Docker: easiest path
`docker pull TODO`

Alternatively, build the image from `Dockerfile`:
```
docker build . -t forces
```
The build is based on `debian:bullseye` and it requires the `forces.tgz` source tarball.

### Compilation: 
#### In directory src/xds_calculation:
```bash
g++ -O3 -o fasta_xds calculate_xds.cpp 
g++ -O3 -o window_scan scan_window.cpp
```
Install the built executables into `PATH`
### In directory src/python
```
python3 setup.py build && python3 setup.py install
```

## Usage:

* Complexity calculation:
```bash
python get_complexity.py "YOUR SEQUENCE HERE"
```

* Double stranded force calculation in scanning window (of given window_length, sliding from start position by shift_size until end position is reached on a given  chromosome. The human genome assembly has to be provided as a fasta file. Output is provided per line for each window): 
```bash
window_scan fasta_file_with_genome chromosome start stop shiftsize window_length outputfilename
```

* Double stranded force calculation per each individual sequence provided in a fastafile:
```bash
fasta_xds fastafilename
```

* `CpG` force calculation in scanning window:
```
compute_sliding_window_force.py [-h] [-L WINDOW] [-c CONTIG] [-d DIMER] [-e END] [-s START] fasta_infile
```
Arguments `CONTIG` to use, `START`, `END` and `WINDOW` length are optional. `DIMER` defaults to `CG`.