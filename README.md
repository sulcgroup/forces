# RNA motif forces
A set of tools to calculate the "force" acting on the maximum length of complementary segments in a given transcript, dinucleotide motif, and sequence complexity

Compilation: In directory src/xds_calculation:

g++ -O3 -o fasta_xds calculate_xds.cpp g++ -O3 -o window_scan scan_window.cpp

Usage:

Complexity calculation: python get_complexity.py "YOUR SEQUENCE HERE"

Double stranded force calculation in scanning window or in a file with fasta sequences respectively

./window_scan fasta_file_with_genome chromosome start stop shiftsize sequence_length outputfilename ./fasta_xds fastafilename

A set of tools to calculate the "force" acting on the maximum length of complementary segments in a given transcript, dinucleotide motif, and sequence complexity
