# BARCODE ERROR CORRECTION PROGRAMS
*Artur Rego-Costa, January 2023*

This repo contains two C++ programs to perform string error correction as usually
necessary when dealing with barcode-amplicon sequencing. One program takes a barcode 
reference list, while the other does not. Programs were written by Alex Nguyen Ba 
and José Rojas Echenique and should be cited as indicated below. 

Both programs can parallelize if run in a multiple-CPU node. In SLURM clusters,
you do that by passing the -c flag: e.g. "#SBATCH -c 6" for 6 CPUs.

Both programs operate based on the principle that two strings A and B that differ by a
Levenshtein distance d will get merged if d <= err and the ratio between the counts
of either string is nA/nB > err_ratio (in which case B strings get converted to A) or
nB/nA > err_ratio (the other way around). The programs do this efficiently by using a
trie-based algorithm. More can be read on (1).

Both programs take as an input a **headerless** tab-separated table, where one of the columns
needs to be error-corrected.

## Dependencies
error_correct_from_dict uses the robin_map class, which is a header-only library maintained
in [github.com/Tessil/robin-map](https://github.com/Tessil/robin-map). For ease of use, I include it in /include/tsl. It would
be common to copy the /include/tsl directory into a local ~/local/include directory.
If you do otherwise, edit the compilation line appropriately.

For Harvad Odyssey cluster users, gcc can be loaded using command:
	module load gcc/8.2.0-fasrc01

## error_correct_from_dict.cpp
Error correct based on a dictionary of known barcodes. This program does not have an
err_ratio criterion: it'll simply correct any string that is similar enough (d <= err) 
to a string in the dictionary, and discard any string for which d > err.

### Compilation
	g++ error_correct_from_dict.cpp -o error_correct_from_dict -O3 -std=c++0x -I$HOME/local/include -lpthread -fopenmp -D_GLIBCXX_PARALLEL

### Usage
	error_correct_from_dict [-cpu i] [-err i] -bc m -dict dict_file input_file > output_file
	cat input_file | error_correct_from_dict [-cpu i] [-err i] -bc -dict dict_file > output_file

### Arguments
	-cpu i				integer number of CPUs per node available for parallelization [default: 1]
	-err i 				integer Levenshtein distance criterion for error correction.
						i.e. two strings get corrected if their Levenshtein distance d <= err.
						[default: 2]
	-bc m				integer column which to error correct (1-indexed, meaning that m == 1 
						for the first column)
	-dict dict_file		dict_file is a single column file with a list of reference barcodes.
	input file			headerless tab-separated input table

### Output
A headerless tab-separated table of the same format as the input, but where the m-th column
is now error-corrected. 
Input lines get removed if the m-th column contains a barcode that is not a perfect match to
neither barcode in the reference list, neither it could not be corrected into any barcode in
the reference list.

### Citation
This code is exactly as used in (2), so cite that.

## fb_trie_rethread_5.cpp
Error correct without a reference list. This program will discard any string which is 
seen less than a number min_count of times, and that cannot be corrected into another
higher-count string.

### Compilation
	g++ fb_trie_rethread_5.cpp -o fb_trie_rethread_5 -O3 -std=c++0x -lpthread -fopenmp -D_GLIBCXX_PARALLEL

### Usage
	fb_trie_rethread_5 [-cpu i] [-err i] [-err_ratio i] [-min_count i] -pop n -bc m file > output_file
	cat file | fb_trie_rethread_5 [-cpu i] [-err i] [-err_ratio i] [-min_count i] -pop n -bc m > output_file

#### Arguments
	-cpu i				number of CPUs per node available for parallelization
	-err i 				Levenshtein distance criterion for error correction.
						i.e. two strings get corrected if their Levenshtein distance d <= err.
						[default: 2]
	-err_ratio i		float ratio nA/nB of number of times strings A and B are seen such that
						(i) if nA/nB > err_ratio, then B gets corrected to A, or (ii) nB/nA, then
						A gets corrected into B. [default: 32]
	-min_count i		integer minimum number of at and below which a string will be dicarded,
						unless it can be corrected into another higher-count string seen in the 
						data. [default: 10]
	-bc m				column which to error correct (1-indexed, meaning that m == 1 for the
						first column)
	input file			headerless tab-separated input table

### Output
A headerless tab-separated table of the same format as the input, but where the m-th column
is now error-corrected. 
Input lines get removed if the m-th column contains a barcode that is seen less than min_count
times and is never error-corrected into another higher-count barcode.

### Citation
This code is very much similar to that of (1), with the exception of parameters err, err_ratio
and min_count being now passed as command line arguments. Cite (1).

## REFERENCES
>(1) Nguyen Ba, Alex N., Ivana Cvijović, José I. Rojas Echenique, Katherine R. Lawrence, Artur
Rego-Costa, Xianan Liu, Sasha F. Levy, and Michael M. Desai. “High-Resolution Lineage Tracking
Reveals Travelling Wave of Adaptation in Laboratory Yeast.” Nature 575, no. 7783 (2019): 494–99.
https://doi.org/10.1038/s41586-019-1749-3.
>
>(2) Nguyen Ba, Alex N., Katherine R. Lawrence, Artur Rego-Costa, Shreyas Gopalakrishnan, Daniel
Temko, Franziska Michor, and Michael M. Desai. “Barcoded Bulk QTL Mapping Reveals Highly Polygenic
and Epistatic Architecture of Complex Traits in Yeast.” ELife 11 (February 1, 2022).
https://doi.org/10.7554/ELIFE.73983.
