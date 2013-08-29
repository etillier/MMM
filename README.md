MMMvII
===

MatrixMatchMaker version 2. Elisabeth Tillier, Robert L. Charleblois, Alexander Rodionov, Alexander Bezginov, Jonathan Rose.
University of Toronto and University Health Network.

This is free software with no guarantees.

please cite:

Rodionov A., Bezginov A., Rose, J.  and  Tillier, ERM( 2011).  A New, Fast Algorithm for Detecting Protein Coevolution using Maximum Compatible Cliques. Algorithms for Molecular Biology 6:17  

Tillier ERM and R.L. Charlebois (2009) The Human Protein Coevolution Network. Genome Research. 19: 1861-1871.




Usage 1: ./mmmWithL -d1 distfile1 -d2 distfile2 -o outfile -a allowance [-brief] [-table [-top]] [-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] [-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb] [-v] [-restore stateFileName]
...or...
Usage 2: ./mmmWithL -b batchfile -o outfile -a allowance [-brief] [-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] [-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb] [-v]

Required parameters:

-d1 distfile1 (if -b batchfile is not specified)
	Specifies the file name of the first distance matrix.
-d2 distfile2 (if -b batchfile is not specified)
	Specifies the file name of the second distance matrix.
-b batchfile (if -d1 distfile1 -d2 distfile2 are not specified)
	Specifies the name of a file containing tab-delimited pairs of distance matrix file names, one pair per line.
-o outfile
	The desired name of the output file.
-a allowance
	Specifies a value between 0.0 and 1.0, with smaller values requiring more precise matching, and larger values tolerating looser matching.

Optional parameters:

-brief
	Outputs a single tab-delimited line of output for the pair of distance matrices.
-table (valid only when -b batchfile is not specified)
	Produces a tab-delimited matrix that lists the fit-weighted proportion of matches for each entry from the two distance matrices, over all common submatrices. For example, if seq1 matched seq2 in two of ten discovered submatrices, with a score of 0.8 and 0.9, respectively, then the entry for seq1_seq2 would be (0.8 + 0.9) / 10 = 0.17
-top (valid only when -table is specified)
	Restricts the -table output to just the best-scoring submatrix.
-u
	Only attempts to match entries from the pair of distance matrices that have the same tag in their name. A tag ends with a pipe ("|") character, e.g. Ecoli in Ecoli|DnaA, or Hsapiens in Hsapiens|myc
-req reqTaxName
	Requires that a submatrix include at least one entry with the specified reqTaxName, e.g. Hsapiens.
-sz strictSize
	Returns common submatrices of size strictSize, if any. Does not attempt to find larger submatrices.
-p distfileNamesPrefix
	The prefix is appended to the path of all distance matrices.
-z
	Forces the output of zeros for combinations that have not reached the minSize (to be used in 'brief' mode).
-minsize minSize
	Returns common submatrices of at least minSize, if any. Will not return smaller submatricies. Use this option to save memory with large matrices.
-maxtrees maxTrees
	Limits the number of returned trees. Saves memory.
-rangeStart start -rangeEnd end
	For distributed computing, launch this program multiple times with non-overlapping ranges, [start, end), start >= 0 and end <= distfile1's matrix size (auto-adjusted when -b is specified). Incompatible with -rndOrder.
-quickmp
	If specified, the quickly but naively estimated match potential is included in batch output.
-v
	Outputs the size of the largest common submatrix found so far, to stdout.
-onlymp
	Does not perform MMM analysis, only prints MatchPotential (brief and MMML modes only).
-restore stateFileName
	Resumes work from a previously interrupted MMM session
-exacttime
	Do not round running time to nearest second, in brief mode
-upper cutoffUpper
	All distances above cutoffUpper are set to zero.
-lower cutoffLower
	All distances below cutoffLower are set to zero.
-log
	Replaces the distances with their logs: dst=ln(dst)+1 if dst > 1.
-noheader
	Omits the header from the output file (brief and MMML modes only).
-onlyheader
	Prints the header to the output file and exits (brief and MMML modes only).
-ignoreerrors
	Does not abort the run even in the case of program errors (intended for debugging only).
-aln2pmb
	Enable recognition of .aln files as FASTA alignments and processing them into PMB matrices.
	Simply specify the alignment filename in place of distance matrix file. This works in both command-line and batch mode.
	A named alignment slice can be defined as comma-separated ranges of columns added after the filename.
	When using multiple ranges to define a slice, the ranges must go in increasing order, and must not overlap.
	An "inverse"slice (remove the slice itself and keep all other columns) can be specified by adding 'r' after the range. Do not use with multiple ranges.
	Format: <filename>|[range name]:<start>-<end>[r][,<start>-<end> ...]
	Examples:
		alignment1.aln
		alignment2.aln|pfam12820:373-546
		alignment3.aln|disorder:16-90,150-150,430-519
		alignment4.aln|:300-400r
MMML extensions:
	-Lt rootedTreeFileName (required)
		File name for a Newick style tree file encompassing all of the taxa to be encountered in the MMM run.
	-LtopCov (optional)
		top-coverage style output.
	-Li ignoreTaxaFileName (optional)
		file listing taxon names to ignore in the MMML analysis.
	-Lb minBootstrapSupport (optional, default=1)
		an integer between 0 and 100 representing the threshold below which to collapse branches in parsing the treefile.
	-List (optional)
		if specified, lists the taxa in each clade in the output.
	-Lv (optional)
		if specified, outputs the clade sets.
