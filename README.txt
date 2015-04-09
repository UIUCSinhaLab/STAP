STAP v2: Sequence To Binding Affinity Prediction (version 2)
Saurabh Sinha’s Lab   <sinhas@illinois.edu>
Initially created by Xin He
Last modified on Oct. 17, 2013

Description

The program makes use of a biophysically motivated computational model to
analyze transcription factor (TF)-DNA binding data, such as ChIP-chip or
ChIP-SEQ data. The program assumes that the measured affinity of a sequence to
a TF (TF_exp) in some ChIP-chip or ChIP-SEQ experiment is determined by: 1)
the number and strength of binding sites of TF_exp in this sequence; 2) the
presence of other sites that may influence the binding affinities of the
neighboring sites to TF_exp. Specifically, it takes as input a set of DNA
sequences, their binding affinities to some TF as measured by experiments
(TF_exp), and the position weight matrices (PWMs) of a set of TFs, including
TF_exp. It will learn the relevant parameters of the biophysical model of
TF-DNA interaction, including those of TF-DNA interaction and those of TF-TF
cooperative interactions. The program can be used for several purposes:

(1) Test if a given TF binding motif can predict the binding affinities of the
sequences . It
predicts the binding sites based on this motif, and computes the theoretical
values of the binding affinities of the sequences. The predicted values will
be compared with the observations to judge the success of the model.

(2) When multiple motifs are given as inputs, the program assumes the first
motif is the one of the experimetnal TF (TF_exp), and the rest are the motifs
that may directly or indirectly interact with TF_exp (meaning that the
neigboring sites of other factors can influence DNA binding of TF_exp). The
program will learn which motifs are likely to interact with TF_exp and further
learn which secondary-motif’s influence is likely to be (a) through long-range
versus short-range interactions with the primary motif, (b) through
synergistic or antagonistic interactions, and (c) through modulation of local
DNA accessibility or direct interactions between TFs.

(3) Once a biophysical model is learned, it can be applied to predict
affinities of other sequences. This would be useful, for example, for
analyzing sequences in a different organism.

Installation

The program needs GNU Scientific Library (GSL). If it is not installed in your
system, go to:
http://www.gnu.org/software/gsl/

Note that after installing GSL, you need to change the start-up script of your
shell, e.g., .bash_profile or .profile at your home directory, which depends
on your OS. Suppose the GSL installation directory is /raid/apps/gsl-1.8/lib
(i.e., my_gsl_dir=/raid/apps/gsl-1.8/lib):

LD_LIBRARY_PATH=/raid/apps/gsl-1.8/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

After extracting the program, change the GSL directory in src/Makefile, e.g.:
GSL_DIR = my_gsl_dir

Then simply type:
make

Running the program

Usage in general:
./seq2binding -s <seqFile> -d <dataFile> -m <motifFile> The program takes
three arguments as input:
seqFile: the FASTA format file of sequences. See examples/Nanog_top_500.fa.

>chr1:136351629-136351631 136351630 -250 +250
>gtggtgatgcccaaccacagaattattttgttgctactttataactgtaattttgatcct
>chr3:122137593-122137593 122137593 -250 +250
atttctagttccagtgactgggagactgaaacaagagagtcacttgagtacaggagtgca

dataFile: the binding data of all sequences in the seqFile. The first column
is the sequence id (must be the same as those used in seqFile, and in the same
order), and the second column is the measured strength of binding. See
examples/Nanog_top_500.txt.


chr1:136351629-136351631        312 chr3:122137593-122137593        307

motifFile: the motifs of the TFs. It could contain multiple motifs. The header
line consists of motif name, length and pseudocount (0.5 should be OK for most
motifs). The first motif should be the one of TF_exp, and the rest are
putative TFs that interact cooperatively with TF_exp. See
examples/Nanog_Oct4.wtmx and examples/Nanog.wtmx.

>Nanog  9	0.5	
20      225	46	209
70      0	19	411
50      66	381	3
434     45	0	21
55      5	66	374
17      32	222	229
74      18	325	83
8       243	146	103
48      145	6	301
<		

Output:
(1) Estimated parameters: binding parameter (how strongly the TF binds with
its binding site); the interaction weight parameters between any pair of TFs
(the order of motifs in the matrix follows the order defined in motifFile):
greater than 1 if cooperative interaction, less than 1 antagonistic, 1 if no
interaction.
(2) Pearson correlation between predicted binding and observed binding.

Examples:
(1) Run with a single factor: test if the provided Nanog motif explains the
binding of top 500 Nanog sequences in experiments: (under examples/ directory)
../src/seq2binding -m Nanog.wtmx -s Nanog_top_500.fa -d Nanog_top_500.txt

(2) Run with multiple factors (TF_exp, and other factors): test if Nanog
interacts cooperatively with Oct4 in the top 500 Nanog sequences: (under
examples/ directory)
../src/seq2binding -m Nanog_Oct4.wtmx -s Nanog_top_500.fa -d Nanog_top_500.txt

Advanced options

-ts <testSeqFile> -td <testDataFile>: test the trained model in additional
testing data. The format of testSeqFile and testDatafile is the same as
seqFile and dataFile.

-n <nExps>: the number of experiments being analyzed. The default value is 1,
i.e. only one experiment (binding data of one TF) is analyzed. When analyzing
binding data of multiple TFs, set nExps as the number of TFs. In this case, it
is assumed that seqFile and dataFile contain the concantenation of data of
multiple factors (assume the number of records of each TF is equal, thus no
explicit delimiter is needed between data of different TFs).


-co <coopOption>: the option of cooperative binding. 0 - no cooperativity at
all; 1 - no self- cooperativity, but hetero-cooperativity; 2 - allow all
cooperativities

-io <interactionOption>: the option for modeling factor-factor interaction. 0
- binary; 1 - linear; 2 - periodic; 3 - normal distribution.

-dt <d_max>: the maximum site-site distance if possible interaction happens
(beyond which there will be no interaction)

-cc <interactionWeightRange> : possible DNA-specific transcription
factor-factor intercation. 1 - possible cooperative factor-factor interation;
-1 - possible antagonistic factor-factor interaction; 11 - both. Different
options cause the interaction weight parameters cast into different ranges.

-oe <r_orientationEffect>: orientation effect if two interacting sites happen
on different strands, suitable for all interaction modeling options. If two
sides are in the same strands, there is no orientation effect, i.e.,
r_orientationEffect =1>

-r <distanceRange>: maximum valid distance range in linear interaction mode.
Outside the range, the site-site interaction starts to decline. By default,
distance range is equal to 50.

-se <spacingEffect>: spacing effect in the periodic interaction. The spacing
effect is a constant term due to DNA looping (exponential of the average free
energy of DNA looping). The other parameters associated with the periodic
interaction mode are as follows:
–	T <period>. 10 is by default.
–	phi <r_phi>: phase angle in periodic model if two sites are in the
same strand
–	phi0 <r_phi0>>: phase angle in periodic model if two sites are in the
different strands

-iv <PPInteractionDecayVariance> : Decay variance of DNA-specific TF-TF
binding affinities in the normal distribution model. A protein-protein
physical interaction represents the biological macromolecule as an elastic
mass-and-spring network. Its force constant exhibites an exponential decay
over the distance between a pair of interacting atoms. The force constant is
defined as an exponential decay function of distance. By default, it is 25. 

-et <energyThr> : general energy threshold for each motif site (by default, it
is 10). Or
-llrt <strMotifSiteThr>: motif site thresholds separated by a punctuation mark
(,). Each motif site threshold is set by the motif’s llr score when pval is
closest to 0.05.  Refer to utils/getMotifsLLRPValue.sh for calculating the
motif site LLR threshold.

-rr <n_RandStarts>: the number of simulation runs, each of which starts with
random parameter inputs.

-pm <paramModelOutFile>: print the learned parameter values in the file
paramModelOutFile. The file format is as follows:
maxBindingWts = 1.02727 1.45616
numOfInFactors=1
inFactorIntMat=
1
numOfOutFactors=1
outFactorIntMat=
1.147017
expRatios = 1.000000

-fstr <featureTrainMsgOutFile>: print the TF-DNA binding features on the
training data set in the file featureTrainMsgOutFile

-fste <featureTestMsgOutFile>: print the TF-DNA binding features on the test
data set in the file featureTestMsgOutFile

-p <trainPredictionOutFile>: print the predicted binding intensities (of the
training sequences in seqFile) in the file trainPredictionOutFile.

-tp <testPredictionOutFile>: print the predicted binding intensities (of the
test sequences in testSeqFile) in the file testPredictionOutFile

-pin <paramModelInFile> [-ot]: give the parameter-value input file which
suggested initial values for model parameters.  –ot is optional.
For example, ../src/seq2binding -m Nanog_Oct4.wtmx -s Nanog_top_500.fa -d
Nanog_top_500.txt -cc 1 -io 0 -co 1 -et 8 -pin paramModelInFile.txt

-ot    ignore the learning or training procedure and predict the binding
affinity just based on  given input parameter values.

-isM1BWFixed  -bw <initBindingWeight>:   force the STAP to fix the primary TF
binding weight as a constant of initBindingWeight, not as a free parameter, in
each simulation iteration. 


There are other parameters that control the factor-factor interaction model.
The file utils/run_pair.sh contains examples of using these parameters. In
most cases, you probably do not need to set these parameters.
 
Utilities

In utils/ directory, some useful scripts are included. However, not all of
them can be ready to execute in your system. Some of them are included only
for the purpose of demonstrating the use of program.

run_pair.sh: demonstrate the use of program (for analyzing cooperative
interactions using binding data of two TFs).

create_null_distr.sh: suppose we want to find cooperative factors of one TF
(TF_exp). First run the program using TF_exp and the test motif, and obtain
the correlation coefficient (CC) of the test motif. Then run this script to
get the null distribution of the CC. The script will sample random motifs from
a specified collection of motif, and calculate CC of the random motifs.
shuffle_wtmx.pl: random shuffling of a motif, used by create_null_distr.sh.

split_wtmx.pl: split a file of many PWMs into multiple files, each of which contains a single motif.

getMotifsLLRPValue.sh: calculate the motif site LLR threshold based on its
Pvalue.

