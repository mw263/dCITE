The command line Python application CIC_msa.py computes the Deflated Cladistic
Information based on Total Entropy (dCITE) statistic which is described in
"dCITE: Measuring Necessary Cladistic Information can help
you Reduce Polytomy Artefacts in Trees"

From that paper's abstract:

"Biologists regularly create phylogenetic trees to better understand the evolutionary
origins of their species of interest, and often use genomes as their data source.
However, as more and more incomplete genomes are published, in many cases 
it may  not be possible to compute genome-based phylogenetic trees due to 
large gaps in the assembled sequences.  In addition, comparison of complete genomes
may not even be desirable due to the presence of horizontally acquired and
homologous genes.  A decision must therefore be made about which gene, or gene
combinations, should be used to compute a tree.

Deflated Cladistic Information based on Total Entropy (dCITE) is proposed
as an easily computed metric for measuring the cladistic information in multiple
sequence alignments  representing a range of taxa, without the need to first compute
the corresponding trees.  dCITE scores can be used to rank candidate genes or
decide whether input sequences provide  insufficient cladistic information, making
artefactual polytomies more likely. The dCITE method can be applied to protein,
nucleotide or encoded phenotypic data, so can be used to select which data-type
is most appropriate, given the choice.  In a series of  experiments the dCITE
method was compared with related measures.  Then, as a practical demonstration,
the ideas developed in the paper were applied to a dataset representing species
from the order Campylobacterales; trees based on sequence combinations, selected
on the basis of their dCITE scores, were compared with a tree constructed to
mimic Multi-Locus Sequence Typing (MLST) combinations of fragments.

We see that the greater the dCITE score the more likely it is that the computed
phylogenetic tree will be free of artefactual polytomies.  Secondly, cladistic
information saturates, beyond which little additional cladistic information can
be obtained by adding additional sequences. Finally, sequences with high cladistic
information produce more consistent trees for the same taxa."

CIC_msa.py 

Running the program

The application is a Unix/Linux command line programs written in Python. To
call it directly you first need to make it callable, e.g.

chmod u+x CIC_msa.py

and then call it in the usual way for a Unix program, e.g.

CIC_msa.py test_MSA.fasta

You can also run a script by typing: python ./script, e.g.
python ./CIC_msa.py test_MSA.fasta

Invoking a script without any arguments will return a summary of the expected
input format and command-line options:

CIC_msa.py

returns

Usage: CIC_msa.py [OPTIONS] <Multi FASTA format Multiple Sequence Alignment file>
The options, in no particular order, are:
        -type [DNA, protein, standard, mixed] - The input data DNA, protein, standard ([0-9]) or mixed (default protein)
        -usegap - Treat gaps as a 21st/5th character state (otherwise ignored)
        -OSasID - Use the OS field, if present in the header, as the sequence ID
        -diversity_cap - Ignore sites where more than 75%% of the alphabet appear (default False)
        -noparsimony - Ignore parsimony assumption on counts (min count can now be 1)(default False)
        -add_deltamin - Add the DeltaMin computation (Wortley & Scotlan, 2006)(default False)
        -o [F|-] - Output file name F or stdout ('-') (default: input basename + .msa_cic)

Dependencies
  - There are no dependencies beyond those available in standard Python 2.7
