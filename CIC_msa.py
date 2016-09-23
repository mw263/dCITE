#!/usr/bin/env python

"""
Compute the total Cladistic Information Content. This is a Python port of
James Cotton's MolNXcharinfo.pl (which implements Cotton & Wilkinson (Taxon, 2008),
based on a formula in Steel and Penny (2005).
DNA and protein alphabets are now properly handled (original favours DNA and miscalls N
as gap).  gap chars not hardwired into the code.

As it makes no difference to the calculation, gaps are treated the same
as missing character - simply fewer real letters. 
Ambiguity codes, and a fortiori X, are here treated as gaps

This CIC is a mapping from a multiple sequence alignment to a score reflecting
the cladistic information inherent in the MSA. It is different to the CIC defined by
Thorley, Wilkinson, Charleston (1998), which works on trees rather than the inputs used
to build the trees

CIC as implemented in MolNXcharinfo.pl tots up a score for each site and then adds
the scores across all the sites. Two issues, duplicate sequences will double the score,
and sites can be linked, in that a subsequence in one sequence is duplicated in other
sequences. Given N characters in one sequence that are also found, in the same sites
in a second sequence, only one character is independent, with the others being dependent,
so the bit score should be discounted.

Total Entropy is like Delta min a score in that identity columns are ignored (score 0) but Delta Min
scores 1 for each character in a column (minus 1), irrespective of how many taxa have that character
at that site. Total Entropy takes account of the proportions.

Optionally, treat gap as another state (rather than X, which conveys no information).

A perpetual issue is that sites are treated as being independent, when in fact they
are linked across sequences forming  quasi haplotypes. The question is how to reduce the impact
of such linkage, which implies that what you learn from one site you also learn from another.
A partial solution is to realise that, ignoring identity sites, the various amino acids
in a site induce a split in the taxa at that site. If another site induces an identical split it can
be ignored as adding new new information. Need a set of sets data structure, but Python can't do that
so use list of lists, or sligthly more efficient cardinality -> tuple list list
"""

import sys, math, re

DNA_ALPHABET = set(['A', 'C', 'G', 'T'])
AA_ALPHABET = set([ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
			 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
STANDARD_ALPHABET = set(['0','1','2','3','4','5','6','7','8','9'])
MIXED_ALPHABET = AA_ALPHABET | STANDARD_ALPHABET

FASTA_PATN = re.compile("\.fasta", re.I)
OS_PATN = re.compile(" *OS= *")
SECOND_FIELD_PATN = re.compile(" *[A-Z][A-Z]= *")
DATA_TYPE_LIST = ["DNA", "protein", "standard", "mixed"]

# The functions can be called by another module ==> seqfile != None
def init(seqfile=None, type="protein", usegap=False) :
  params_dict = {"type":type, "usegap":usegap, "o":None, "OSasID":False, "diversity_cap":False,
		 "noparsimony":False, "add_deltamin":False}
  if seqfile == None :  # call from command line
    args = sys.argv
    if len(args) == 1:
      sys.stderr.write("Usage: %s [OPTIONS] <Multi FASTA format Multiple Sequence Alignment file>\n"\
				%   args[0])
      sys.stderr.write("The options, in no particular order, are:\n")
      sys.stderr.write("\t-type [DNA, protein, standard, mixed] - The input data DNA, protein, standard ([0-9]) or mixed (default protein)\n")
      sys.stderr.write("\t-usegap - Treat gaps as a 21st/5th character state (otherwise ignored)\n")
      sys.stderr.write("\t-OSasID - Use the OS field, if present in the header, as the sequence ID\n")
      sys.stderr.write("\t-diversity_cap - Ignore sites where more than 75%% of the alphabet appear (default False)\n")
      sys.stderr.write("\t-noparsimony - Ignore parsimony assumption on counts (min count can now be 1)(default False)\n")
      sys.stderr.write("\t-add_deltamin - Add the DeltaMin computation (Wortley & Scotlan, 2006)(default False)\n")
      sys.stderr.write("\t-o [F|-] - Output file name F or stdout ('-') (default: input basename + .msa_cic)\n")
      sys.exit(0)

    i = 1
    while i < len(args) :
      if args[i][0] == '-' :
	if args[i] in  ["-usegap", "-OSasID", "-diversity_cap", "-noparsimony", "-add_deltamin"] :
	  params_dict[args[i][1:]] = True
	  i += 1
	elif args[i] == "-o" :
	  params_dict["o"] = args[i+1]   # file name will be replaced with actual file
	  i += 2
	elif args[i] == "-type" :
	  if not args[i+1] in DATA_TYPE_LIST :
	    sys.stderr.write("Mishap. Possible input datatypes are: DNA, protein, standard or mixed\n")
	    sys.exit(1)
	  params_dict[args[i][1:]] = args[i+1]
	  i += 2
	else:
	  sys.stderr.write("Unknown option: %s\n" % args[i])
	  sys.exit(1)
      else :
	break

    if i == len(args) :
      sys.stderr.write("No Multi-FASTA file specified\n")
      sys.exit(1)
    infilename = args[i]
  else:
    infilename = seqfile

  try:
    infile = open(infilename, 'r')
  except IOError:
    sys.stderr.write("%s: Cannot open input file %s\n" % (args[0], args[i]))
    sys.exit(1)

  if params_dict["o"] == '-' :
    params_dict["o"] = sys.stdout
  elif seqfile == None : # A file that needs to be opened for the output because call not from another module
    if params_dict["o"] != None :
      outfile_name = params_dict["o"]
    else:
      outfile_name = FASTA_PATN.split(args[i])[0] + ".msa_cic"
    try:
      params_dict["o"] = open(outfile_name, 'w')
    except IOError:
      sys.stderr.write("Cannot open output file %s\n" % outfile_name)
      sys.exit(1)

  # Rather than using a dictionary, to single ID mapping, which has the effect of deleting duplicates,
  # keep a list of duplicates so they can be reported (the first is the one actually returned)
  seq_dict = {}  # seq_dict: seq -> ID list
  there_are_duplicates = False
  len_ref_seq = None
  seq = ""
  nseqs0 = 0
  ID = None
  for line in infile:
    line = line.strip()
    if line == "" :
      continue

    if line[0] == '>' :
      if seq != "" :
	if len_ref_seq != None :  # same length sanity check
	  if len(seq) != len_ref_seq :
	    sys.stderr.write("The length of two sequences from the input multiple sequence allignment differ, %d vs %d!!\n"\
			% (len(seq), len_ref_seq))
	    sys.exit(1)
	else:
	  len_ref_seq = len(seq)
	if seq_dict.has_key(seq) :
	  seq_dict[seq].append(ID)
	  there_are_duplicates = True
	else:
	  seq_dict[seq] = [ID]
	nseqs0 += 1
      seq = ""
      if params_dict["OSasID"] :
	fields = OS_PATN.split(line)
	ID = SECOND_FIELD_PATN.split(fields[1])[0]
	if ID == "" :
	  ID = line.split()[0][1:]
      else:  # The ID is just the original sequence ID
	ID = line.split()[0][1:]
    else:
      seq = seq + line

  if seq != "" :
    if len_ref_seq != None : 
      if len(seq) != len_ref_seq :
        sys.stderr.write("The length of two sequences from the input MSA differ %d vs %d\n"\
			% (len(seq), len_ref_seq))
        sys.exit(1)
    if seq_dict.has_key(seq) :
      seq_dict[seq].append(ID)
      there_are_duplicates = True
    else:
      seq_dict[seq] = [ID]
    nseqs0 += 1

  if params_dict["type"] == "DNA" :
    alphabet = DNA_ALPHABET
  elif params_dict["type"] == "protein" :
    alphabet = AA_ALPHABET
  elif params_dict["type"] == "standard" :
    alphabet = STANDARD_ALPHABET
  else:
    alphabet = MIXED_ALPHABET
  if params_dict["usegap"] :
    alphabet.add('-')

  return(seq_dict, len_ref_seq, alphabet, nseqs0, params_dict, there_are_duplicates)

def CIC(seq_dict, len_ref_seq, alphabet, params_dict) :
  seq_list = seq_dict.keys()
  ntaxa = len(seq_list)
  if ntaxa == 1 :
    return(0, 0, 0, 0, 0, 0)
  if params_dict["diversity_cap"] :
    diversity_cap = min(int(0.75*ntaxa), int(0.75*len(alphabet)))
  else:
    diversity_cap = 100  # Any large number will do!
  tax_sets_dict = {} # The splits induced by characters
  totalCIC = 0.0
  CITE = 0.0
  n_informative_sites = 0
  dCITE = 0.0
  n_d_informative_sites = 0
  deltamin = 0
  for col in range(len_ref_seq) :
    nletters = 0     # nletters will be the count of seqs if gaps are counted, or will vary dept on gaps
    char_IDset_dict = {}   # char -> ID list (for just this site)
    for t in range(ntaxa) :
      char = seq_list[t][col]
      if char in alphabet :
	nletters += 1 
	if char_IDset_dict.has_key(char) :
	  char_IDset_dict[char].append(seq_dict[seq_list[t]][0])
	else:
	  char_IDset_dict[char] = [seq_dict[seq_list[t]][0]]

    for char in char_IDset_dict.keys() :
      char_IDset_dict[char].sort()

    Sc_1 = 0
    for ch in char_IDset_dict.keys() :
      if len(char_IDset_dict[ch]) > 1 :
        Sc_1 += 1
    deltamin +=  max(Sc_1 -1, 0)

    CIC = 0.0
    if nletters > 3 :
      for j in range(3, nletters - len(char_IDset_dict.keys()) + 2) :
        bj = 0
        for ch in char_IDset_dict.keys() :
          if len(char_IDset_dict[ch]) >= j :
            bj += 1
        firstbit = math.log(2 * j - 3, 2)
        CIC +=  firstbit * (1-bj)
        # print 'col', col, 'j', j, 'bj', bj, 'firstbit', firstbit, 'CIC', CIC, "Char Set", char_IDset_dict.keys()
    totalCIC += CIC

    Entropy = 0.0
    if params_dict["noparsimony"] :
      min_letters = 1
    else:
      min_letters = 3
    if nletters > min_letters :
      proportions = {}
      for ch in char_IDset_dict.keys() :
        if len(char_IDset_dict[ch]) > 1 or params_dict["noparsimony"]:
          proportions[ch] = float(len(char_IDset_dict[ch])) / nletters
	else:
	  del char_IDset_dict[ch]   # splits with only 1 member are deleted from further consideration this site
      if len(proportions) > 1 and len(proportions) <= diversity_cap :
	n_informative_sites += 1
	for letter, p in proportions.items() :
          Entropy += -p*math.log(p, 2)
	  # print letter, p
        CITE += Entropy

      card_of_split = len(char_IDset_dict) 
      list_of_split = char_IDset_dict.values()  # list of ID lists
      list_of_split.sort()  # each level must be separately sorted
      if card_of_split > 1 and len(proportions) <= diversity_cap and not seen_before(list_of_split, card_of_split, tax_sets_dict) :
	n_d_informative_sites += 1
	dCITE += Entropy
	if tax_sets_dict.has_key(card_of_split) :
 	  tax_sets_dict[card_of_split].append(list_of_split)
        else:
	  tax_sets_dict[card_of_split] = [list_of_split]

  return(totalCIC, deltamin, CITE, n_informative_sites, dCITE, n_d_informative_sites)

"""
given a sorted list of ID lists, and a dictionary mapping to sorted lists of ID lists,
is the input list of ID lists found among similar sized lists contained in the dictionary?
"""

def seen_before(list_of_split, card_of_split, tax_sets_dict) :
  if not tax_sets_dict.has_key(card_of_split) :
    return(False)

  list_of_split_list = tax_sets_dict[card_of_split]
  for list_split_i in list_of_split_list :
    same = True
    for j in range(card_of_split):
      if list_of_split[j] != list_split_i[j] :
	same = False
	break
    if same:
      return(True)
  return(False)
    
if __name__ == "__main__" :
  seq_dict, seqlen, alphabet, nseqs0, params_dict, there_are_duplicates = init()
  outfile = params_dict["o"]
  totalCIC, deltamin,  CITE, n_informative_sites, dCITE, n_d_informative_sites = CIC(seq_dict, seqlen, alphabet, params_dict)
  outfile.write("Total CIC %0.2f Length %d  Nseqs %d From %d starting sequences \n"\
		 % (totalCIC, seqlen, len(seq_dict), nseqs0))
  outfile.write("Total Entropy %0.2f from %d informative sites\n" % (CITE, n_informative_sites))
  outfile.write("Deflated Total Entropy %0.2f from %d informative sites\n" % (dCITE, n_d_informative_sites))
  if params_dict["add_deltamin"] :
    outfile.write("Delta min %d\n" % deltamin)
  if there_are_duplicates :
    for seq, IDlist in seq_dict.items() :
      if len(IDlist) > 1 :
	outfile.write("%s has duplicates:\n" % IDlist[0])
	for ID in IDlist[1:] :
	  outfile.write("\t%s\n" % ID)
