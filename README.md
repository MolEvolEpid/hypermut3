## Hypermut 3.0

Web version (also has more details about method): https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html

**Reference for Hypermut 1.0:**
Rose, PP and Korber, BT. 2000. Detecting hypermutations in viral sequences with an emphasis on G -> A hypermutation. Bioinformatics 16(4): 400-401.

**Primary purpose:** Analysis and detection of APOBEC-induced hypermutation. 
See [here](https://www.hiv.lanl.gov/content/sequence/HYPERMUT/Readme.html) for more details on hypermutation. 

**General purpose:** To document the nature and context of nucleotide substitutions in a sequence population relative to a reference sequence.

## Overview

Hypermut 3.0 allows searching for mutations fitting a pattern you specify. 
The positions that match the upstream context pattern, followed by the specified mutation (relative to the reference sequence, 
assumed to be the first entered, and treated as ancestral) followed by the downstream context will be found. 
Likewise, matches to the control pattern will be shown for comparison. 
The context requirements can be enforced on the reference sequence, or on the query sequence (recommended, especially if the reference is distant) or both. 
Normally only the contexts should differ between the control pattern and the test pattern. 
Fisher's exact test is then used to detect any increase of mutation for the specified context compared to the control context.

## Installation

Hypermut 3.0 is written in Python3 and has the following package requirements: sys, getopt, re, scipy.stats, itertools, warnings.

To clone this repo:

```
git clone https://github.com/MolEvolEpid/hypermut
```

## Running Hypermut 3.0

Aligned sequences in fasta file format are required as input.
The first sequence in the alignment will be used as the reference sequence, and each of the other sequences will be used as a query sequence. 
Please choose the reference sequence carefully (see details below). 

To search for hypermutation by APOBEC3G or APOBEC3F using the example fasta file, you can run the command:

```
python mutsearch.py -u example_summary_output.csv -o example_verbose_output.txt 'G,A,G,A,.,RD,.,YN|RC' < example.fasta
```

The input string takes the following form:

```
'mutfrom,mutto,controlfrom,controlto,primaryupstream,primarydownstream,controlupstream,controldownstream'
```

The script can also take the following optional command line arguments:

```
-h, --help         help menu
-s, --start        position at which to start searching for mutations (default: 0)
-f, --finish       position at which to end searching for mutations (default: end of sequence)
-e, --enforce      what sequence to enforce the context on: ancestor (A), descendant (D), or both (B) (default: A)
-m, --multistate   how to treat multistate characters (and - all nucleotides are present, or - one of the nucleotides is present, ignore - don't consider mutation sites with multistate characters) (default: and)
-o, --outfile      verbose output file including potential sites and whether there was the correct mutation at those sites
-u, --summaryfile  summary of mutation counts and potential sites for primary and control contexts
```

## Details

- Context:
  - As in regular expressions, the symbol "|" means "OR". Thus GGT|GAA matches GGT or GAA.
  - () can be used for grouping (i.e., one could also write G(GT|AA).
  - All of the IUPAC codes are supported (e.g., R means G or A, while D means not C and a vertical bar ("|") means "OR".
  - Contexts can be multiple characters, but mutations can only be one character. 
  - For technical reasons, the upstream context pattern must always match a fixed number of nucleotides.
    For example, A|(TC) is not allowed as an upstream pattern because it could have length 1 or 2.
  - The primary and control contexts cannot be overlapping.
- Reference sequence:
  - The first sequence in the fasta file.
  - Can only contain non-multistate characters and gaps (`-`).
  - For an intrapatient set, the reference could be the consensus of all the sequences, assuming that the majority are not hypermutated.
  - For a set of unrelated sequences, the reference should probably be the consensus sequence for the appropriate subtype.
- Query sequence(s):
  - Can contain IUPAC characters and gaps (`-`).
  - If the query sequence contains multistate characters, they can be treated as follows:
    - All multistate characters are considered to be present in the sample ("and"). This makes sense if it is from sequencing of a population.
      In this case, if the correct mutation is present in the multistate character, the match is considered partial (1/n_states_in_multistate_char). 
    - A multistate character means that one of the bases is present in the sample ("or"). This makes sense if it is from sequencing a single clone.
      In this case, if the correct mutation is present in the multistate character, the match is considered a complete match (1).
    - Ignore any contexts where the mutation is a multistate character.
      This makes sense if you have single clones and don't want to overcount potential matches.
      In this case, these sites are not considered potential mutations even if the context is correct.
  - Contexts where the mutation in the query is a gap are ignored and not considered potential mutations.

## Output

There are two outputs:

- Summary output:
  - 1 row for each sequence
  - Columns for:
    - Sequence name
    - Number of actual primary mutations 
    - Number of potential primary mutations 
    - Number of actual control mutations 
    - Number of potential control mutations
    - Rate ratio of primary vs. control
    - Fisher's exact p-value
- Verbose output:
  - Input pattern
  - Regular expression pattern
  - Row for each potential mutation site (with the correct context) including the columns:
    - Sequence number
    - Sequence name
    - Whether the site matches the control pattern/context (1) or primary pattern/context (0)
    - Potential mutation site
    - Whether the expected mutation was present (1) or not (0)

