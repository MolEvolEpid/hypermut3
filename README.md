## Hypermut 3.0

Web version (also has more details about method): https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html

**Reference for Hypermut 1.0:**
Rose, PP and Korber, BT. 2000. Detecting hypermutations in viral sequences with an emphasis on G -> A hypermutation. Bioinformatics 16(4): 400-401.\
https://academic.oup.com/bioinformatics/article/16/4/400/187265

**Primary purpose:** Analysis and detection of APOBEC3F- and APOBEC3G-induced hypermutation. 
See [here](https://www.hiv.lanl.gov/content/sequence/HYPERMUT/Readme.html) for more details on hypermutation. 

**General purpose:** To document the nature and context of nucleotide substitutions in a sequence population relative to a reference sequence.

## Overview

Hypermut 3.0 allows searching for mutations fitting a pattern you specify. 
The positions that match the upstream context pattern, followed by the specified mutation (relative to the reference sequence, 
assumed to be the first entered, and treated as ancestral) followed by the downstream context will be found. 
Matches to the opposite control pattern will be shown for comparison. 
The context requirements can be enforced on the reference sequence, or on the query sequence (recommended, especially if the reference is distant) or both. 
Fisher's exact test is then used to detect any increase of mutation for the specified context compared to the control context.

## Installation

Hypermut 3.0 is written in Python3 and requires the `scipy` package. 

To clone this repo:

```
git clone https://github.com/MolEvolEpid/hypermut
```

The `hypermut_env.yaml` file can be used to create a conda environment with the (very minimal) requried dependencies, if desired. First, install [conda](https://github.com/conda-forge/miniforge). Then, run the following command from the `hypermut` directory:

```
mamba env create -f hypermut_env.yaml
```

This will create a hypermut conda environment that can be activated using:

```
conda activate hypermut
```

## Running Hypermut 3.0

Aligned sequences in fasta file format are required as input.
The first sequence in the alignment will be used as the reference sequence, and each of the other sequences will be used as a query sequence. 
Please choose the reference sequence carefully (see details below). 

To search for hypermutation by APOBEC3G or APOBEC3F using the example fasta file, you can run the command:

```
python hypermut.py example/example.fasta G A . RD -s example_summary_output.csv -p example_positions_output.txt
```

The positional inputs are as follows:

```
  fasta                 Alignment file in fasta format
  mutationfrom          Base in the reference to consider as a site of interest for nucleotide substitution
  mutationto            Base in the query to consider a nucleotide substitution of interest
  upstreamcontext       Upstream nucleotide context of interest
  downstreamcontext     Downstream nucleotide context of interest
```

The optional arguments include:

```
-h, --help            show this help message and exit
--positionsfile POSITIONSFILE, -p POSITIONSFILE
                      Optional file path to output potential mutation sites and
                      whether there was the correct mutation at those sites
--summaryfile SUMMARYFILE, -s SUMMARYFILE
                      Optional file path to output a summary of mutation counts and
                      potential sites for primary and control contexts
--enforce {A,D,B}, -e {A,D,B}
                      What sequence to enforce the context on:
                      ancestor/reference (A), descendant/query (D, default), or both (B)
--match {strict,partial}, -m {strict,partial}
                      Whether to include only complete matches (strict, default),
                      or also include partial matches (not completely overlapping
                      bases between query and context, partial)
--keepgaps, -k        Flag indicating to keep gaps in the alignment when i
                      dentifying pattern matches (default without flag is to remove gaps)
--begin BEGIN, -b BEGIN
                      Position at which to start searching for mutations (default: 0).
                      Note that the context may fall outside of these positions.
--finish FINISH, -f FINISH
                      Position at which to end searching for mutations (default: end of sequence).
                      Note that the context may fall outside of these positions.
```

## Details

- Context:
  - As in regular expressions, the symbol "|" means "OR". Thus GGT|GAA matches GGT or GAA.
  - Unlike Hypermut 2.0, () **CANNOT** be used for grouping (i.e.,  G(GT|AA) is wrong, instead use GGT|GAA).
  - All of the IUPAC codes are supported (e.g., R means G or A, while D means not C) and a vertical bar ("|") means "OR".
  - Contexts can be multiple characters, but mutations can only be one character. 
  - The upstream context patterns must always match a fixed number of nucleotides.
    For example, A|TC is not allowed as a pattern because it could have length 1 or 2.
- Reference sequence:
  - The first sequence in the fasta file.
  - Can only contain non-multistate characters (ACGT) and gaps (`-`).
  - For an intrapatient set, the reference could be the consensus of all the sequences, assuming that the majority are not hypermutated.
    - For more details about consensus making, and a webtool, see [here](https://www.hiv.lanl.gov/content/sequence/CONSENSUS/consensus.html).
  - For a set of unrelated sequences, the reference should probably be the consensus sequence for the appropriate subtype.
    - For pre-made subtype consensus sequences for HIV, see [here](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html). 
  
- Query sequence(s):
  - Can contain [IUPAC nucleotide codes](https://www.bioinformatics.org/sms/iupac.html) (T, not U) and gaps (`-`).
  - Contexts where the mutation in the query is a gap are ignored and not considered potential mutations.
  - If the query sequence contains multistate characters, they can be treated as follows: **ADD FIGURE FROM MANUSCRIPT ONCE COMPLETE**
    - **Strict** (default): Only completely inclusive matches containing multistate characters are considered (for the mutation and the context). 
     - For a mutation site, the entire site is not considered if there is a partial match, e.g. if the context is correct but the primary mutation is `A` and the query mutation is `R`. 
      - For the context, if the primary downstream context is `DT`, then `RT` would be considered the correct context. However, `NT` would not be considered the correct context. 
      - This makes sense if the sequencing is from single clones and you don't want to consider ambiguous matches.
    - **Partial**: Partially overlapping matches (for the mutation and the context) are considered using the equation: **ADD EQUATION FROM PAPER ONCE COMPLETE**.  
       - For a mutation site, if the primary mutation is `A` and the query mutation is `R`, then this would be considered a 50% match. 
      - For the context, if the primary downstream context is `DT`, then a query `NT` context would be split between primary (75%) and control (25%) patterns. 
      - This makes sense if the sequence is derived from a population.
 

## Output

There are two outputs:

- Summary output:
  - 1 row for each sequence
  - Columns for:
    - Sequence name
    - Number of actual primary mutations 
    - Number of potential primary mutations (correct context)
    - Number of actual control mutations 
    - Number of potential control mutations (correct context)
    - Rate ratio of primary vs. control
    - Fisher's exact p-value
- Verbose output:
  - Input pattern
  - Regular expression pattern
  - Row for each potential mutation site (with the correct context) including the columns:
    - Sequence number
    - Sequence name
    - Whether the site matches the control pattern/context (1) or primary pattern/context (0)
    - Proportion of the site that matches the control or primary pattern/context
    - Potential mutation site
    - Whether the expected mutation was present or not

## Example code for cumulative plot

Sometimes it is useful to look at the plot of cumulative number of potential match sites vs. cumulative number of actual matches. Here is R code that you can use to create this plot:

```
# load library
library(tidyverse)

# read in positions file
positions <- read_csv('example_positions_output.csv', comment = '#')

# cumulative plot (for primary)
positions %>% 
  filter(control == 0) %>% 
  arrange(potential_mut_site) %>% 
  group_by(seq_name) %>% 
  mutate(cum_potential = cumsum(prop_control),
         cum_match = cumsum(mut_match)) %>% 
  ggplot(aes(x = cum_potential, y = cum_match, col = seq_name)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Cumulative number of potential sites', y = 'Cumulative number of matches', col = '')
```

Here is a comparison of the example data run in strict and partial modes using code similar to the above (see the `example` folder and corresponding `README` for more details on exactly how this was generated):

![image](example/example.pdf)


## Tests

To run the unit tests for the functions used in `hypermut.py`, you need `pytest` (which is included in the conda environment). Then you can run the command:

```
pytest test_hypermut.py
```

To also get the code coverage, run:
```
pytest test_hypermut.py --cov --cov-report=html
```

You can open 'htmlcov/index.html' to browse the code coverage. 
