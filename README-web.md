# Hypermut 3

**Primary purpose:** Analysis and detection of APOBEC3F- and APOBEC3G-induced hypermutation. 
See [here](https://www.hiv.lanl.gov/content/sequence/HYPERMUT/Readme.html) for more details on hypermutation (the "What is hypermutation" section). 

**General purpose:** To document the nature and context of nucleotide substitutions in a sequence population relative to a reference sequence.

**Note:** This page describes the inputs and outputs for the Hypermut webtool. A command line version of the tool can be downloaded from [GitHub](https://github.com/MolEvolEpid/hypermut/tree/main). 

## Overview

Hypermut 3 allows searching for mutations fitting a pattern you specify. 
The positions that match the upstream context pattern, followed by the specified mutation (relative to the reference sequence, 
assumed to be the first entered, and treated as ancestral) followed by the downstream context will be found. 
Matches to the opposite control pattern will be shown for comparison. 
The context requirements can be enforced on the reference sequence, or on the query sequence (recommended, especially if the reference is distant) or both. 
Fisher's exact test is then used to detect any increase of mutation for the specified context compared to the control context.

## Inputs

To run the program successfully:

1. Select the [alignment format](https://www.hiv.lanl.gov/content/sequence/HelpDocs/SEQsamples.html) you are using to present your sequences. 
1. Input your sequence alignment file. **The program designates the first sequence in the file as the reference sequence and considers all other sequences as queries to compare to the reference.** Please choose the reference sequence carefully (see details below). 
1. Choose whether to view the complete sequence (leave boxes blank) or a subregion. If you choose to view a subregion, enter the range of the desired subregion in the boxes.
1. Define the to and from mutations, and the upstream and downstream nucelotide contexts (see below for details). For typical analyses of APOBEC3G- and APOBEC3F-induced hypermutation, these options should be left in their default settings. For more detailed analysis, you can edit the mutation pattern in the provided boxes to search for any desired pattern. Please see the [Hypermut 2.0](https://www.hiv.lanl.gov/content/sequence/HYPERMUT/background.html) description for more information about specific context patterns and what they mean.
1. Decide whether to enforce the context on the reference sequence, query sequence, or both (default: query, recommended). 
1. Decide whether to include only complete matches (strict, default), or also include partial matches (not completely overlapping bases between query and context, partial). Strict matching should be used for sequences that come from single genomes or clones, partial matching should be used for population sequences (see below for details). 
1. Decide whether to skip (default, recommended) or keep gaps in the alignment when identifying pattern matches. 

### Input details

**Mutations:**

- The mutations to and from must only be one nucleotide. 

**Context:**

- As in regular expressions, the symbol "|" means "OR". Thus GGT|GAA matches GGT or GAA.
- Unlike Hypermut 2.0, () **CANNOT** be used for grouping (i.e.,  G(GT|AA) is wrong, instead use GGT|GAA).
- All of the [IUPAC codes](https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html) are supported 
  (e.g., R means G or A, while D means not C).
- Contexts can be multiple characters, but mutations can only be one character. 
- The upstream and downstream context patterns must always match a fixed number of nucleotides.
  For example, A|TC is not allowed as a pattern because it could have length 1 or 2.

**Reference sequence:**

- The first sequence in the fasta file.
- In strict matching mode (see below), can contain IUPAC characters and gaps (`-`). 
- In partial matching mode (see below), can only contain non-multistate characters (ACGT) and gaps (`-`).
- For an intrapatient set, the reference could be the strict/majority consensus of all the sequences, assuming that the majority are not hypermutated.
  - For more details about strict/majority consensus making, and a webtool, see [here](https://www.hiv.lanl.gov/content/sequence/CONSENSUS/consensus.html).
- For a set of unrelated sequences, the reference should probably be the strict/majority consensus sequence for the appropriate subtype.
  - For pre-made subtype strict/majority consensus sequences for HIV, see [here](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html).

**Query sequence(s):**

- Can contain [IUPAC nucleotide codes](https://www.bioinformatics.org/sms/iupac.html) (T, not U) and gaps (`-`).
- Contexts where the mutation in the query is a gap are ignored and not considered potential mutations.
- If the query sequence contains multistate characters, they can be treated as follows: 
  - **Strict** (default): Only completely inclusive matches containing multistate characters are considered (for the mutation and the context). 
    - For a mutation site, the entire site is not considered if there is a partial match, e.g. if the context is correct but the primary mutation is `A` and the query mutation is `R`. 
    - For the context, if the primary downstream context is `DT`, then `RT` would be considered the correct context. However, `NT` would not be considered the correct context. 
    - This makes sense if the sequencing is from single clones and you don't want to consider ambiguous matches.
  - **Partial**: Partially overlapping matches (for the mutation and the context) are considered.  
    - For a mutation site, if the primary mutation is `A` and the query mutation is `R`, then this would be considered a 50% match. 
    - For the context, if the primary downstream context is `DT`, then a query `NT` context would be split between primary (75%) and control (25%) patterns. 
    - This makes sense if the sequence is derived from a population.
 
## Outputs

Hypermut 3 outputs: 

1. A list of input arguments defined.
1. A summary of the sequences and their statistics. The Fisher Exact P-value may be a useful way to determine if a specific sequence is a hypermutant. 
1. A graphical output allows you to view the cumulative number of contexts and mutations across the sequences.

The underlying data corresponding to each of these outputs can be downloaded. For more details about the downloaded files, please refer to the Output section of the Hypermut 3 [GitHub](https://github.com/MolEvolEpid/hypermut/tree/main) page. 
