#!/usr/bin/env python
###!/sw/bin/python3

# Finds matches and prints (1) summary to stdout, and (2) match locations to a file
# Note main argument now has 8 fields.
   
import argparse
import re
from scipy.stats import fisher_exact
from itertools import product
import warnings
from math import prod

def check_width(regexpstring, culprit):
    try:
        tmp=re.compile("(?<="+regexpstring + ")")
    except:
        raise ValueError(f"{culprit} pattern must correspond to a fixed length expression (example: T|GC is not allowed) Please try again using a {culprit} pattern that has only one possible width.")

def check_chars(chars, good_chars, error_message):
    bad_chars = [x for x in set(chars) if x not in good_chars]
    if len(bad_chars):
        raise ValueError(f"{error_message}. Yours contains: {bad_chars}")
    
def check_input_patterns(mutfrom, mutto, primaryupstream, primarydownstream, iupac_dict):
    # check that all patterns consist of IUPAC characters
    check_chars(primaryupstream, list(iupac_dict.keys())+['.','|'],
                "The upstream pattern must include only IUPAC characters, '.', and '|'.")
    check_chars(primarydownstream, list(iupac_dict.keys())+['.','|'],
                "The downstream pattern must include only IUPAC characters, '.', and '|'.")
    # stop with error if any pattern isn't fixed width
    check_width(primaryupstream, "upstream")
    check_width(primarydownstream, "downstream")
    # require mutation to be only one base
    try:
        mutfrom = iupac_dict[mutfrom]
        mutto = iupac_dict[mutto]
    except:
        raise ValueError("Mutation from and to must each be a single IUPAC character.")        
    # check for (undesired) overlap in mutation from and to
    if sum([1 for x in mutfrom if x in mutto]) != 0:
        warnings.warn("Mutation from and to have overlapping bases.")
    # check for (undesired) overlap in context
    primary_pattern = re.sub('\\.', "", '|'.join([y + x for x in primarydownstream.split('|') for y in primaryupstream.split('|')]))
    base_info_primary = [[iupac_dict[y] for y in list(x)] for x in primary_pattern.split('|')]
    contexts_primary = [''.join(list(y)) for x in base_info_primary for y in product(*x)]
    if len(contexts_primary) != len(set(contexts_primary)):
        raise ValueError("Context is redundant. Please provide non-redundant patterns.")

def compute_context_prop(refseq, seq, context, enforce, iupac_dict):
    prop = 1
    if context[0] != '.':
      prop = 0
      if enforce == 'D':
        for u in context:
            prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(seq, list(u))])
      elif enforce == 'A':
        for u in context:
            prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(refseq, list(u))])
      elif enforce == 'B':
        # this is only allowed in the strict context, so it should ultimately be 0 or 1
        for u in context:
            prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(seq, list(u))])
            prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(refseq, list(u))])
        prop = prop/2
    return prop

def slice_seq(sequence, start, context_len, context_type, keep_gaps):
    if context_type == 'upstream':
        if keep_gaps:
            seq_sliced = list(sequence[:start][-context_len:])
        else:
            seq_sliced = list(re.sub('-', '', sequence[:start])[-context_len:])
    if context_type == 'downstream':
        if keep_gaps:
            seq_sliced = list(sequence[start+1:][:context_len])
        else:
            seq_sliced = list(re.sub('-', '', sequence[start+1:])[:context_len])
    return seq_sliced

def find_match_weight(refseq, sequence, start, end, 
                      mutto, upstream_context, downstream_context, 
                      enforce, iupac_dict, match, keep_gaps):
    base = sequence[start:end]
    site_primary = matchval_primary = site_control = matchval_control = 0
    if base != '-':
      upstream_ref = upstream_seq = downstream_ref = downstream_seq = []
      up_same_len = down_same_len = True
      if enforce == 'D':
        if upstream_context[0] != '.':
            upstream_seq = slice_seq(sequence, start, len(upstream_context[0]), 'upstream', keep_gaps) 
        if downstream_context[0] != '.':
            downstream_seq = slice_seq(sequence, start, len(downstream_context[0]), 'downstream', keep_gaps) 
            down_same_len = len(downstream_context[0]) == len(downstream_seq)
      elif enforce == 'A':
        if upstream_context[0] != '.':
            upstream_ref = slice_seq(refseq, start, len(upstream_context[0]), 'upstream', keep_gaps) 
            up_same_len = len(upstream_context[0]) == len(upstream_ref)
        if downstream_context[0] != '.':
            downstream_ref = slice_seq(refseq, start, len(downstream_context[0]), 'downstream', keep_gaps) 
            down_same_len = len(downstream_context[0]) == len(downstream_ref)
      elif enforce == 'B':
          if upstream_context[0] != '.':
            upstream_ref = slice_seq(refseq, start, len(upstream_context[0]), 'upstream', keep_gaps) 
            upstream_seq = slice_seq(sequence, start, len(upstream_context[0]), 'upstream', keep_gaps) 
            up_same_len = len(upstream_context[0]) == len(upstream_ref) and len(upstream_context[0]) == len(upstream_seq)
          if downstream_context[0] != '.':
            downstream_ref = slice_seq(refseq, start, len(downstream_context[0]), 'downstream', keep_gaps) 
            downstream_seq = slice_seq(sequence, start, len(downstream_context[0]), 'downstream', keep_gaps) 
            down_same_len = len(downstream_context[0]) == len(downstream_ref) and len(downstream_context[0]) == len(downstream_seq)
      if up_same_len and down_same_len:
        upstream_seq_prop_primary = compute_context_prop(upstream_ref, upstream_seq, upstream_context, enforce, iupac_dict)
        downstream_seq_prop_primary = compute_context_prop(downstream_ref, downstream_seq, downstream_context, enforce, iupac_dict)
        # 1 means complete primary or control site, fraction means partial primary or control site
        site_primary = upstream_seq_prop_primary * downstream_seq_prop_primary
        site_control = 1-site_primary
        chars_seq = iupac_dict[base]
        chars_mut = iupac_dict[mutto]
        correct_mut = sum([x in chars_mut for x in chars_seq])/len(chars_seq)
        matchval_primary=correct_mut*site_primary
        matchval_control=correct_mut*site_control
      if keep_gaps:
        if '-' in upstream_ref+upstream_seq+downstream_ref+downstream_seq:
            site_primary = matchval_primary = site_control = matchval_control = 0
      if match == "strict" and (int(site_primary) != site_primary or int(matchval_primary) != matchval_primary): # ignore if multistate strict mode and not complete overlap
        site_primary = matchval_primary = site_control = matchval_control = 0
    return site_primary,matchval_primary,site_control,matchval_control

def summarize_matches(refseq, queryseq, start, finish, 
                      potentialre, enforce,
                      mutto_orig, upstream_orig, downstream_orig, 
                      iupac_dict, match, keep_gaps, seqs, name, positionsfile):
    sites_primary = matches_primary = sites_control = matches_control = 0  # sites= potentials also passing secondre
    if finish is None:
        potentials=potentialre.finditer(refseq,start)
    else:
        potentials=potentialre.finditer(refseq,start,finish)
    for mymatch in potentials:
        site_primary, matchval_primary, site_control, matchval_control = \
            find_match_weight(refseq, queryseq, mymatch.start(), mymatch.end(), 
                              mutto_orig, upstream_orig, downstream_orig, 
                              enforce, iupac_dict, match, keep_gaps)
        sites_primary+=site_primary
        matches_primary+=matchval_primary
        sites_control+=site_control
        matches_control+=matchval_control
        if(positionsfile is not None):
            if site_primary != 0:
                positionsfile.write(str(seqs) + "," + name + ",0," + str(mymatch.start()+1) + "," + str(matchval_primary)+'\n')
            if site_control != 0:
                positionsfile.write(str(seqs) + "," + name + ",1," + str(mymatch.start()+1) + "," + str(matchval_control)+'\n')
    return sites_primary,matches_primary,sites_control,matches_control

def calc_fisher(primarysites, primaries, controlsites, controls):
    return fisher_exact([[primaries, primarysites-primaries],[controls,controlsites-controls]], alternative = 'greater')

def check_positive(value):
    try:
        ival = int(value)
    except:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if ival != float(value) or ival < 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return(ival)

iupac_dict = {"A": list("A"), "C": list("C"), "G": list("G"), "T": list("T"),
                  "R": list("AG"), "Y": list("CT"), "S": list("GC"), "W": list("AT"), "K": list("GT"), "M": list("AC"),
                  "B": list("CGT"), "D": list("AGT"), "H": list("ACT"), "V": list("ACG"), "N": list("ACGT"), "-": list("-")}

# main starts here #
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Hypermut 3.0', description='Identify mutations in a user-defined context')
    parser.add_argument('fasta', type=str,
                        help="Alignment file in fasta format")
    parser.add_argument('mutationfrom', type=str, 
                        help="Base in the reference to consider as a site of interest for nucleotide substitution")
    parser.add_argument('mutationto', type=str,
                        help="Base in the query to consider a nucleotide substitution of interest")
    parser.add_argument('upstreamcontext', type=str,
                        help="Upstream nucleotide context of interest")
    parser.add_argument('downstreamcontext', type=str,
                        help="Downstream nucleotide context of interest")
    parser.add_argument('--positionsfile', '-p', type=str, 
                        help="Optional file path to output potential mutation sites and whether there was the correct mutation at those sites")
    parser.add_argument('--summaryfile', '-s', type=str,
                        help="Optional file path to output a summary of mutation counts and potential sites for primary and control contexts")
    parser.add_argument('--enforce', '-e', type=str, 
                        choices=['A','D','B'], default='D',
                        help="What sequence to enforce the context on: ancestor/reference (A), descendant/query (D, default), or both (B)")
    parser.add_argument('--match', '-m', type=str, 
                        choices=['strict','partial'], default='strict',
                        help="Whether to include only complete matches (strict, default), or also include partial matches (not completely overlapping bases between query and context, partial)")
    parser.add_argument('--keepgaps', '-k', action='store_true',
                        help="Flag indicating to keep gaps in the alignment when identifying pattern matches (default without flag is to remove gaps)")
    # also check that this is positive
    parser.add_argument('--begin', '-b', type=check_positive, default=0,
                        help="Position at which to start searching for mutations (default: 0)")
    # also check that this is positive
    parser.add_argument('--finish', '-f', type=check_positive, 
                        help="Position at which to end searching for mutations (default: end of sequence)")
    args=parser.parse_args()

    # only allow partial matches when context is enforced on query sequence only
    if args.match == 'partial' and args.enforce != 'D':
        raise ValueError("When match is partial, enforce must be D.")
    
    # make sure all input bases are uppercase
    mutationfrom=args.mutationfrom.upper() 
    mutationto=args.mutationto.upper() 
    upstreamcontext=args.upstreamcontext.upper() 
    downstreamcontext=args.downstreamcontext.upper() 
    
    check_input_patterns(mutationfrom, mutationto, upstreamcontext, downstreamcontext, iupac_dict)

    # convert different context options into a list
    primaryupstream_orig = upstreamcontext.split('|')
    primarydownstream_orig = downstreamcontext.split('|')

    # prep pattern for allowing gaps and multistate characters in the regular expression
    # primaryfrom will only match ACGT regardless because no multistate characters are allowed in the reference sequence
    primaryfromre = re.compile('['+''.join(iupac_dict[mutationfrom])+']',re.I)

    # open fasta file for reading
    fa = open(args.fasta, 'r')

    # prep for writing summary file
    sf = None
    if args.summaryfile is not None:
        sf = open(args.summaryfile, 'w')
        sf.write("Sequence,Primary Matches,Out of (Potential Primary Sites),Control Matches,Out of (Potential Controls),Rate Ratio(A/B)/(C/D),Fisher Exact P-value\n")
    # prep for writing positions file
    pf = None
    if args.positionsfile is not None:
        pf = open(args.positionsfile, 'w')
        pf.write("#regexps="+'from '+mutationfrom+',to '+mutationto+',up '+upstreamcontext+',down '+downstreamcontext+'\n')
        pf.write("seq_num,seq_name,control,potential_mut_site,mut_match\n")    

    # start reading in fasta file
    line=fa.readline()
    if line[0] != '>':
        raise ValueError('Input alignment must be in FASTA format.')
    refname=line[1:].strip()
    refseq=""
    line=fa.readline()
    while line and line[0]!=">":
        refseq+=line
        line=fa.readline()
    refseq=refseq.replace("\n","").upper()
    # check reference sequence
    check_chars(refseq, list('ACGT-'), "The reference sequence must contain only the following characters: ACGT-.")

    seqs=0
    while line:
        name=line[1:].strip()
        sequence=""
        line=fa.readline()
        while line and line[0]!=">":
            sequence=sequence+line
            line=fa.readline()
        sequence=sequence.replace("\n","").upper()
        check_chars(sequence, list(iupac_dict.keys())+['-'], 'Query sequences must contain only IUPAC characters or - (for gap).')

        seqs+=1
        primarysites, primaries, controlsites, controls = \
            summarize_matches(refseq, sequence, args.begin, args.finish, 
                              primaryfromre, args.enforce,
                              mutationto, primaryupstream_orig, primarydownstream_orig, 
                              iupac_dict, args.match, args.keepgaps, seqs, name, pf)

        odds_ratio, pval = calc_fisher(primarysites, primaries, controlsites, controls)
        try:
            ratio = "%0.2f" %(primaries*controlsites/(1.0*primarysites*controls))
        except:
            if primaries*controlsites > 0:
                ratio="inf" 
            else:
                ratio="undef" 
    
        if args.summaryfile is not None:
            sf.write(name+","+str(primaries)+","+str(primarysites)+","+str(controls) +","+str(controlsites) +","+ratio +",%.6g" %float(pval)+'\n')

    # if args.summaryfile is not None:
    #     sf.close()
    # if args.positionsfile is not None:
    #     pf.close()