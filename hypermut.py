#!/usr/bin/env python
###!/sw/bin/python3

# Finds matches and prints (1) summary to stdout, and (2) match locations to a file
# Note main argument now has 8 fields.
   
import sys, getopt
import re
from scipy.stats import fisher_exact
from itertools import product
import warnings
from math import prod

usage="usage: hypermut.py [-s start] [-f finish] [-h] [-e (A|B|D)] [-i (True/False)] [-m (strict|partial)] [-o outfile] [-u summaryfile] 'primaryfrom,primaryto,primaryupstream,primarydownstream' < inputseqs.fasta"

def check_width(regexpstring, culprit):
    try:
        tmp=re.compile("(?<="+regexpstring + ")")
    except:
        raise ValueError(f"{culprit} pattern must correspond to a fixed length expression (example: T|GC is not allowed) Please try again using a {culprit} pattern that has only one possible width.")

def check_chars(chars, good_chars, error_message):
    bad_chars = [x for x in set(chars) if x not in good_chars]
    if len(bad_chars):
        raise ValueError(f"{error_message}. Yours contains: {bad_chars}")
    
def check_input_patterns(patterns, iupac_dict):
    # check that all patterns consist of IUPAC characters
    check_chars(patterns, list(iupac_dict.keys())+[',','.','|'],
                "All patterns and mutations must include only IUPAC characters, '.', and '|'.")

    (primaryfrom, primaryto, primaryupstream, primarydownstream)=str.split(patterns, ",")

    # stop with error if any pattern isn't fixed width
    check_width(primaryupstream, "upstream")
    check_width(primarydownstream, "downstream")

    try:
        mutfrom = iupac_dict[primaryfrom]
        mutto = iupac_dict[primaryto]
    except:
        # require mutation to be only one base
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

def slice_seq(sequence, start, context_len, context_type, ignore_gaps):
    if context_type == 'upstream':
        if ignore_gaps:
            seq_sliced = list(re.sub('-', '', sequence[:start])[-context_len:])
        else:
            seq_sliced = list(sequence[:start][-context_len:])
    if context_type == 'downstream':
        if ignore_gaps:
            seq_sliced = list(re.sub('-', '', sequence[start+1:])[:context_len])
        else:
            seq_sliced = list(sequence[start+1:][:context_len])
    return seq_sliced

def find_match_weight(refseq, sequence, start, end, 
                      mutto, upstream_context, downstream_context, 
                      enforce, iupac_dict, multistate, ignore_gaps):
    base = sequence[start:end]
    site_primary = matchval_primary = site_control = matchval_control = 0
    if base != '-':
      upstream_ref = upstream_seq = downstream_ref = downstream_seq = []
      up_same_len = down_same_len = True
      if enforce == 'D':
        if upstream_context[0] != '.':
            upstream_seq = slice_seq(sequence, start, len(upstream_context[0]), 'upstream', ignore_gaps) 
        if downstream_context[0] != '.':
            downstream_seq = slice_seq(sequence, start, len(downstream_context[0]), 'downstream', ignore_gaps) 
            down_same_len = len(downstream_context[0]) == len(downstream_seq)
      elif enforce == 'A':
        if upstream_context[0] != '.':
            upstream_ref = slice_seq(refseq, start, len(upstream_context[0]), 'upstream', ignore_gaps) 
            up_same_len = len(upstream_context[0]) == len(upstream_ref)
        if downstream_context[0] != '.':
            downstream_ref = slice_seq(refseq, start, len(downstream_context[0]), 'downstream', ignore_gaps) 
            down_same_len = len(downstream_context[0]) == len(downstream_ref)
      elif enforce == 'B':
          if upstream_context[0] != '.':
            upstream_ref = slice_seq(refseq, start, len(upstream_context[0]), 'upstream', ignore_gaps) 
            upstream_seq = slice_seq(sequence, start, len(upstream_context[0]), 'upstream', ignore_gaps) 
            up_same_len = len(upstream_context[0]) == len(upstream_ref) and len(upstream_context[0]) == len(upstream_seq)
          if downstream_context[0] != '.':
            downstream_ref = slice_seq(refseq, start, len(downstream_context[0]), 'downstream', ignore_gaps) 
            downstream_seq = slice_seq(sequence, start, len(downstream_context[0]), 'downstream', ignore_gaps) 
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
      if not ignore_gaps:
        if '-' in upstream_ref+upstream_seq+downstream_ref+downstream_seq:
            site_primary = matchval_primary = site_control = matchval_control = 0
      if multistate == "strict" and (int(site_primary) != site_primary or int(matchval_primary) != matchval_primary): # ignore if multistate strict mode and not complete overlap
        site_primary = matchval_primary = site_control = matchval_control = 0
    return site_primary,matchval_primary,site_control,matchval_control

def summarize_matches(refseq, queryseq, start, finish, 
                      potentialre, enforce,
                      mutto_orig, upstream_orig, downstream_orig, 
                      iupac_dict, multistate, ignore_gaps, seqs, name):
    sites_primary = matches_primary = sites_control = matches_control = 0  # sites= potentials also passing secondre
    if finish is None:
        potentials=potentialre.finditer(refseq,start)
    else:
        potentials=potentialre.finditer(refseq,start,finish)
    for mymatch in potentials:
        site_primary, matchval_primary, site_control, matchval_control = \
            find_match_weight(refseq, queryseq, mymatch.start(), mymatch.end(), 
                              mutto_orig, upstream_orig, downstream_orig, 
                              enforce, iupac_dict, multistate, ignore_gaps)
        sites_primary+=site_primary
        matches_primary+=matchval_primary
        sites_control+=site_control
        matches_control+=matchval_control
        if(sys.stdout):
            if site_primary != 0:
                print(str(seqs) + "," + name + ",0," + str(mymatch.start()+1) + "," + str(matchval_primary))
            if site_control != 0:
                print(str(seqs) + "," + name + ",1," + str(mymatch.start()+1) + "," + str(matchval_control))
    return sites_primary,matches_primary,sites_control,matches_control

def calc_fisher(primarysites, primaries, controlsites, controls):
    return fisher_exact([[primaries, primarysites-primaries],[controls,controlsites-controls]], alternative = 'greater')
 
# main starts here #
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:s:f:e:m:o:u:i:", ["help", "start=","finish=","enforce=","multistate=","outfile=", "summaryfile=", "ignoregaps="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))
        print(usage)
        sys.exit(2)

    # defaults
    start=0
    finish=outfile=sumfile=None
    enforce="D"  # first letter of Ancestor, Descendant, or Both
    multistate = "strict" # strict or partial
    ignore_gaps = True # True or False

    for o, a in opts:
        if o in ("-o", "--outfile"):
            outfile=open(a,'w')
        if o in ("-u", "--summaryfile"):
            sumfile=open(a,'w')
        if o in ("-h", "--help"):
            print(usage)
            sys.exit()
        if o in ("-s", "--start"):
            try:
                start = int(a)
            except ValueError:
                print("Start must be integer")
            if start < 0:
                print("Don't understand negative start")
        if o in ("-f", "--finish"):
            try:
                finish = int(a)
            except ValueError:
                print("Finish must be integer")
            if finish < 0:
                print("Can't understand negative finish")
        if o in ("-e", "--enforce"):
            enforce = a[0].upper()
            if enforce not in ['A', 'D', 'B']:
                raise ValueError("enforce must be A, D, or B")
        if o in ("-m", "--multistate"):
            multistate=a.lower()
            if multistate not in ['strict', 'partial']:
                raise ValueError("multistate must be strict or partial")
        if o in ("-i", "--ignoregaps"):
            if a[0].upper() not in ['T', 'F']:
                raise ValueError('ignoregaps must be True or False')
            ignore_gaps = a[0].upper() == 'T'

    if not (len(args)==1):
        print(usage)
        sys.exit()

    if sumfile:
        sumfile.write("Sequence,Primary Matches,Out of (Potential Primary Sites),Control Matches,Out of (Potential Controls),Rate Ratio(A/B)/(C/D),Fisher Exact P-value\n")

    origstdout=sys.stdout
    sys.stdout=outfile

    arg=args[0]
    if(sys.stdout):
        print("#arg=", arg)

    # only allow partial matches when context is enforced on query sequence only
    if multistate == 'partial' and enforce != 'D':
        raise ValueError("When multistate is partial, enforce must be D.")

    arg=arg.upper() # make sure all input bases are uppercase
    iupac_dict = {"A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
                  "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
                  "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
                  "N": ["A", "C", "G", "T"], "-": ["-"]}
    
    check_input_patterns(arg, iupac_dict)

    # original patterns
    (primaryfrom_orig, primaryto_orig, primaryupstream_orig, primarydownstream_orig)=str.split(arg, ",")
    # convert different context options into a list
    primaryupstream_orig = primaryupstream_orig.split('|')
    primarydownstream_orig = primarydownstream_orig.split('|')

    # prep pattern for allowing gaps and multistate characters in the regular expression
    # primaryfrom will only match ACGT regardless because no multistate characters are allowed in the reference sequence
    primaryfromre = re.compile('['+''.join(iupac_dict[primaryfrom_orig])+']',re.I)

    if(sys.stdout):
        print("#regexps=",arg)
        print("seq_num,seq_name,control,potential_mut_site,mut_match")    

    # start reading in fasta file
    line=sys.stdin.readline()
    refname=line[1:].strip()
    refseq=""
    line=sys.stdin.readline()
    while line and line[0]!=">":
        refseq=refseq+line
        line=sys.stdin.readline()
    refseq=refseq.replace("\n","").upper()
    # check reference sequence
    check_chars(refseq, list('ACGT-'), "The reference sequence must contain only the following characters: ACGT-.")

    seqs=0
    while line:
        name=line[1:].strip()
        sequence=""
        line=sys.stdin.readline()
        while line and line[0]!=">":
            sequence=sequence+line
            line=sys.stdin.readline()
        sequence=sequence.replace("\n","").upper()
        check_chars(sequence, list(iupac_dict.keys())+['-'], 'Query sequences must contain only IUPAC characters or - (for gap).')

        seqs+=1
        primarysites, primaries, controlsites, controls = \
            summarize_matches(refseq, sequence, start, finish, 
                              primaryfromre, enforce,
                              primaryto_orig, primaryupstream_orig, primarydownstream_orig, 
                              iupac_dict, multistate, ignore_gaps, seqs, name)

        sys.stdout=origstdout
        odds_ratio, pval = calc_fisher(primarysites, primaries, controlsites, controls)

        try:
            ratio = "%0.2f" %(primaries*controlsites/(1.0*primarysites*controls))
        except:
            if primaries*controlsites > 0:
                ratio="inf" 
            else:
                ratio="undef" 
    
        if sumfile:
            s = name+","+str(primaries)+","+str(primarysites)+","+str(controls) +","+str(controlsites) +","+ratio +",%.6g" %float(pval)+'\n' 
            sumfile.write(s)

        sys.stdout.flush()
        sys.stdout=outfile

    sys.stdout=origstdout

    if sumfile:
        sumfile.close()
    if outfile:
        outfile.close()
    origstdout.close()