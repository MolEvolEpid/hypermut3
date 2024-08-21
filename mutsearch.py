#!/usr/bin/env python
###!/sw/bin/python2.4

# this version is just supposed to find matches and print them to a file
# and print a summary to stdout
# Note main argument now has 8 fields.

# Print summary to stdout, locations to outfile
   
import sys, getopt
import re
from scipy.stats import fisher_exact
import itertools
import warnings
from math import prod
usage="usage: mutsearch.py [-s start][-f finish][-h][-e (A|B|D)][-m (strict|partial)] [-o outfile] [-u summaryfile] 'mutfrom,mutto,controlfrom,controlto,mutupstream,mutdownstream,controlupstream,controldownstream' < inputseqs.fasta"

def isfixedwidth(regexpstring):
    try:
        tmp=re.compile("(?<="+regexpstring + ")")
    except:
        return(False)
    return(True)

def die_widthnotfixed(culprit):
    try:
        sys.stdout=origstdout
    except:
        pass
    print("error: ", culprit," pattern must correspond to a fixed length expression")
    print("(example: T|GC is not allowed).")
    print("Please try again using a ", culprit, "pattern that has only one possible width.")
    sys.exit(1)


def widthonly(regexpstring):
    if not isfixedwidth(regexpstring):
        raise ValueError("not fixed width")
    dotstring=""
    while (not isfixedwidth(regexpstring+"|"+dotstring)) and len(dotstring)<500:
        dotstring = dotstring + "."
    if len(dotstring)==500:
        raise ValueError("regexp too long (greater than 500)")
    return dotstring

def check_chars(chars, good_chars, error_message):
    bad_chars = [x for x in set(chars) if x not in good_chars]
    if len(bad_chars):
        raise ValueError(f"{error_message}. Yours contains: {bad_chars}")

def check_width(context):
    context=re.sub('\\.', "", context)
    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(context, ",")
    # stop with error if any pattern isn't fixed width
    if (not isfixedwidth(mutupstream)) or (not isfixedwidth(controlupstream)):
        die_widthnotfixed("upstream")
    if (not isfixedwidth(mutdownstream)) or (not isfixedwidth(controldownstream)):
        die_widthnotfixed("dowwnstream")
    if (not isfixedwidth(mutfrom)) or (not isfixedwidth(mutto)):
        die_widthnotfixed("primary mutation")
    if (not isfixedwidth(controlfrom)) or (not isfixedwidth(controlto)):
        die_widthnotfixed("control mutation")
    if len(widthonly(mutupstream)) != len(widthonly(controlupstream)):
        raise ValueError("Upstream primary and control patterns must be the same length")
    if len(widthonly(mutdownstream)) != len(widthonly(controldownstream)):
        raise ValueError("Downstream primary and control patterns must be the same length")
    # return length of upstream+downstream context
    return(len(widthonly(mutupstream))+len(widthonly(mutdownstream))) 

def check_context(patterns):
    width = check_width(patterns)
    patterns=re.sub('A'," A ",patterns)
    patterns=re.sub('C'," C ",patterns)
    patterns=re.sub('G'," G ",patterns)
    patterns=re.sub('T'," T ",patterns)
    # standard IUPAC codes translated to regexps
    patterns=re.sub('N'," ACGT ",patterns) 
    patterns=re.sub('R'," AG ",patterns)
    patterns=re.sub('Y'," CT ",patterns)
    patterns=re.sub('B'," CGT ",patterns)
    patterns=re.sub('D'," AGT ",patterns)
    patterns=re.sub('H'," ACT ",patterns)
    patterns=re.sub('V'," ACG ",patterns)
    patterns=re.sub('M'," AC ",patterns)
    patterns=re.sub('K'," GT ",patterns)
    patterns=re.sub('S'," CG ",patterns)
    patterns=re.sub('W'," AT ",patterns) 
    patterns=re.sub('\\.', "", patterns)

    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(patterns, ",")    

    mut_pattern = '|'.join([y + x for x in mutdownstream.split('|') for y in mutupstream.split('|')])
    control_pattern = '|'.join([y + x for x in controldownstream.split('|') for y in controlupstream.split('|')])

    base_info_mut = [[list(y) for y in x.split()] for x in mut_pattern.split('|')]
    base_info_control = [[list(y) for y in x.split()] for x in control_pattern.split('|')]

    contexts_mut = [''.join(list(y)) for x in base_info_mut for y in itertools.product(*x)]
    contexts_control = [''.join(list(y)) for x in base_info_control for y in itertools.product(*x)]

    overlap = set(contexts_mut) & set(contexts_control)
    n_patterns = len(set(contexts_mut)) + len(set(contexts_control))
    n_patterns_expected = 4**width
    if len(overlap) != n_patterns_expected:
        if len(overlap):
            raise ValueError(f"Partially overlapping primary and control pattern contexts: {overlap}. This means that positions cannot be categorized uniquely into primary and control groups.")
        if n_patterns != n_patterns_expected:
            warnings.warn(f"Primary and control patterns (n= {str(n_patterns)}) do not create the full complement of possible patterns (n={str(n_patterns_expected)})", stacklevel=2)

def check_input_patterns(patterns, iupac_codes):
    check_chars(patterns, iupac_codes+[',','.','|'],
                "All patterns and mutations must include only IUPAC characters, '.', and '|'.")
    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(patterns, ",")    
    # require mutation to be only one base
    check_chars(mutfrom, iupac_codes, "Primary mutation from must be a single IUPAC character.")
    check_chars(mutto, iupac_codes, "Primary mutation to must be a single IUPAC character.")
    check_chars(controlfrom, iupac_codes, "Control mutation from must be a single IUPAC character.")
    check_chars(controlto, iupac_codes, "Control mutation to must be a single IUPAC character.")
    # check for consistent context width and (undesired) overlap between primary and control contexts
    check_context(patterns)

def allow_gaps_multistate(pattern, multistate):
    # allow for gaps
    pattern=re.sub('A',"[-]*[A][-]*",pattern)
    pattern=re.sub('C',"[-]*[C][-]*",pattern)
    pattern=re.sub('G',"[-]*[G][-]*",pattern)
    pattern=re.sub('T',"[-]*[T][-]*",pattern)
    # standard IUPAC codes translated to regexps, complete matches only
    pattern=re.sub('R',"[-]*[AGR][-]*",pattern)
    pattern=re.sub('Y',"[-]*[CTY][-]*",pattern)
    pattern=re.sub('M',"[-]*[ACM][-]*",pattern)
    pattern=re.sub('K',"[-]*[GTK][-]*",pattern)
    pattern=re.sub('S',"[-]*[CGS][-]*",pattern)
    pattern=re.sub('W',"[-]*[ATW][-]*",pattern)
    pattern=re.sub('B',"[-]*[CGTYBSK][-]*",pattern)
    pattern=re.sub('D',"[-]*[AGTRDWK][-]*",pattern)
    pattern=re.sub('H',"[-]*[ACTHYWM][-]*",pattern)
    pattern=re.sub('V',"[-]*[ACGRVSM][-]*",pattern)
    pattern=re.sub('N',"[-]*[ACGTRYBDHVNWSKM][-]*",pattern)
    pattern=re.sub('\\.', "", pattern)
     # add in partial matches
    if multistate == 'partial': 
        pattern=re.sub(re.escape('[-]*[A][-]*'),"[-]*[ARDHVNWM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*[C][-]*'),"[-]*[CYBHVNSM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*[G][-]*'),"[-]*[GRBDVNSK][-]*",pattern)
        pattern=re.sub(re.escape('[-]*[T][-]*'),"[-]*[TYBDHNWK][-]*",pattern)
        pattern=re.sub(re.escape('[-]*[AGR][-]*'),"[-]*[AGRBDHVNWSKM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*[CTY][-]*'),"[-]*[CTYBDHVNWSKM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*ACM[-]*'),"[-]*[ACMRYBDHVNWS][-]*",pattern)
        pattern=re.sub(re.escape('[-]*GTK[-]*'),"[-]*[GTKRYBDHVNWS][-]*",pattern)
        pattern=re.sub(re.escape('[-]*CGS[-]*'),"[-]*[CGSRYBDHVNKM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*ATW[-]*'),"[-]*[ATWRYBDHVNKM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*CGTYBSK[-]*'),"[-]*[CGTYBSKRDHVNWM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*AGTRDWK[-]*'),"[-]*[AGTRDWKYBHVNSM][-]*",pattern)
        pattern=re.sub(re.escape('[-]*ACTHYWM[-]*'),"[-]*[ACTHYWMRBDVNSK][-]*",pattern)
        pattern=re.sub(re.escape('[-]*ACGRVSM[-]*'),"[-]*[ACGRVSMYBDHNWK][-]*",pattern)

    return(pattern)

def compile_regex(mutation, upstream=None, downstream=None):
    if upstream == None and downstream == None:
        regex = re.compile(mutation,re.I)
    elif upstream != None and downstream != None:
        # use ?= "lookahead" and ?<= lookbehind so we can find overlapping patterns
        regex = re.compile("(?<="+ upstream +")[-]*("+ mutation +")[-]*(?="+ downstream + ")",re.I)
    else:
        ValueError('upstream and downstream must both be None or neither be None.')
    return regex

def get_potential_matches(refseq, start, finish, regex):
    if(finish):
        potentialmatches=regex.finditer(refseq,start,finish)
    else:
        potentialmatches=regex.finditer(refseq,start)
    return(potentialmatches)

def find_match_weight(sequence, start, end, mutto, upstream_context, downstream_context, iupac_dict, multistate):
    base = sequence[start:end]
    if base == '-':
      site = 0
      match_val = 0
    else:
      upstream_seq = None
      if upstream_context[0] != '.':
        upstream_seq = sequence[start-len(upstream_context[0]):start]
      downstream_seq = None
      if downstream_context[0] != '.':
        downstream_seq = sequence[start+1:end+len(downstream_context[0])]
      chars_seq = iupac_dict[base]
      chars_mut = iupac_dict[mutto]
      upstream_seq_prop = 1
      if upstream_context[0] != '.':
        upstream_seq_prop = 0
        for u in upstream_context:
          upstream_seq_prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(list(re.sub('-', '', upstream_seq)), list(u))])
      downstream_seq_prop = 1
      if downstream_context[0] != '.':
        downstream_seq_prop = 0
        for d in downstream_context:
          downstream_seq_prop += prod([len([x for x in iupac_dict[s] if x in iupac_dict[c]])/len(iupac_dict[s]) for s,c in zip(list(re.sub('-', '', downstream_seq)), list(d))])

      # 1 means complete primary or control site, fraction means partial primary or control site
      site = upstream_seq_prop * downstream_seq_prop
      match_val = 0
      if site != 0:
        match_val=sum([x in chars_mut for x in chars_seq])/len(chars_seq)*site
      if multistate == "strict" and (int(site) != site or int(match_val) != match_val): # ignore if multistate strict mode and not complete overlap
        site=0
        match_val=0

    return site,match_val

def summarize_matches(refseq, queryseq, start, finish, 
                      potentialre, secondre, 
                      mutto_orig, upstream_orig, downstream_orig, 
                      iupac_dict, multistate, seqs, name):
    (sites,matches)=(0,0)  # sites= potentials also passing secondre
    potentials=get_potential_matches(refseq,start,finish,potentialre)
    for mymatch in potentials:
        if secondre.match(queryseq,mymatch.start()): #optional match arg
            site, match_val = find_match_weight(queryseq, mymatch.start(), mymatch.end(), 
                                                mutto_orig, upstream_orig, downstream_orig, 
                                                iupac_dict, multistate)
            sites+=site
            matches+=match_val
            if(sys.stdout):
                print(str(seqs) + "," + name + ",0," + str(mymatch.start()+1) + "," + str(match_val))
    return sites,matches

def calc_fisher(mutsites, muts, controlsites, controls):
    return fisher_exact([[muts, mutsites-muts],[controls,controlsites-controls]], alternative = 'greater')
 
# main starts here #

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:f:e:m:o:u:", ["help", "start=","finish=","enforce=","multistate=","outfile=", "summaryfile="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))
        print(usage)
        sys.exit(2)

    # defaults
    start=0
    finish=None
    outfile=None
    sumfile=None
    enforce="D"  # first letter of Ancestor, Descendant, or Both
    multistate = "strict" # strict or partial

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

    if not (len(args)==1):
        print(usage)
        sys.exit()

    if sumfile:
        sumfile.write("Sequence,Muts(Match Sites),Out of(Potential Mut Sites),Controls(Control Muts),Out of(Potential Controls),Rate Ratio(A/B)/(C/D),Fisher Exact P-value\n")

    origstdout=sys.stdout
    sys.stdout=outfile

    arg=args[0]
    if(sys.stdout):
        print("#arg=", arg)

    arg=arg.upper()
    iupac_codes = list('ACGTRYBDHVNWSKM')
    iupac_dict = {"A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
                "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
                "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
                "N": ["A", "C", "G", "T"]}

    check_input_patterns(arg, iupac_codes)

    # original patterns
    (mutfrom_orig, mutto_orig, controlfrom_orig, controlto_orig, mutupstream_orig, mutdownstream_orig, controlupstream_orig, controldownstream_orig)=str.split(arg, ",")
    # convert different context options into a list
    mutupstream_orig = mutupstream_orig.split('|')
    mutdownstream_orig = mutdownstream_orig.split('|')
    controlupstream_orig = controlupstream_orig.split('|')
    controldownstream_orig = controldownstream_orig.split('|')

    # prep pattern for allowing gaps and multistate characters in the regular expression
    # mutfrom and controlfrom will only match ACGT regardless because no multistate characters are allowed in the reference sequence
    arg = allow_gaps_multistate(arg, multistate)
    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(arg, ",")

    if(sys.stdout):
        print("#regexps=",arg)
        print("seq_num,seq_name,control,potential_mut_site,mut_match")    

    # Note that the regexps will be applied to different strings depending on the value of enforce.
    # potential always checked on ancestor
    # second confirms potential in descendant
    any_base = '[ACGTRYBDHVNWSKM]'
    if enforce=="D": # descendant
        potentialmutre=compile_regex(mutfrom) 
        potentialcontrolre=compile_regex(controlfrom) 
        secondmutre=compile_regex(any_base, mutupstream, mutdownstream) 
        secondcontrolre=compile_regex(any_base, controlupstream, controldownstream)
    else: # both or ancestor
        potentialmutre=compile_regex(mutfrom, mutupstream, mutdownstream) 
        potentialcontrolre=compile_regex(controlfrom, controlupstream, controldownstream)
        if enforce=="A": # ancestor
            secondmutre=compile_regex(any_base) 
            secondcontrolre=compile_regex(any_base) 
        elif enforce=="B": # both
            secondmutre=compile_regex(any_base, mutupstream, mutdownstream) 
            secondcontrolre=compile_regex(any_base, controlupstream, controldownstream) 
        else:
            raise ValueError("unknown enforce value")

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
        check_chars(sequence, iupac_codes+['-'], 'Query sequences must contain only IUPAC characters or - (for gap).')

        seqs+=1

        mutsites, muts = summarize_matches(refseq, sequence, start, finish, 
                                           potentialmutre, secondmutre, 
                                           mutto_orig, mutupstream_orig, mutdownstream_orig, 
                                           iupac_dict, multistate, seqs, name)
        
        controlsites, controls = summarize_matches(refseq, sequence, start, finish, 
                                                   potentialcontrolre, secondcontrolre, 
                                                   mutto_orig, controlupstream_orig, controldownstream_orig, 
                                                   iupac_dict, multistate, seqs, name)

        sys.stdout=origstdout
        odds_ratio, pval = calc_fisher(mutsites, muts, controlsites, controls)

        try:
            ratio = "%0.2f" %(muts*controlsites/(1.0*mutsites*controls))
        except:
            if muts*controlsites > 0:
                ratio="inf" 
            else:
                ratio="undef" 
    
        if sumfile:
            s = name+","+str(muts)+","+str(mutsites)+","+str(controls) +","+str(controlsites) +","+ratio +",%.6g" %float(pval)+'\n' 
            sumfile.write(s)

        sys.stdout.flush()
        sys.stdout=outfile

    sys.stdout=origstdout

    if sumfile:
        sumfile.close()
    if outfile:
        outfile.close()
    origstdout.close()
