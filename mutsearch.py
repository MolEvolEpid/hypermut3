#!/usr/bin/env python
###!/sw/bin/python2.4

# this version is just supposed to find matches and print them to a file
# and print a summary to stdout
# Note main argument now has 8 fields.

# Print summary to stdout, locations to outfile
   
import sys, getopt
import re
import string
from subprocess import *
from scipy.stats import fisher_exact
import itertools
import warnings
usage="usage: mutsearch.py [-s start][-f finish][-h][-e (A|B|D)][-m multiplier] [-o outfile] [-u summaryfile] 'mutfrom,mutto,controlfrom,controlto,mutupstream,mutdownstream,controlupstream,controldownstream' < inputseqs.tbl"

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

def check_overlap(mut_pattern, control_pattern, width, loc):

    base_info_mut = [[list(y) for y in x.split()] for x in mut_pattern.split('|')]
    base_info_control = [[list(y) for y in x.split()] for x in control_pattern.split('|')]

    contexts_mut = [''.join(list(y)) for x in base_info_mut for y in itertools.product(*x)]
    contexts_control = [''.join(list(y)) for x in base_info_control for y in itertools.product(*x)]

    overlap = set(contexts_mut) & set(contexts_control)
    n_patterns = len(set(contexts_mut)) + len(set(contexts_control))
    n_patterns_expected = 4**width
    if len(overlap) != n_patterns_expected:
        if len(overlap):
            raise ValueError(f"Partially overlapping {loc} pattern contexts: {overlap}. This means that positions cannot be categorized uniquely into primary and control groups.")
        if n_patterns != n_patterns_expected:
            warnings.warn(f"{loc} primary and control patterns (n= {str(n_patterns)}) do not create the full complement of possible patterns (n={str(n_patterns_expected)})", stacklevel=2)

def check_width(context):
    # standard IUPAC codes translated to regexps
    context=re.sub('R',"[AGR]",context)
    context=re.sub('Y',"[CTY]",context)
    context=re.sub('M',"[ACM]",context)
    context=re.sub('K',"[GTK]",context)
    context=re.sub('S',"[CGS]",context)
    context=re.sub('W',"[ATW]",context)
    context=re.sub('B',"[CGTYBSK]",context)
    context=re.sub('D',"[AGTRDWK]",context)
    context=re.sub('H',"[ACTYHWM]",context)
    context=re.sub('V',"[ACGRVSM]",context)
    context=re.sub('N',"[ACGTRYBDHVNWSKM]",context)
    context=re.sub('\\.', "", context)

    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(context, ",")

    if (not isfixedwidth(mutupstream)) or (not isfixedwidth(controlupstream)):
        die_widthnotfixed("upstream")
    if (not isfixedwidth(mutfrom)) or (not isfixedwidth(mutto)):
        die_widthnotfixed("mutation")
    if (not isfixedwidth(controlfrom)) or (not isfixedwidth(controlto)):
        die_widthnotfixed("control mutation")
    if len(widthonly(mutupstream)) != len(widthonly(controlupstream)):
        raise ValueError("Upstream primary and control patterns must be the same length")
    if len(widthonly(mutdownstream)) != len(widthonly(controldownstream)):
        raise ValueError("Downstream primary and control patterns must be the same length")
    
    return(len(widthonly(mutupstream))+len(widthonly(mutdownstream))) 

def check_context(context):
    width = check_width(context)

    context=re.sub('A'," A ",context)
    context=re.sub('C'," C ",context)
    context=re.sub('G'," G ",context)
    context=re.sub('T'," T ",context)
    context=re.sub('N'," ACGT ",context) # standard IUPAC codes translated to regexps
    context=re.sub('R'," AG ",context)
    context=re.sub('Y'," CT ",context)
    context=re.sub('B'," CGT ",context)
    context=re.sub('D'," AGT ",context)
    context=re.sub('H'," ACT ",context)
    context=re.sub('V'," ACG ",context)
    context=re.sub('M'," AC ",context)
    context=re.sub('K'," GT ",context)
    context=re.sub('S'," CG ",context)
    context=re.sub('W'," AT ",context) 
    context=re.sub('\\.', "", context)
    
    (mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(context, ",")    

    mut_pattern = '|'.join([y + x for x in mutdownstream.split('|') for y in mutupstream.split('|')])
    control_pattern = '|'.join([y + x for x in controldownstream.split('|') for y in controlupstream.split('|')])

    base_info_mut = [[list(y) for y in x.split()] for x in mut_pattern.split('|')]
    base_info_control = [[list(y) for y in x.split()] for x in control_pattern.split('|')]

    contexts_mut = [''.join(list(y)) for x in base_info_mut for y in itertools.product(*x)]
    contexts_control = [''.join(list(y)) for x in base_info_control for y in itertools.product(*x)]

    overlap = set(contexts_mut) & set(contexts_control)
    n_patterns = len(set(contexts_mut)) + len(set(contexts_control))
    n_patterns_expected = 4**(width)
    if len(overlap) != n_patterns_expected:
        if len(overlap):
            raise ValueError(f"Partially overlapping upstream-downstream pattern contexts: {overlap}. This means that positions cannot be categorized uniquely into primary and control groups.")
        if n_patterns != n_patterns_expected:
            warnings.warn(f"{loc} primary and control patterns (n= {str(n_patterns)}) do not create the full complement of possible patterns (n={str(n_patterns_expected)})", stacklevel=2)
 
# main starts here #

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:f:e:m:o:u:", ["help", "start=","finish=","enforce=","multiplier=","outfile=", "summaryfile="])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))
    print(usage)
    sys.exit(2)

start=0
finish=None
outfile=None
sumfile=None
enforce="D"  # first letter of Ancestor, Descendant, or Both
multiplier = 1.0

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
    if o in ("-m", "--multiplier"):
        multiplier=a


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
arg_context=arg
iupac = list('ACGTRYBDHVNWSKM')

bad_chars = [x for x in set(arg) if x not in iupac+[',','.','|']]
if len(bad_chars):
    raise ValueError(f"All patterns and mutations must include only IUPAC characters, '.', and '|'. You included non-IUPAC characters: {bad_chars}")

check_context(arg_context)

# allow for gaps
arg=re.sub('A',"[A][-]*",arg)
arg=re.sub('C',"[C][-]*",arg)
arg=re.sub('G',"[G][-]*",arg)
arg=re.sub('T',"[T][-]*",arg)
# standard IUPAC codes translated to regexps
arg=re.sub('R',"[AGR][-]*",arg)
arg=re.sub('Y',"[CTY][-]*",arg)
arg=re.sub('M',"[ACM][-]*",arg)
arg=re.sub('K',"[GTK][-]*",arg)
arg=re.sub('S',"[CGS][-]*",arg)
arg=re.sub('W',"[ATW][-]*",arg)
arg=re.sub('B',"[CGTYBSK][-]*",arg)
arg=re.sub('D',"[AGTRDWK][-]*",arg)
arg=re.sub('H',"[ACTYHWM][-]*",arg)
arg=re.sub('V',"[ACGRVSM][-]*",arg)
arg=re.sub('N',"[ACGTRYBDHVNWSKM][-]*",arg)
arg=re.sub('\\.', "", arg)
if(sys.stdout):
   print("#regexps=",arg)
   print("seq_num,seq_name,control,potential_mut_site,mut_match")

(mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(arg, ",")

mutfrom = re.sub('\\[|\\]|-|\\*', '', mutfrom)
mutto = re.sub('\\[|\\]|-|\\*', '', mutto)
controlfrom = re.sub('\\[|\\]|-|\\*', '', controlfrom)
controlto = re.sub('\\[|\\]|-|\\*', '', controlto)

# require mutation to be only one base
if mutfrom not in iupac or mutto not in iupac or controlfrom not in iupac or controlto not in iupac:
    raise ValueError("Mutations must be single IUPAC characters.") 

# use ?= "lookahead" and ?<= lookbehind so we can find overlapping patterns
# Note that the regexps will be applied to different
# strings depending on the value of enforce.

# potential always checked on ancestor
# second confirms potential in descendant
# actual means mutation fits pattern

if enforce=="D": # descendant
    potentialmutre=re.compile(mutfrom,re.I)
    secondmutre=re.compile("(?<="+ mutupstream +")[-]*("+widthonly(mutto) +")[-]*(?="+ mutdownstream+ ")",re.I)
    potentialcontrolre=re.compile(controlfrom,re.I)
    secondcontrolre=re.compile("(?<="+ controlupstream +")[-]*("+widthonly(controlto) +")[-]*(?="+ controldownstream+ ")",re.I)
    actualmutre=re.compile(mutto,re.I)
    actualcontrolre=re.compile(controlto,re.I)
else:       # both or ancestor
    potentialmutre=re.compile("(?<="+ mutupstream +")[-]*("+mutfrom +")[-]*(?="+ mutdownstream+ ")",re.I)
    potentialcontrolre=re.compile("(?<="+ controlupstream +")[-]*("+controlfrom +")[-]*(?="+ controldownstream+ ")",re.I)
    if enforce=="A": # ancestor
        actualmutre=re.compile(mutto,re.I)
        secondmutre=re.compile(widthonly(mutto),re.I)  # nullstring might work too?
        actualcontrolre=re.compile(controlto,re.I)
        secondcontrolre=re.compile(widthonly(controlto),re.I)
    elif enforce=="B": # both
        actualmutre=re.compile("(?<="+ mutupstream +")[-]*("+mutto +")[-]*(?="+ mutdownstream+ ")",re.I)
        secondmutre=re.compile("(?<="+ mutupstream +")[-]*("+widthonly(mutto) +")[-]*(?="+ mutdownstream+ ")",re.I)
        actualcontrolre=re.compile("(?<="+ controlupstream +")[-]*("+controlto +")[-]*(?="+ controldownstream+ ")",re.I)
        secondcontrolre=re.compile("(?<="+ controlupstream +")[-]*("+widthonly(controlto) +")[-]*(?="+ controldownstream+ ")",re.I)
    else:
        raise ValueError("unknown enforce value")

#print('potentialmutre', potentialmutre)
#print('potentialcontrolre', potentialcontrolre)
#print('actualmutre', actualmutre)
#print('secondmutre', secondmutre)
#print('actualcontrolre', actualcontrolre)
#print('secondcontrolre', secondcontrolre)

seqs=0

line=sys.stdin.readline()
refname=line[1:].strip()
refseq=""
line=sys.stdin.readline()
while line and line[0]!=">":
    refseq=refseq+line
    line=sys.stdin.readline()
refseq=refseq.replace("\n","").upper()

bad_chars = [x for x in set(refseq) if x not in ['A', 'C', 'G', 'T', '-']]
if len(bad_chars):
    raise ValueError(f'The reference sequence must contain only the following characters: ACGT-. Yours contains: {bad_chars}')

seqs=0
while line:
    (mutsites,muts)=(0,0)  # mutsites= potentialmuts also passing secondre
    (controlsites,controls)=(0,0)
    name=line[1:].strip()
    sequence=""
    line=sys.stdin.readline()
    while line and line[0]!=">":
        sequence=sequence+line
        line=sys.stdin.readline()
    sequence=sequence.replace("\n","").upper()

    bad_chars = [x for x in set(sequence) if x not in iupac+['-']]
    if len(bad_chars):
        raise ValueError(f'The reference sequence must contain only IUPAC characters or - (for gap). Yours contains: {bad_chars}')

    seqs+=1

    if(finish):
        potentialmuts=potentialmutre.finditer(refseq,start,finish)
    else:
        potentialmuts=potentialmutre.finditer(refseq,start)

    for mymatch in potentialmuts:
        if secondmutre.match(sequence,mymatch.start()): #optional match arg
            base = sequence[mymatch.start():mymatch.end()]
            if base == '-':
                continue
            yval=0  # 0 is potential match, 1 will be actual match
            mutsites+=1
            if actualmutre.match(sequence,mymatch.start()):
                yval = 1
                muts += 1
            if(sys.stdout):
                print(str(seqs) + "," + name + ",0," + str(mymatch.start()+1) + "," + str(yval))

    if(finish):
        potentialcontrols=potentialcontrolre.finditer(refseq,start,finish)
    else:
        potentialcontrols=potentialcontrolre.finditer(refseq,start)

    for mymatch in potentialcontrols:
        if secondcontrolre.match(sequence,mymatch.start()): #optional match arg
            base = sequence[mymatch.start():mymatch.end()]
            if base == '-':
                continue
            yval=0  # 0 is potential match, 1 will be actual match
            controlsites+=1
            if actualcontrolre.match(sequence,mymatch.start()):
                yval = 1
                controls+=1
            if(sys.stdout):
               print(str(seqs) + "," + name + ",1," + str(mymatch.start()+1) + "," + str(yval))

    sys.stdout=origstdout
    odds_ratio, pval = fisher_exact([[muts, mutsites-muts],[controls,controlsites-controls]], alternative = 'greater')

    try:
        ratio = "%0.2f" %(muts*controlsites/(1.0*mutsites*controls))
    except:
        if muts*controlsites > 0:
            ratio ="inf" 
        else:
            ratio ="undef" 
   
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
