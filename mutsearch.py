#!/usr/bin/env python
###!/sw/bin/python2.4

# this version is just supposed to find matches and print them to a file
# and print a summary to stdout
# Note main argument now has 8 fields.

# now supports both fasta and tabular!



# Print summary to stdout, locations to outfile
   
import sys, getopt
import re
import string
from subprocess import *
usage="usage: mutsearch.py [-s start][-f finish][-h][-v][-e (A|B|D)][-m multiplier] [-o outfile] [-u summaryfile] 'mutfrom,mutto,controlfrom,controlto,mutupstream,mutdownstream,controlupstream,controldownstream' < inputseqs.tbl"

fishertest="perl -wT fisher_exact_test.pl"

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
    if verbose:
        print("</table>")
    print("error: ", culprit," pattern must correspond to a fixed length expression")
    print("(example: T|GC is not allowed).")
    if verbose:
        print("<br><b>")
    print("Please try again using a ", culprit, "pattern that has only one possible width.")
    if verbose:
        print("</body></html><font size=1 color=\"silver\">")
    sys.exit(1)


def widthonly(regexpstring):
    if not isfixedwidth(regexpstring):
        raise error("not fixed width")
    dotstring=""
    while (not isfixedwidth(regexpstring+"|"+dotstring)) and len(dotstring)<500:
        dotstring = dotstring + "."
    if len(dotstring)==500:
        raise error("regexp too long (greater than 500)")
    return dotstring
        

    

# main starts here #

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:f:ve:m:o:u:", ["help", "start=","finish=","verbose","enforce=","multiplier=","outfile=", "summaryfile="])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))
    print(usage)
    sys.exit(2)

verbose = False
start=0

#finish=sys.maxint
# actually, there's a 'bug' in python that breaks if
# we do that on 64 bit machines.  The following is less
# correct but fine for HIV.
# finish=200000000
# actually we can fix it when finish is called
finish=None
#outfile=sys.stdout don't want this if no file specified
outfile=None
sumfile=None

enforce="D"  # first letter of Ancestor, Descendant, or Both
multiplier = 1.0

for o, a in opts:
    if o in ("-v", "--verbose"):
        verbose = True  # this means use html output
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
        enforce = a[0].upper

    if o in ("-m", "--multiplier"):
        multiplier=a


if not (len(args)==1 ):
    print(usage)
    print('aa')
    sys.exit()

if sumfile:
    sumfile.write("Sequence,Muts(Match Sites),Out of(Potential Mut Sites),Controls(Control Muts),Out of(Potential Controls),Rate Ratio(A/B)/(C/D),Fisher Exact P-value\n")

if verbose:
    print("<table>")
    print("<tr><th> <th><u>Sequence:</u>&nbsp;<br><font size=-2>(Select for&nbsp;<br>graphing)</font><th><u>Muts:</u> <br><font size=-2>(Match<br> Sites)</font> <th> <u>Out of:</u> <br><font size=-2>(Potential <br> Mut Sites)</font>")
    print("            <th> &nbsp;<th><u>Controls:</u><br><font size=-2>(Control <br> Muts)</font> <th> <u>Out of:</u> <br><font size=\"-2\">(Potential <br> Controls)</font>")
    print(" <th><u>Rate&nbsp;Ratio:</u><br><font size=\"-2\">(A/B)/(C/D)<br>&nbsp;                      <th colspan=2> <u>Fisher Exact P-value:</u><br>&nbsp;<font size=-2>(=P(Muts,Poten.Muts-Muts,<br>&nbsp;&nbsp;Cntrls,Poten.Cntrls-Cntrls))</font></tr>")
    sys.stdout.flush()

origstdout=sys.stdout
sys.stdout=outfile

arg=args[0]
if(sys.stdout):
   print("#arg=", arg)

arg=arg.upper()
arg=re.sub('N',".",arg) # standard IUPAC codes translated to regexps
arg=re.sub('R',"[AG]",arg)
arg=re.sub('Y',"[CT]",arg)
arg=re.sub('B',"[^A]",arg)
arg=re.sub('D',"[^C]",arg)
arg=re.sub('H',"[^G]",arg)
arg=re.sub('V',"[^T]",arg)
arg=re.sub('M',"[AC]",arg)
arg=re.sub('K',"[GT]",arg)
arg=re.sub('S',"[CG]",arg)
arg=re.sub('W',"[AT]",arg)
if(sys.stdout):
   print("#regexps=",arg)

(mutfrom, mutto, controlfrom, controlto, mutupstream, mutdownstream, controlupstream, controldownstream)=str.split(arg, ",")


# use ?= "lookahead" and ?<= lookbehind so we can find overlapping patterns
# Note that the regexps will be applied to different
# strings depending on the value of enforce.

# potential always checked on ancestor
# second confirms potential in descendant
# actual means mutation fits pattern

if (not isfixedwidth(mutupstream)) or (not isfixedwidth(controlupstream)):
    die_widthnotfixed("upstream")
if (not isfixedwidth(mutfrom)) or (not isfixedwidth(mutto)):
    die_widthnotfixed("mutation")
if (not isfixedwidth(controlfrom)) or (not isfixedwidth(controlto)):
    die_widthnotfixed("control mutation")

    
if enforce=="D": # descendant
    potentialmutre=re.compile(mutfrom,re.I)
    secondmutre=re.compile("(?<="+ mutupstream +")("+widthonly(mutto) +")(?="+ mutdownstream+ ")",re.I)
    potentialcontrolre=re.compile(controlfrom,re.I)
    secondcontrolre=re.compile("(?<="+ controlupstream +")("+widthonly(controlto) +")(?="+ controldownstream+ ")",re.I)
    actualmutre=re.compile(mutto,re.I)
    actualcontrolre=re.compile(controlto,re.I)
else:       # both or ancestor
    potentialmutre=re.compile("(?<="+ mutupstream +")("+mutfrom +")(?="+ mutdownstream+ ")",re.I)
    potentialcontrolre=re.compile("(?<="+ controlupstream +")("+controlfrom +")(?="+ controldownstream+ ")",re.I)
    if enforce=="A": # ancestor
        actualmutre=re.compile(mutto,re.I)
        secondmutre=re.compile(widthonly(mutto),re.I)  # nullstring might work too?
        actualcontrolre=re.compile(controlto,re.I)
        secondcontrolre=re.compile(widthonly(controlto),re.I)
    elif enforce=="B": # both
        actualmutre=re.compile("(?<="+ mutupstream +")("+mutto +")(?="+ mutdownstream+ ")",re.I)
        secondmutre=re.compile("(?<="+ mutupstream +")("+widthonly(mutto) +")(?="+ mutdownstream+ ")",re.I)
        actualcontrolre=re.compile("(?<="+ controlupstream +")("+controlto +")(?="+ controldownstream+ ")",re.I)
        secondcontrolre=re.compile("(?<="+ controlupstream +")("+widthonly(controlto) +")(?="+ controldownstream+ ")",re.I)
    else:
        raise error("unknown enforce value")

seqs=0
seqlist=[]
fasta=False

line=sys.stdin.readline()
try:
    (refname, refseq) = line.split("\t")
except:
    if line[0]!=">":
        (refname, refseq) = line.split()
    else:
        fasta=True
        refname=line[1:]
        refseq=""
        line=sys.stdin.readline()
        while line and line[0]!=">":
            refseq=refseq+line
            line=sys.stdin.readline()
        refseq=refseq.replace("\n","")

if not fasta:
    line=sys.stdin.readline()


seqs=0
while line:
    (mutsites,muts)=(0,0)  # mutsites= potentialmuts also passing secondre
    (controlsites,controls)=(0,0)
    if not fasta:
        (name, sequence) = line.split()
    else:
        name=line[1:]
        sequence=""
        line=sys.stdin.readline()
        while line and line[0]!=">":
            sequence=sequence+line
            line=sys.stdin.readline()
        sequence=sequence.replace("\n","")

        
    seqs+=1
    if(sys.stdout):
        print("#seqno= "+ str(seqs))
        print("#name= "+ name)

    if(finish):
        potentialmuts=potentialmutre.finditer(refseq,start,finish)
    else:
        potentialmuts=potentialmutre.finditer(refseq,start)

    for mymatch in potentialmuts:
        if secondmutre.match(sequence,mymatch.start()): #optional match arg
            yval=0  # 0 is potential match, 1 will be actual match
            mutsites+=1
            if actualmutre.match(sequence,mymatch.start()):
                yval = 1
                muts += 1
            if(sys.stdout):
                print(mymatch.start()+1, yval)
    if(sys.stdout):
      print("")
      print("# "+ name +" controls")

    if(finish):
        potentialcontrols=potentialcontrolre.finditer(refseq,start,finish)
    else:
        potentialcontrols=potentialcontrolre.finditer(refseq,start)

    for mymatch in potentialcontrols:
        if secondcontrolre.match(sequence,mymatch.start()): #optional match arg
            yval=0  # 0 is potential match, 1 will be actual match
            controlsites+=1
            if actualcontrolre.match(sequence,mymatch.start()):
                yval = 1
                controls+=1
            if(sys.stdout):
               print(mymatch.start()+1, yval)


    if(sys.stdout):            
      print("")
    sys.stdout=origstdout
    mycmd=fishertest+" "+str(muts)+" "+str(mutsites-muts)+" "+str(controls)+" "+str(controlsites-controls)
    mypipe=Popen(mycmd,shell=True,stdout=PIPE)
    pval=mypipe.stdout.read()

    try:
        ratio = "%0.2f" %(muts*controlsites/(1.0*mutsites*controls))
    except:
        if muts*controlsites > 0:
            ratio ="inf" 
        else:
            ratio ="undef" 
       
    if verbose:
        print("<tr><td>  <INPUT type=\"checkbox\" name=\"Seq"+str(seqs)+"\" value=\""+name+"\"checked><td>"+name)
        print("<td align=right>",muts,"<td align=right>",mutsites,"<td><td align=right>",controls,"<td align=right>",controlsites, "<td align=center>", ratio) 
        #try:
        #    print "<td align=center> %0.2f" %( muts*controlsites/(1.0*mutsites*controls))
        #except:
        #    if muts*controlsites > 0:
        #        print "<td align=center> inf"
        #    else:
        #        print "<td align=center> undef"
        print("<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<td> ",pval,"</tr>")
    else:
        print(muts,mutsites,controls,controlsites,ratio ,pval)
   
    if sumfile:
        s = name.strip()+","+str(muts)+","+str(mutsites)+","+str(controls) +","+str(controlsites) +","+ratio +","+str(float(pval))+'\n' 
        print(s)
        sumfile.write(s)

    sys.stdout.flush()
    sys.stdout=outfile

    # get next line
    if not fasta:
        line=sys.stdin.readline()

sys.stdout=origstdout
if verbose:
    print("</table><br>")
    print('<INPUT TYPE="HIDDEN" NAME="Sequences" VALUE="',seqs,'">')

if sumfile:
   sumfile.close()
if outfile:
   outfile.close()
origstdout.close()
