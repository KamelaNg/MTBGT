#!/usr/bin/env python
import sys
import re
import argparse
import TBGT_IO as tbgtio
import TBGT_Tests as tbgtt
"""
Author: Kamela Ng, Conor Meehan <Feb 2019>

This script converts M. tuberculosis whole genome sequence data into the most likely output for different rifampicin resistance detection molecular tests

Input:
A folder containing all the vcf or tab files that are to be processed
A map file that states the relationship between codons where mutations are known to be detected by the RR-TB RDTs (based on these studies:Ng KC, Meehan CJ, Torrea G, Goeminne L, Diels M, Rigouts L, de Jong BC, and Andre E. 2018a. Potential Application of Digitally Linked Tuberculosis Diagnostics for Real-Time Surveillance of Drug-Resistant Tuberculosis Transmission: Validation and Analysis of Test Results. JMIR Med Inform 6:e12. 10.2196/medinform.9309 and Ng KCS, van Deun A, Meehan CJ, Torrea G, Driesen M, Gabriels S, Rigouts L, Andre E, and de Jong BC. 2018c. Xpert Ultra Can Unambiguously Identify Specific Rifampin Resistance-Conferring Mutations. J Clin Microbiol 56. 10.1128/JCM.00686-18) and positions in the genome
A filename prefix for the outputs
The tests to run
Whether input is tab or vcf
If overview tables are to be made

The map file is optional. If not supplied, it is assumed the standard H37Rv NC000962.3 was used and that the file default_mapfile.txt supplied with the script is used
An example of this file is (tab separated columns):
Codon	Ref	Gen	Pos
170	760314	
428	761088
430	761094
431	761097
432	761100
434	761106
435	761109
437	761115
441	761127
445	761139
446	761142
450	761154
452	761160
491	761277

NOTE:
All tests except Sanger use the above codons except 170 and 491, which are not detceted by commercial RDTs. 
The Sanger module uses the above default codon list and displays any change in amino acid at these positions. To add more codons, they must be added to the dictionary at the top of the geno_to_Probe_Tests.py module, including the wildtype nucloetides for the 3 codon positions.


The filename prefix for the output is optional. If not supplied, output will be genomeToProbe (e.g. genomeToProbe_Xpert.txt etc. )

The tests to run are supplied as letter inputs. To run all the tests the --tests flag is not needed or can supply --tests all. For the individual tests the following codes are used (lower or upper case)
X	Xpert
U	Ultra
H	Hain
N	Nipro
S	Sanger

Thus, using --tests xu will run classic Xpert and Ultra only. The letters can be supplied in any order (i.e. --tests xhn and --tests HnX are the same)

For vcf files use --type vcf and for MTBseq tab output use --type tab (upper or lower case are both acceptable)

If an overview table of the results from each requested test is to be made the --summary y options should be used. This is on by default


Output:
The output files contain, for each input tab file:
Sample name	RIF Resistance	Codon number	Codon	xx
Where xx is the results of the probe presence/absence

If the summary option is selected, a tab delimited file with the frequency and total count of each probe (Xpert, Ultra, Haina nd Nipro) or mutation (Sanger) along with the sensitive is created.
Tests are in the order Xpert, Ultra, Hain, Nipro, Sanger (or a subset as requested)

Usage:
python TBGT.py --folder <folderName> --type tab/vcf --map <codonMapFile> (optional) --out <outputFilename> (optional) --tables y/n (optional)
"""

#parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument('--folder', required=True, help='Folder containing the vcf/tab files to be processed')
parser.add_argument('--map', required=False, default="default_mapfile.txt", help='File listing the relationship between codon and genome position (default is "default_mapfile.txt" supplied with script)')
parser.add_argument('--out', required=False, default="genomeToProbe", help='Filename prefix for output (default is "genomeToProbe")')
parser.add_argument('--tests', required=False, default="all", help='The different tests to run (see helpfile)')
parser.add_argument('--type', required=True, help='Whether input is vcf or tab')
parser.add_argument('--summary', required=False, default="Y", help='Option to create an tab delimited overview tables file')


args = parser.parse_args()

#read in map file 
try:
	mapf=open(args.map, 'r')
except IOError:
	print("\n map file not found.")
	sys.exit()

#read in map and create a dictionary of the codons to genome positions
#each key is a position and the value is the codon_position (e.g. position 1 in the 450 codon is 450_1)
mapdict={}
genomePos=[]
while 1:
	line=mapf.readline()
	if not line:
		break
	line=line.rstrip()
	if re.match("Codon",line):#on an info line so skip
		continue
	sections=line.split("\t")
	pos=sections[1]
	mapdict[pos]=sections[0]+"_0"
	mapdict[str(int(pos)+1)]=sections[0]+"_1"
	mapdict[str(int(pos)+2)]=sections[0]+"_2"
	genomePos.append(pos)
	genomePos.append(str(int(pos)+1))
	genomePos.append(str(int(pos)+2))
	codonpos=mapdict.keys() #get a list of all the genome positions, grouped as triplets
mapf.close()



#read in the tab or vcf files, keeping the SNP that is present at each of the positions in the genomePos
#create a dictionary where the filename (before the underscore of the filename) is the key and a dictionary of all genome positions (keys) and their SNPs (if any) (values) is returned
if args.type.upper()=='TAB':
	sampleGenomePos= tbgtio.tabInput(genomePos,args.folder)
elif args.type.upper()=='VCF':
	sampleGenomePos= tbgtio.vcfInput(genomePos,args.folder)
else:
	print("\n type input not recognised. Please use tab or vcf")
	sys.exit()

#run the modules as requested by the user

tests={"H":0, "U":0, "X":0, "N":0, "S":0,}

if args.tests.upper()=='ALL':
	tests={"H":1, "U":1, "X":1, "N":1, "S":1,}
else:
	testLetters=list(args.tests.upper())
	for letter in testLetters:
		tests[letter]=1

if tests["X"]==1:	
	print("Creating classic Xpert outputs")
	tbgtt.xpert(mapdict, sampleGenomePos, args.out)
if tests["U"]==1:	
	print("Creating Xpert Ultra outputs")
	tbgtt.ultra(mapdict, sampleGenomePos, args.out)
if tests["H"]==1:	
	print("Creating Hain LPA outputs")
	tbgtt.hain(mapdict, sampleGenomePos, args.out)
if tests["N"]==1:	
	print("Creating Nipro LPA outputs")
	tbgtt.nipro(mapdict, sampleGenomePos, args.out)
if tests["S"]==1:	
	print("Creating rpoB Sanger sequencing outputs")
	tbgtt.sanger(mapdict, sampleGenomePos, args.out)

print("RDT test conversion complete")


#Create tables as requested by the user
if args.summary.upper()=="Y":
	print("Creating summary tables")
	tbgtio.tables(tests, args.out)

sys.exit()
