import glob 
import re
import sys
from collections import defaultdict
'''
Function to read in a folder of tab delimited variant files output from MTBseq.
Input is the list of genome positions to extract and the folder name containing the tab files
First line begins with an # so skip this.
Genome position is in the 1st column, variant type is 4th column, and the SNP is in the 5th column so extract these
A dictionary where the filename (before the underscore of the filename) is the key and a dictionary of all genome positions (keys) and their SNPs (if any) (values) is returned
'''
def tabInput(genomePos,folder):
	filesGenomePos={}
	OpenDir = glob.glob(folder + "/*")	
	print('Reading in tab files')									  
	for File in OpenDir: 
		if "variants" in File and File.endswith(".tab"): 
			try:
				tabF=open(File,'r')
			except IOError:
				print(File+" file not found.")
				sys.exit()
			sampleName=re.sub("_.*","",File)
			sampleName=re.sub(".*/","",sampleName)
			print(sampleName)
			genomeSNPmap=dict.fromkeys(genomePos)
			while 1:
				line=tabF.readline()
				if not line:
					break
				line=line.rstrip()
				if re.match("#",line):	#on an info line so skip
					continue
				sections=line.split("\t")
				mutpos=sections[0]
				vartype=sections[3]
				altbase=sections[4]
				if vartype.upper() == "SNP" and mutpos in genomePos:
					genomeSNPmap[mutpos]= altbase
			filesGenomePos[sampleName]=genomeSNPmap
			tabF.close()
	return(filesGenomePos)
	
'''
Function to read in a folder of vcfs.
Input is the list of genome positions to extract and the folder name containing the vcfs
Comment lines begins with an # so skip these.
Genome position is in the 2nd column and the SNP is in the 4th column so extract these
A dictionary where the filename (before the underscore of the filename) is the key and a dictionary of all genome positions (keys) and their SNPs (if any) (values) is returned
'''
def vcfInput(genomePos,folder):
	filesGenomePos={}
	OpenDir = glob.glob(folder + "/*")	
	print('Reading in vcf files')									  
	for File in OpenDir:
		if File.endswith(".vcf"): 			
			try:
				vcfF=open(File,'r')
			except IOError:
				print(File+" file not found.")
				sys.exit()
			sampleName=re.sub("_.*","",File)
			sampleName=re.sub("\.vcf","",sampleName)
			sampleName=re.sub(".*/","",sampleName)
			genomeSNPmap=dict.fromkeys(genomePos)
			while 1:
				line=vcfF.readline()
				if not line:
					break
				line=line.rstrip()
				if re.match("#",line):	#on an info line so skip
					continue
				sections=line.split("\t")
				mutpos=sections[1]
				altbase=sections[4]
				info=sections[7]
				if "TYPE=SNP" in info.upper() and mutpos in genomePos:
					genomeSNPmap[mutpos]= altbase
			filesGenomePos[sampleName]=genomeSNPmap
			vcfF.close()
	return(filesGenomePos)	

'''
Function to write count and proportion tables to file.
Input is the test name, the save handle, the count dictionary and the proprotion dictionary
Output is the header per table and the table tab separated
'''
def writeTables(tsF,testName,cd,pd):
	tsF.write(testName+" counts\n")
	countkeys=list(cd.keys())
	countkeys.sort()
	for ck in countkeys:
		tsF.write(ck+"\t"+str(cd[ck])+"\n")
	tsF.write("\n"+testName+" proportions\n")
	propkeys=list(pd.keys())
	propkeys.sort()
	for pk in propkeys:
		tsF.write(pk+"\t"+str(round(pd[pk],2))+"\n")
	tsF.write("\n")

'''
Function to create table summaries of Xpert classic and Ultra results
'''
def xpertTables(F,tsF,testName):
	#go into the Xpert output file and grab all the prob listings and count them up. Combination of probes are listed separate from probes by themselves. 'NOT DETECTED' is counted as sensitive
	countdict=defaultdict(int)
	total=0
	
	while 1:
		line=F.readline()
		if not line:
			break
		line=line.rstrip()	
		sections=line.split("\t")
		if sections[0]=="Sample name":#header line
			continue	
		if sections[1]=="NOT DETECTED":
			countdict["Rif Sensitive"]+=1
			total+=1
		else:
			countdict[sections[4]]+=1	
			total+=1
	#create the proportion count of each entry
	propdict={}	
	for entry in countdict:
		propdict[entry]=float(countdict[entry])/total
	#save the two tables to file
	writeTables(tsF,testName,countdict,propdict)	

'''
Function to create table summaries of LPA tests (Hain and Nipro)
The wt and mut probes are combined for the summaries
'''
def lpaTables(F,tsF,testName):
	#go into the LPA output file and grab all the prob listings and count them up. Combination of probes are listed separate from probes by themselves. 'NOT DETECTED' is counted as sensitive
	countdict=defaultdict(int)
	total=0
	
	while 1:
		line=F.readline()
		if not line:
			break
		line=line.rstrip()	
		sections=line.split("\t")
		if sections[0]=="Sample name":#header line
			continue	
		if sections[1]=="NOT DETECTED":
			countdict["Rif Sensitive"]+=1
			total+=1
		else:
			if sections[5]!="":
				probes=sections[4]+","+sections[5]
			else:
				probes=sections[4]	
			countdict[probes]+=1	
			total+=1
	#create the proportion count of each entry
	propdict={}	
	for entry in countdict:
		propdict[entry]=float(countdict[entry])/total
	#save the two tables to file
	writeTables(tsF,testName,countdict,propdict)	
		
	
'''
Function to create an excel file of overview tables per test
Input is the list of tests performed and the prefix used so that the relevant output files can be read in

'''
def tables(tests,outPrefix):	
	#Create output files
	try:
		tablesSave = open(outPrefix+"_summaryTables.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	
	if tests["X"]==1:
		try:
			F=open(outPrefix+"_XpertClassic.txt",'r')
		except IOError:
			print(File+" file not found.")
			sys.exit()	
		xpertTables(F,tablesSave,"Xpert classic")
		F.close()
		
	if tests["U"]==1:	
		try:
			F=open(outPrefix+"_XpertUltra.txt",'r')
		except IOError:
			print(File+" file not found.")
			sys.exit()	
		xpertTables(F,tablesSave,"Xpert Ultra")
		F.close()
		
	if tests["H"]==1:	
		try:
			F=open(outPrefix+"_HainLPA.txt",'r')
		except IOError:
			print(File+" file not found.")
			sys.exit()
		lpaTables(F,tablesSave,"Hain LPA")
		F.close()	
		
	if tests["N"]==1:	
		try:
			F=open(outPrefix+"_NiproLPA.txt",'r')
		except IOError:
			print(File+" file not found.")
			sys.exit()
		lpaTables(F,tablesSave,"Nipro LPA")
		F.close()	
		
	if tests["S"]==1:#sanger input is same	format as the Xpert so can use that function
		try:
			F=open(outPrefix+"_Sanger.txt",'r')
		except IOError:
			print(File+" file not found.")
			sys.exit()	
		xpertTables(F,tablesSave,"Sanger")
		F.close()
