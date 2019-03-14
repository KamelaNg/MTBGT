import sys
import re
import copy

#a dictionary of the wildtype codons where the number (MTB numbering, not E. coli) is the key and the codon as a list is the value
codon_nuc = {
			'170' : ['G', 'T', 'C'],
			'428' : ['A', 'G', 'C'], 
			'430' : ['C', 'T', 'G'], 
			'431' : ['A', 'G', 'C'], 
			'432' : ['C', 'A', 'A'], 
			'434' : ['A', 'T', 'G'], 
			'435' : ['G', 'A', 'C'], 
			'437' : ['A', 'A', 'C'], 
			'441' : ['T', 'C', 'G'], 
			'445' : ['C', 'A', 'C'], 
			'446' : ['A', 'A', 'G'],
			'450' : ['T', 'C', 'G'], 
			'452' : ['C', 'T', 'G'],
			'491' : ['A', 'T', 'C']
	}

#a dictionary of all the codon to amino acid 3 letter codes
codon_aa3l = {
'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys', 'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys', 'TTA': 'Leu', 'TCA': 'Ser', 'TAA': 'STOP', 'TGA': 'STOP', 'TTG': 'Leu', 'TCG': 'Ser', 'TAG': 'STOP', 'TGG': 'Trp',
'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg', 'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', 'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', 'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser', 'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', 'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', 'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly', 'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', 'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'	
}

'''
Function to create classic Xpert outputs
Input:
	The map of genome position to codons (mapdict)
	The SNPs at each of those positiosn in each sample (sampleGenomePos)
	The prefix to put on output file (outPrefix)
Output:
	The final file showing the probes and mutations as found by Xpert	
'''
def xpert(mapdict, sampleGenomePos, outPrefix):#Hardcode the codon to wildtype nucleotides
	

	#hardcode the codon positions with the mutated codons and the associated probe change
	classicXpertmut = {
		'428' : {'AGG':['ProbeA']},
		'430' : {'CCG':['ProbeA']},
		'431' : {'GGC':['ProbeA','ProbeB']}, 
		'432' : {'GAA':['ProbeB']},
		'434' : {'GTG':['ProbeB'],'ATA':['ProbeB'],'ACG':['ProbeB']}, 
		'435' : {'GTC':['ProbeB'],'TAC':['ProbeB_delayed'],'TTC':['ProbeB_delayed'],'GGC':['ProbeB']}, 
		'437' : {'GAC':['ProbeC']},
		'441' : {'CAG':['ProbeC'],'TTG':['ProbeC']},
		'445' : {'TAC':['ProbeD'],'GAC':['ProbeD'],'CTC':['ProbeD'],'AAC':['ProbeD'],'CGC':['ProbeD_delayed'],'GGC':['ProbeD'],'AGC':['ProbeD'],'TCC':['ProbeD'],'CAG':['ProbeD'],'CAA':['ProbeD'],'ACC':['ProbeD']},
		'446' : {'CAG':['ProbeD']},
		'450' : {'TTG':['ProbeE'],'TGG':['ProbeE'],'TTC':['ProbeE']},
		'452' : {'CCG':['ProbeE_delayed']},
	}
	
	#output probe pattern
	outputList=["1","1","1","1","1","1","1","1"]
	outputKeys=["ProbeA","ProbeB","ProbeB_delayed","ProbeC","ProbeD","ProbeD_delayed","ProbeE","ProbeE_delayed"]

	#Create output files
	try:
		xpertSave = open(outPrefix+"_XpertClassic.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	xpertSave.write("Sample name"+"\t"+"RIF Resistance"+"\t"+"Codon number"+"\t"+"Codon"+"\t"+"Capturing Probe"+"\t"+"ProbeA"+"\t"+"ProbeB"+"\t"+"ProbeBdelayed"+"\t"+"ProbeC"+"\t"+"ProbeD"+"\t"+"ProbeDdelayed"+"\t"+"ProbeE"+"\t"+"ProbeEdelayed\n")


	#for each file make a copy of the WT codon patterns, mutate it based on the tab/vcf file and then record the associated result as a string per codon
	#then print out the probe presence absence (0/1) and the mutations
	codonList=list(codon_nuc.keys())
	codonList.sort()
	
	for sample in sampleGenomePos:
		sampleCodons=copy.deepcopy(codon_nuc)
		samplePos=sampleGenomePos[sample]

		for pos in samplePos:#check each position for a SNP and if found, modify the codon
			if samplePos[pos]==None: #if WT then skip
				continue
			#continue only if a SNP is found	
			codon=re.sub("_.*","",mapdict[pos])
			codonPos=int(re.sub(".*_","",mapdict[pos]))
			sampleCodons[codon][codonPos]=samplePos[pos]

		#take the codons from the sample and covert to probes for outputs	
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		capprobeSet=set() #use set as multiple mutations may cause the same probe change so don't count twice
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		capprobeList=list(capprobeSet)
		capprobeList.sort()
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos] and codonPos in classicXpertmut.keys(): #the codons are not the same the codon is part of the Xpert set
				codonCheck="".join(sampleCodons[codonPos])
				codnucres=classicXpertmut[codonPos]
				if codonCheck in codnucres: #if the mutated codon is a known one, get the probe change and change in the output					
					mutCodonValList.append(codonCheck)
					mutCodonList.append(codonPos)					 
					for capprobe in (codnucres[codonCheck]):						
						sampleOutputList[outputKeys.index(capprobe)]="0"
						if capprobe not in capprobeList:
							capprobeList.append(capprobe)
							capprobeSet.add(capprobe)
		
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			xpertSave.write(sample+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(capprobeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			xpertSave.write(sample+"\t"+"NOT DETECTED\n")
	xpertSave.close()	
	
'''
Function to create Hain LPA outputs
Input:
	The map of genome position to codons (mapdict)
	The SNPs at each of those positiosn in each sample (sampleGenomePos)
	The prefix to put on output file (outPrefix)
Output:
	The final file showing the probes and mutations as found by Hain LPA	
'''
def hain(mapdict, sampleGenomePos, outPrefix):#Hardcode the codon to wildtype nucleotides
	

	#hardcode the codon positions with the mutated codons and the associated probe change
	Hainplusmut = {
		'428' : {'AGG':[['WT1'],['']]},
		'430' : {'CCG':[['WT2'],['']]},
		'431' : {'GGC':[['WT2'],['']]}, 
		'432' : {'GAA':[['WT2','WT3'],['']]},
		'434' : {'GTG':[['WT3'],['']],'ATA':[['WT3'],['']],'ACG':[['WT3'],['']]}, 
		'435' : {'GTC':[['WT3','WT4'],['MUT1']],'TAC':[['WT3','WT4'],['']],'TTC':[['WT3','WT4'],['']],'GGC':[['WT3','WT4'],['']]}, 
		'437' : {'GAC':[['WT4'],['']]},
		'441' : {'CAG':[['WT5','WT6'],['']],'TTG':[['WT5','WT6'],['']]},
		'445' : {'TAC':[['WT7'],['MUT2A']],'GAC':[['WT7'],['MUT2B']],'CTC':[['WT7'],['']],'AAC':[['WT7'],['']],'CGC':[['WT7'],['']],'GGC':[['WT7'],['']],'AGC':[['WT7'],['']],'TCC':[['WT7'],['']],'CAG':[['WT7'],['']],'CAA':[['WT7'],['']],'ACC':[['WT7'],['']]},
		'446' : {'CAG':[['WT7'],['']]},
		'450' : {'TTG':[['WT8'],['MUT3']],'TGG':[['WT8'],['']],'TTC':[['WT8'],['']]},
		'452' : {'CCG':[['WT8'],['']]},
	}
	#output probe pattern
	outputList=["1","1","1","1","1","1","1","1","0","0","0","0"]
	outputKeys=["WT1","WT2","WT3","WT4","WT5","WT6","WT7","WT8","MUT1","MUT2A","MUT2B","MUT3"]


	#Create output files
	try:
		hainSave = open(outPrefix+"_HainLPA.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	hainSave.write("Sample name"+"\t"+"RIF Resistance"+"\t"+"Codon number"+"\t"+"Codon"+"\t"+"Absent Probe"+"\t"+"Developing Probe"+"\t"+"WT1"+"\t"+"WT2"+"\t"+"WT3"+"\t"+"WT4"+"\t"+"WT5"+"\t"+"WT6"+"\t"+"WT7"+"\t"+"WT8"+"\t"+"MUT1"+"\t"+"MUT2A"+"\t"+"MUT2B"+"\t"+"MUT3\n")
  

	#for each file make a copy of the WT codon patterns, mutate it based on the tab/vcf file and then record the associated result as a string per codon
	#then print out the probe presence absence (0/1) and the mutations
	codonList=list(codon_nuc.keys())
	codonList.sort()
	
	for sample in sampleGenomePos:
		sampleCodons=copy.deepcopy(codon_nuc)
		samplePos=sampleGenomePos[sample]
		for pos in samplePos:#check each position for a SNP and if found, modify the codon
			if samplePos[pos]==None: #if WT then skip
				continue
				
			#continue only if a SNP is found	
			codon=re.sub("_.*","",mapdict[pos])
			codonPos=int(re.sub(".*_","",mapdict[pos]))
			sampleCodons[codon][codonPos]=samplePos[pos]
		
		#take the codons from the sample and covert to probes for outputs	
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		wtprobeList=[]
		mutprobeList=[]
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos] and codonPos in Hainplusmut.keys(): #the codons are not the same and the codon is part of the Hain set
				codonCheck="".join(sampleCodons[codonPos])
				codnucres=Hainplusmut[codonPos]				
				if codonCheck in codnucres:
					mutCodonValList.append(codonCheck)
					mutCodonList.append(codonPos)
					for wtprobe in (codnucres[codonCheck][0]):
						sampleOutputList[outputKeys.index(wtprobe)]="0"
						if wtprobe not in wtprobeList:
							wtprobeList.append(wtprobe)
					for mutprobe in (codnucres[codonCheck][1]):	
						if mutprobe in outputKeys:
							sampleOutputList[outputKeys.index(mutprobe)]="1"						
							mutprobeList.append(mutprobe)
		
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			hainSave.write(sample+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(wtprobeList)+"\t"+",".join(mutprobeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			hainSave.write(sample+"\t"+"NOT DETECTED\n")
	hainSave.close()
	
'''
Function to create Nipro LPA outputs
Input:
	The map of genome position to codons (mapdict)
	The SNPs at each of those positiosn in each sample (sampleGenomePos)
	The prefix to put on output file (outPrefix)
Output:
	The final file showing the probes and mutations as found by the Nipro LPA
'''
def nipro(mapdict, sampleGenomePos, outPrefix):#Hardcode the codon to wildtype nucleotides
	
	#hardcode the codon positions with the mutated codons and the associated capturing probe
	Nipromut = {
			'428' : {'AGG':[['S1'],['']]},
			'430' : {'CCG':[['S1'],['']]}, 
			'431' : {'GGC':[['S1'],['']]}, 
			'432' : {'GAA':[['S1'],['']]},
			'434' : {'GTG':[['S2'],['']],'ATA':[['S2'],['']],'ACG':[['S2'],['']]}, 
			'435' : {'GTC':[['S2'],['R2']],'TAC':[['S2'],['']],'TTC':[['S2'],['']],'GGC':[['S2'],['']]}, 
			'437' : {'GAC':[['S2'],['']]},
			'441' : {'CAG':[['S3'],['']],'TTG':[['S3'],['']]},
			'445' : {'TAC':[['S4'],['R4a']],'GAC':[['S4'],['R4b']],'CTC':[['S4'],['']],'AAC':[['S4'],['']],'CGC':[['S4'],['']],'GGC':[['S4'],['']],'AGC':[['S4'],['']],'TCC':[['S4'],['']],'CAG':[['S4'],['']], 'CAA':[['S4'],['']], 'ACC':[['S4'],['']]},
			'446' : {'CAG':[['S4'],['']]},
			'450' : {'TTG':[['S5'],['R5']],'TGG':[['S5'],['']],'TTC':[['S5'],['']]}, 
			'452' : {'CCG':[['S5'],['']]},
	}
	#output probe pattern
	outputList=["1","1","1","1","1","0","0","0","0"]
	outputKeys=["S1","S2","S3","S4","S5","R2","R4a","R4b","R5"]



	#Create output files
	try:
		niproSave = open(outPrefix+"_NiproLPA.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	niproSave.write("Sample name"+"\t"+"RIF Resistance"+"\t"+"Codon number"+"\t"+"Codon"+"\t"+"Absent Probe"+"\t"+"Developing Probe"+"\t"+"S1"+"\t"+"S2"+"\t"+"S3"+"\t"+"S4"+"\t"+"S5"+"\t"+"R2"+"\t"+"R4a"+"\t"+"R4b"+"\t"+"R5\n")


	#for each file make a copy of the WT codon patterns, mutate it based on the tab/vcf file and then record the associated result as a string per codon
	#then print out the probe presence absence (0/1) and the mutations
	codonList=list(codon_nuc.keys())
	codonList.sort()
	
	for sample in sampleGenomePos:
		sampleCodons=copy.deepcopy(codon_nuc)
		samplePos=sampleGenomePos[sample]
		for pos in samplePos:#check each position for a SNP and if found, modify the codon
			if samplePos[pos]==None: #if WT then skip
				continue
				
			#continue only if a SNP is found	
			codon=re.sub("_.*","",mapdict[pos])
			codonPos=int(re.sub(".*_","",mapdict[pos]))
			sampleCodons[codon][codonPos]=samplePos[pos]
		
		#take the codons from the sample and covert to probes for outputs	
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		wtprobeList=[]
		mutprobeList=[]
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos] and codonPos in Nipromut.keys(): #the codons are not the same and the codon is part of the Hain set 
				codonCheck="".join(sampleCodons[codonPos])
				codnucres=Nipromut[codonPos]
				if codonCheck in codnucres:
					mutCodonValList.append(codonCheck)
					mutCodonList.append(codonPos)
					for wtprobe in (codnucres[codonCheck][0]):
						sampleOutputList[outputKeys.index(wtprobe)]="0"
						if wtprobe not in wtprobeList:
							wtprobeList.append(wtprobe)
					for mutprobe in (codnucres[codonCheck][1]):	
						if mutprobe in outputKeys:
							sampleOutputList[outputKeys.index(mutprobe)]="1"						
							mutprobeList.append(mutprobe)
		
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			niproSave.write(sample+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(wtprobeList)+"\t"+",".join(mutprobeList)+"\t"+"\t".join(sampleOutputList)+"\n")
		else: #is WT so inform
			niproSave.write(sample+"\t"+"NOT DETECTED\n")
	niproSave.close()

'''
Function to create Xpert Ultra outputs
Input:
	The map of genome position to codons (mapdict)
	The SNPs at each of those positiosn in each sample (sampleGenomePos)
	The prefix to put on output file (outPrefix)
Output:
	The final file showing the probes and mutations as found by Xpert Ultra	
'''
def ultra(mapdict, sampleGenomePos, outPrefix):#Hardcode the codon to wildtype nucleotides
	
	Ultramut = {
			'428' : {'AGG':[['rpoB1'],['3.5']]},
			'430' : {'CCG':[['rpoB1'],['5.9-6.3']]}, 
			'431' : {'GGC':[['rpoB1'],['2.9']]}, 
			'432' : {'GAA':[['rpoB1'],['3.4']]},
			'434' : {'GTG':[['rpoB1'],['2.5']],'ATA':[['rpoB2'],['3.2']],'ACG':[['rpoB2'],['3.3']]}, 
			'435' : {'GGC':[['rpoB2'],['3.3']],'TAC':[['rpoB2'],['4.0-4.4']],'TTC':[['rpoB2'],['5.3']],'GTC':[['rpoB2'],['3.5-3.7']]}, 
			'441' : {'CAG':[['rpoB2','rpoB3'],['4.7','2.3']],'TTG':[['rpoB2','rpoB3'],['3.0','2.3']]},
			'445' : {'TAC':[['rpoB3'],['3.2-3.3']],'GAC':[['rpoB3'],['3.7-3.9']],'CTC':[['rpoB3'],['3.5-3.6']],'AAC':[['rpoB3'],['3.5-3.6']],'CGC':[['rpoB3'],['3.5-3.6']],'GGC':[['rpoB3'],['3.5-3.6']],'AGC':[['rpoB3'],['3.5-3.6']],'TCC':[['rpoB3'],['3.5-3.6']],'CAG':[['rpoB3'],['3.5-3.6']],'CAA':[['rpoB3'],['3.5-3.6']],'ACC':[['rpoB3'],['4.9']]}, 
			'446' : {'CAG':[['rpoB4B'],['5.0']]},
			'450' : {'TTG':[['rpoB3','rpoB4A'],['2.5-2.9','6.0-6.5']],'TGG':[['rpoB3','rpoB4A'],['2.3-2.7','3.3-3.7']],'TTC':[['rpoB3'],['4.0']]}, 
			'452' : {'CCG':[['rpoB4B'],['5.7-6.1']]},
	}
	#output probe pattern
	outputList=["1","1","1","1","1"]
	outputKeys=["rpoB1","rpoB2","rpoB3","rpoB4A","rpoB4B"]




	#Create output files
	try:
		ultraSave = open(outPrefix+"_XpertUltra.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	ultraSave.write("Sample name"+"\t"+"RIF Resistance"+"\t"+"Codon number"+"\t"+"Codon"+"\t"+"Capturing Probe"+"\t"+"rpoB1"+"\t"+"rpoB2"+"\t"+"rpoB3"+"\t"+"rpoB4A"+"\t"+"rpoB4B"+"\t"+"deltaTm"+"\n")


	#for each file make a copy of the WT codon patterns, mutate it based on the tab/vcf file and then record the associated result as a string per codon
	#then print out the probe presence absence (0/1) and the mutations
	codonList=list(codon_nuc.keys())
	codonList.sort()
	
	for sample in sampleGenomePos:
		sampleCodons=copy.deepcopy(codon_nuc)
		samplePos=sampleGenomePos[sample]
		for pos in samplePos:#check each position for a SNP and if found, modify the codon
			if samplePos[pos]==None: #if WT then skip
				continue
				
			#continue only if a SNP is found	
			codon=re.sub("_.*","",mapdict[pos])
			codonPos=int(re.sub(".*_","",mapdict[pos]))
			sampleCodons[codon][codonPos]=samplePos[pos]
		
		#take the codons from the sample and covert to probes for outputs	
		sampleOutputList=copy.deepcopy(outputList)
		mutCodonList=[]
		mutCodonValList=[]
		capprobeList=[]
		deltaTmList=[]
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos] and codonPos in Ultramut.keys(): #the codons are not the same and the codon is part of the Hain set 
				codonCheck="".join(sampleCodons[codonPos])
				codnucres=Ultramut[codonPos]				
				if codonCheck in codnucres:					
					mutCodonValList.append(codonCheck)
					mutCodonList.append(codonPos)
					for capprobe in (codnucres[codonCheck][0]):
						sampleOutputList[outputKeys.index(capprobe)]="0"
						deltaTmList.append(codnucres[codonCheck][1])
						if capprobe not in capprobeList:
							capprobeList.append(capprobe)
		
		#if the sample output differs from the WT output, write relevant information
		if sampleOutputList!=outputList:
			ultraSave.write(sample+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(capprobeList)+"\t"+"\t".join(sampleOutputList)+"\t"+",".join((deltaTmList.pop()))+"\n")
		else: #is WT so inform
			ultraSave.write(sample+"\t"+"NOT DETECTED\n")
	ultraSave.close()	

'''
Function to create Sanger sequence outputs
Input:
	The map of genome position to codons (mapdict)
	The SNPs at each of those positions in each sample (sampleGenomePos)
	The prefix to put on output file (outPrefix)
Output:
	The final file showing the mutations as found by Sanger sequencing at the given codon positions (based on dictionary at top of this module)
'''
def sanger(mapdict, sampleGenomePos, outPrefix):

	#Create output files
	try:
		sangerSave = open(outPrefix+"_Sanger.txt", "w")  
	except IOError:
		print('no room for save file')				
		sys.exit()
	sangerSave.write("Sample name"+"\t"+"RIF Resistance"+"\t"+"Codon number"+"\t"+"Codon"+"\t"+"RR mutation"+"\n")


	#for each file make a copy of the WT codon patterns, mutate it based on the tab/vcf file and then record the associated result as a string per codon
	#then print out the mutations
	codonList=list(codon_nuc.keys())
	codonList.sort()
	
	for sample in sampleGenomePos:
		sampleCodons=copy.deepcopy(codon_nuc)
		samplePos=sampleGenomePos[sample]
		for pos in samplePos:#check each position for a SNP and if found, modify the codon
			if samplePos[pos]==None: #if WT then skip
				continue
				
			#continue only if a SNP is found	
			codon=re.sub("_.*","",mapdict[pos])
			codonPos=int(re.sub(".*_","",mapdict[pos]))
			sampleCodons[codon][codonPos]=samplePos[pos]
		
		mutCodonList=[]
		mutCodonValList=[]
		RRconferringmut=[]
		
		#compare the wt and the sample codon lists and where they differ, put a 0 in the appropriate probe
		for codonPos in codonList:
			if codon_nuc[codonPos]!=sampleCodons[codonPos] : #the codons are not the same
				wtCodon="".join(codon_nuc[codonPos])
				mutCodon="".join(sampleCodons[codonPos])
				wtAA=codon_aa3l[wtCodon]
				mutAA=codon_aa3l[mutCodon]
				RRmut=wtAA+codonPos+mutAA
				mutCodonValList.append(mutCodon)
				mutCodonList.append(codonPos)
				RRconferringmut.append(RRmut)
		
		#if the sample output differs from the WT output, write relevant information
		if sampleCodons!=codon_nuc:
			sangerSave.write(sample+"\t"+"DETECTED"+"\t"+",".join(mutCodonList)+"\t"+",".join(mutCodonValList)+"\t"+",".join(RRconferringmut)+"\n")
		else: #is WT so inform
			sangerSave.write(sample+"\t"+"NOT DETECTED\n")
			
	sangerSave.close()
