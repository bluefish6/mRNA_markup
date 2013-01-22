#!/usr/bin/env python
#-*- coding: utf-8 -*-
import Bio.SeqIO
import Bio.Alphabet
import Bio.Blast.NCBIWWW
import Bio.Blast.Applications
import Bio.Blast.NCBIXML
import sys
import subprocess

if len(sys.argv)<8:
	print "Not enough arguments provided."
	print "Usage: "+sys.argv[0]+" sequence_file.fasta evalue bacterial_db vector_db domain_db refprot_db protein_db"
	print "For example: "+sys.argv[0]+" mRNA_in.fasta 1e-20 bacteriaDB uniVec CDD atpeptair10 uniprot90"
	print "Warning: all databases should be in this directory."
	sys.exit()
	
input = sys.argv[1] 
e_value = sys.argv[2]
baza_bakteryjna = sys.argv[3]
baza_wektorow = sys.argv[4]
baza_domen = sys.argv[5]
baza_referencyjna_bialek = sys.argv[6]
baza_wszystkich_bialek = sys.argv[7]


sequences = list(Bio.SeqIO.parse(input,"fasta", Bio.Alphabet.generic_rna))

#krok 1 - eliminacja zanieczyszczenia wektorami
blastn_cline = Bio.Blast.Applications.NcbiblastnCommandline(query=input, db=baza_wektorow, evalue=e_value, outfmt=5, out="step1.xml") #przygotowanie blasta
stdout, stderr = blastn_cline() #uruchomienie blasta
blastXMLhandle1=open("step1.xml")

wyniki = Bio.Blast.NCBIXML.parse(blastXMLhandle1) #analiza wyprodukowanego xmla przez blasta
step1Good = []
step1Bad = []
i=0
for wynik in wyniki:
	if len(wynik.descriptions)>0: #sprawdzamy czy byly jakies hity
		#sekwencja zanieczyszczona wektorami
		step1Bad.append(sequences[i])
	else:
		#sekwencja niezanieczyszczona
		step1Good.append(sequences[i])
	i+=1
Bio.SeqIO.write(step1Bad,"VC-mRNA.fas", "fasta") #wypisywanie do pliku zanieczyszczonych sekwencji
Bio.SeqIO.write(step1Good,"not_VC-mRNA.fas", "fasta") #wypisywanie do pliku niezanieczyszczonych sekwencji

#krok 2 - eliminacja zanieczyszczenia bakteriami
blastn_cline = Bio.Blast.Applications.NcbiblastnCommandline(query="not_VC-mRNA.fas", db=baza_bakteryjna, evalue=e_value, outfmt=5, out="step2.xml") #przygotowanie blasta
stdout, stderr = blastn_cline() #uruchomienie blasta
blastXMLhandle2=open("step2.xml")

wyniki = Bio.Blast.NCBIXML.parse(blastXMLhandle2) #analiza wyprodukowanego xmla przez blasta
step2Good = []
step2Bad = []
i=0
for wynik in wyniki:
	if len(wynik.descriptions)>0: #sprawdzamy czy byly jakies hity
		#sekwencja zanieczyszczona bakteriami
		step2Bad.append(step1Good[i])
	else:
		step2Good.append(step1Good[i]) #sekwencja niezanieczyszczona
	i+=1
Bio.SeqIO.write(step2Bad,"BC-mRNA.fas", "fasta") #wypisywanie do pliku zanieczyszczonych sekwencji
Bio.SeqIO.write(step2Good,"not_BC-mRNA.fas", "fasta") #wypisywanie do pliku niezanieczyszczonych sekwencji

#krok 3 - porownywanie z baza bialkowa
step3Good = []
step3Bad = []
blastx_cline = Bio.Blast.Applications.NcbiblastxCommandline(query="not_BC-mRNA.fas", db=baza_referencyjna_bialek, evalue=e_value, outfmt=0, out="step3.txt") #przygotowanie blasta do museqbxa
stdout, stderr = blastx_cline() #uruchomienie blasta
blastx_cline = Bio.Blast.Applications.NcbiblastxCommandline(query="not_BC-mRNA.fas", db=baza_referencyjna_bialek, evalue=e_value, outfmt=5, out="step3.xml") #przygotowanie blasta do pozostalych analiz
stdout, stderr = blastx_cline() #uruchomienie blasta
subprocess.call("./MuSeqBox -q -o RA.msb -l 40 -p 4 -infname step3.txt", shell=True) #uruchomienie museqboxa

blastXMLhandle3=open("step3.xml") 
wyniki = Bio.Blast.NCBIXML.parse(blastXMLhandle3) #analiza wyprodukowanego xmla przez blasta
i=0
for wynik in wyniki:
	if len(wynik.descriptions)>0: #sprawdzamy czy byly jakies hity
		#sekwencja ma odpowiedniki w bazie bialkowej
		step3Good.append(step2Good[i])
	else:
		#sekwencja nie ma odpowiednikow w bazie bialkowej
		step3Bad.append(step2Good[i])
	i+=1
Bio.SeqIO.write(step3Bad,"not_RA-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji bez hitow do bazy bialkowej
Bio.SeqIO.write(step3Good,"RA-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji z hitami do bazy bialkowej


#krok 3.1 - identyfikacja potencjalnych pelnych sekwencji kodujacych
subprocess.call("./MuSeqBox -q -o output.txt -l 40 -p 4 -n 1 -s 1 -F 5 5 10 10 90 60 fullcds.txt -infname step3.txt", shell=True) #uruchomienie museqboxa

museq31=open("fullcds.txt") #analiza wynikow museqboxa
set31 = set()
set31seq = set()
set31seqBad = set()
record_flag = False #czy przebrnelismy przez naglowek i czy nie dotarlismy do konca pliku
for line in museq31:
	if "---------------------" in line: #oddziela naglowek i stopke od tresci
		record_flag = not(record_flag)
	elif (record_flag):
		set31.add(line.split(' ')[0].strip()) #identyfikator jest na poczatku
for seq in step3Good:
	if (seq.name in set31) or (seq.name.split("|")[-1] in set31): #identyfikator z museqboxa moze nie miec identycznej struktury jak w pliku oryginalnym
		set31seq.add(seq)
	else:
		set31seqBad.add(seq)
		
Bio.SeqIO.write(set31seq,"FL-mRNA.fas", "fasta") #zapisanie do pliku pelnych sekwencji
Bio.SeqIO.write(set31seqBad,"not_FL-mRNA.fas", "fasta") #zapisanie niepelnych sekwencji

#krok 3.2 - identyfikacja potencjalnych sekwencji chimerycznych
subprocess.call("./MuSeqBox -M 10 50 -l 40 -p 4 -i RA.msb -o PC.msb", shell=True)

museq32=open("PC.msb") #analiza wynikow museqboxa z poprzedniego kroku
set32 = set()
set32seq = set()
set32seqBad = set()
for line in museq32:
  if "!Potential" in line: #jest taki napis przy rekordach ktore sa potencjalnymi sekwencjami chimerycznymi
	set32.add(line.split(":")[1].split("matches")[0].strip().partition(' ')[2]) #tak uzyskujemy identyfikator
for seq in set31seqBad:
	if seq.name in set32 or seq.name.split("|")[-1] in set32: #identyfikator z museqboxa moze nie miec identycznej struktury jak w pliku oryginalnym
		set32seq.add(seq)
	else:
		set32seqBad.add(seq)
		
Bio.SeqIO.write(set32seq,"PC-mRNA.fas", "fasta") #wypisujemy potencjalne sekwencje chimeryczne
Bio.SeqIO.write(set32seqBad,"not_PC-mRNA.fas", "fasta") #wypisujemy sekwencje niechimeryczne

#krok 4 - porownanie z odpowiednia baza bialkowa
step4Good = []
step4Bad = []
blastx_cline = Bio.Blast.Applications.NcbiblastxCommandline(query="not_RA-mRNA.fas", db=baza_wszystkich_bialek, evalue=e_value, outfmt=5, out="step4.xml") #przygotowanie blasta
stdout, stderr = blastx_cline() #uruchomienie blasta
blastXMLhandle4=open("step4.xml")
wyniki = Bio.Blast.NCBIXML.parse(blastXMLhandle4) #analiza wynikow
i=0
for wynik in wyniki:
	if len(wynik.descriptions)>0:
		#sekwencja ma odpowiedniki w bazie
		step4Good.append(step3Bad[i])
	else:
		#sekwencja nie ma odpowiednikow w bazie
		step4Bad.append(step3Bad[i])
	i+=1
Bio.SeqIO.write(step4Bad,"not_AA-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji bez hitow do bazy bialkowej
Bio.SeqIO.write(step4Good,"AA-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji z hitami do bazy bialkowej

#krok 5 - porownanie do bazy domen bialkowych
step5Good = []
step5Bad = []
rpstblastn_cline = Bio.Blast.Applications.NcbirpstblastnCommandline(query="not_AA-mRNA.fas", db=baza_domen, evalue=e_value, outfmt=5, out="step5.xml") #przygotowanie blasta
stdout, stderr = rpstblastn_cline() #uruchomienie blasta
blastXMLhandle5=open("step5.xml")
wyniki = Bio.Blast.NCBIXML.parse(blastXMLhandle5) #analiza wynikow
i=0
for wynik in wyniki:
	if len(wynik.descriptions)>0:
		#sekwencja ma odpowiedniki w bazie
		step5Good.append(step4Bad[i])
	else:
		#sekwencja nie ma odpowiednikow w bazie
		step5Bad.append(step4Bad[i])
	i+=1
Bio.SeqIO.write(step5Bad,"not_CD-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji bez hitow do bazy domenowej
Bio.SeqIO.write(step5Good,"CD-mRNA.fas", "fasta") #wypisywanie do pliku sekwencji z hitami do bazy domenowej

summary = open("summary.txt", "w") #zapisanie podsumowania do pliku
summary.write("mRNAmarkup Report\n  Number of input sequences:\t\t\t\t"+str(len(sequences))+"\n  Number of potential vector-contaminated sequences:\t"+str(len(step1Bad))+" (file: VC_mRNA.fas)\n  Number of potential bacterial-contaminated sequences:\t"+str(len(step2Bad))+" (file: BC_mRNA.fas)\n  Number of sequences matching the ReferenceDB:\t\t"+str(len(step3Good))+" (file: RA_mRNA.fas)\n    Number of potential full-length coding sequences:\t\t"+str(len(set31seq))+" (file: FL_mRNA.fas)\n    Non-qualifying sequences:\t\t\t\t\t"+str(len(set31seqBad))+" (file: not_FL_mRNA.fas)\n      Number of potential chimeric sequences:\t\t\t  "+str(len(set32seq))+" (file: PC_mRNA.fas)\n      Non-qualifying sequences:\t\t\t\t\t  "+str(len(set32seqBad))+" (file: not_PC_mRNA.fas)\n  Number of sequences matching the AllProteinDB:\t"+str(len(step4Good))+" (file: AA_mRNA.fas)\n  Number of sequences matching the ProteinDomainDB:\t"+str(len(step5Good))+" (file: CD_mRNA.fas)\n  Number of remaining sequences:\t\t\t"+str(len(sequences)-len(step1Bad)-len(step2Bad)-len(step3Good)-len(step4Good)-len(step5Good))+" (file: remaining-mRNA.fas)\n")
