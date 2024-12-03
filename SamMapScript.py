#!/usr/bin/python3
# -*- coding : utf-8 -*-

#licence auteur et version

__authors__ = ("Loik Galtier")

import os
import re
import csv
import sys
import pandas as pd


######## Setting the file ########
# check the selected file exist and have the .sam
def findFile():
	if len(sys.argv) == 1:
		noFileFound = True
		while noFileFound:
			name = input("Quel est le chemin d'accés au fichier sam : ?")  # Ask the path
			if os.path.exists(name) and re.search(r".sam", name):  # if the path exist AND the .sam exist
				noFileFound = False  # quit the while
				print("Il y as un fichier sam")  # continue -> Ask what question to respond (1, 2, 3, 4)
				return name  # return the file open
			else:
				print("Il n'y as pas de fichier sam, try again \n")  # Ask again for a file.
	else:
		if os.path.exists(sys.argv[1]) and re.search(r".sam", sys.argv[1]):  # if the path exist AND the .sam exist
			print("Il y as un fichier sam")  # continue -> Ask what question to respond (1, 2, 3, 4)
			return sys.argv[1]  # return the file open
		else:
			print("Erreur : Aucun fichier .sam n'a était trouvé \n")
			sys.exit()


#### Define quality ####
def askQuality():
	if len(sys.argv) == 1:
		noAnswer = True
		while noAnswer:
			needQualityMin = input(
				"Voulez vous appliqué un filtre de qualité minimal pour la suite du traitement : (Y/N)").upper()  # Ask the path
			if needQualityMin == "Y":
				QualityMin = input("Quel valeur minimal voulais vous appliquer : [0-255] \n")
				if QualityMin.isdigit():
					print(
						"Toute les prochaines opérations se feront sur les reads ayant un qualité supérieur à : " + QualityMin + "\n")
					noAnswer = False
					return int(QualityMin)
				else:
					print("Cela n'est pas une valeur correcte \n")
			elif needQualityMin == "N":
				noAnswer = False
				print("\n")
				return 0
			else:
				print("Valeur non correcte \n")
				
	elif len(sys.argv) >= 3 :
		if sys.argv[2].isdigit():
			print("Toute les prochaines opérations se feront sur les reads ayant un qualité supérieur à : " + sys.argv[
				2] + "\n")
			return int(sys.argv[2])
		elif sys.argv[2].isalpha():
			print("Erreur : Cela n'est pas une valeur correcte pour la qualité")
			sys.exit()
		else:
			return 0
	
	else:
		return 0




######## Extract part and use header ########
def sam_reading(sam_file_path):
	# Ouvrir le fichier SAM
	with open(sam_file_path, "r") as sam_file:
		
		# Parcourir chaque ligne du fichier SAM
		sequenceNames = []
		
		flags = []
		rnames = []
		mapqs = []
		cigars = []
		
		numberOfReadTotal = 0
		
		# affichage info header
		hdTable = {}
		sqTable = {}
		pgTable = {}
		
		# Extraire les informations du header
		for line in sam_file:
			
			if line.startswith("@"):
				if line.startswith("@HD"):
					vN = re.search(r"VN:([^\t|\n]+)", line)
					sO = re.search(r"SO:([^\t|\n]+)", line)
					if vN:
						hdTable[vN.group(1)] = ""
						if sO:
							hdTable[vN.group(1)]= [sO.group(1)]
				
				elif line.startswith("@PG"):
					idnamepg = re.search(r"ID:([^\t|\n]+)", line)
					versionpg = re.search(r"VN:([^\t|\n]+)", line)
					if idnamepg:
						pgTable[idnamepg.group(1)] = ""
						if versionpg:
							pgTable[idnamepg.group(1)] = versionpg.group(1)
						
				elif line.startswith("@SQ"):
					sN = re.search(r"SN:([^\t|\n]+)", line)
					lN = re.search(r"LN:([^\t|\n]+)", line)
					if sN and lN:
						sequenceNames.append(sN.group(1))
						sqTable[sN.group(1)] = lN.group(1)
						
				elif line.startswith("@RG"):
					idgroup = re.search(r"ID:([^\t|\n]+)", line)
					if idgroup:
						print("Un groupe de lecture a était formé, son ID unique est : " + idgroup.group(1) + "\n")
		sam_file.seek(0)  # pour arreter la boucle
		
		print("_________________ \n Entête :\n")
		
		if (hdTable):
			print("Les données de ce fichier sont :")
			data = [(f"{cle} |", f"{val} |") for cle, val in hdTable.items()]
			t_data = pd.DataFrame(data, columns=["Format Version |", "Ordre |"])
			print(t_data, "\n")
		
		if (sqTable):
			print("Des séquences références en était utilisé :")
			data = [(f"{cle} |", f"{val} |") for cle, val in sqTable.items()]
			t_data = pd.DataFrame(data, columns=["séquence name |", "Length |"])
			print(t_data, "\n")
		
		if (pgTable):
			print("Des programes en était utilisé :")
			
			data = [(f"{cle} |", f"{val} |") for cle, val in pgTable.items()]
			t_data = pd.DataFrame(data, columns=["Programme id |", "Version |"])
			print(t_data, "\n")
		# Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation
		sam_reader = csv.reader(sam_file, delimiter='\t')
		for row in sam_reader:  # Pour chaque ligne de mon fichier SAM
			if row[0].startswith("@"):  # Ignorer les lignes d'en-tête qui commencent par '@'
				continue
			if (int(row[4]) >= QualityMin and int(row[4]) < 255):  # si la ligne a une qualité supérieur à la qualité demandé et inférieur a 255
				
				# Accéder aux informations de chaque colonne
				flag = int(row[1])  # Flag
				flags.append(flag)  # rajouter le flag de chaque ligne dans la liste
				
				rname = row[2]  # Nom de la séquence de référence
				rnames.append(rname)
				
				mapq = row[4]
				mapqs.append(mapq)
				
				cigar = row[5]  # Chaîne CIGAR
				cigars.append(cigar)
				
				numberOfReadTotal += 1
		# print(f"QNAME: {qname}, FLAG: {flag}, RNAME: {rname}, POS: {pos}, CIGAR: {cigar}, SEQ: {seq}")
	return sequenceNames, flags, rnames,mapqs , cigars, numberOfReadTotal


######## Convert to binary ########
#fonction flags_to_binary pour convertir le flag en binaire
def flags_to_binary(flags):
	for i in range(len(flags)):  # Boucle qui va parcourir de 0 à la taille des flags
		flags[i] = bin(flags[i])  # Convertir en binaire
		flags[i] = flags[i][2:].zfill(
			12)  # Éliminer les 2 premiers caractères du flag (0 et b) & remplir de 0 sur la gauche jusqu'a 12
	# print(flags[i])
	return flags  # Retour des flags en binaire


######## Question 1 ########
# fonction analysis_of_mapped_reads, prend en paramètre les flags en binaire
def analysis_of_mapped_reads(binary_flags):
	nbMappedTable = {}
	nbSemmappedTable = {}
	nbUnmapped = 0
	
	# mapped part
	for sqname in sequenceRefName:
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if (rname[i] == sqname):  # verifier si le read est mappé sur la sequence de reference
				flag = binary_flags[i]  # Mettre un flag en binaire dans la variable flag
				if (flag[-3] == "0"):  # le bit 3 code pour l'info "read non mappé". Si "0" = mappé
					if (sqname in nbMappedTable):
						nbMappedTable[sqname] += 1  # Rajouter 1 pour compter le nombre de read
					else:
						nbMappedTable[sqname] = 1  # Sinon initialisé
	
	# semimapped part
	for sqname in sequenceRefName:
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if (rname[i] == sqname):  # verifier si le read est mappé sur la sequence de reference
				flag = binary_flags[i]  # Mettre un flag en binaire dans la variable flag
				if flag[-2] == "1":  # le bit 2 code pour l'info "aligné". Si "1" = aligné.
					if cigars[i] != "100M":  # si 100M, alors il n'y as pas de problème d'alignement, en ignore
						if (sqname in nbSemmappedTable):  # Si sqname existe dans la table des semi mappé
							nbSemmappedTable[sqname] += 1  # Rajouter 1 pour compter le nombre de read
						else:
							nbSemmappedTable[sqname] = 1  # Sinon initialise
		
		if sqname in nbSemmappedTable and sqname in nbMappedTable:  # Si la séquence se trouve dans ma table semi mappé ET mappé,
			nbMappedTable[sqname] -= nbSemmappedTable[
				sqname]  # enlevé les occurences semi mappé de la table mappé complétement
	
	# unmapped part
	for i in range(len(binary_flags)):
		flag = binary_flags[i]  # Mettre un flag en binaire dans la variable flag
		if flag[-3] == "1" or rname[i] == "*" or rname[i] not in sequenceRefName:  # le bit 3 code pour l'info "read non mappé". Si "1" = non mappé
			nbUnmapped += 1  # Rajouter 1 pour compter le nombre de read
	
	print("_________________ \n Alignement :\n")
	print("Parmis", totalNumberOfRead, "reads, les reads sont :")
	
	if (nbMappedTable):
		data1 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" Mappé |") for cle, val in
		         nbMappedTable.items()]  #
		t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |", "Quantité |", " 100% |", " Alligné ?|"])
		t_data = pd.concat([t_data1], axis=0)
	
	if (nbSemmappedTable and nbMappedTable):
		data2 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" Semi-mappé |") for
		         cle, val in nbSemmappedTable.items()]  #
		t_data2 = pd.DataFrame(data2, columns=[" Sequence reference |", "Quantité |", " 100% |", f" Alligné ?|"])
		t_data = pd.concat([t_data1, t_data2], axis=0)
	
	if (nbUnmapped > 0):
		data3 = [(" Non aligné |", f"{nbUnmapped} |", f" {round(nbUnmapped / totalNumberOfRead, 4) * 100} |",
		          f" non mappé |")]  #
		t_data3 = pd.DataFrame(data3, columns=[" Sequence reference |", "Quantité |", " 100% |", f" Alligné ?|"])
		if (nbSemmappedTable and nbMappedTable):
			t_data = pd.concat([t_data1, t_data2, t_data3], axis=0)
		elif (nbMappedTable):
			t_data = pd.concat([t_data1, t_data3], axis=0)
		else:
			t_data = pd.concat([t_data3], axis=0)
	
	# t_data = pd.concat([t_data1, t_data2, t_data3], axis=0)
	print(t_data, "\n")


#######  ########
def count_flag_number():

	print("_________________ \n Quantité des flags :\n")
	
	flags_dico_description = {
		"Flags2048": "supplementary alignment",
		"Flags1024": "PCR or optical duplicate",
		"Flags512": "not passing filters",
		"Flags256": "secondary alignment",
		"Flags128": "the last segment in the template",
		"Flags64": "the first segment in the template",
		"Flags32": "SEQ of the next segment reverse complemented",
		"Flags16": "SEQ being reverse complemented",
		"Flags8": "next segment in the template unmapped",
		"Flags4": "segment unmapped",
		"Flags2": "each segment properly aligned",
		"Flags1": "template having multiple segments sequencing"
	}
	
	flags_dico = {
		"Flags2048": 0,
		"Flags1024": 0,
		"Flags512": 0,
		"Flags256": 0,
		"Flags128": 0,
		"Flags64": 0,
		"Flags32": 0,
		"Flags16": 0,
		"Flags8": 0,
		"Flags4": 0,
		"Flags2": 0,
		"Flags1": 0
	}
	
	flags_index = list(flags_dico.keys())
	pd.set_option('display.max_colwidth', 40)
	for i in range(len(binary_flags)):
		# Mettre un flag en binaire dans la variable flag
		flag = binary_flags[i]
		
		for j in range(12):
			if flag[j] == "1":
				flags_dico[flags_index[j]] += 1
	
	data1 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |") for cle, val in
	         flags_dico.items()]  #
	t_data1 = pd.DataFrame(data1, columns=[" Flags |", "Quantité |", " 100% |"])
	
	data2 = [(f"{val} |") for cle, val in flags_dico_description.items()]  #
	t_data2 = pd.DataFrame(data2, columns=[" Description |"])
	t_data = pd.concat([t_data1, t_data2], axis=1)
	print(t_data)

### ####
def pair_Analysis():
	dico_Pair = {
		# Pos : Pnext        actuel position : position d'une possible pair
	    }
	
	with open(sam_file_path, "r") as sam_file:
		sam_reader = csv.reader(sam_file, delimiter='\t')
		for row in sam_reader:  # Pour chaque ligne de mon fichier SAM
			if row[0].startswith("@"):  # Ignorer les lignes d'en-tête qui commencent par '@'
				continue
			if (int(row[4]) >= QualityMin and int(row[4]) < 255):  # si la ligne a une qualité supérieur à la qualité demandé et inférieur a 255
				if((row[6] == row[2] and row[6] != "*") or row[6] == "="):
					dico_Pair[row[3]] = row [7]
					
	
				
#### ####
def quality_analysis():
	
	
	qual_dico_sequence={}
	
	for sqname in sequenceRefName:
		qual_dico = {
			'0':0,
			'1-30':0,
			'31-60':0,
			'61 -':0
		}
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if (rname[i] == sqname):        # verifier si le read est mappé sur la sequence de reference
				if(mapqs[i].isdigit() and mapqs[i] != '0'):     #verifier si mapQ est bien un chiffre
					if(int(mapqs[i]) > 0 and int(mapqs[i]) <= 30):
						qual_dico['1-30'] += 1  # Rajouter 1 pour compter le nombre de read
					elif(int(mapqs[i]) > 30 and int(mapqs[i]) <= 60):
						qual_dico['31-60'] += 1  # Rajouter 1 pour compter le nombre de read
					else:
						qual_dico['61 -'] += 1  # Rajouter 1 pour compter le nombre de read
				else:
					qual_dico['0'] += 1  # Rajouter 1 pour compter le nombre de read
		
		qual_dico_sequence[sqname] = qual_dico  # Stock
		print("\n_________________ \n Tranche de qualité :\n")
		
		for upKey, upValue in qual_dico_sequence.items():
			data1 = [(f"{upKey} |", f"{key} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |") for key, val in upValue.items()]
			t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |"," Qualité |", " Quantité |", " 100% |"])
			t_data = pd.concat([t_data1], axis=1)
			print(t_data)
			

####  #####
def cigar_analysis():
	

	print("\n_________________ \n Quantité des cigars :\n")
	
	longueur_théorique = 0
	longueur_constaté = 0
	
	cigar_dico_description = {
		"M": "",
		"I": "",
		"D": "",
		"N": "",
		"S": "",
		"H": "",
		"P": "",
		"=": "",
		"X": ""
	}
	
	cigar_dico = {
		"M" : 0,
		"I" : 0,
		"D" : 0,
		"N" : 0,
		"S" : 0,
		"H" : 0,
		"P" : 0,
		"=" : 0,
		"X" : 0
	}
	
	for sqname in sequenceRefName:
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if (rname[i] == sqname):        # verifier si le read est mappé sur la sequence de reference
				if(re.search(r"(\*)", cigars[i])):
					continue
				else:
					M = re.search(r"(\d*)M", cigars[i])
					if M:
						cigar_dico["M"] += int(M.group(1))
						longueur_constaté += int(M.group(1))
					
					I = re.search(r"(\d*)I", cigars[i])
					if I:
						cigar_dico["I"] += int(I.group(1))
						longueur_constaté += int(I.group(1))
						
					D = re.search(r"(\d*)D", cigars[i])
					if D:
						cigar_dico["D"] += int(D.group(1))
						
					N = re.search(r"(\d*)N", cigars[i])
					if N:
						cigar_dico["N"] += int(N.group(1))
						
					S = re.search(r"(\d*)S", cigars[i])
					if S:
						cigar_dico["S"] += int(S.group(1))
						longueur_constaté += int(S.group(1))
						
					H = re.search(r"(\d*)H", cigars[i])
					if H:
						cigar_dico["H"] += int(H.group(1))
					
					P = re.search(r"(\d*)P", cigars[i])
					if P:
						cigar_dico["P"] += int(P.group(1))
						
					E = re.search(r"(\d*)=", cigars[i])
					if E:
						cigar_dico["="] += int(E.group(1))
						longueur_constaté += int(E.group(1))

					X = re.search(r"(\d*)X", cigars[i])
					if X:
						cigar_dico["X"] += int(X.group(1))
						longueur_constaté += int(X.group(1))
					
						
	
	data1 = [(f"{cle} |", f"{val} |", f"{round(val/longueur_constaté, 5) * 100} |") for cle, val in cigar_dico.items()]  #
	t_data1 = pd.DataFrame(data1, columns=[" Code |", " Valeur |", " 100% |"])
	t_data = pd.concat([t_data1], axis=1)
	print(t_data)
	

######## Start ########
if len(sys.argv) == 1:
	print(
		"Vous n'avez pas rentré d'argument, vous pouvez répondre au question suivante ou utilisez -h pour plus de détails ou automatiser le procéssus \n")
	
elif "-h" in sys.argv or "--help" in sys.argv:
	print("Pour faire fonctionner ce script :")
	print("SamMapScript.py [file...]")
	print("SamMapScript.py [file...] [Quality Min]")
	print("SamMapScript.py [file...] [Quality Min] [Alignement]\n")
	
	print("Quality Min : Number between 0 - 255\n")
	print("Alignement : \tPrendre seulement les séquences :")
	print(" -ALL \t Tout prendre")
	print(" -M \t Seulement les mappés")
	print(" -SM \t Seulement les mappés et semimappé\n")
	
	sys.exit()


#### Open file and take quality ####

sam_file_path = findFile()
QualityMin = askQuality()
# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals
sequenceRefName, flags, rname, mapqs, cigars, totalNumberOfRead = sam_reading(sam_file_path)
binary_flags = flags_to_binary(flags)  # Conversion des flags en binaire sur 12 bits

#### Number of reads ####
analysis_of_mapped_reads(binary_flags)
count_flag_number()
#pair_Analysis()
quality_analysis()
cigar_analysis()

print("\n")
