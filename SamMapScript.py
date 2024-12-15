#!/usr/bin/python3
# -*- coding : utf-8 -*-

#licence auteur et version

# trop de dur, interpretation plus clair et informative, # finir selection alignement
# demander granularité pour analyse qualité et limite qualité #desactivé cigar mais regarder les indels
#verification partie analyse flags entre pair (8/16 32/64 128/256) et présenté si erreur (128/256 = 50 ou incoérence entre valeur 218/256 avec

__authors__ = "Loik Galtier"
__Version__ = "1.0.0"
__Licence__ = ("This program is free software: you can redistribute it and/or modify"
               " it under the terms of the GNU General Public License as published by)"
               " the Free Software Foundation, either version 3 of the License, or "
               "(at your option) any later version. "
               "This program is distributed in the hope that it will be useful,"
               " but WITHOUT ANY WARRANTY; without even the implied warranty of"
               " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
               " GNU General Public License for more details. "
               "You should have received a copy of the GNU General Public License"
               " along with this program. If not, see <https://www.gnu.org/licenses/>.")

import os
import re
import csv
import sys

import pandas as pd

####### Help #######
def show_help():
	print("Pour faire fonctionner ce script :")
	print("SamMapScript \t posera quelque question pour lancer l'analyse")
	print("SamMapScript -h \t presentera cette aide")
	print("SamMapScript.py [file] \t lancera l'analyse du fichier ")
	print("SamMapScript.py [file] [option...] \t lancera l'analyse du fichier avec les options adapters: \n \n")
	print("--[Option]-- \n")
	print("qual[int] \t définie la qualité minimale des reads a analysé à 'int'")
	print("step_qual[int] \t définie le pas pour les tranches de qualité de 'int'")
	print("specific_sequence[seq] \t permet l'anayse seulement sur les séquences de référence seq, pour prendre plusieurs séquences : specific_sequence[seq1,seq2] (attention au espaces) \n")
	
	print("analyse_map")
	print("count_flag")
	print("analyse_quality")
	print("analyse_cigar")

######## Setting the file ########
# check the selected file exist and have the .sam
def find_file():
	if len(sys.argv) == 1:
		no_file_found = True
		while no_file_found:
			name = input("Quel est le chemin d'accés au fichier sam : ?")  # Ask the path
			if os.path.exists(name) and re.search(r".sam", name):  # if the path exist AND the .sam exist
				no_file_found = False  # quit the while
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
def ask_quality():
	if len(sys.argv) == 1:
		no_answer = True
		while no_answer:
			need_quality_min = input(
				"Voulez vous appliqué un filtre de qualité minimal pour la suite du traitement : (Y/N)").upper()  # Ask the path
			if need_quality_min == "Y":
				quality_min = input("Quel valeur minimal voulais vous appliquer : [0-255] \n")
				if quality_min.isdigit():
					print(
						"Toute les prochaines opérations se feront sur les reads ayant un qualité supérieur à : " + quality_min + "\n")
					no_answer = False
					return int(quality_min)
				else:
					print("Cela n'est pas une valeur correcte \n")
			elif need_quality_min == "N":
				no_answer = False
				print("\n")
				return 0
			else:
				print("Valeur non correcte \n")
	else:
		return 0

######## Extract part and use header ########
def sam_reading(sam_file_path):
	# Ouvrir le fichier SAM
	with open(sam_file_path, "r") as sam_file:
		
		# Parcourir chaque ligne du fichier SAM
		sequence_names = []
		
		flags = []
		rnames = []
		mapqs = []
		cigars = []
		
		number_of_read_total = 0
		
		# affichage info header
		hd_table = {}
		sq_table = {}
		pg_table = {}
		
		# Extraire les informations du header
		for line in sam_file:
			
			if line.startswith("@"):
				if line.startswith("@HD"):
					vN = re.search(r"VN:([^\t|\n]+)", line)
					sO = re.search(r"SO:([^\t|\n]+)", line)
					if vN:
						hd_table[vN.group(1)] = ""
						if sO:
							hd_table[vN.group(1)]= [sO.group(1)]
				
				elif line.startswith("@PG"):
					idnamepg = re.search(r"ID:([^\t|\n]+)", line)
					versionpg = re.search(r"VN:([^\t|\n]+)", line)
					if idnamepg:
						pg_table[idnamepg.group(1)] = ""
						if versionpg:
							pg_table[idnamepg.group(1)] = versionpg.group(1)
						
				elif line.startswith("@SQ"):
					sN = re.search(r"SN:([^\t|\n]+)", line)
					lN = re.search(r"LN:([^\t|\n]+)", line)
					if sN and lN:
						sequence_names.append(sN.group(1))
						sq_table[sN.group(1)] = lN.group(1)
						
				elif line.startswith("@RG"):
					idgroup = re.search(r"ID:([^\t|\n]+)", line)
					if idgroup:
						print("Un groupe de lecture a était formé, son ID unique est : " + idgroup.group(1) + "\n")
		sam_file.seek(0)  # pour arreter la boucle
		
		print("_________________ \n Entête :\n")
		
		if hd_table:
			print("Les données de ce fichier sont :")
			data = [(f"{cle} |", f"{val} |") for cle, val in hd_table.items()]
			t_data = pd.DataFrame(data, columns=["Format Version |", "Ordre |"])
			print(t_data, "\n")
		
		if sq_table:
			print("Des séquences références en était trouvé :")
			data = [(f"{cle} |", f"{val} |") for cle, val in sq_table.items()]
			t_data = pd.DataFrame(data, columns=["séquence name |", "Length |"])
			print(t_data, "\n")
		
		if pg_table:
			print("Des programmes en était utilisé :")
			
			data = [(f"{cle} |", f"{val} |") for cle, val in pg_table.items()]
			t_data = pd.DataFrame(data, columns=["Programme id |", "Version |"])
			print(t_data, "\n")
		# Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation
		sam_reader = csv.reader(sam_file, delimiter='\t')
		for row in sam_reader:  # Pour chaque ligne de mon fichier SAM
			if row[0].startswith("@"):  # Ignorer les lignes d'en-tête qui commencent par '@'
				continue
			if quality_min <= int(row[4]) < 255:  # si la ligne a une qualité supérieur à la qualité demandé et inférieur a 255
				if (len(specific_sequence) > 0 and row[2] in specific_sequence) or len(specific_sequence) == 0:
					# Accéder aux informations de chaque colonne
					flag = int(row[1])  # Flag
					flags.append(flag)  # rajouter le flag de chaque ligne dans la liste
					
					rname = row[2]  # Nom de la séquence de référence
					rnames.append(rname)
					
					mapq = row[4]
					mapqs.append(mapq)
					
					cigar = row[5]  # Chaîne CIGAR
					cigars.append(cigar)
					
					number_of_read_total += 1
		# print(f"QNAME: {qname}, FLAG: {flag}, RNAME: {rname}, POS: {pos}, CIGAR: {cigar}, SEQ: {seq}")
	if len(specific_sequence) > 0 :
		return specific_sequence, flags, rnames, mapqs, cigars, number_of_read_total
	else:
		return sequence_names, flags, rnames,mapqs , cigars, number_of_read_total

######## Convert to binary ########
#fonction flags_to_binary pour convertir le flag en binaire
def flags_to_binary(flags):
	for i in range(len(flags)):  # Boucle qui va parcourir de 0 à la taille des flags
		flags[i] = bin(flags[i])  # Convertir en binaire
		flags[i] = flags[i][2:].zfill(12)  # Éliminer les 2 premiers caractères du flag (0 et b) & remplir de 0 sur la gauche jusqu'a 12 #si plus de 12, passera avec le plus
	# print(flags[i])
	return flags  # Retour des flags en binaire

######## Question 1 ########
# fonction analysis_of_mapped_reads, prend en paramètre les flags en binaire
def analysis_of_mapped_reads(binary_flags):
	nb_mapped_table = {}
	nb_semmapped_table = {}
	nb_unmapped = 0

	# mapped and semimapped part
	for sqname in sequenceRefName:
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if rname[i] == sqname:  # verifier si le read est mappé sur la sequence de reference
				flag = binary_flags[i]  # Mettre un flag en binaire dans la variable flag
				if flag[-2] == "1":  # le bit 2 code pour l'info "aligné". Si "1" = aligné.
					if not re.fullmatch(r"(\d*M+)+", cigars[i]):  # si 100M, alors il n'y as pas de problème d'alignement, en ignore # regex integer M seulement pour
						if sqname in nb_semmapped_table:  # Si sqname existe dans la table des semi mappé
							nb_semmapped_table[sqname] += 1  # Rajouter 1 pour compter le nombre de read
						else:
							nb_semmapped_table[sqname] = 1  # Sinon initialise
					else:
						if sqname in nb_mapped_table:
							nb_mapped_table[sqname] += 1  # Rajouter 1 pour compter le nombre de read
						else:
							nb_mapped_table[sqname] = 1  # Sinon initialisé
						

	
	# unmapped part
	for i in range(len(binary_flags)):
		flag = binary_flags[i]  # Mettre un flag en binaire dans la variable flag
		if flag[-3] == "1" or rname[i] == "*" or rname[i] not in sequenceRefName:  # le bit 3 code pour l'info "read non mappé". Si "1" = non mappé
			nb_unmapped += 1  # Rajouter 1 pour compter le nombre de read
	
	print("_________________ \n Statistique d'Alignement :\n")
	print("Parmis", totalNumberOfRead, "reads, les reads sont :")
	
	if totalNumberOfRead != 0:
		data1 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" Mappé |") for cle, val in nb_mapped_table.items()]  #
		t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |", "Quantité |", " % |", " Aligné ?|"])
		t_data = pd.concat([t_data1], axis=0)
		
		data2 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" Semi-mappé |") for cle, val in nb_semmapped_table.items()]  #
		t_data2 = pd.DataFrame(data2, columns=[" Sequence reference |", "Quantité |", " % |", f" Aligné ?|"])
		t_data = pd.concat([t_data1, t_data2], axis=0)
		
		if nb_unmapped > 0:
			data3 = [(" Non aligné |", f"{nb_unmapped} |", f" {round(nb_unmapped / totalNumberOfRead, 4) * 100} |", f" non mappé |")]  #
			t_data3 = pd.DataFrame(data3, columns=[" Sequence reference |", "Quantité |", " % |", f" Aligné ?|"])
			if nb_semmapped_table:
				t_data = pd.concat([t_data1, t_data2, t_data3], axis=0)
			elif nb_mapped_table:
				t_data = pd.concat([t_data1, t_data3], axis=0)
		print(t_data.to_string(index=False), "\n")

def count_flag():
	
	flags_counts = {}
	
	for i in range(len(binary_flags)):
		# Mettre un flag en binaire dans la variable flag
		flag = binary_flags[i]
		
		for j in range(len(flag)):
			if flag[-j - 1] == "1":
				if (str(2 ** j) + " bits") not in flags_counts:
					flags_counts[str(2 ** j) + " bits"] = 1
				else:
					flags_counts[str(2 ** j) + " bits"] += 1
			else:
				if (str(2 ** j) + " bits") not in flags_counts:
					flags_counts[str(2 ** j) + " bits"] = 0
	return flags_counts
	
def analysis_flag():
	print("_________________ \n analyse des flags :\n")
	for i in range(len(flags_counts)):
		if i == 0:
			print(f"Il y as {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) des reads sont appariés et donc {100 - round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}% sont indépendants. \n")
		if i == 2:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i+1)) + " bits"]:
				print(f"{flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) reads et leurs reads associés ne sont pas alignés sur une séquence. \n")
			else:
				print(f"Des reads sont aligné mais pas leur complémentaire. \n Il y as {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) reads non alignés et {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (1+i)) + ' bits'] / totalNumberOfRead * 100, 2)}%) reads associés non aligné. \n")
		
		if i == 4:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i+1)) + " bits"]:
				print("Tous les reads et leurs associés sont dans la même orientation.")
				print(f"Il y a {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) dans le sens inverse complémentaire et autant dans le sens non inverse.\n")
			else:
				print("Il y a plus de reads dans un sens que dans l'autre !")
				print(f"Il y a {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) dans le sens inverse complémentaire mais {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] / totalNumberOfRead * 100, 2)}%) dans le sens non inverse.\n")
		if i == 6:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i + 1)) + " bits"]:
				print("Tous les reads sont en pair.")
				print(f"Il y a {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) en premier segment et autant comme seconds.\n")
			else:
				print("Il y a un problème dans le fichier.")
				print(f"Il y a {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) en premier segment et {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] / totalNumberOfRead * 100, 2)}%) en seconds.\n")


####### Optionnel : se  fait que si parametre "detail_flag" est appellé, certaine valeur son en dur ici, mais cette def n'est pas compris dans l'utilisation du projet initial ########
def count_flag_number():

	print("_________________ \n Quantité de flags :\n")
	
	pd.set_option('display.max_colwidth', 40)
	
	if totalNumberOfRead != 0:
		data1 = [(f"{cle} |",f"{val} |", f" {round(val / totalNumberOfRead * 100, 2)} |") for cle, val in flags_counts.items()]  #
		t_data1 = pd.DataFrame(data1, columns=["Bits |","Quantité |", " % |"])
		
		t_data = pd.concat([t_data1], axis=1)
		
		print(t_data.to_string(index=False))
		
#### ####
def quality_analysis():
	
	qual_dico_sequence={}
	
	for sqname in sequenceRefName:
		qual_dico = {
			'0 <= quality':0,
		}
		
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if rname[i] == sqname:  # verifier si le read est mappé sur la sequence de reference
				if mapqs[i].isdigit() and mapqs[i] != '0':  # verifier si mapQ est bien un chiffre
					j = 0
					while j < 156:
						tmp = str(j) + " < quality <= " + str(j + quality_step)
						if j < int(mapqs[i]) <= j + quality_step:
							if tmp in qual_dico:
								qual_dico[tmp] += 1  # Rajouter 1 pour compter le nombre de read
							else:
								qual_dico[tmp] = 1
							j = 156 #permet d'arreter la boucle
						else:
							if tmp not in qual_dico:
								qual_dico[tmp] = 0  #permet d'initialiser toute les valeurs, jusqu'a la plus grande
						j += quality_step
				else:
					qual_dico['0 <= quality'] += 1  # Rajouter 1 pour compter le nombre de read
		
		qual_dico_sequence[sqname] = qual_dico  # Stock
		

		print(f"\n_________________ \n Tranche de qualité pour la séquence : {sqname}, avec une qualité minimal de {quality_min}, et par pas de {quality_step}\n")
		
		data1 = [(f" |", f" |", f" |", f"  |")]  # remplacer par %
		t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |", " Qualité |", " Quantité |", " % |"])
		t_data = pd.concat([t_data1], axis=1)
		
		for upKey, upValue in qual_dico_sequence.items():
			if totalNumberOfRead != 0:
				data1 = [(f"{upKey} |", f"{key} |", f"{val} |", f" {round(val / totalNumberOfRead * 100, 4)} |") for key, val in upValue.items()] #remplacer par %
				t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |"," Qualité |", " Quantité |", " % |"])
				t_data = pd.concat([t_data1], axis=0)
				#t_data.reset_index(drop=True, inplace=True)
		print(t_data)
		
####  #####
def cigar_analysis():
	# pas obligatoire pas dur

	print("\n_________________ \n Quantité des cigars :\n")
	
	longueur_theorique = 0
	longueur_constate = 0
	
	
	cigar_dico = {}
	
	for sqname in sequenceRefName:
		for i in range(len(binary_flags)):  # Boucle qui parcours de 0 à la quantité de reads
			if rname[i] == sqname:        # verifier si le read est mappé sur la sequence de reference
				if re.search(r"(\*)", cigars[i]):
					continue
				else:
					tmps = re.findall(fr"(\d+)([A-Za-z])", cigars[i])
					for tmp in tmps:
						if tmp[0]:
							number = int(tmp[0])
						else:
							number= 1  # Par défaut, 1 si le nombre est absent
						letter = tmp[1]
						
						if letter in cigar_dico:
							cigar_dico[letter] += number  # Rajouter nombre trouverplus tot pour compter le nombre de read
						else:
							cigar_dico[letter] = number
						longueur_theorique += number
						
	if longueur_theorique != 0:
		data1 = [(f"{cle} |", f"{val} |", f"{round(val / longueur_theorique, 5) * 100} |") for cle, val in cigar_dico.items()]  #
		t_data1 = pd.DataFrame(data1, columns=[" Code |", " Valeur |", " % |"])
		t_data = pd.concat([t_data1], axis=1)
		print(t_data)
	
####### Start ########

only_step = False
quality_set = False

if len(sys.argv) == 1:
	print("Vous n'avez pas rentré d'argument, vous pouvez répondre au question suivante ou utilisez -h pour plus de détails ou automatiser le procéssus \n")

for arg in enumerate(sys.argv[1:]):
	if arg[1] in ("-h","--help"): #changer pour liste d'action pour lancer la qual puis le find pui le read puis binary puis les actions contextuels ave l'affichage de l'entete en plus  b
		show_help()

#### Open file and take quality ####

sam_file_path = find_file()
quality_min = 0
quality_step = 10
specific_sequence = []
size_flag = 12

flags_counts = {}



for arg in enumerate(sys.argv[1:]):
	if arg[1].startswith("qual["): #changer pour liste d'action pour lancer la qual puis le find pui le read puis binary puis les actions contextuels ave l'affichage de l'entete en plus  b
		qual = re.search(r"qual\[(\d*)\]",arg[1])
		quality_min = int(qual.group(1))
		print(f"La qualité minimale sera de {quality_min}")
		quality_set = True
	if arg[1].startswith("step_qual["):
		qualstep = re.search(r"step_qual\[(\d*)\]",arg[1])
		quality_step = int(qualstep.group(1))
		print(f"L'analyse de la répartition des qualités se fera par pas de {quality_step}")
	if arg[1].startswith("specific_sequence["):
		seq = re.search(r"specific_sequence\[(.*)\]",arg[1])
		if seq:
			specific_sequence = seq.group(1).split(',')
			print(f"Seul les reads alligné sur la séquence : {specific_sequence} seront utilisé.")

if not quality_set:
	ask_quality()

# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals
sequenceRefName, flags, rname, mapqs, cigars, totalNumberOfRead = sam_reading(sam_file_path)
binary_flags = flags_to_binary(flags)  # Conversion des flags en binaire sur 12 bits

for arg in enumerate(sys.argv[1:]):
	if arg[1] == "analyse_map":
		analysis_of_mapped_reads(binary_flags)
		only_step = True
	if arg[1] == "count_flag" or arg[1] == "analyse_flag":
		flags_counts = count_flag()
	if arg[1] == "analyse_flag":
		analysis_flag()
		only_step = True
	if arg[1] == "count_flag":
		count_flag_number()
		only_step = True
	if arg[1] == "analyse_quality":
		quality_analysis()
		only_step = True
	if arg[1] == "analyse_cigar":
		cigar_analysis()
		only_step = True
		
#### Number of reads ####
if not only_step:
	flags_counts = count_flag()
	analysis_of_mapped_reads(binary_flags)
	analysis_flag()
	#count_flag_number()
	#pair_Analysis()
	quality_analysis()
	#cigar_analysis()

print("\n")
