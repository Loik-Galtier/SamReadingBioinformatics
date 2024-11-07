#!/usr/bin/python3
#-*- coding : utf-8 -*-
import os
import re
import csv
import sys

#### Setting the file####
#check the selected file exist and have the .sam
def findFile ():
    noFileFound = True
    while noFileFound :
        name = input("Path and Name of the file : ?") #Ask the path
        if os.path.exists(name) and re.search(r".sam", name): # if the path exist AND the .sam exist
            noFileFound = False #quit the while
            print("There is a file") # continue -> Ask what question to respond (1, 2, 3, 4)
            return name #return the file open
        else:
            print("Don't have a file") # Ask again for a file.


def sam_reading(sam_file_path):
    # Ouvrir le fichier SAM
    with open(sam_file_path, "r") as sam_file:
        # Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation car le fichier SAM est séparé par des tabulations
        sam_reader = csv.reader(sam_file, delimiter='\t')

        # Parcourir chaque ligne du fichier SAM
        flags = []
        quals = []
        coverage = {}
        cigars = []
        # Pour chaque ligne de mon fichier SAM
        for row in sam_reader:
            # Ignorer les lignes d'en-tête qui commencent par '@'
            if row[0].startswith("@"):
                continue

            # Accéder aux informations de chaque colonne
            qname = row[0]  # Nom du read
            flag = int(row[1])  # Flag
            flags.append(flag)  # rajouter le flag de chaque ligne dans la liste
            rname = row[2]  # Nom de la séquence de référence
            start_pos = int(row[3])  # Position de début de l'alignement
            mapq = int(row[4])  # Qualité de l'alignementt
            cigar = row[5]  # Chaîne CIGAR
            cigars.append(cigar)
            rnext = row[6]  # Référence du read suivant dans le cas des paires
            pnext = int(row[7])  # Position du read suivant dans le cas des paires
            tlen = int(row[8])  # Longueur du fragment pour les paires
            seq = row[9]  # Séquence de l'ADN
            seq_length = len(seq)
            if len(row) > 10: #voir si existe collone 10
                qual = row[10]  # Qualité de chaque base dans la séquence
                quals.append(qual)  # rajouter la qual de chaque ligne dans la liste

            # pour la question 3 localisation des reads sur la séquence

            for pos in range(start_pos, start_pos + seq_length):
                if pos in coverage:
                    coverage[pos] += 1
                else:
                    coverage[pos] = 1

            # Afficher les informations du read
            # print(f"QNAME: {qname}, FLAG: {flag}, RNAME: {rname}, POS: {pos}, CIGAR: {cigar}, SEQ: {seq}")
    return flags, quals, coverage, cigars


# Je définie une fonction nommée flags_to_binary pour convertir le flag en binaire parceque le flag contient les infos en bit
def flags_to_binary(flag_size, flags):
    # Boucle qui va parcourir de 0 à la taille des flags-1
    for i in range(len(flags)):
        flags[i]=bin(flags[i])
        # Éliminer les 2 premiers caractères du flag (0 et b) & les remplir de 0 sur la gauche -> atteindre la taille souhaitée
        flags[i]=flags[i][2:].zfill(flag_size)
        # Faire retourner une liste de flag transformée en binaire
        #print(flags[i])
    return flags

#### Question 1 ####
# Je définie une fonction nommée number_of_mapped_reads qui prend en paramètre la taille des flags et les flags en binaire
def number_of_mapped_reads (flag_size, binary_flags):
    nbr = 0
    for i in range(len(binary_flags)):
        # Mettre le flag en binaire de la ligne en question (i) dans la variable flag
        flag=binary_flags[i]
        # le bit 3 code pour l'info "read mappé ou pas", donc on soustrait le chiffre 3 de la taille du flag. Si "0" = mappé, sans prendre en compte si il est mal mappé
        if (flag[-3] == "0"):
            # Rajouter 1 pour compter le nombre de read
            nbr = nbr + 1
    print ("J'ai", nbr, "read mappés")


def number_of_unmapped_reads(flag_size, binary_flags):
    nbr = 0
    for i in range(len(binary_flags)):
        # Mettre le flag en binaire de la ligne en question (i) dans la variable flag
        flag=binary_flags[i]
        # le bit 3 code pour l'info "read mappé ou pas", donc on soustrait le chiffre 3 de la taille du flag. Si "1" = non mappé
        if flag[-3] == "1":
        # Rajouter 1 pour compter le nombre de read
            nbr = nbr + 1
    print ("J'ai", nbr, "read non mappés")

def number_of_semimapped_reads(flag_size, binary_flags, cigars):
    nbr = 0
    for i in range(len(binary_flags)):
        # Mettre le flag en binaire de la ligne en question (i) dans la variable flag
        flag=binary_flags[i]
        # le bit 2 code pour l'info "alignement", donc on soustrait le chiffre 2 de la taille du flag. Si "1" = alors il y as une information
        if flag[-2] == "1":
            if cigars[i] != "100M": #si 100M, alors il n'y as pas de problème
                # Rajouter 1 pour compter le nombre de read
                nbr = nbr + 1
    print ("J'ai", nbr, "read semi mappés")


### Start ###
sam_file_path = findFile()
# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals
flags, quals, coverage, cigars = sam_reading(sam_file_path)

flag_size = 12
binary_flags = flags_to_binary(flag_size ,flags) #Conversion des flags en binaire sur 12 bits
number_of_mapped_reads(flag_size, binary_flags)
number_of_unmapped_reads(flag_size, binary_flags)
number_of_semimapped_reads(flag_size, binary_flags, cigars)