#!/usr/bin/python3
#-*- coding : utf-8 -*-
import os
import re
import csv
import sys


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
            print(f"Qual : {qual}")
    return flags, quals, coverage






### Start ###
sam_file_path = findFile()
flag_size = 12
# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals
flags, quals, coverage = sam_reading(sam_file_path)