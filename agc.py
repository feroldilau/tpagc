#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "FEROLDI Laurent et MEUNIER Valentin"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["FEROLDI Laurent et MEUNIER Valentin"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "FEROLDI Laurent et MEUNIER Valentin"
__email__ = " feroldilau@eisti.eu et meunierval@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()



# Dé-duplication en séquence “complète” 

#prend deux arguments correspondant au fichier fasta.gz 
# et à la longueur minimale des séquences 
# et retourne un générateur de séquences de longueur l >= minseqlen: yield sequence.

def read_fasta(amplicon_file, minseqlen):
       list1 = []
       n = 0
    with gzip.open(amplicon_file, "rt") as file1:
        seq = ""
        for i in file1:
            line1 = i.strip("\n")
            if line1[:1]=="T" or line1[:1]=="G" or line1[:1]=="C" or line1[:1]=="A":
                seq = seq + line1
            else:
                list1.append(seq.strip("\n"))
                seq = ""
        for j in list1:
            n = len(j)
            if n>=minseqlen:
                yield j

#Prend trois arguments correspondant au fichier fasta,  
# la longueur minimale des séquences et leur comptage minimum.
#  Elle fait appel au générateur fourni par read_fasta et 
# retourne un générateur des séquences uniques
#  ayant une occurrence O>=mincount ainsi que leur occurrence. 
# Les séquences seront retournées par ordre décroissant d’occurrence: yield [sequence, count]

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dictampli = {} 
    listseq = read_fasta(amplicon_file, minseqlen)
    for i in listseq:
        if i in dictampli:
            dictampli[i] += 1
        else :
            dictampli[i] = 1
    for j, k in sorted(dictampli.items(), key=lambda item: item[1], reverse = True): 
        if k >= mincount:
            yield [j, k]


# Recherche de séquences chimériques par approche “de novo”

#prend une séquence et une longueur de segment l: chunk_size et 
# retourne une liste de sous-séquences de taille l non chevauchantes.
#  A minima 4 segments doivent être obtenus par séquence.

def get_chunks(sequence, chunk_size):
    liste1 = []
    for i in range(0, len(sequence), chunk_size):
        if i+chunk_size<=len(sequence):
            liste1.append(sequence[i:i+chunk_size])
    return liste1

#
def get_unique(ids):
    return {}.fromkeys(ids).keys()

#
def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

#prend une séquence et une longueur de k-mer et 
# retourne un générateur 
# de tous les mots de longueur k présents dans cette séquence,
#   yield kmer

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]


#prend un alignement (sous forme de liste) et
#  calcule le pourcentage d’identité entre les deux séquences

def get_identity(alignment_list):
    a = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            a += 1
    a = ( a/len(alignment_list[0])) * 100
    return a

#prend un dictionnaire ayant pour clé un dictionnaire de kmer (vide au départ) et
#  pour valeur une liste d’identifiant des séquences dont ils proviennent.
#  Et retourne un dictionnaire de kmer 
# contenant les kmers uniques présents dans chaque séquence 
# pour une longueur de donnée de kmer.

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    liste1 = list(cut_kmer(sequence, kmer_size))
    for i in liste1:
        if i in kmer_dict:
            kmer_dict["kmer"].append(id_seq)
        else:
            kmer_dict["kmer"] = [id_seq]
    return kmer_dict


#prend un dictionnaire ayant pour clé un index de kmer et
#  pour valeur une liste d’identifiant des séquences dont ils proviennent,
#  une séquence et une longueur de kmer: kmer_size. 
#Vous utiliserez Counter de la librairie Collections et
#  sa fonction most_common pour identifier
#  les 8 séquences les plus similaires à notre séquence entrée.

def search_mates(kmer_dict, sequence, kmer_size):
    liste1 = []
    liste2 = []
    liste3 = list(cut_kmer(sequence, kmer_size))
    for i in liste3:
        if i in kmer_dict:
            for j in kmer_dict[i]:
                liste1.append(j)
    for k in Counter(liste1).most_common(8):
        liste2.append(k[0])
    return liste2


#prend une matrice donnant par segment
#  le taux d’identité entre la séquence candidate et deux séquences parentes et 
# retourne un booléen 
# indiquant si la séquence candidate est une chimère (True) ou ne l’est pas (False). 

def detect_chimera(perc_identity_matrix):
    a = False
    b = False
    liste1 = []
    for i in perc_identity_matrix:
        liste1.append(statistics.stdev(i))
        if i[0] > i[1]:
            a = True
        if i[0] < i[1]:
            b = True
        if statistics.mean(liste1) > 5 and a and b:
            return True
        return False


#Fait appel au générateur fourni par dereplication_fulllength 
# et retourne un générateur des séquences non chimérique 
# au format: yield [sequence, count]

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    dict1 = {}
    liste1 = []
    a = False
    idk = 0
    for i, seq in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):
        if i%100 == 0:
            print(i)
        Lchunk = get_chunks(seq[0], chunk_size)
        list2 = []
        for chunk in Lchunk:
            list2.append(search_mates(dict1, chunk, kmer_size))
        Lcommon = []
        for j in range(len(list2)):
            Lcommon = common(Lcommon, list2[j])
        if len(Lcommon) > 1:
            for c in Lcommon[0:2]:
                chunkref = get_chunks(liste1[c], chunk_size)
                IMa = [[]*4]
                for k in range(len(Lchunk)):
                    ali = nw.global_align(Lchunk[k], chunkref[k],
                        gap_open=-1, gap_extend=-1, matrix="agc/MATCH")
                    IMa[k].append(get_identity(ali))
            a = detect_chimera(IMa)
        if not a:
            dict1 = get_unique_kmer(dict1, seq[0], idk, kmer_size)
            liste1.append(seq[0])
            idk += 1
            yield seq



#Pour chaque séquence non-chimérique,
#  vous devez s’il existe une séquence ayant une abondance plus importante
#   qui lui soit similaire à >97%.
#Si oui, cette séquence n’est pas une OTU et sera ignorée.
#Si non, cette séquence est une OTU et rejoindra notre liste d’OTU finale.
#fait appel à chimera_removal et
#  réalise également des mesures d’identité à l’aide de get_identity.
#  Elle retourne une liste d’OTU, 
# cette liste indiquera pour chaque séquence son occurrence (count)

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    liste1 = []
    liste2 = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    for i in range(len(liste2)):
        a = True
        b, c = liste2[i]
        for j in range(i+1, len(liste2)):
            b2, c2 = liste2[j]
            if b != b2 and get_identity([b, b2])>97 and c2>c:
                a = False
                break
        if a:
            liste1.append([b, c])
    return liste1


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

#Prend une liste d’OTU et le chemin vers un fichier de sortie 
# et affiche les OTU au format:
#>OTU_{numéro partant de 1} occurrence:{nombre d’occurrence à la déréplication}
#{séquence au format fasta}

def write_OTU(otu_list, output_file):
    a = 0
    with open(output_file, "w") as file1:
        for i, j in otu_list:
            file1.write(">OTU_"+str(a+1)+" occurrence:"+ str(j) + "\n")
            file1.write(fill(str(i)))
            file1.write("\n")
            a += 1


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    output_file = args.output_file
    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    write_OTU(OTU_list, output_file)

if __name__ == '__main__':
    main()