#Programme pour passer d'un format de positions de jmol (.xyz) a un un
#format de siesta


nom_in = raw_input("Nom du fichier format siesta?\n")
nom_out = raw_input("Nom du fichier de sortie?\n")
nb_types = input("Nombre de types d'atomes?\n")

types = []
for k in range(nb_types):
    tp = raw_input("Symbole du %ie type?\n"%(k+1))
    types.append(tp)
def lecture_listes_reels(nomFichier):
#On ouvre le fichier
   fichier_lu = open(nomFichier,'r')
#On le lit
   lignes = fichier_lu.readlines()
#On va chercher le nombre de colonnes
   longueur = len(lignes[2].split())
#On va chercher le nombre de lignes
   nbLignes = len(lignes)
#On cree une liste qui contiendra une liste par colonne
   liste = []
#On initialise les sous-listes
   for i in range(longueur):
       liste.append([])
    #On separe chaque ligne en ses elements
   for i in range(0,nbLignes):
       donnees = lignes[i].split()
       #Et on ajoute a chaque sous-liste l'element adequoit
       for j in range(longueur):
           liste[j].append((donnees[j]))
           #On ferme le fichier
   fichier_lu.close()
           #En sortie: la liste de liste
   return liste

#Methode d'ecriture
def ecriture_listes_reels(nomFichier,liste):
    #On ouvre le fichier
    sortie = open(nomFichier,'w')
    #On va chercher le nombre de lignes a ecrire
    nbLignes = len(liste[1])
    #Et le nombre de colonnes
    nbColonnes = len(liste)
    #On cree une liste qui va contenir lignes qui seront ecrites
    lignes = range(nbLignes+2)
    lignes[0] = "%i\n"%nbLignes
    lignes[1] = "\n"   
    for i in range(nbLignes):
        lignes[i+2] = "" #On commence par creer une chaine
        for k in range(nb_types):
            if liste[3][i] == "%i"%(k+1):
                lignes[i+2] = lignes[i+2] + "%i"%(k+1)
        for j in range(0,3):
            lignes[i+2] = lignes[i+2] + " %s"%liste[j][i]	#On y ajoute les points
        lignes[i+2] = lignes[i+2] + "\n" #Puis un retour de chariot
        lignes[i+2] = lignes[i+2].lstrip() #Et on enleve l'espace du debut.
            #On ecrit ca dans le fichier
    sortie.writelines(lignes)
    sortie.close()

liste = lecture_listes_reels(nom_in)
ecriture_listes_reels(nom_out,liste)
