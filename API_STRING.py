#!/usr/bin/env python3.8
##################################################################
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
from Bio import Entrez
import re 
import requests ## python -m pip install requests
##################################################################

#####rs_ID = input("Inserte rs: ")
file1 = open(input("Lista de SNP: "), 'r')
file2 = open('Genes_afectados.txt', 'w+')
file2.write("Rs\tTipo de Polimorfismo\tGenes afectados\tGene_Id\tGene_name\tGene_Id_2\tGene_name_2\n")
Lines = file1.readlines()
Lines = [w.replace('\n', '') for w in Lines]
print(str(len(Lines)) + " Variantes geneticas a buscar")
count = 0
for line in Lines:
    ### a rs10953105 tiene 2 genes uno codificante y el otro no codificante
    count+=1
    print(str(count/len(Lines)*100)+"% de Variantes buscadas ("+ str(count) +" de " + str(len(Lines)) + ")")
    Entrez.email = "adolfo.rh@postqyf.uchile.cl" 
    handle = Entrez.efetch(db="snp", id=line, retmode="text")
    Variant_Type = "N.A."
    NAME = list()
    GENE_ID = list()    
    for element in handle.readline().strip().split("<"):
        if re.match("SNP_CLASS>\w+", element):
            Variant_Type = element.replace("SNP_CLASS>", "")
        if re.match("NAME>\w+", element):
            NAME1 = element.replace("NAME>", "")
            NAME.append(NAME1)
        if re.match("GENE_ID>\w+", element):
            GENE_ID1 = element.replace("GENE_ID>", "")
            GENE_ID.append(GENE_ID1)
    if len(GENE_ID) == 2:
        file2.write("\t".join([line, Variant_Type, str(len(GENE_ID)), GENE_ID[0], NAME[0], GENE_ID[1], NAME[1]]) + '\n')
    elif len(GENE_ID) == 1:
        file2.write("\t".join([line, Variant_Type, str(len(GENE_ID)), GENE_ID[0], NAME[0]]) + '\n')
    elif len(GENE_ID) > 2:
        file2.write(line + "\t" + Variant_Type + "\t" + str(len(GENE_ID)) + "\t" + "Revisar, aumentar capacidad de busquenda" + '\n')
    elif len(GENE_ID) == 0:
        file2.write(line + "\t" + "¿Intergenico?, no presenta un gen asociado en dbSNP" + '\n')
    handle.close() 

##############################  STRING  #########################################

    def STRING():
   
        string_api_url = "https://string-db.org/api"
        output_format = "tsv"
        method = "interaction_partners" 

        request_url = "/".join([string_api_url, output_format, method]) ## Construct URL

        my_genes = list()
        params = {}
        my_genes = NAME[0]
        params = {
        "identifiers" : my_genes, # your protein
        "species" : 9606, # species NCBI identifier 
        "required_score" : 700
        }
        if len(NAME) == 2:
            my_genes = [NAME[0], NAME[1]]
            params = {
            "identifiers" : "%0d".join(my_genes), # your protein
            "species" : 9606, # species NCBI identifier 
            "required_score" : 700
            }
                   
        print(my_genes)
        print(len(NAME))
 
        response = requests.post(request_url, data=params) ## Call STRING
    
        Interactors= list()

        for line in response.text.strip().split("\n"):
            try:
                l = line.strip().split("\t")
                query_ensp = l[0]
                query_name = l[2]
                partner_ensp = l[1]
                partner_name = l[3]
                combined_score = l[5]
                print("\t".join([query_ensp, query_name, partner_name, combined_score]))
                if partner_name != None:
                    Interactors.append(partner_name)
            except:
                pass
                print("Not found in String")
    #STRING()




    