import os
import re
os.chdir('/Users/James/Desktop/Dissertation/16S_Vs_WGS/UrbanRuralChina-master/16SrRNA/inputData/qiime/')

#Borrowed code for multiple string replacement
def replaceMultiple(mainString, toBeReplaces, newString):
    # Iterate over the strings to be replaced
    for elem in toBeReplaces :
        # Check if string is in the main string
        if elem in mainString :
            # Replace the string
            mainString = mainString.replace(elem, newString)
    return  mainString

taxa_levels=("kingdom","phylum","class","order","family","genus")
taxa_letter=("k","p","c","o","f","g")

for k in range(2,7):
    #solving conflict between zero based and one based indexing
    k_index=k-1
    
    filename="otu_table_L"+str(k)+"_taxa_as_col.txt"
    infile=open(filename,'r')
    
    header=infile.readline()
    header=header.strip("\n")
    columns=header.split("\t")
    
    taxa_names=["Sample_ID"]
    garbage_count=0
    for i in columns[1:]:
        temp=i.split(";")
        temp_taxa_name=temp[k_index].replace(taxa_letter[k_index]+"__","")
        if(temp_taxa_name==""):
            parent_taxa=temp[k_index-1].replace(taxa_letter[k_index-1]+"__","")
            if(parent_taxa==""):
                parent_taxa="Qiime_Garbage_Count_"+str(garbage_count)
                garbage_count=garbage_count+1
            temp_taxa_name="Unclassified_"+parent_taxa
        taxa_names.append(temp_taxa_name)
        print(taxa_names)
    
    output_name="Qiime_china_ruralurban_"+taxa_levels[k_index]+".tsv"
    outfile=open(output_name,'w')
    outfile.writelines("\t".join(taxa_names)+"\n")
    for line in infile:
        outfile.writelines(line)
    outfile.close()
    infile.close()
    