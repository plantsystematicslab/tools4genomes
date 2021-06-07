import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq


############## Apresentação do splitbyscaffold ################
print("splitbyscaffold ver 1 beta  Jun 2021")
print("A quick and efficient procedure for grouping scaffolds.")
print("by Vitor Miranda & Saura R. Silva")
print("São Paulo State University - Jaboticabal, Brazil")
print("e-mail: vitor.miranda@unesp.br")
print()
###############################################################


########## Caso não tenha sido informato o INFILE #######
if len(sys.argv)==1:
       print("Usage:  python splitbyscaffold INFILE ")
       print("INFILE = file in FASTA")
       print()
       sys.exit()       
#########################################################
       

File_In=sys.argv[1]      #Arquivo de entrada FASTA

set_scaff=[]             #lista de nomes de scaffolds presentes no infile

records_In = list(SeqIO.parse(File_In, "fasta"))

print("*** Reading",len(records_In),"FASTA sequences...")
print()

############ Percorre sequência por sequência para criar ######
############ lista de IDs de scaffolds sem repetição ##########
i=0
while (i<len(records_In)):
    #print(records_In[i].id)

    ###### cria lista de nomes de scaffolds #######
    pos_ = records_In[i].id.find(":", 0, len(records_In[i].id)) 
    
    if records_In[i].id[:pos_] not in set_scaff:  
        set_scaff += [records_In[i].id[:pos_]]

    #print(len(set_scaff),set_scaff)
    ###############################################

    i+=1
    
###############################################################
    
       
######## cria arquivos de saida com nomes curtos de scaffolds ########
j=0
while (j<len(set_scaff)):
    File_Out = open(set_scaff[j] + ".out.fas", "w")
    j+=1
######################################################################


################ cria lista com nomes de arquivos Out  ###############
set_scaff_names_fileout=[]
i=0
while (i<len(set_scaff)):
    set_scaff_names_fileout += [set_scaff[i] + ".out.fas"]
    i+=1
######################################################################



############ Percorre sequência por sequência do arquivo In ######
############ para coletar e salvar scaffold em arquivo Out #######
j=0
while (j<len(set_scaff)):

    ###### Percorre sequência por sequência do arquivo In ######
    i=0
    while (i<len(records_In)):
        #print(records_In[i].id)

        ###### Grava apenas a scaffold correspondente ao arquivo #####
        if set_scaff[j]==records_In[i].id[:records_In[i].id.find(':',0,len(records_In[i].id))] : 
            with open(set_scaff_names_fileout[j], 'a') as arquivo:
                arquivo.write('>'+str(records_In[i].id)+'\n')
                arquivo.write(str(records_In[i].seq)+'\n')
        ##############################################################
 
        i+=1

    arquivo.close()
    j+=1
    ###########################################################
        
###############################################################


###################### SUMMARY REPORT #########################
print('The job is done!')
print()
print('***Summary report:')

if i>1:
    print(i,'sequences were analyzed.')
else:
    print(i,'sequence was analyzed.')

if j>1:
    print(j, 'scaffolds were found and grouped in your own files.')
else:
    print('Only', j, 'scaffold was found.')


###############################################################



