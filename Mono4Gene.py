
import os.path
import glob
import time
import sys
import csv
import pandas as pd
import numpy as np

from ete3 import Tree


############## Apresentação do splitbyscaffold ################
print("Mono4Gene ver 1 beta  Jun 2021")
print("A procudere for testing the monophyly of a putative group over different gene trees.")
print("by Vitor Miranda & Saura R. Silva")
print("São Paulo State University - Jaboticabal, Brazil")
print("e-mail: vitor.miranda@unesp.br")
print()
###############################################################


dir_gene_trees      = "gene_trees"          # Name of the directory containing the gene tree files
ext_tree            = [".nwk", ".newick", ".treefile", ".tre", ".txt"]      # Newick tree file extensions

#_tested_group_nwk   = '(Mguttatus,Usimulans, Ufoliosa);'           

_report_CSV          ='report.csv'   # lista com resultados



########## Caso não tenha sido informato o INFILE=tested group #######
if len(sys.argv)==1:
       print("Usage:  python Mono4Gene.py INFILE ")
       print("\nINFILE = File (newick format) with putative monophyletic group to be test.")
       print()
       sys.exit()
######################################################################
       

_tested_group_nwk = sys.argv[1]         # Newick file with group to be tested


###########################################################
############### Input trees from file #####################

if not os.path.isdir(dir_gene_trees):
    print("Tree files should be included in the directory: " + dir_gene_trees + ".")
    sys.exit() 
	
start_time = time.time()

tree_files = []
filenames = os.path.join(dir_gene_trees, '*')

for file in glob.glob(filenames):
    root, ext = os.path.splitext(file)
    if ext in ext_tree:
        tree_files.append(file)
        #print(file)
        
        
if len(tree_files) == 0:
    print("Not found Tree (newick) files in the directory: " + dir_gene_trees + ".")
    sys.exit() 

###########################################################
###########################################################


###########################################################
##################### GROUP TO BE TESTED ##################
t = Tree(_tested_group_nwk)

tested_group = []

for leaf in t:
    #print (leaf.name)
    tested_group += [leaf.name]

###########################################################
###########################################################



################### DATAFRAME CREATION #####################

df = pd.DataFrame([], index=range(0,len(tree_files)), columns=["Tree file(s)", "Group status", 
                                                                "# Group(s) avoiding monophyly", 
                                                                "Terminal(s) not present in tree",
                                                                "% Absence",
                                                                "Index S"])

#df = pd.DataFrame([], index=range(0,len(tree_files)), columns=["col1", "col2", "col3", "col4","col5"])
############################################################



################################################################
######## TEST IF GROUP IS MONOPHYLETIC IN THE TREE_GENE ########

'''
print('Group to be tested in the tree(s): ', '{: <20}'.format(str(tested_group)[1:-1]))

print('')
print('------------------------------------------------------------------------------------------------------------------------------------------------------------')
print('                                                                          # Group(s)          Terminal(s) not                          ')
print('Tree file(s)                                                Group status  avoiding monophyly  present in tree')
print('------------------------------------------------------------------------------------------------------------------------------------------------------------')

'''

print('Checking file(s) and computing...')

_i=-1

for files_n in tree_files:
    _i+=1
    Gene_tree = Tree(files_n, format=1)
    #print(Gene_tree)

    X = Gene_tree.check_monophyly(values=tested_group, target_attr="name", ignore_missing=True, unrooted=True)
    if X[0]: status_group = 'monophyletic'
    else: status_group = 'polyphyletic'


    ###### Test if each element of tested group is present in each tree ######
    Absent_elements = []
    for element in tested_group:
        if element not in Gene_tree: Absent_elements.append(element) 
    ##########################################################################


    if tested_group==Absent_elements: status_group = "N/A"   # case none taxa of tested group is present in the tested tree
    
   
    '''
    ########## print line in screen #########         
    print('{:-<59}'.format(files_n[len(dir_gene_trees)+1:]) + ' '
        +'{: <14}'.format(status_group)
        +'{: <19}'.format(str(len(X[2])))
        +'{: <30}'.format(str(Absent_elements)[1:-1]))
    #########################################
    '''

    N_tips_avoinding_monophyly = len(X[2])
    N_Absent_elements = len(Absent_elements)
    _PP_Absence = N_Absent_elements / len(tested_group) 
    _PP_Absence = round(_PP_Absence,4)


    ####### record lines in report file SCV ######         
    l=[]
    
    l.append(files_n[len(dir_gene_trees)+1:])
    l.append(status_group) 
    l.append(str(len(X[2])))
    l.append(str(Absent_elements)[1:-1])
    l.append(_PP_Absence)
    
    

    ###### Compute the Index S #######
    
    positive_points = 0
    negative_points = 0

    if status_group == 'monophyletic': positive_points = 100
    elif status_group == 'polyphyletic': positive_points = 50


    negative_points = N_tips_avoinding_monophyly

    __index = (positive_points - negative_points) * (1- _PP_Absence)
    
    if status_group == "N/A": __index = -9


    '''

    if status_group == 'monophyletic':
        if N_tips_avoinding_monophyly==0 : __index = 1
        else: __index = 1/N_tips_avoinding_monophyly

    elif status_group == 'polyphyletic':
        __index = 1/ (N_tips_avoinding_monophyly+N_Absent_elements)

    

    if status_group == "N/A": __index = -9


    '''


    l.append(round(__index,3))
    ###################################


    df.loc[_i] = l
    ##############################################



print('\nJob is done!!!\n')

sorted_df = df.sort_values(["Index S"], ascending=False)

print('Indexed genes by Index S')
print('----------------------------------------------------------------------------------------------------------------------------------------------------------------')
print (sorted_df)

print()


df.to_csv (_report_CSV, index = None, header=True) 
sorted_df.to_csv ('indexed_'+_report_CSV, index = None, header=True)


##################### Print elapsed time #######################
end_time = time.time()

print('The putative group was tested in ', len(tree_files), 'tree(s): ', '{: <20}'.format(str(tested_group)[1:-1]))
print()
print("Report created in CSV file        : "+ _report_CSV)
print("Indexed report created in CSV file: indexed_"+_report_CSV)
print("\nElapsed time: {0} sec".format(end_time - start_time))
print()
################################################################

