#Conor Elling
#DSP439 Final Exam
#May 10 2023



from itertools import product
import pandas as pd

def find_all_sizes(gene = pd.Series(['ATTTGGATT'])):
    
    #The goal of this code is to find the k values necessary for a gene sequence
    
    max_size = []
    for i in range(0, len(gene)):
        max_size.append(len(gene[i]))
    k_list = []
   
    #k_list lists possible gene sequence sizes
    
    max_list = [] 
    gene_list = []
    
    #using the 'for' function, the code will loop through each gene length and return the k_list
    
    for i in range(0, len(max_size)):
        for j in range(0, max_size[i]):
            k_list.append(j+1)
            max_list.append(max_size[i])
            gene_list.append(gene[i])
    return k_list

def all_combos(letters = ['A', 'G', 'C', 'T'], k = 2):
     #This line of code is intended to find all possible combinations of size "k"
     #Replacement is used in this situation to make sure a value of k is not used twice
    combos_list = product(letters, repeat = k)
    combos = pd.DataFrame(combos_list).agg(''.join, axis=1)
    #link together the tuples so the return function produces a series
    return combos

def create_all_size_combos(k_list = [1,2,3], gene = pd.Series(['ATTTGGATT', 'ATGTCTGTCTGTA']), letters = ['A', 'G', 'C', 'T']):
    #This line of code is intended to find all combinations for possibles sizes of "k"
    combos_df_list = []
    for i in range(0, len(k_list)):
        combos = all_combos(letters = letters, k = k_list[i])
        combos_df_list.append(combos)
    combos_df = pd.concat(combos_df_list)
    return combos_df
    #The return function produces all combinations of letters for all possible sizes of "k"

def ifexists_and_count(combos = ['AT', 'GC'], gene = pd.Series(['ATGTCTGTCTGTA'])):
    #Firstly, make the series into a list
    gene_series_tolist = gene.to_list()
    #Creating empty lists will allow for loops to run though 
    combo_list = []
    count_list = []
    combo_found_list = []
    gene_list = []

   
   #Each string in the gene list is looped through 
    for j in range (0, len(gene_series_tolist)):
        gene_value = gene_series_tolist[j]
        #Find how many times the pattern occurrs in the string
        for i in range(0, len(combos)):
            count = gene_value.count(combos[i])
            count_list.append(count)
            combo_list.append(combos[i])
            gene_list.append(gene_value)
            if count > 0:
                combo_found_list.append(1)
            else:
                combo_found_list.append(0)

   #A Dataframe is created to return information from the function about the patteren
    df = pd.DataFrame()
    df['GeneString'] = gene_list
    df['SearchFor'] = combo_list
    df['Count'] = count_list
    df['IfExists'] = combo_found_list
    search_len_list = []
    possible_strings = []
    for i in range(0, len(df)):
        search_len_list.append(len(combo_list[i]))
        if len(gene_list[i]) == len(combo_list[i]):
            if combo_found_list[i] == 1:
                possible_strings.append(1)
            else:
                possible_strings.append(0)
        else:
            possible_strings.append(1)
    df['SearchLen'] = search_len_list
    df['Possible'] = possible_strings
    return df


#Count the number of possible substrings of any size, specified as k

def possible(k, glen = 9): 
    #4^k is less than the length of g
    if 4**k < glen:
        #if k equals 1 then 4 is assigned to p
        if k == 1:
            p = 4
        #If k is not 1 and 4^k is less than the length of g, then p is length of the gene sequence
        else:
            x = glen^k
            p = x
    # If 4^k is greater than length of gene sequence, p is calculated
    else:
        x = glen - (k-1)
        p = x
    return p

def main_genes():
    #This function will compute linguistic complexity as well as the amount of possible substrings and the number of observed substrings in this gene sequence
    #Here the user will enter the genetic sequence
    gene_input = input("Enter a sequence: ")
    gene = pd.Series([gene_input])

    #This calculates the counts for the gene sequence of all possible combinations
    letters = ['A', 'G', 'C', 'T']
    k_list = find_all_sizes(gene=gene)
    combos_df = create_all_size_combos(k_list = k_list, gene = gene)
    df = ifexists_and_count(combos = combos_df.to_list(),
                            gene = gene)
    
    #Data is grouped to find out how many subtrings were found
    final = df.groupby(by=['GeneString', 'SearchLen']).sum(numeric_only=True).reset_index()
    genestringlist = final['GeneString'].to_list()
    searchlenlist = final['SearchLen'].to_list()
    z_possible = []
    for i in range (0, len(final)):
        z = possible(searchlenlist[i], len(genestringlist[i]))
        z_possible.append(z)
    final['Possible'] = z_possible

    #Data is subsetted, columns renamed
    observedtable = final[['GeneString', 'SearchLen', 'IfExists', 'Possible']]
    observedtable.columns = ['Gene', 'k', 'Observed', 'Possible']

    #Linguistic Complexity 
    final = final[['GeneString', 'IfExists', 'Possible']]
    dffinal = final.groupby(by=['GeneString']).sum(numeric_only=True).reset_index()
    dffinal['LinguisticComplexity'] = dffinal['IfExists']/dffinal['Possible']
    dffinal.columns = ['Genetic Sequence', 'Observed Substrings', 'Possible Substrings', 'Linguistic Complexity']   
    

    vargene = dffinal.iloc[0, 0]
    varobserved = dffinal.iloc[0, 1]
    varpossible = dffinal.iloc[0, 2]
    varcomplexity = dffinal.iloc[0, 3]

    #Findings are printec
    print("Genetic Sequence: ", vargene)
    viewobserved = observedtable[['k', 'Observed', 'Possible']]
    print(viewobserved)
    print("Observed Substrings: ", varobserved)
    print("Possible Substrings: ", varpossible)
    print("Linquistic Complexity: ", varcomplexity)
    return dffinal, observedtable
#dffinal, observedtable = main_genes()







    
