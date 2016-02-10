import numpy as np
import pandas as pd 
import collections
import pickle

def genes_indexing_data(annotation_filename,data_filename):
	#Returns dictionary of gene symbols linking to a dictionary of data indexed by condition
	#Requires that you type 'probe' above the probe column in the datafile
	ann_df = pd.read_csv(annotation_filename)
	ann_df = ann_df[~pd.isnull(ann_df["Gene Symbol"])]
	data_df = pd.read_csv(data_filename)
	probeset_to_gene_symbol = dict(zip(ann_df["Probeset ID"], ann_df["Gene Symbol"]))
	gene_symbol_to_row = {}
	for _, row in data_df.iterrows():
		if row.probe in probeset_to_gene_symbol:
			gene_symbol = probeset_to_gene_symbol[row.probe]
			#appends NDArray, after removing probe data, which will serve as output
			gene_symbol_to_row[gene_symbol] = row.values[1:]
	return (gene_symbol_to_row,data_df.columns)

#$import 
def genes_indexing_promoters(promoter_filename,filter_alternatives=True):
	promoter_df = pd.read_csv(promoter_filename,header=None)
	genes_to_promoters = {}
	gene_names = []
	new_seq = ''
	for row in promoter_df.iterrows():
		if row[1][0][0] == '>':
			#This code changes the promoter sequence into an array with 4 positions for every 1 bp possibility
			new_seq_updated = np.zeros((len(new_seq))*4)
			for i in range(0,len(new_seq)-1):
				seq_letter = new_seq[i]
				if seq_letter == 'A':
					new_seq_updated[i*4] = 1
				elif seq_letter == 'G':
					new_seq_updated[i*4+1] = 1
				elif seq_letter == 'C':
					new_seq_updated[i*4+2] = 1
				elif seq_letter == 'T':
					new_seq_updated[i*4+3] = 1
			try: 
				genes_to_promoters[gene_name] = new_seq_updated
			except:
				pass
			#change the below to keep the alternative promoter information 
			#This is hacky redo this part
			gene_name = row[1][0].split(' ',2)[1].split('_',1)[0]
			gene_names.append(gene_name)
			new_seq = ''
		else:
			new_seq += row[1][0]
	if filter_alternatives:
		duplicate_promoters = [item for item, count in collections.Counter(gene_names).items() if count > 1]
		for i in duplicate_promoters:
			genes_to_promoters.pop(i)
	return genes_to_promoters


def return_nn_input_from_promoters_and_data(promoter_filename,annotation_filename,data_filename,filter_alternatives=True):
	#currently takes around 72 seconds to run
	gene_symbol_to_row = genes_indexing_data(annotation_filename,data_filename)[0]
	genes_to_promoters = genes_indexing_promoters(promoter_filename,filter_alternatives=True)
	output = []
	genes_not_indexed_by_promoters = []
	for i in gene_symbol_to_row.keys():
		try:
			output.append((genes_to_promoters[i],gene_symbol_to_row[i]))
		except:
			genes_not_indexed_by_promoters.append(i)
	return output

hi = return_nn_input_from_promoters_and_data("Mouse_promoter_neg_500_to_0.csv","1_STGene_Mouse_Affy_Annotations.csv","new_microarray2.csv")

pickle.dump(hi,open('promoters_microarray.p','wb'))
