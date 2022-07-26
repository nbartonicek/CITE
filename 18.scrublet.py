import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

input_dir = '/share/ScratchGeneral/nenbar/projects/CITE/scripts'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/genes.tsv'))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.08,sim_doublet_ratio=10)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
	min_cells=3, 
	min_gene_variability_pctl=85, 
	n_prin_comps=100)

np.savetxt("doublet_scores.txt",doublet_scores)
np.savetxt("predicted_doublets.txt",predicted_doublets)
