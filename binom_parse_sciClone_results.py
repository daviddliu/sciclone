#!/bin/python

import csv
import pandas as pd
import scipy
from sklearn.metrics import adjusted_rand_score
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import ipdb

plt.switch_backend('PS') 
import matplotlib
matplotlib.use('PS')


def binom_cluster_analysis(params_dir, clusters_file_name, true_file_name, output_name, coverage=500):

    with open(clusters_file_name, "r+") as clusters_file, open(true_file_name, 'r+') as true_file:
        """
        Read in true data
        """
        content = true_file.readlines()
        n_samples = int(content[2])
        # Index 3 has the number of SNVs
        n_SNVs = int(content[3])
        
        # Now get the true clusters. Reverse the content array.

        clusters = map(lambda cluster_str: map(int, cluster_str.split(";")), content[::-1][1].split(" "))
        # Convert to a list where each SNV has a cluster from 0 to n-1.
        true_cluster_assignments = [0] * n_SNVs 
        for i, cluster in enumerate(clusters):
            for SNV in cluster:
                true_cluster_assignments[SNV] = i
        
        true_VAFs = []

        # Cross reference to get the cluster frequencies
        # Offset 3
        for i in range(1+3, n_samples+1+3):
            sample_VAFs = map(float, content[i].split(" "))
            cluster_VAFs = set()
            for i, VAF in enumerate(sample_VAFs):
                cluster_VAFs.add((VAF, true_cluster_assignments[i]))

            true_sample_VAFs = [0] * len(clusters)
            for cluster_VAF in cluster_VAFs:
                true_sample_VAFs[cluster_VAF[1]] = cluster_VAF[0]

            true_VAFs.append(true_sample_VAFs)

        a_file_name = params_dir + 'binom_params_as'
        b_file_name = params_dir + 'binom_params_bs' 

        with open(a_file_name) as a_file, open(b_file_name) as b_file:
            a_params = pd.DataFrame.from_csv(a_file, sep= "\t")
            b_params = pd.DataFrame.from_csv(b_file, sep= "\t")

        """
        Read in putative clusters
        """

        VAFs = pd.DataFrame.from_csv(clusters_file, sep = "\t")

        # Drop useless columns
        column_names = list(VAFs.columns.values)
        for name in column_names:
            if ("var" not in name.split(".")) and ("ref" not in name.split(".")) and (name != "cluster"):
                del VAFs[name]

        # Cluster assignments
        putative_cluster_assignments = VAFs['cluster'].tolist()

        """
        Make some calculations.
        """

        # Pool reads for each cluster for each sample
        num_samples = (len(list(VAFs.columns.values)) - 1)/2
        num_clusters = len(set(putative_cluster_assignments))

        # Calculate ARI
        rand_score = adjusted_rand_score(true_cluster_assignments, putative_cluster_assignments)
        
        # Calculate SNV parameter accuracy
        param_error_vector = []

        """
        Plot: real VAFs and cluster versus these inferred ones.
        """
        with PdfPages(output_name) as pdf:
            fig = plt.figure(figsize=(8, max(6, 2*num_samples)))

            # Assign each cluster a color
            individual_VAF_colors = cm.rainbow(scipy.linspace(0,1, len(set(true_cluster_assignments))))

            # Plot for each cluster 
            for sample_index in range(num_samples):
                # Pool reads for the clusters of this sample
                # -> {0: [var, ref], ... ,p: [var, ref]]
                cluster_params = []
                for i in range(num_clusters):
                    cluster_params.append([0,0])
                sample_var_reads = VAFs['Sample %d.var' % (sample_index+1)].tolist()
                sample_ref_reads = VAFs['Sample %d.ref' % (sample_index+1)].tolist()
                for i, (var, ref) in enumerate(zip(sample_var_reads, sample_ref_reads)):
                    cluster_params[putative_cluster_assignments[i]-1][0] += int(var)
                    cluster_params[putative_cluster_assignments[i]-1][1] += int(ref)
                sample_pooled_cluster_reads = cluster_params
                for true_VAF in true_VAFs[sample_index]:
                    param_error_vector.append(min([abs((sr[0]/float(sr[0]+sr[1]) - true_VAF)) for sr in sample_pooled_cluster_reads]))

                # Plot the actual VAFs
                ax1 = fig.add_subplot(num_samples, 1, sample_index + 1)
                ax1.set_title("Sample %d" % (sample_index + 1), y=0.80)
                ax1.set_xlabel("VAF")
                ax1.set_ylabel("PDF")
                if coverage > 500:
                    ylim = coverage/2 
                else:
                    ylim = coverage
                ax1.set_ylim(0, ylim)

                # Individual stuff
                for i, (var_read, ref_read) in enumerate(zip(sample_var_reads, sample_ref_reads)):
                    a = var_read
                    b = ref_read
                    x = scipy.linspace(scipy.stats.beta.ppf(0.01, a, b), scipy.stats.beta.ppf(0.99, a, b), 100)
                    # Get the colors based on the cluster membership.
                    color = individual_VAF_colors[true_cluster_assignments[i]]
                    ax1.plot(x, scipy.stats.beta.pdf(x, a, b), color=color, linewidth=0.5)

                # True clusters
                plotted = dict() 
                for i, VAF in enumerate(true_VAFs[sample_index]):
                    color = individual_VAF_colors[i]
                    if VAF in plotted:
                        plotted[VAF] -= 0.3
                    else:
                        plotted[VAF] = 2.0
                    ax1.plot((VAF, VAF), (0, ylim), color=color, linewidth=plotted[VAF])

                # Overlay the clusters
                sample_as = a_params['Sample %d.var' % (sample_index+1)]
                sample_bs = b_params['Sample %d.var' % (sample_index+1)]
                iteration = 0
                for (params, a_val, b_val) in zip(sample_pooled_cluster_reads, sample_as.iteritems(), sample_bs.iteritems()):
                    a =  params[0] + a_val[1]#+ iteration * 10# var reads + a
                    b =  params[1] + b_val[1]#+ iteration * 10# ref reads + b
                    x = scipy.linspace(scipy.stats.beta.ppf(0.01, a, b), scipy.stats.beta.ppf(0.99, a, b), 100)
                    # x1 = scipy.linspace(scipy.stats.beta.ppf(0.01, params[0], params[1]), scipy.stats.beta.ppf(0.99, params[0], params[1]), 100)
                    ax1.plot(x, scipy.stats.beta.pdf(x, a, b), 'black')
                    # ax1.plot(x1, scipy.stats.beta.pdf(x, params[0], params[1]), 'brown')
                    iteration += 1
                
                ax1.set_xlim(-0.005, ax1.get_xlim()[1])
                
                if sample_index == 0:
                    true_patch = mpatches.Patch(color='white', label='True cluster frequencies and reads')
                    black_patch = mpatches.Patch(color='black', label='Cluster posterior distributions')
                    # brown_patch = mpatches.Patch(color='brown', label='Pooled reads beta distribution')

                    plt.legend(handles=[true_patch, black_patch], prop={'size':6})

            # Show the plot
            parameter_error = scipy.mean(param_error_vector)
            fig.suptitle("Binomial clustering | True clusters: %d, found %d | ARI: %f | Parameter error: %f" % (len(set(true_cluster_assignments)), num_clusters, rand_score, parameter_error))
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            return rand_score, parameter_error


binom_cluster_analysis('data/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_0/', 'data/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_0/binom_clusters', 'data/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_0.true', 'results.pdf', coverage=1000)
