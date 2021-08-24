import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as plt
import matplotlib.pyplot as plt
#import pyreadr
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--geneexpfile', required = True,
                    help = 'path to the geneexpfile file with cohort participant id as row names')
parser.add_argument('-c', '--clusterfile', required = True,
                    help = 'path to file with cohort_participantid and cluster number')
parser.add_argument('-t', '--clustertype', required = True,
                    help = 'column name with cluster type in cluster file')
parser.add_argument('-o', '--output_prefix', required = True,
                    help = "Both top and bottom features plots will start with same prefix")
parser.add_argument('-v', '--vtest_outname', required = True,
                    help = "Output name with vtest results")
args = parser.parse_args()

# Reading cluster file and renaming column name
clusteringtype = args.clustertype
clusters = pd.read_csv(args.clusterfile, sep="\t")
clusters = clusters.rename(columns={
    "Kids_First_Biospecimen_ID": "Kids_First_Biospecimen_ID"})[["Kids_First_Biospecimen_ID", clusteringtype]]

# Reading in geneexpfile file
gene_features = pd.read_table(args.geneexpfile, sep="\t").T


#print (gene_features)

gene_features = gene_features.reset_index().rename(columns={
    "index": "Kids_First_Biospecimen_ID"}).set_index("Kids_First_Biospecimen_ID")

#print (gene_features)

## Merging cluster labels and geneexpfile scores
gene_with_cluster = clusters.merge(gene_features, on='Kids_First_Biospecimen_ID')
gene_with_cluster = gene_with_cluster.iloc[:,1:] # Removing sample name
#gene_with_cluster = gene_with_cluster.iloc[:,1:] # Removing sample name

#print (gene_with_cluster)

gene_with_cluster = gene_with_cluster.astype(float)
gene_with_cluster[clusteringtype] = gene_with_cluster[clusteringtype].astype(int)

## Mean and var(sigma sq) for all features from geneexpfile score datafram
gene_means = gene_with_cluster.iloc[:, 1:].mean(axis = 0)
gene_var = gene_with_cluster.iloc[:, 1:].var()


## This part of the code calculates the formula from slide 12 in these
###. slides - http://eric.univ-lyon2.fr/~ricco/cours/slides/en/classif_interpretation.pdf
gene_groupedbycluster = gene_with_cluster.groupby(clusteringtype)
vt_sign_values = pd.DataFrame({})

# iterate over each group
for group_name, df_group in gene_groupedbycluster:
    gene_df = df_group.iloc[:, 1:]
    num_samples = gene_df.shape[0]
    numerator = gene_df.mean() - gene_means
    denominator = np.sqrt(((len(gene_with_cluster) - num_samples)/(len(gene_with_cluster) -1)) * (gene_var/num_samples))
    vt_sign_values[group_name] = numerator/denominator


# Renaming and choosing improtant features based on v.test scores
vt_sign_values = vt_sign_values.reset_index().rename(columns={"index": "feature"})
vt_sign_values_plotinput_reformatted = vt_sign_values.set_index('feature').stack().reset_index(name='vt_score').rename(columns={'level_1':'cluster'})

#print (vt_sign_values)
#print (vt_sign_values_plotinput_reformatted)

vt_sign_values_plotinput_reformatted.to_csv(args.vtest_outname, sep="\t", index=None)

# Removing index from vt_sign_values_plotinput_top5
vt_sign_values_plotinput_top5 = vt_sign_values_plotinput_reformatted.groupby(["cluster"]).apply(
    lambda x: x.sort_values(["vt_score"], ascending = False).head(10))
vt_sign_values_plotinput_top5.index = range(len(vt_sign_values_plotinput_top5))


vt_sign_values_plotinput_bottom5 = vt_sign_values_plotinput_reformatted.groupby(["cluster"]).apply(
    lambda x: x.sort_values(["vt_score"]).head(10))
vt_sign_values_plotinput_bottom5.index = range(len(vt_sign_values_plotinput_bottom5))

# Plotting top 10 features for every cluster
width_in_inches = 6
height_in_inches = 10
colors = ["rosybrown", "darkgray", "brown", "peru", "coral", "mediumpurple", "silver", "honeydew", "wheat",
         "lightcyan", "olive", "darkseagreen", "slategray", "teal", "skyblue", "lavender", "plum", "goldenrod",
         "khaki", "lightpink", "palevioletred", "y", "lightblue", "lemonchiffon", "seagreen"]

customPalette = sns.set_palette(sns.color_palette(colors))
plt.figure(figsize=(width_in_inches, height_in_inches))
sns.stripplot(x="cluster", y="vt_score",hue="feature", data=vt_sign_values_plotinput_top5,
              size=8, palette="bright",edgecolor='black', linewidth=0.5,jitter=0.3)

gfg = sns.stripplot(x="cluster", y="vt_score",hue="feature", data=vt_sign_values_plotinput_top5,
                            size=8, palette="bright",edgecolor='black', linewidth=0.5,jitter=0.3)

plt.legend(bbox_to_anchor=(1, 1))
plt.setp(gfg.get_legend().get_texts(), fontsize='8')
plt.title("Top 10 splicing events for every cluster")
plt.savefig(args.output_prefix+"_top10features.png", bbox_inches = "tight")

# Plotting bottom 5 features for every cluster
width_in_inches = 6
height_in_inches = 10
colors = ["rosybrown", "seagreen", "brown", "peru", "coral", "mediumpurple", "silver", "honeydew", "wheat",
         "lightcyan", "olive", "darkseagreen", "slategray", "teal", "skyblue", "lavender", "plum", "goldenrod",
         "khaki", "lightpink", "palevioletred", "y", "lightblue", "lemonchiffon"]

customPalette = sns.set_palette(sns.color_palette(colors))
plt.figure(figsize=(width_in_inches, height_in_inches))
sns.stripplot(x="cluster", y="vt_score",hue="feature", data=vt_sign_values_plotinput_bottom5,
              size=8, palette="bright", edgecolor='black', linewidth=0.5,jitter=0.3)

gfg = sns.stripplot(x="cluster", y="vt_score",hue="feature", data=vt_sign_values_plotinput_top5,
                            size=8, palette="bright",edgecolor='black', linewidth=0.5,jitter=0.3)

plt.setp(gfg.get_legend().get_texts(), fontsize='8')
plt.legend(bbox_to_anchor=(1, 1))
plt.title("Bottom 10 splicing events for every cluster")
plt.savefig(args.output_prefix+"_bottom10features.png", bbox_inches = "tight")
