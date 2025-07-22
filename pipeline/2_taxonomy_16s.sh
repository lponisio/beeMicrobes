#!/usr/bin/env bash

## *****************************************************************************
##                                   16s 
## *****************************************************************************
# 1. ASSIGN TAXONOMY 16s
## *****************************************************************************

# The aim of this script is to assign classifications to the representative sequences (rep-seqsX.qza)
# that were clustered with denoising algorithms in the seq_prep.sh script

# navigate to the dropbox folder containing saved data files
# cd ~/Dropbox/beeMicrobes_saved/beeMicrobes_pipeline_output

# Use a information-rich database that is clustered at 99% sequence similarity at least 
# In our case, using Silva for 16s, and NCBI AND RDP for rbcl
# we have to "train" the classifier dataset just once.

## IF YOU ALREADY HAVE THE CLASSIFIERS FOLDER, SKIP TO 1c

# 1a. Download the newest silva 132 database into a new working directory from https://www.arb-silva.de/download/archive/qiime

# use 99_otus_16S.fasta and  consensus_taxonomy_7_levels.txt to create the training set.
# mkdir 16s-trainingclassifier
# cd 16s-trainingclassifier
# wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
# unzip Silva_132_release.zip
# rm SILVA_132_release.zip

#downloaded from: https://docs.qiime2.org/2023.9/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier
#silva 136 SSURef N99 seq and tax files

# 1b. Train feature classifier
# import reference sequences from silva data as a qiime2 artifact

## go back into qiime
docker run -it \
  -v /Volumes/bombus/ncullen/University\ of\ Oregon\ Dropbox/Nevin\ Cullen/beeMicrobes_saved/beeMicrobes_pipeline_output:/mnt/beeMicrobes_pipeline_output \
  qiime2/core:2019.1

#updated silva classifier is in 16s_classifier folder
cd ../../mnt/beeMicrobes_pipeline_output/training_classifiers/16s

## 99 is 99% match between our seq and the database
# qiime tools import \
# --type 'FeatureData[Sequence]' \
# --input-path silva-138-99-seqs.qza \
# --output-path 99_otus_16S.qza
# 
# # import taxonomy strings. Check and see if your taxonomy file is a tab-seperated file without a header.
# # if it doesnt have a header, specify "headerlessTSVTaxonomyFormat" since the default source formats usually have headers
# 
# qiime tools import \
# --type 'FeatureData[Taxonomy]' \
# --input-format TSVTaxonomyFormat \
# --input-path taxonomy/16S_only/99/majority_taxonomy_7_levels.txt \
# --output-path 99_otus_16S_taxonomy.qza

# We now have two Qiime2 artifacts, 99_otus_16s.qza (reference sequences) and 99_otus_16s_taxonomy.qza (taxonomic names). 
# trim silva to my region using my sequencing primers. We tell the algorithm our genomic primer forward and reverse sequences
# we do this because taxonomic classification is more accurate when a naive bayes classifier is trained on the region
# of the 16s sequence that we sequenced (Werner et al. 2012).

# NOTE: we will deal with 2 similarly named files, ref-seqs16s.qza (reference sequences) and
# rep_seqs16s.qza (representative sequences). These are different files and are not interchangeable.

qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CMGGATTAGATACCCKGG \
  --p-r-primer AGGGTTGCGCTCGTTG \
  --o-reads ref-seqs16s.qza

# visualize:
qiime feature-table tabulate-seqs \
  --i-data ref-seqs16s.qza \
  --o-visualization ref-seqs16s.qzv
#qiime tools view ref-seqs16s.qzv # (but this is big and took forever to load...).

## may need to clean up docker memory usage
#docker system prune

#Train the classifier:
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs16s.qza \
  --i-reference-taxonomy silva-138-99-tax.qza  \
  --o-classifier classifier16s.qza

# 1c. classify rep seqs and get the resulting taxonomic ids 

# may need to install  scikit learn
# cd ../ until you get into your main folder of computer/wherever you install stuff
# pip install -U scikit-learn==0.19.1 #OR WHATEVER VERSION IS COMPATIBLE WITH THE CONTAINER! Might be an old version

qiime feature-classifier classify-sklearn \
  --i-classifier classifier16s.qza \
  --i-reads  ../../merged/16s/rep-seqs16s.qza \
  --o-classification  ../../merged/16s/taxonomy16s.qza

# switch to the newest version of qiime
### Why do we need to switch between versions of qiime for these two tasks?
exit

docker run -it \
  -v /Volumes/bombus/ncullen/University\ of\ Oregon\ Dropbox/Nevin\ Cullen/beeMicrobes_saved/beeMicrobes_pipeline_output:/mnt/beeMicrobes_pipeline_output \
  qiime2/core

# visualize. navigate back to where you have your taxonomy files
cd ../../mnt/beeMicrobes_pipeline_output/training_classifiers/16s

qiime metadata tabulate \
  --m-input-file taxonomy16s.qza \
  --o-visualization taxonomy16s.qzv

### ************************************************************************
# 2. FILTERING STEPS (PART A)
### ************************************************************************

# Go through and filter out unwanted reads on the merged tables and repseqs.
# We can't do all the filtering steps on the merged files, 
# but we can do some, then we make a phylogenetic tree, then go back and filter round-specific issues

cd ../../mnt/beeMicrobes_pipeline_output/merged/16s

#2b. filter 1: out the chloroplast and mitochondria reads 

qiime taxa filter-table \
  --i-table table16s.qza \
  --i-taxonomy taxonomy16s.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table tablefilt1.qza
  
qiime feature-table summarize \
  --i-table tablefilt1_16s.qza \
  --o-visualization tablefilt1_16s.qzv 

#2c. filter 2: remove sequences only found in one sample

qiime feature-table filter-features \
  --i-table tablefilt1.qza \
  --p-min-samples 2 \
  --o-filtered-table tablefilt2.qza
  
qiime feature-table summarize \
  --i-table tablefilt2.qza \
  --o-visualization tablefilt2.qzv

#2d: Filter our rep seqs file so that you can see what samples you have left after filtering and subsampling

qiime feature-table filter-seqs \
  --i-data rep-seqs-16s.qza \
  --i-table tablefilt2.qza \
  --o-filtered-data rep-seqs-16s-filtered.qza

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-16s-filtered.qza \
  --o-visualization rep-seqs-16s-filtered.qzv
#view rep-seqs-16s-filtered.qzv in qiime2 viewer in browser

## here is a step we are changing for 2023 -- since we are using the trees to determine microbial phylogenetic distance, we need to 
## filter out unwanted sequences before we generate the trees.
## also, we will not split sequences into runs and instead will filter out contaminants across all runs

# 2e. Filter out large taxonomic groups of bacteria that are known to be contaminants

# there are known bacteria are salt-loving and often found in buffers. 
# they include  Halomonas, Shewanella, Oceanospirillales, and the acne bacteria Propionibacterium. 

qiime taxa filter-table \
  --i-table table16s.qza \
  --i-taxonomy taxonomy16s.qza \
  --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales \
  --o-filtered-table tablef1.qza

# 2f. Let's specifically look at sequences in our controls and remove the bacteria that are in them
qiime taxa barplot \
  --i-table tablef1.qza \
  --i-taxonomy taxonomy16s.qza \
  --m-metadata-file ../../maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
  --o-visualization taxa-bar-plots-f1.qzv

# now filter out taxa that are in the controls!
# manually selected from the table16s file which ASVs were found in controls and wrote them down here

########################
### Manual Filtering ###
########################

# use the following rationale for selecting bacterial sequences to remove:
# look at what bacteria show up our control samples. we then, for each one of these
# bacteria, look up if it is in all the samples or just that control and made a call about whether or not to remove it
# if bacteria is present in just the control sample, you want to remove it
# however, some are also present in a lot of other samples! we don't want to lose data. if bacteria is in the control AND in more than 30 samples just leave it and don't filter it out (10% of samples)
# you may want to make an exception if the bacterial contaminant is obviously present in one just one plate, because even though its in a lot of samples its likely a contaminant
# another exception is if you have the contaminant in a lot of samples BUT also in a lot of the controls, get rid of it

qiime taxa filter-table \
  --i-table tablef1.qza \
  --i-taxonomy taxonomy16s.qza \
  --p-mode exact \
  --p-exclude "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Microbacteriaceae;g__Galbitalea",\
  "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae",\
  "d__Bacteria;p__Tenericutes;c__Mollicutes;o__Entomoplasmatales;f__Spiroplasmataceae;g__Spiroplasma",\
  "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Kineosporiales;f__Kineosporiaceae;g__Kineosporia;s__uncultured;bacterium",\
  "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__Segetibacter",\
  "d__Bacteria",\
  "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;__",\
  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Comamonadaceae",\
  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Acidovorax",\
  "d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium",\
  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Ralstonia",\
  "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus",\
  "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Phenylobacterium",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Burkholderiaceae;__;__",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Neisseriaceae;D_5__Snodgrassella;D_6__Snodgrassella alvi",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Orbales;D_4__Orbaceae;D_5__Gilliamella;Ambiguous_taxa" \
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Lactococcus",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Negativicutes;D_3__Selenomonadales;D_4__Veillonellaceae;D_5__Dialister;Ambiguous_taxa",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Acetobacterales;D_4__Acetobacteraceae;D_5__Saccharibacter;Ambiguous_taxa" \
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Carnobacteriaceae;D_5__Granulicatella",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Lactococcus",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Cardiobacteriales;D_4__Cardiobacteriaceae;D_5__Cardiobacterium;__" \
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;D_5__Bacillus",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Planococcaceae;D_5__Planomicrobium",\
  "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Actinomycetales;D_4__Actinomycetaceae;D_5__Actinomyces",\
  "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Corynebacteriales;D_4__Corynebacteriaceae;D_5__Corynebacterium;D_6__Corynebacterium durum",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Aerococcaceae;D_5__Abiotrophia;D_6__unidentified",\
  "D_0__Bacteria;D_1__Fusobacteria;D_2__Fusobacteriia;D_3__Fusobacteriales;D_4__Fusobacteriaceae;D_5__Fusobacterium",\
  "D_0__Bacteria;D_1__Fusobacteria;D_2__Fusobacteriia;D_3__Fusobacteriales;D_4__Leptotrichiaceae;D_5__Leptotrichia;D_6__Leptotrichia sp. oral taxon 212",\
  "D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Betaproteobacteriales;D_4__Neisseriaceae;D_5__Neisseria",\
  "D_0__Bacteria;D_1__Tenericutes;D_2__Mollicutes;D_3__Entomoplasmatales;D_4__Spiroplasmataceae;D_5__Spiroplasma" \
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Carnobacteriaceae;D_5__Granulicatella",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Enterococcaceae;D_5__Enterococcus",\
  "D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Lactococcus",\
  "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Corynebacteriales;D_4__Corynebacteriaceae;D_5__Corynebacterium 1"
  "Unassigned"
  --o-filtered-table tablef2.qza

# qiime taxa barplot \
#   --i-table tablef2.qza \
#   --i-taxonomy taxonomy16s.qza \
#   --m-metadata-file maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
#   --o-visualization taxa-bar-plots-f2.qzv

#Filter out controls by making new map file excluding controls
qiime feature-table filter-samples \
  --i-table tablef2.qza \
  --m-metadata-file ../../maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
  --o-filtered-table tablef3.qza

#now need to filter the rep-seqs FeatureTable[sequence] to the updated filtered FeatureTable[frequency]

qiime feature-table filter-seqs \
--i-data rep-seqs-16s-filtered.qza \
--i-table tablef3.qza \
--o-filtered-data rep-seqs-16s-complete-filter.qza


### ************************************************************************
## 3. #GENERATE TREE FOR PHYLOGENETIC DIVERSITY ANALYSES
### ************************************************************************

# generate a tree with the merged data

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-16s-final-filter.qza \
  --o-alignment aligned_repseqs16s.qza \
  --o-masked-alignment masked_aligned_repseqs16s.qza \
  --o-tree unrooted_tree16s.qza \
  --o-rooted-tree rooted-tree16s.qza

## NOTE: need to double check if the error below is still occurring
## now to fix the error: The table does not appear to be completely represented by the phylogeny.
#need to drop from the tree features not in the table anymore after filtering

# qiime fragment-insertion filter-features \
# --i-table tablef3.qza \
# --i-tree rooted-tree16s.qza \
# --o-filtered-table tablef4.qza \
# --o-removed-table removed-tablef4.qza
### *************************************************************************
#4. DETERMINE SUBSAMPLE DEPTH
### *************************************************************************

#We want to make a rarefaction curve to see how subsampling depths influence our alpha diversity metrics.
#Open visualization in Qiime2 for tablef5. View interactive sample detail and look at visualized table to see impact of different depths
# record number of reads, number of samples at specified sampling depth, and how many features retained

#For subsampling:  840 reads.  715 (90.62%) samples at the specifed sampling depth.Retained 600,600 (5.71%) features

qiime diversity alpha-rarefaction \
  --i-table tablef5.qza \
  --i-phylogeny rooted-tree16s.qza \
  --p-max-depth 10000 \
  --m-metadata-file ../../maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
  --o-visualization alphararefact16s.qzv

qiime feature-table summarize --i-table tablef5.qza --o-visualization tablef5.qzv

### *************************************************************************
# 5. ALPHA AND BETA DIVERSITY
### *************************************************************************

#Generate core metrics folder. This command makes the new directory and outputs a bunch of files in it
# make sured youre in "merged/16s"

## NOTE: NEED TO RENAME THIS FOLDER
mkdir final

## NOTE: Need to find map files and generate combined map file.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree16s.qza \
  --i-table tablef5.qza \
  --p-sampling-depth 840 \
  --m-metadata-file ../../maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
  --output-dir final/core_metrics16s \
  --verbose

#now we have lots of files in that core_metrics directory. gotta export rarefied_table.qza so we can get a qzv and convert to a csv
#using combined full map
qiime taxa barplot \
  --i-table final/core_metrics16s/rarefied_table.qza \
  --i-taxonomy taxonomy16s.qza \
  --m-metadata-file ../../maps/16s/beeMicrobes_combined_map_no-ctrls.txt \
  --o-visualization final/core_metrics16s/rarefiedtable16s.qzv

#drag the new qzv into qiime 2 view, set taxonomic level to 7, and download a csv! put this into "final asv tables" folder

# end of script, proceed to X_taxonomy_RBCL.sh