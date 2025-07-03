#!/usr/bin/env bash

## *****************************************************************************
##                                   16s 
## *****************************************************************************
# 1. ASSIGN TAXONOMY 16s
## *****************************************************************************
#navigate to Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw

# We need taxonomic classifications 

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
  -v ~/Dropbox/BMMA2025:/mnt/BMMA2025 \
  qiime2/core:2019.1

#updated silva classifier is in 16s_classifier folder
cd ../mnt/BMMA2025/data/qiime2_data/training_classifiers/16s

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
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CMGGATTAGATACCCKGG \
  --p-r-primer AGGGTTGCGCTCGTTG \
  --o-reads ref-seqs16s.qza

# visualize:
qiime feature-table tabulate-seqs \
  --i-data ../../merged/16s/ref-seqs16s.qza \
  --o-visualization ../../merged/16s/ref-seqs16s.qzv
#qiime tools view ref-seqs16s.qzv # (but this is big and took forever to load...).

## may need to clean up docker memory usage
#docker system prune

#Train the classifier:
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ../../merged/16s/ref-seqs16s.qza \
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
exit

docker run -it \
  -v ~/Dropbox/BMMA2025:/mnt/BMMA2025 \
  qiime2/core

# visualize. navigate back to where you have your taxonomy files
cd ../mnt/BMMA2025/data/qiime2_data/training_classifiers/16s_classifier

qiime metadata tabulate \
  --m-input-file ../../merged/16s/taxonomy16s.qza \
  --o-visualization ../../merged/16s/taxonomy16s.qzv


#cant do this step without a master map, which we don't have for the merged files. 
#skip for now. we still wanna do the other steps on the merged files
# cd ../../
# qiime taxa barplot --i-table merged/16s/table16s.qza --i-taxonomy merged/16s/taxonomy16s.qza --m-metadata-file merged/16s/maps/sky2018map16s.txt --o-visualization merged/16s/taxa-bar-plots.qzv


### ************************************************************************
# 2. FILTERING STEPS (PART A)
### ************************************************************************

# Go through and filter out unwanted reads on the merged tables and repseqs.
# We can't do all the filtering steps on the merged files, 
# but we can do some, then we make a phylogenetic tree, then go back and filter round-specific issues

# 2a. Prefilter step: The silva database does not have Lactobacillus michnerii, so our reads for Lactobacillus kunkeei are likely wrong. 

# first visualize taxonomy16s.qzv. from the visualizer, download taxonomy16s.tsv. 
# make a copy of taxonomy16s.tsv and rename it taxonomy16sfixed.tsv - we will make corrections here
# then visualize repseqs16sfiltered.qzv in qiime2view

## unmerged
# qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv


qiime feature-table tabulate-seqs --i-data rep-seqs-16s.qza --o-visualization rep-seqs-16s.qzv


# using these two documents, select the sequence corresponding with the feature ids for L. kunkeei

# Make corrections to taxonomy16sfixed.tsv. convert it to txt manually on your computer.
# Convert this text file back into a qza object using qiime

# we didnt have to end making these modifications in sky islands or SI 2019 because we didnt have L. kunkeei or michinerii. but left this for future projects
# if you have to make these changes,  make sure to change the paths below appropraitely. we kept everything labeled taxonomy16s.qza, NOT fixed


qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy16s.qza/taxonomy16sfixed.txt \
--output-path taxonomy16s.qza/taxonomy16sfixed.qza


#2b. filter 1: out the chloroplast and mitochondria reads 

qiime taxa filter-table --i-table table16s.qza --i-taxonomy taxonomy16s.qza --p-exclude mitochondria,chloroplast --o-filtered-table tablefilt1.qza
qiime feature-table summarize --i-table tablefilt1.qza --o-visualization tablefilt1.qzv 


#2c. filter 2: remove sequences only found in one sample 

qiime feature-table filter-features --i-table tablefilt1.qza --p-min-samples 2 --o-filtered-table tablefilt2.qza
qiime feature-table summarize --i-table tablefilt2.qza --o-visualization tablefilt2.qzv


#2d: Filter our rep seqs file so that you can see what samples you have left after filtering and subsampling


qiime feature-table filter-seqs --i-data rep-seqs-16s.qza --i-table tablefilt2.qza --o-filtered-data rep-seqs-16s-filtered.qza
qiime feature-table tabulate-seqs --i-data rep-seqs-16s-filtered.qza --o-visualization rep-seqs-16s-filtered.qzv
#view rep-seqs-16s-filtered.qzv in qiime2 viewer in browser


## here is a step we are changing for 2023 -- since we are using the trees to determine microbial phylogenetic distance, we need to 
## filter out unwanted sequences before we generate the trees.
## also, we will not split sequences into runs and instead will filter out contaminants across all runs



# 2e. Filter out large taxonomic groups of bacteria that are known to be contaminants

# there are known bacteria are salt-loving and often found in buffers. 
# they include  Halomonas, Shewanella, Oceanospirillales, and the acne bacteria Propionibacterium. 


qiime taxa filter-table --i-table table16s.qza --i-taxonomy taxonomy16s.qza --p-exclude halomonas,lautropia,rothia,haemophilus,desulfuromonadales,pseudomonas,devosia,sphingomonas,solirubrobacterales,escherichia-shigella,staphylococcus,prevotellaceae,blastococcus,aeromonas,streptococcus,shewanella,oceanospirillales,propionibacteriales --o-filtered-table tablef1.qza

# 2f. Let's specifically look at sequences in our controls and remove the bacteria that are in them


qiime taxa barplot --i-table tablef1.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/combined-map-2018-2021.txt --o-visualization taxa-bar-plots-f1.qzv

# now filter out taxa that are in the controls!
# manually selected from the table16s file which ASVs were found in controls and wrote them down here

#2018 controls
qiime taxa filter-table --i-table tablef1.qza --i-taxonomy taxonomy16s.qza --p-mode exact --p-exclude "d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Microbacteriaceae;g__Galbitalea",\
"d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae",\
"d__Bacteria;p__Tenericutes;c__Mollicutes;o__Entomoplasmatales;f__Spiroplasmataceae;g__Spiroplasma",\
"d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Kineosporiales;f__Kineosporiaceae;g__Kineosporia;s__uncultured;bacterium",\
"d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Chitinophagales;f__Chitinophagaceae;g__Segetibacter",\
"d__Bacteria",\
"Unassigned" --o-filtered-table tablef2.qza

qiime taxa barplot --i-table tablef2.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/combined-map-2018-2021.txt --o-visualization taxa-bar-plots-f2.qzv

#2021 controls
qiime taxa filter-table --i-table tablef2.qza --i-taxonomy taxonomy16s.qza --p-mode exact --p-exclude "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium;__",\
"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Comamonadaceae",\
"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Acidovorax",\
"d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium",\
"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Ralstonia",\
"d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus",\
"d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Caulobacterales;f__Caulobacteraceae;g__Phenylobacterium",\
"Unassigned" --o-filtered-table tablef3.qza

qiime taxa barplot --i-table tablef3.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/combined-map-2018-2021.txt --o-visualization taxa-bar-plots-f3.qzv

#Filter out controls by making new map file excluding controls
qiime feature-table filter-samples --i-table tablef3.qza --m-metadata-file maps/combined-map-2018-2021-noCtrl.txt --o-filtered-table tablef4.qza

#now need to filter the rep-seqs FeatureTable[sequence] to the updated fitered FeatureTable[frequency]

qiime feature-table filter-seqs \
--i-data rep-seqs-16s-filtered.qza \
--i-table tablef4.qza \
--o-filtered-data rep-seqs-16s-final-filter.qza
  
  
### ************************************************************************
## 3. #GENERATE TREE FOR PHYLOGENETIC DIVERSITY ANALYSES
### ************************************************************************

# generate a tree with the merged data

#Now alignment of the reads using MAFFT.
qiime alignment mafft --i-sequences rep-seqs-16s-final-filter.qza --o-alignment aligned_repseqs16s.qza

#Then filter out the unconserved, highly gapped columns from alignment
qiime alignment mask --i-alignment aligned_repseqs16s.qza --o-masked-alignment masked_aligned_repseqs16s.qza

#And make a tree with Fasttree which creates a maximum likelihood phylogenetic trees from aligned sequences
#for more information on fastree, http://www.microbesonline.org/fasttree/
qiime phylogeny fasttree --i-alignment masked_aligned_repseqs16s.qza --o-tree unrooted_tree16s.qza

#now root it
qiime phylogeny midpoint-root --i-tree unrooted_tree16s.qza --o-rooted-tree rooted-tree16s.qza


qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-16s-final-filter.qza --o-alignment aligned_repseqs16s.qza --o-masked-alignment masked_aligned_repseqs16s.qza --o-tree unrooted_tree16s.qza --o-rooted-tree rooted-tree16s.qza


## now to fix the error: The table does not appear to be completely represented by the phylogeny.
#need to drop from the tree faetures not in the table anymore after filtering

qiime fragment-insertion filter-features \
--i-table tablef4.qza \
--i-tree rooted-tree16s.qza \
--o-filtered-table tablef5.qza \
--o-removed-table removed-tablef5.qza
### *************************************************************************
#4. DETERMINE SUBSAMPLE DEPTH
### *************************************************************************

#We want to make a rarefaction curve to see how subsampling depths influence our alpha diversity metrics.
#Open visualization in Qiime2 for tablef5. View interactive sample detail and look at visualized table to see impact of different depths
# record number of reads, number of samples at specified sampling depth, and how many features retained

#For subsampling:  840 reads.  715 (90.62%) samples at the specifed sampling depth.Retained 600,600 (5.71%) features

qiime diversity alpha-rarefaction --i-table tablef5.qza --i-phylogeny rooted-tree16s.qza --p-max-depth 10000 --m-metadata-file maps/combined-map-2018-2021.txt --o-visualization alphararefact16s.qzv

qiime feature-table summarize --i-table tablef5.qza --o-visualization tablef5.qzv





### *************************************************************************
# 5. ALPHA AND BETA DIVERSITY
### *************************************************************************

#Generate core metrics folder. This command makes the new directory and outputs a bunch of files in it
# make sured youre in "merged/16s"

mkdir final

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree16s.qza --i-table tablef5.qza --p-sampling-depth 840 --m-metadata-file maps/combined-map-2018-2021.txt --output-dir final/core_metrics16s --verbose


#if you want to specify the number of parallel cores running, add this to the code
--p-n-jobs 8

#now we have lots of files in that core_metrics directory. gotta export rarefied_table.qza so we can get a qzv and convert to a csv
#using combined full map
qiime taxa barplot --i-table final/core_metrics16s/rarefied_table.qza --i-taxonomy taxonomy16s.qza --m-metadata-file maps/combined-map-2018-2021-noCtrl.txt --o-visualization final/core_metrics16s/rarefiedtable16s.qzv


#drag the new qzv into qiime 2 view, set taxonomic level to 7, and download a csv! put this into "final asv tables" folder

