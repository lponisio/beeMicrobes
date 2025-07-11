#!/usr/bin/env bash 

# Bioinformatics pipeline code for processing 16s and RBCL reads from fastq files. 

# If you have the demultiplexed sequences (prefix_demux.qza) start from step X

##########################################
### 1. Run docker pull for QIIME 1 & 2 ###
##########################################

#1: Download docker on your computer, allows you to run qiime and
# other software without downloading them.

# On whatever computer will be running the tasks, innitiate a docker
# container for qiime1

#only need to do this docker pull step once per user
### FOR QIIME1 ###
docker pull mbari/qiime1

### FOR QIIME2 ###
docker pull quay.io/qiime2/core:2024.10

#Helpful notes:

## 1. if running on the shared lab computer, you may get errors if
## your zipped sequence results are not synched to your dropbox
## account. Double clicking the file will usually make it appear. If
## you look at the filesize and it is 0, this is likely because the
## file is not synched to the local.

## 2. If you get an error about memory or space issues, check that the
##     other users do not have their dropbox files on the local
##     computer, you can fix this by logging in and making all dropbox
##     files online only. This should clear up enough space to fix
##     this error.

## Use docker to open a container for qiime1
## Mount the ffar pipeline folder into the qiime1 container

#######################################
### 2. Mount directory to container ###
#######################################

docker run -it \
  -v ~/Dropbox/beeMicrobes_saved/beeMicrobes_pipeline_output:/mnt/beeMicrobes_pipeline_output \
  mbari/qiime1

# move your directory into the mounted folder
cd ./mnt/beeMicrobes_pipeline_output

# Check that all of the appropriate files/folders are in the mounted directory
ls

# boot up qiime1 within the container
source activate qiime1 

## All of the barcode parsing and demuxing was done in the indivudal
## project pipelines so we are starting after the demuxing. 

##############################
### 6. Demultiplex samples ###
##############################

cd ./data/qiime2_data/demux_files

#9: Demultiplex 16s reads first. Only works in version Qiime2 2019.1
#note: this step takes ~2 hours!

### DEMULTIPLEX SKY ISLANDS SAMPLES ###

## 16S DEMUX - SKY ISLANDS ##
# !Just one example presented here. Modify code as necessary if
# !needing to demultiplex other sequences!

# Run 1
qiime demux emp-paired --i-seqs seqs.qza \
                       --m-barcodes-file maps/sky2018map16s.txt \
                       --m-barcodes-column barcodesequence \
                       --o-per-sample-sequences demux16s.qza 

## Visualize demux seqs

### 16s QZV
ls -lh SF_R0_demux16s.qza
## SF R0
qiime demux summarize --i-data /mnt/BMMA2025/data/qiime2_data/demux_files/SF_R0_demux16s.qza \
                      --o-visualization ../qzv_files/SF_R0_demux16s.qzv
qiime tools view SF_R0_demux16s.qzv

## SF R1
qiime demux summarize --i-data SF_R1_demux16s.qza \
                      --o-visualization ../qzv_files/SF_R1_demux16s.qzv
qiime tools view SF_R1_demux16s.qzv

## SF R2
qiime demux summarize --i-data SF_R2_demux16s.qza \
                      --o-visualization ../qzv_files/SF_R2_demux16s.qzv
qiime tools view SF_R2_demux16s.qzv

## SF R3
qiime demux summarize --i-data SF_R3_demux16s.qza \
                      --o-visualization ../qzv_files/SF_R3_demux16s.qzv
qiime tools view SF_R3_demux16s.qzv

## SF R4
qiime demux summarize --i-data SF_R4_demux16s.qza \
                      --o-visualization SF_R4_demux16s.qzv
qiime tools view SF_R4_demux16s.qzv

## RBCL QZV

## SF R0
qiime demux summarize --i-data SF_R0_demuxRBCL.qza \
                      --o-visualization ../qzv_files/SF_R0_demuxRBCL.qzv
qiime tools view SF_R0_demuxRBCL.qzv

## SF R1
qiime demux summarize --i-data SF_R1_demuxRBCL.qza \
                      --o-visualization ../qzv_files/SF_R1_demuxRBCL.qzv
qiime tools view SF_R1_demuxRBCL.qzv

## SF R2
qiime demux summarize --i-data SF_R2_demuxRBCL.qza \
                      --o-visualization ../qzv_files/SF_R2_demuxRBCL.qzv
qiime tools view SF_R2_demuxRBCL.qzv

## SF R3
qiime demux summarize --i-data SF_R3_demuxRBCL.qza \
                      --o-visualization ../qzv_files/SF_R3_demuxRBCL.qzv
qiime tools view SF_R3_demuxRBCL.qzv

## SF R4
qiime demux summarize --i-data SF_R4_demuxRBCL.qza \
                      --o-visualization SF_R4_demuxRBCL.qzv
qiime tools view SF_R4_demuxRBCL.qzv


##################################
### 7. Cluster reads into ASVs ###
##################################

## click the tab Interactive Quality from the main visualization tab

## when interpretating the quality boxes, you can use the bottom of
## the black box as a conservative measure for the phred score (not the
## whiskers and not the middle of the box). Quality score around 30 is
## acceptable, lower not so good. Select the number before the
## quality scores begin to decline. 

# Note: this step takes hours!

### 16S Denoising (ASV Assignment) ###

#############################
### Sky Islands Denoising ###
#############################

## 16s Denoising

## Run 1 2018
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2018_demux16s.qza  \
  --p-trunc-len-f 180  \
  --p-trunc-len-r 220  \
  --p-trim-left-f 0  \
  --p-n-threads 0  \
  --output-dir dada2-16s  \
  --o-representative-sequences ../rep_seqs/SI_R2018_rep-seqs16s.qza  \
  --o-table ../feature_tables/SI_R2018_table16s.qza \
  --o-denoising-stats ../denoising_stats/SI_R2018_denoising-stats16s.qza


qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## Run 2 2023 lane1
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2023_L1_demux16s.qza  \
  --p-trunc-len-f 205  \
  --p-trunc-len-r 237  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences ../rep_seqs/SI_R2023_L1_rep-seqs16s.qza  \
  --o-table ../feature_tables/SI_R2023_L1_table16s.qza \
  --o-denoising-stats ../denoising_stats/SI_R2023_L1_denoising-stats16s.qza
  
qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## Run 3 2023 lane2
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2023_L2_demux16s.qza  \
  --p-trunc-len-f 147  \
  --p-trunc-len-r 192  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SI_R2023_L2_rep-seqs16s.qza  \
  --o-table feature_tables/SI_R2023_L2_table16s.qza \
  --o-denoising-stats denoising_stats/SI_R2023_L2_denoising-stats16s.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## repeat the steps for RBCL

## 2018 RBCL
cd R2023/lane2
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2018_demuxRBCL.qza  \
  --p-trunc-len-f 180  \
  --p-trunc-len-r 218  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SI_R2018_rep-seqsRBCL.qza  \
  --o-table feature_tables/SI_R2018_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SI_R2018_denoising-statsRBCL.qza
 
## 2023 Lane 1 RBCL
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2023_L1_demuxRBCL.qza  \
  --p-trunc-len-f 232  \
  --p-trunc-len-r 225  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SI_R2023_L1_rep-seqsRBCL.qza  \
  --o-table feature_tables/SI_R2023_L1_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SI_R2023_L1_denoising-statsRBCL.qza
 
## 2023 Lane 2 RBCL
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SI_R2023_L1_demuxRBCL.qza  \
  --p-trunc-len-f 138  \
  --p-trunc-len-r 192  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SI_R2023_L2_rep-seqsRBCL.qza  \
  --o-table feature_tables/SI_R2023_L2_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SI_R2023_L1_denoising-statsRBCL.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

###########################
### Sunflower Denoising ###
###########################

## SF 16S Run 0
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R0_demux16s.qza  \
  --p-trunc-len-f 256  \
  --p-trunc-len-r 261  \
  --p-trim-left-f 17  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SF_R0_rep-seqs16s.qza  \
  --o-table feature_tables/SF_R0_table16s.qza \
  --o-denoising-stats denoising_stats/SF_R0_denoising-stats16s.qza

## SF 16S Run 1
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R1_demux16s.qza  \
  --p-trunc-len-f 150  \
  --p-trunc-len-r 225  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SF_R1_rep-seqs16s.qza  \
  --o-table feature_tables/SF_R1_table16s.qza \
  --o-denoising-stats denoising_stats/SF_R1_denoising-stats16s.qza
 
 ## SF 16S Run 2
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demux16s.qza  \
  --p-trunc-len-f 152  \
  --p-trunc-len-r 238  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SF_R2_rep-seqs16s.qza  \
  --o-table feature_tables/SF_R2_table16s.qza \
  --o-denoising-stats denoising_stats/SF_R2_denoising-stats16s.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## SF 16S Run 3
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R3_demux16s.qza  \
  --p-trunc-len-f 207  \
  --p-trunc-len-r 226  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SF_R3_rep-seqs16s.qza  \
  --o-table feature_tables/SF_R3_table16s.qza \
  --o-denoising-stats denoising_stats/SF_R3_denoising-stats16s.qza
 
## SF 16S Run 4
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R4_demux16s.qza  \
  --p-trunc-len-f 233  \
  --p-trunc-len-r 254  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/SF_R4_rep-seqs16s.qza  \
  --o-table feature_tables/SF_R4_table16s.qza \
  --o-denoising-stats denoising_stats/SF_R4_denoising-stats16s.qza

qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv

## Repeat for RBCL

## SF RBCL Run 0
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R0_demuxRBCL.qza  \
  --p-trunc-len-f 134  \
  --p-trunc-len-r 176  \
  --p-trim-left-f 5  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/SF_R0_rep-seqsRBCL.qza  \
  --o-table feature_tables/SF_R0_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SF_R0_denoising-statsRBCL.qza

## SF RBCL Run 1
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demuxRBCL.qza  \
  --p-trunc-len-f 178  \
  --p-trunc-len-r 178  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/SF_R1_rep-seqsRBCL.qza  \
  --o-table feature_tables/SF_R1_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SF_R1_denoising-statsRBCL.qza
 
## SF RBCL Run 2
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demuxRBCL.qza  \
  --p-trunc-len-f 150  \
  --p-trunc-len-r 180  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/SF_R2_rep-seqsRBCL.qza  \
  --o-table feature_tables/SF_R2_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SF_R2_denoising-statsRBCL.qza

qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCL.qzv
qiime feature-table summarize --i-table dada2-RBCL/table16s.qza --o-visualization dada2-RBCL/table16s.qzv

## RBCL Run 3
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demuxRBCL.qza  \
  --p-trunc-len-f 178  \
  --p-trunc-len-r 178  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/SF_R3_rep-seqsRBCL.qza  \
  --o-table feature_tables/SF_R3_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SF_R3_denoising-statsRBCL.qza
 
## RBCL Run 4
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demuxRBCL.qza  \
  --p-trunc-len-f 178  \
  --p-trunc-len-r 178  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/SF_R3_rep-seqsRBCL.qza  \
  --o-table feature_tables/SF_R3_tableRBCL.qza \
  --o-denoising-stats denoising_stats/SF_R4_denoising-statsRBCL.qza

qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCL.qzv
qiime feature-table summarize --i-table dada2-RBCL/table16s.qza --o-visualization dada2-RBCL/table16s.qzv

###########################
### OR Survey Denoising ###
###########################

## 16s Denoising

## OR 16S Run 2023
# modify appropriatley
cd theCorrectFilePath
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demux16s.qza  \
  --p-trunc-len-f 156  \
  --p-trunc-len-r 213  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/OR_R2023_rep-seqs16s.qza  \
  --o-table feature_tables/OR_R2023_table16s.qza \
  --o-denoising-stats denoising_stats/OR_R2023_denoising-stats16s.qza

## OR 16S Run 2024
# modify appropriatley
cd theCorrectFilePath
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demux16s.qza  \
  --p-trunc-len-f numberHere  \
  --p-trunc-len-r numberHere  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-16s  \
  --o-representative-sequences rep_seqs/OR_R2024_rep-seqs16s.qza  \
  --o-table feature_tables/OR_R2024_table16s.qza \
  --o-denoising-stats denoising_stats/OR_R2024_denoising-stats16s.qza
 
## Repeat for RBCL

## RBCL Run 1
# modify appropriatley
cd ~/Dropbox/BMMA2025/data/qiime2_data
qiime dada2 denoise-paired  \
  --i-demultiplexed-seqs SF_R2_demuxRBCL.qza  \
  --p-trunc-len-f numberHere  \
  --p-trunc-len-r numberHere  \
  --p-trim-left-f 0  \
  --p-n-threads 2  \
  --output-dir dada2-RBCL  \
  --o-representative-sequences rep_seqs/OR_R2023_rep-seqsRBCL.qza  \
  --o-table feature_tables/OR_R2023_tableRBCL.qza \
  --o-denoising-stats denoising_stats/OR_R2023_denoising-statsRBCL.qza

# check outputs to make sure you didn't lose too many samples. 
# We found being more conservative and doing shorter truncations gives you the same number of sequences, 
# but sorted into fewer features, likely cause trimmed off poor quality reads

## *****************************************************************************
##       MERGE files from runs ONCE THEY EXIST
## *****************************************************************************
# time to merge the files from your different runs.
# NOTE: much of this comes from https://john-quensen.com/tutorials/merging-dada2-results-in-qiime2/ 

# cd back to your main folder (in this case the qiime2_data folder)

## make a new merge directory, and separate 16s and RBCL directories within it
cd ../
mkdir merged
cd merged
mkdir 16s
mkdir RBCL
cd ../

#### 1. 16s ###
# 1a. first merge the table files. do this from your main SI_pipeline folder
qiime feature-table merge \
      --i-tables feature_tables/SI_R2018_table16s.qza \
      --i-tables feature_tables/SI_R2023_L1_table16s.qza \
      --i-tables feature_tables/SI_R2023_L2_table16s.qza \
      --i-tables feature_tables/SF_R0_table16s.qza \
      --i-tables feature_tables/SF_R1_table16s.qza \
      --i-tables feature_tables/SF_R2_table16s.qza \
      --i-tables feature_tables/SF_R3_table16s.qza \
      --i-tables feature_tables/SF_R4_table16s.qza \
      --i-tables feature_tables/OR_R2023_table16s.qza \
      --o-merged-table merged_data/16s/table16s.qza

# 1b: next merge the rep-seqs
qiime feature-table merge \
      --i-tables rep_seqs/SI_R2018_rep-seqs16s.qza \
      --i-tables rep_seqs/SI_R2023_L1_rep-seqs16s.qza \
      --i-tables rep_seqs/SI_R2023_L2_rep-seqs16s.qza \
      --i-tables rep_seqs/SF_R0_rep-seqs16s.qza \
      --i-tables rep_seqs/SF_R1_rep-seqs16s.qza \
      --i-tables rep_seqs/SF_R2_rep-seqs16s.qza \
      --i-tables rep_seqs/SF_R3_rep-seqs16s.qza \
      --i-tables rep_seqs/SF_R4_rep-seqs16s.qza \
      --i-tables rep_seqs/OR_R2023_rep-seqs16s.qza \
      --o-merged-table merged/16s/rep-seqs16s.qza

# 2a: merge the feature tables
qiime feature-table merge \
      --i-tables feature_tables/SI_R2018_tableRBCL.qza \
      --i-tables feature_tables/SI_R2023_L1_tableRBCL.qza \
      --i-tables feature_tables/SI_R2023_L2_tableRBCL.qza \
      --i-tables feature_tables/SF_R0_tableRBCL.qza \
      --i-tables feature_tables/SF_R1_tableRBCL.qza \
      --i-tables feature_tables/SF_R2_tableRBCL.qza \
      --i-tables feature_tables/SF_R3_tableRBCL.qza \
      --i-tables feature_tables/SF_R4_tableRBCL.qza \
      --i-tables feature_tables/OR_R2023_tableRBCL.qza \
      --o-merged-table merged/RBCL/tableRBCL.qza

# 2b: next merge the rep-seqs
qiime feature-table merge \
      --i-tables rep_seqs/SI_R2018_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SI_R2023_L1_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SI_R2023_L2_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SF_R0_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SF_R1_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SF_R2_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SF_R3_rep-seqsRBCL.qza \
      --i-tables rep_seqs/SF_R4_rep-seqsRBCL.qza \
      --i-tables rep_seqs/OR_R2023_rep-seqsRBCL.qza \
      --o-merged-table merged/RBCL/rep-seqsRBCL.qza

# End of script. Move to 2_0_pipeline_16s.sh to proceed with data prep
