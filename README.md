# Genomic_shock_LTF_BTF
Chapter 5:

The raw transcriptomic data for this chapter can be found in SharePoint folder, available upon request.

Our reference genome was assembled from Hi-C and PacBio HiFi reads following the VGP pipelines for these kinds of data in Galaxy Australia (https://usegalaxy.org.au). The curated assembly is found as "LTFr_final_renamed_hap1.fa". The gene annotation for the reference genome used the RNA from two tissues and the zebra finch proteome, which were uploaded to Galaxy. The resulting annotation is named "LTFr_annotation_braker3.gtf".

The TE annotation of the reference genome was carried out with EarlGrey and using the script "RepeatM_LTF.sh". The output was sent to Oliver Tam from the Molly Gale Hammell Laboratory (https://www.mghlab.org) to produce TE annotation files compatible with TEtranscripts and TElocal, "PoHecki.filtered_TE.gtf" and "PoHecki.filtered_TE.gtf.locInd", respectively.

Read alignment allowing for multimapping was carried out using the "Prep1_STAR_index_genome.sh" and "Prep2_STAR_alignment.sh" scripts.  
 
TEtranscripts and TElocal were run for the aligned BAM files using the "02_TEcount.sh" using their respective gene and TE annotation files (i.e. "LTFr_annotation_braker3.gtf" and "PoHecki.filtered_TE.gtf" for TEtranscripts, or "PoHecki.filtered_TE.gtf.locInd" for TElocal).

The resulting count files were collected for all samples in "Counts_final.csv" and "Counts_final_TElocal.csv"

The differential expression analysis and downstream analysis, including the production of the figures in the chapter for TEtranscripts was done with "TEtranscripts_script.R". The differential expression analysis and downstream analysis, including the production of figures for TElocal was done with "TElocal.R".

A complete repository of the PhD thesis has been created elsewhere since file size limitations make it hard to upload all data, and can be found here: https://drive.google.com/file/d/1W5KDuDYTOoX9cbLOTBuyc-lM4j5ymqMH/view?usp=sharing

For any additional questions: sergiogonmoll[at]gmail.com or h.l.dugdale[at]rug.nl
