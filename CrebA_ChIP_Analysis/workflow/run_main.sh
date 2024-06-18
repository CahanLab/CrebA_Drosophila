
##### this is the standard format for the chip-seq data #####
# chrom - Name of the chromosome (or contig, scaffold, etc.).
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# name - Name given to a region (preferably unique). Use "." if no name is assigned.
# score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were "'0"' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
# strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
# signalValue - Measurement of overall (usually, average) enrichment for the region.
# pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
# qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
# peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

##### let's build the dictionary for the flybase id ######
python build_flybase_dictionary.py

##### running FIMO on all the binding regions -- CrebA motifs #####
# curate the chen format to meme format 
python chen_to_meme.py 

##### get the promoter regions for all the genes ######
python get_promoter.py 

##### get the intersecting region for all 3 tracks ######
mkdir ../output/intersect_region_CrebA
bedtools intersect -a ../input/GEO_processed/GSM4213092_oregon_CrebA_peaks.narrowPeak \
                    -b ../input/GEO_processed/GSM4213094_fkh_CrebA_peaks.narrowPeak \
                    ../input/GEO_processed/GSM4213096_sage_CrebA_peaks.narrowPeak > ../output/intersect_region_CrebA/intersect_regions.narrowPeak

##### get the nearest gene that is bounded #####
python match_nearest_region_to_gene.py 

##### get the raw sequences ######
sed 's/^chr//; s/%$//' ../output/pruned_intersect_peaks/CrebA_fkh_sage_intersect_250.narrowPeak | sed 's/_CP007111v1_random//g' > ../output/pruned_intersect_peaks/chrom_removed.narrowPeak
bedtools getfasta -fi ../input/reference_genome/dmel-all-chromosome-r6.33.fasta \
                  -bed ../output/pruned_intersect_peaks/chrom_removed.narrowPeak \
                  > ../output/pruned_intersect_peaks/peak_250_sequence_regions.fasta

##### get the FIMO analysis of these regions ######

