
######################################################################################## Differential Expression Analysis using edgeR package ########################################################################################


######################################## refine each file containing only the circRNA ID and raw count and save to a new folder called DE_analysis_edgeR ########################################
names = list.files(path = "/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/circRNA_outputs/refined_bed_files/", pattern = "*bed", all.files = FALSE, full.names = FALSE);
names;

j = 1;
while ( j <= 52)
{
  setwd('/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/circRNA_outputs/refined_bed_files/');
  data = read.table(file.path(names[j]), header = TRUE);
  
  data = data[, c(5,9)];
  
  colnames(data) = c("ID", "counts");
  data$counts = as.numeric(data$counts);
  
  setwd('/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/circRNA_outputs/DE_analysis_edgeR/');
  write.table(data, names[j] , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE);
  j = j + 1;
}


######################################## read files into edgeR using readDGE ######################################## 
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR);

setwd('/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/circRNA_outputs/DE_analysis_edgeR/');
names = list.files(path = "/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/circRNA_outputs/DE_analysis_edgeR/", pattern = "*bed", all.files = FALSE, full.names = FALSE);
names;



files = names;
rawdata = readDGE(files);
y = DGEList(rawdata);

names = gsub('.{17}$','',names);
rownames(y$samples) = names;
colnames(y$counts) = names;

y$samples$group = c(rep("ASD",13), rep("CTRL",39));
View(y$samples);
View(y$counts);

write.table(y$counts, file = "merged_circRNAs_counts", append = FALSE, quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE);
write.table(y$samples, "total_number_of_reads_for_each_sample", col.names = TRUE, row.names = TRUE, quote = FALSE);

######################################## filter out the circRNAs with total read number in all samples less than 16 and found in less than 3 samples (24,417 -> 294)##################

# keep = rowSums(y$counts > 2) >= 4; # 2427 remained, none DE circRNAs
# keep = rowSums(y$counts) >= 10;  # 5105 remained, 43 down regulated, 0 up regulated
# keep = rowSums(y$counts) >= 50;  # 1214 remained,2 down, 0 up;
# keep = rowSums(y$counts) >= 20; # 2735 remained, 8 down, 0 up;
# keep = rowSums(y$counts) >= 15; #3473 remained, 12 down, 0 up;
keep = rowSums(cpm(y)>1) >= 3;  # 46108 to 15874 remained, 376 down, 158 up; this is the used results
table(keep)

y = y[keep, , keep.lib.sizes = FALSE];

write.table(y$counts, file = "merged_circRNAs_counts_filtered_by_cpm", append = FALSE, quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE);
write.table(y$samples, "total_number_of_reads_for_each_sample_filtered_by_cpm", col.names = TRUE, row.names = TRUE, quote = FALSE);

######################################## TMM normalization (trimmed mean of M values) to account for compositional difference between the libraries ###############################

y = calcNormFactors(y);
y$samples

y$samples$lib.size*y$samples$norm.factors;
################################## examine the library relationship ##################################

jpeg("samples_distribution")
plotMDS(y);
dev.off();
########################### MD plot for each library, using first 1 as an example here ###########################

jpeg("Sample_distribution_MD")
plotMD(y, column = 1);
abline(h = 0, col = "red", lty = 2, lwd = 2);
dev.off();

################################## Design matrix: tell R which sample belongs to which group, only the order matters #################################

design = model.matrix(~0+group, data = y$samples);
design;
rownames(design) = colnames(y);

################################# estimating the dispersion: the variability between biological replicates: #################################
install.packages("statmod");
library(statmod);

y = estimateDisp(y, design,robust = TRUE);
y$common.dispersion;
jpeg("Sample_Variations_BCV");
plotBCV(y);
dev.off();


################################################################# differential expression analysis #################################################################
################################ GLM test method : 
fit = glmFit(y, design);
lrt = glmLRT(fit, contrast = c(1,-1));
result_lrt = topTags(lrt, n = 20);
View(result_lrt);
summary(decideTests(lrt));
View(lrt$table);

plotMD(lrt);
abline( h = c(-1,1), col = "blue");

################################ QL F-test method:
fit = glmQLFit(y, design, robust = TRUE);
plotQLDisp(fit);
qlf = glmQLFTest(fit, contrast = c(1,-1));
result_qlf = topTags(qlf, n = 534);
View(result_qlf);
summary(decideTests(qlf));

jpeg('MD_plot_qlf_test');
plotMD(qlf);
abline( h = c(-1, 1), col = "blue");
dev.off();

################################ save data to the working directory "/media/jun/circRNAs/RNAseq/DE_analysis_edgeR" : 
write.table(result_qlf, file = "DE_result_qlf_test_filter_cpm1_3", append = FALSE, quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE);
write.table(result_lrt, file = "DE_result_lrt_test", col.names = TRUE, row.names = TRUE, quote = FALSE);
