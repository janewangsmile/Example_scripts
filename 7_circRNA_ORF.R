#### read the .fa file into the desiable format:
library(Biostrings)
faFile <- readDNAStringSet("Desktop/m6A_model_validataion_data/circRNA/circRNADb/circRNA_fasta")
seq_name <- names(faFile);
sequence <- paste(faFile);
circRNAs <- data.frame(seq_name, sequence);
#### install package for ORF predictor:
source("https://bioconductor.org/biocLite.R")
biocLite("systemPipeR")
library(systemPipeR)
#### predict the ORF for each circRNA:
circRNA_ORF=predORF(faFile, n = 1, type = "df", strand = "sense")   #### circRNAs with no ORFs are not included in circRNA_ORF 
write.csv(circRNA_ORF,"Desktop/m6A_model_validataion_data/circRNA/circRNA_ORF.csv",row.names = FALSE)
circRNA_ORF = read.csv("Desktop/m6A_model_validataion_data/circRNA/circRNA_ORF.csv",header = TRUE)
#### find all the DRACH sites in each circRNA
circRNAs$DRACH_all = NA;
i = 1;
while ( i <= nrow(circRNAs))
{
  circRNAs$DRACH_all[i] = gregexpr("[AGT][AG]AC[ACT]", circRNAs$sequence[i]);
  i = i + 1;
}
#### put each DRACH site in different rows:
circRNAsDRACHsites <- circRNAs;
circRNAsDRACHsites <- circRNAsDRACHsites[,-3];
circRNAsDRACHsites$DRACH_site = NA;

i = 1; # circRNA row
n = 1; # DRACH site row
while ( i <= nrow(circRNAs))
{
  j = 1; # DRACH in each lncRNA
  all_DRACH = unlist(circRNAs$DRACH_all[i]);
  while ( j <= length(all_DRACH))
  {
    circRNAsDRACHsites[n,c(1,2)] = circRNAs[i,c(1,2)];
    circRNAsDRACHsites[n,3] = all_DRACH[j];
    n = n + 1;
    j = j + 1;
  }
  i = i + 1;
}
#### find the Adenosine site 
circRNAsDRACHsites$A_site = NA;
circRNAsDRACHsites$A_site = circRNAsDRACHsites$DRACH_site+2;
write.csv(circRNAsDRACHsites,"Desktop/m6A_model_validataion_data/circRNA/all_DRACH_sites_in_all_circRNAs_nofilterd.csv",row.names = FALSE)
#### feature engineering:
#### remove the A site with flanking sequence <100-bp
circRNAsDRACHsites$sequence = as.character(circRNAsDRACHsites$sequence);
circRNAsDRACHsites$seq_length = nchar(circRNAsDRACHsites$sequence)
circRNAsDRACHsites = circRNAsDRACHsites[circRNAsDRACHsites$A_site>=101, ]
circRNAsDRACHsites = circRNAsDRACHsites[circRNAsDRACHsites$A_site <= (circRNAsDRACHsites$seq_length -100), ]
#### extract the 100nt*2 +1 sequences
circRNAsDRACHsites$seq_201 = subseq(circRNAsDRACHsites$sequence, circRNAsDRACHsites$A_site-100, circRNAsDRACHsites$A_site + 100)
min(nchar(circRNAsDRACHsites$seq_201))
#### check all the seq is ATGC no N:
"N" %in% circRNAsDRACHsites$seq_201
#### one-hot encoding:
circRNAsDRACHsites$encoded_seq = as.character("X");

i = 1;
while ( i <= nrow(circRNAsDRACHsites))
{
  code_seq_A = gsub("A", "1000", circRNAsDRACHsites$seq_201[i]);
  code_seq_T = gsub("T", "0100", code_seq_A);
  code_seq_G = gsub("G", "0010", code_seq_T);
  code_seq_C = gsub("C", "0001", code_seq_G);
  circRNAsDRACHsites$encoded_seq[i] = as.character(code_seq_C);
  i = i + 1;
}

min(nchar(circRNAsDRACHsites$encoded_seq))
write.csv(circRNAsDRACHsites,"Desktop/m6A_model_validataion_data/circRNA/filtered_circRNAsDRACHsites_one_hot_encoded.csv",row.names = FALSE)
circRNAsDRACHsites = read.csv("Desktop/m6A_model_validataion_data/circRNA/filtered_circRNAsDRACHsites_one_hot_encoded.csv",header = TRUE)

#### converted to 3D array for deep learning prediction:
circRNA_array = array(data = 0, dim = c(4,201, nrow(circRNAsDRACHsites)));
dim(circRNA_array)
circRNAsDRACHsites$encoded_seq = as.character(circRNAsDRACHsites$encoded_seq)

s = 1;
r = 1;
c = 1;
j = 1;
while (s <= nrow(circRNAsDRACHsites) )
{
  c = 1;
  while ( c <= 201)
  {
    r = 1;
    while (r <= 4)
    {
      circRNA_array[r,c,s] = as.integer(substr(circRNAsDRACHsites$encoded_seq[s],j,j))
      r = r + 1;
      j = j + 1;
    }
    c = c + 1;
  }
  s = s+1;
  j = 1;
}

dim(circRNA_array)
any(is.na(circRNA_array))
#### change the index order to the ones keras package can use: 
circRNA_array_new = aperm(circRNA_array, c(3,2,1));
dim(circRNA_array_new)
any(is.na(circRNA_array_new))
saveRDS(circRNA_array_new,"Desktop/m6A_model_validataion_data/circRNA/circRNA_array_data.rds")
#### predict with deep learning model in python...
#### get the predicted results
circRNAs_result = read.table("Desktop/m6A_model_validataion_data/circRNA/predict_result.txt", sep = "\n")
circRNAs_result = read.csv("Desktop/m6A_model_validataion_data/circRNA/predict_result_prob.csv",header = TRUE);

#### combine with the circRNAs all feature dataframe
circRNAsDRACHsites$predict_result = NA;
circRNAsDRACHsites$predict_result = circRNAs_result;

circRNAsDRACHsites = cbind(circRNAsDRACHsites,circRNAs_result)

#### remove useless column:
circRNAsDRACHsites = circRNAsDRACHsites[,-c(6,7)]

#### combine with ORF information 
circRNAsDRACHsites$ORF_start = NA;
circRNAsDRACHsites$ORF_end = NA;
circRNAsDRACHsites$ORF_length= NA;

circRNAsDRACHsites$seq_name = as.character(circRNAsDRACHsites$seq_name)
circRNA_ORF$seqnames = as.character(circRNA_ORF$seqnames)

i = 1; # circRNAs_A_sites row
while ( i <= nrow(circRNAsDRACHsites))
{
  r = NA;
  r = which(circRNA_ORF$seqnames == as.character(circRNAsDRACHsites$seq_name[i]))
  if( length(r) > 0)
  {
    circRNAsDRACHsites$ORF_start[i] = circRNA_ORF$start[r];
    circRNAsDRACHsites$ORF_end[i] = circRNA_ORF$end[r];
    circRNAsDRACHsites$ORF_length[i] = circRNA_ORF$width[r]
  }
  i = i + 1
}

#### set a filter for ORF length: 150nt
#### remove rows with no ORFs:
filtered_circRNAsDRACHsites = circRNAsDRACHsites[!is.na(circRNAsDRACHsites$ORF_length), ]
write.csv(filtered_circRNAsDRACHsites,"Desktop/m6A_model_validataion_data/circRNA/filtered_circRNAsDRACHsites_withORFs.csv",row.names = FALSE);

circRNAs_150_ORF_positive = filtered_circRNAsDRACHsites[(filtered_circRNAsDRACHsites$ORF_length >= 150 & filtered_circRNAsDRACHsites[,6] == 1), ]
circRNAs_150_ORF_negative = filtered_circRNAsDRACHsites[(filtered_circRNAsDRACHsites$ORF_length >= 150 & filtered_circRNAsDRACHsites[,6] == 0), ]

circRNAs_150_ORF_positive = filtered_circRNAsDRACHsites[(filtered_circRNAsDRACHsites$ORF_length >= 150 & filtered_circRNAsDRACHsites[,7] >= 0.8), ]

####select the positive m6A before ORF_start site, then circRNAs with translatability
circRNAs_150_ORF_positive = filtered_circRNAsDRACHsites[(filtered_circRNAsDRACHsites$ORF_length >= 150 & filtered_circRNAsDRACHsites[,7] > 0.5), ]
circRNAs_150_ORF_positive_before_ORF = circRNAs_150_ORF_positive[circRNAs_150_ORF_positive$A_site < circRNAs_150_ORF_positive$ORF_start, ]

circRNAs_150_ORF_positive_before_ORF_100nt = circRNAs_150_ORF_positive_before_ORF[circRNAs_150_ORF_positive_before_ORF$A_site > (circRNAs_150_ORF_positive_before_ORF$ORF_start-100), ]
h <- hist(circRNAs_150_ORF_positive_before_ORF_100nt$V2)
plot(h,xlim=c(0.5,1.0),xlab = "Threshold of Predicted m6A",main = "circRNAs m6A of Potential Translatability", lwd=2)
circRNAs_high_translatable = circRNAs_150_ORF_positive_before_ORF_100nt[circRNAs_150_ORF_positive_before_ORF_100nt$V2>=0.9,]

circRNAs_annotation = read.table("Desktop/m6A_model_validataion_data/circRNA/circRNADb/circRNA_dataset.txt",sep="\t")
circRNAs_annotation$V1 = as.character(circRNAs_annotation$V1)
circRNAs_high_translatable$seq_name = as.character(circRNAs_high_translatable$seq_name)
circRNAs_annotation = circRNAs_annotation[circRNAs_annotation[,1] %in% circRNAs_high_translatable$seq_name, ]

ribosomal_circRNAs = read.csv("Desktop/m6A_model_validataion_data/circRNA/ribosome_associated_circRNAs.csv",header = TRUE)
ribosomal_circRNAs$NM_id = as.character(ribosomal_circRNAs$NM_id)
circRNAs_annotation$V8 = as.character(circRNAs_annotation$V8)
ribosomal_circRNAs = ribosomal_circRNAs[(ribosomal_circRNAs$NM_id %in% circRNAs_annotation$V8) & (ribosomal_circRNAs$start %in% circRNAs_annotation$V3)&(ribosomal_circRNAs$end %in% circRNAs_annotation$V4), ]


#### calculate the relative position of positive m6A with ORFs
circRNAs_150_ORF_positive = circRNAs_150_ORF_positive[,-c(2,3,6)]
circRNAs_150_ORF_negative = circRNAs_150_ORF_negative[,-c(2,3,6)]

circRNAs_150_ORF_positive = circRNAs_150_ORF_positive[,-c(2,3)]
write.csv(circRNAs_150_ORF_positive,"Desktop/m6A_model_validataion_data/circRNA/150_ORF_positive_m6A.csv",row.names = FALSE)
write.csv(circRNAs_150_ORF_positive,"Desktop/m6A_model_validataion_data/circRNA/150_ORF_positive_m6A_prob_0.8.csv",row.names = FALSE)
#### m6A site distribution in circRNAs
par(mfrow=c(2,1))
###positive
circRNAs_150_ORF_positive$A_relative_position = circRNAs_150_ORF_positive$A_site/circRNAs_150_ORF_positive$seq_length
d <- density(circRNAs_150_ORF_positive$A_relative_position)

circRNAs_150_ORF_positive$ORF_Start_relative_position = circRNAs_150_ORF_positive$ORF_start/circRNAs_150_ORF_positive$seq_length
d2 <- density(circRNAs_150_ORF_positive$ORF_Start_relative_position)

circRNAs_150_ORF_positive$ORF_End_relative_position = circRNAs_150_ORF_positive$ORF_end/circRNAs_150_ORF_positive$seq_length
d3 <- density(circRNAs_150_ORF_positive$ORF_End_relative_position)

plot(range(d$x,d2$x,d3$x),range(d$y,d2$y,d3$y), type= "n", xlab = "Relative circRNA Position", ylab="Density",
     main = "circRNAs: Predicted m6A and ORF_start/_end Distribution");
lines(d, col='red', lwd =2)
lines(d2, col ='black', lwd =2 )
lines(d3, col = 'blue',lwd = 2)
#lines(d3, col = 'blue', lwd = 2)
legend(0.5, 3.0, legend = c("m6A_site","ORF_start_site","ORF_end_site"), col = c('red','black','blue'),
       lty = 1:1:1, lwd = 2:2:2, cex = 1)

####negative
circRNAs_150_ORF_negative$A_relative_position = circRNAs_150_ORF_negative$A_site/circRNAs_150_ORF_negative$seq_length
d3 <- density(circRNAs_150_ORF_negative$A_relative_position)

circRNAs_150_ORF_negative$ORF_Start_relative_position = circRNAs_150_ORF_negative$ORF_start/circRNAs_150_ORF_negative$seq_length
d4 <- density(circRNAs_150_ORF_negative$ORF_Start_relative_position)

plot(range(d3$x,d4$x),range(d3$y,d4$y), type= "n", xlab = "Relative Position in circRNAs", ylab="Density",
     main = "Negative_Distribution of predicted m6A sites and ORF start sites");
lines(d3, col='red', lwd=2)
lines(d4, col ='black', lwd=2 )
legend(0.5, 3.0, legend = c("m6A_site","ORF_start_site"), col = c('red','black'),
       lty = 1:1, lwd = 2:2, cex = 1, box.lty = 0)

##### 
circRNAs_150_ORF_positive$relative_position = NA;
circRNAs_150_ORF_positive$relative_position = (circRNAs_150_ORF_positive$A_site - circRNAs_150_ORF_positive$ORF_start + 1)/circRNAs_150_ORF_positive$ORF_length

d = density(circRNAs_150_ORF_positive$relative_position)
plot(d, xlim=c(-6,6),xlab = "Relative ORF Position",main = "circRNAs: Predicted m6A Distribution against ORF", lwd=2)
abline(v=c(0,1),col=c("blue","blue"),lty=c(2,1),lwd=c(2,2))
median(circRNAs_150_ORF_positive$relative_position)
abline(v=0.9532164,col="red",lty=2,lwd=2)
legend(2,0.2, c("ORF_start = 0","median = 0.95","ORF_end = 1"),col = c("blue","red","blue"),lty = c(2,2,1),lwd = c(2,2,2),cex = 1)
