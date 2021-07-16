###### Part I: Import database and library #####
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(httr)
library(jsonlite)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg19
source('./find_effect_variation.R')

###### Part II: VCF file cleaning ######

### Import vcf file as a dataframe (Header is not included)
filename <- "Challenge_data_(1).vcf"
vcf <- read.table(filename, stringsAsFactors = FALSE, comment.char="#")

### Read field name from vcf file
field_names <- c()
con <- file(filename, "r")
while (TRUE) {
  field_names <- readLines(con, 1)
  if (length(field_names) == 0) break
  else if (grepl("^#", field_names) & !(grepl("^#{2}", field_names))) break
}
close(con)
if (length(field_names) == 0) stop("Column names were not found in VCF file.")
field_names <- sub("#", "", field_names)
field_names <- strsplit(field_names, "\t", fixed = TRUE)[[1]]
if (length(field_names) != ncol(vcf)) stop("One or more field names were missing.")
colnames(vcf) <- field_names

### Transform INFO into data.frame
vcf_INFO_list <- strsplit(vcf$INFO, ";", fixed = T)
if (length(unique(lengths(vcf_INFO_list))) != 1L) stop("One or more records were not complete.")
vcf_INFO_list_m <- lapply(vcf_INFO_list, function(x) {matrix(unlist(strsplit(x, "=", fixed = T)), 
                                                      nrow = length(x), byrow = T)}) 
                                                      #Separate key and value
vcf_INFO <- as.data.frame(t(sapply(vcf_INFO_list_m, function(x) x[,2])))
colnames(vcf_INFO) <- vcf_INFO_list_m[[1]][,1]

### Tranform FORMAT into data.frame
n <- which(field_names == "FORMAT")
ns <- length(field_names[-(1:n)]) #Determine the number of samples
FORMAT_key <- strsplit(vcf[1,n], ":", fixed = TRUE)[[1]]
vcf_FORMAT <- list()
for (i in 1:ns){
  vcf_FORMAT_list <- strsplit(vcf[,n+i], ":", fixed = TRUE)
  if (length(unique(lengths(vcf_FORMAT_list))) != 1L) stop("One or more records were not complete.")
  vcf_FORMAT[[i]] <- as.data.frame(t(sapply(vcf_FORMAT_list, function(x) x)))
  colnames(vcf_FORMAT[[i]]) <- FORMAT_key
  names(vcf_FORMAT)[i] <- field_names[n+i]
}

### Combine all fields from INFO and FORMAT -> vcf_extend -> save as csv file
vcf$RECORD_ID <- paste0("line", 1:nrow(vcf))
vcf_extend <- vcf[, c(ncol(vcf), 1:(which(field_names == "INFO")-1))]
vcf_extend <- cbind(vcf_extend, vcf_INFO, vcf_FORMAT)
write.csv(vcf_extend, "vcf_extend.csv", row.names = F)


###### Part III: VCF file annotation ######

### 0. Initialize annotated dataframe
vcf_annotation <- vcf_extend[,c('RECORD_ID','CHROM','POS','REF','ALT','LEN')]

### 1. Type of variation and their effect
vcf_annotation$Variation_Type <- vcf_extend$TYPE
effect_severity <- c(0, 1, 2, 2, 2, 3, 3, 3, 3)  ## score the severity  of variation effect
names(effect_severity) <- c("Silent mutation",   ## 0:no risk
                            "Intergenic/Intragenic variant",  ## 1:low risk
                            "Missense mutation", "Inframe insertion", "Inframe deletion", ## 2:moderate
                            "Nonsense mutation", "Stop codon lost", "Start codon lost", "Frameshift") ## 3:high risk
#annotate variation effect
vcf_annotation$Variation_Effect <- character(length = nrow(vcf_annotation))
for (i in 1772:nrow(vcf_annotation)){
  print(i)
  vcf_type <- strsplit(as.character(vcf_annotation$Variation_Type[i]), ',')[[1]]
  if (length(vcf_type) == 0) stop("Variation type is missing.")
  alt <- strsplit(vcf_annotation$ALT[i], ',', fixed = TRUE)[[1]]
  effect_temp <- c()
  for (j in 1:length(vcf_type)){
    effect_temp <- c(effect_temp, switch(vcf_type[j], 
                                         "snp"=find_effect_snp(chr = paste0('chr',vcf_annotation$CHROM[i]), 
                                                               pos = vcf_annotation$POS[i], 
                                                               ref = vcf_annotation$REF[i],
                                                               alt = alt[j]),
                                         "ins"=find_effect_insdel(chr = paste0('chr',vcf_annotation$CHROM[i]),
                                                                  pos = vcf_annotation$POS[i], 
                                                                  ref = vcf_annotation$REF[i],
                                                                  alt = alt[j]),
                                         "del"=find_effect_insdel(chr = paste0('chr',vcf_annotation$CHROM[i]),
                                                                  pos = vcf_annotation$POS[i], 
                                                                  ref = vcf_annotation$REF[i],
                                                                  alt = alt[j]),
                                         "complex"=find_effect_complex(chr = paste0('chr',vcf_annotation$CHROM[i]), 
                                                                       pos = vcf_annotation$POS[i], 
                                                                       ref = vcf_annotation$REF[i],
                                                                       alt = alt[j]),
                                         "mnp"=find_effect_mnp(chr = paste0('chr',vcf_annotation$CHROM[i]), 
                                                               pos = vcf_annotation$POS[i], 
                                                               alt = alt[j])))
  }
  ## find the most deleterious effect
  vcf_annotation$Variation_Effect[i] <- effect_temp[which.max(effect_severity[effect_temp])]
}

### 2. Depth of sequence coverage at the site of variation
vcf_annotation$Depth <- as.numeric(as.character(vcf_extend$DP))

### 3. Number of reads supporting the variant
AO <- as.character(vcf_extend$AO)
AO_split <- strsplit(AO, ",", fixed = T)
AO_split <- sapply(AO_split, function(x) sum(as.numeric(x)))
vcf_annotation$Variant_reads <- AO_split

### 4. Percentage of reads supporting the variant versus those supporting reference reads
vcf_annotation$Variant_Reference_ratio <- vcf_annotation$Variant_reads/as.numeric(as.character(vcf_extend$RO))

### 5. Allele frequency of variant from ExAC API
vcf_annotation$Frequency_ExAC <- numeric(length = nrow(vcf_annotation))
ExAC_call <- paste0(vcf_annotation$CHROM,'-',vcf_annotation$POS,'-',vcf_annotation$REF,'-',vcf_annotation$ALT)
ExAC_query <- POST(url="http://exac.hms.harvard.edu/rest/bulk/variant/variant", body=toJSON(ExAC_call), encode = "json")
ExAC_content <- content(ExAC_query)
for (i in 1:nrow(vcf_annotation)){
  vcf_annotation$Frequency_ExAC[i] <- if(!(is.null(ExAC_content[[ExAC_call[i]]]$allele_freq))) ExAC_content[[ExAC_call[i]]]$allele_freq else NA
}

###### Part IV: Save VCF annotation files ######
write.csv(vcf_annotation, 'vcf_annotation.csv', row.names = F)

