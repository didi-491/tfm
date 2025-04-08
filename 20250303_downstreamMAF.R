#Jenifer Brea
#This script has been designed to analyze the downstream analysis of mutect2 results for FIS lung immunotherapy project

#DEPENDENCIES
library(maftools)
library(ggplot2)
library(forcats)
library(xlsx)
library(dplyr)
library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg19)
library(NMF)
library(pheatmap)
library(data.table)
# library(enrichR)
library(paletteer)

library(vidger)
library(ggpubr)
library(scales)
library(tidyverse)
library(gdata)

# Setting working directory
setwd("C://01_MariaGallardo/9040_ABACUS-PE/")

# Clinical data
load("C://01_MariaGallardo/clinical_data/20241016
     _MIBC_clinicaldata.RData")
clinical <- bladderDB


# Loading variants
sampleList <- clinical$WGS_pre_sample_name 
# c('17-2243_A2', '17_29658_F1', '17-3027_A3', '17_31535_A', '17_3209_A2', '17_5108', '17-61_A', '17-7736_B2', '16-26860', '16-35576', '18_32128','18_35034','18_36196', '18_38523', '19_17715', '19_30173', '20_08518', '20_16062')
sampleListM <- sampleList
# c('17-2243-A2', '17-29658-F1', '17-3027-A3', '17-31535-A', '17-32096-A2', '17-5108', '17-61-A', '17-7736-B2', '16-26860-B2', '16-35576-A', '18-32128','18-35034', '18-36196', '18-38523','19-17715', '19-30173','20-08518', '20-16062')

#'17-2243_A2' and 19-35018 are from the same patient
# '19_21244' have no clinical data

# Arranging clinical
# clinical <- clinical[clinical$ID.Muestra %in% sampleListM,]
# clinical <- clinical[match(sampleListM, clinical$ID.Muestra),]
clinical$sampleID <- sampleList
clinical$Tumor_Sample_Barcode <- sampleList
clinical$Response_firstLine <- clinical$response
clinical$qualResponse_firstLine <- clinical$response_type

clinicalResp <- clinical[clinical$Response_firstLine == 'Responder',]
clinicalNon <- clinical[clinical$Response_firstLine == 'NonResponder',]


# Loading variants

# See: https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl (lines 1014:1039)
vc_nonsyn = c("Targeted_Region", 
              "Splice_Site", 
              "Nonsense_Mutation", 
              "Frame_Shift_Del", 
              "Frame_Shift_Ins", 
              "Nonstop_Mutation",
              "Translation_Start_Site", 
              "In_Frame_Ins", 
              "In_Frame_Del", 
              "Missense_Mutation", 
              "Intron", 
              "Splice_Region", 
              "Silent",
              "RNA", "5'UTR", "3'UTR", "IGR", "5'Flank", "3'Flank")


for(i in sampleList){
  print(i)
  nam <- paste(i, '_maf', sep="")
  assign(nam, read.maf(sprintf('Data/%s/%s', i, list.files((sprintf('Data/%s',i)), pattern = "\\.PASS.maf$")), vc_nonSyn = vc_nonsyn)) #"\\pattern$" indicates the pattern as a suffix
}

# Merging all variants
sampleMafList <- c(paste0(sampleList, '_maf'))
all <- merge_mafs(mget(sampleMafList), vc_nonSyn = vc_nonsyn)

all_samples <- all@data


all_samples$variantID <- paste0(all_samples$Chromosome, ':', all_samples$Start_Position, '-', all_samples$End_Position, '_', all_samples$Reference_Allele, '>', all_samples$Tumor_Seq_Allele2, '_', all_samples$Strand, '|', all_samples$Hugo_Symbol, '_', all_samples$HGVSc, '_', all_samples$HGVSp, '_', all_samples$Tumor_Sample_Barcode)
all_samples$variantID_noSample <- paste0(all_samples$Chromosome, ':', all_samples$Start_Position, '-', all_samples$End_Position, '_', all_samples$Reference_Allele, '>', all_samples$Tumor_Seq_Allele2, '_', all_samples$Strand, '|', all_samples$Hugo_Symbol, '_', all_samples$HGVSc, '_', all_samples$HGVSp)


# Loading conflicting genome regions

#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=918997873_vqLdb739jXJdJWbLxqNQYHQSShaQ&clade=mammal&org=Human&db=hg19&hgta_group=varRep&hgta_track=dbVar_conflict&hgta_table=0&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=primaryTable&hgta_outFileName=GIAB_call_conflict.bed
dbVarConflict <- read.delim('Databases/dbVar_Conflict_SV.bed', header = F)
dbVarConflict <- dbVarConflict[,1:3]
#Sonia or Martin has pass me the file in some moment, but I am not sure from were did they download it
encodeBlacklist <- read.delim('Databases/ENCODE_blacklist_hg19.bed', header = F)
encodeBlacklist <- encodeBlacklist[,1:3]
#https://webs.iiitd.edu.in/raghava/humcfs/download.html
fragileSites <- read.delim('Databases/fragile_sites.bed', header = F)
fragileSites <- fragileSites[,1:3]
#https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=problematic
giabCallConflict <- read.delim('Databases/GIAB_call_conflict.bed', header = F)
giabCallConflict <- giabCallConflict[,1:3]

conflictRegions <- as.data.frame(rbind(dbVarConflict, encodeBlacklist, fragileSites, giabCallConflict))
conflictRegionsRanges <- GRanges(seqnames=conflictRegions$V1,
                                 ranges=paste(conflictRegions$V2, conflictRegions$V3, sep='-'))

# Removing conflicting regions

allGranges <- GRanges(seqnames=all_samples$Chromosome,
                      ranges=paste(all_samples$Start_Position, all_samples$End_Position, sep='-'))
values(allGranges) <- as.data.frame(all_samples$variantID)

allNonConflictive <- subsetByOverlaps(allGranges, conflictRegionsRanges, invert = T)
allNonConflictive_data0 <- as.data.frame(mcols(allNonConflictive))

allNonConflictive_data <- all_samples[all_samples$variantID %in% allNonConflictive_data0[,1],]
allNonConflictive_data$shortVariantID <- gsub('_.\\|.*$', '',allNonConflictive_data$variantID)
allNonConflictive_data <- as.data.frame(allNonConflictive_data)

#Loading probableArtifacts, obtained with: /mnt/netapp1/genomes_tmp/9040_ABACUS-PE/6_snv/scripts/downstreamMutect_manualFiltering/awkDownstreamFilter.sh
probableArtifacts <- read.delim('Databases/ProbableArtifacts.tsv', header=F, sep=' ')

artifactsindata <- allNonConflictive_data$shortVariantID %in% probableArtifacts$V1

allNonConflictive_data <- allNonConflictive_data[!(artifactsindata),]

write.table(allNonConflictive_data, 'Results/allVariants.maf', sep='\t', col.names=T, row.names=F, quote=F)


################################################################################################
############################ 1. All samples
################################################################################################

###########################################
##########################################
########################################
## ALL VARIANTS
######################################
####################################

all <- read.maf('allVariants.maf', vc_nonSyn = vc_nonsyn, clinicalData = clinical)
all_samples <- all@data
all_samples <- merge(all_samples, clinical, by='Tumor_Sample_Barcode', )

lp <- read.maf('all_variants_potentiallyPathogenic.maf', vc_nonSyn = vc_nonsyn, clinicalData = clinical)
lp_samples <- lp@data
lp_samples <- merge(lp_samples, clinical, by='Tumor_Sample_Barcode', )

write.table(lp_samples, "20241216_somaticVariants_MIBC_WGShg19_likelyPathogenic.txt", sep = "\t", dec = ".", quote = F, row.names = F)


nbMut <- all_samples %>% group_by(Tumor_Sample_Barcode) %>% count()
nbMut <- merge(nbMut, clinical, by=c('Tumor_Sample_Barcode'), all.x=T)

ggplot(data=nbMut, aes(x=nbTD2_New150pb_1.0.2, y=n, fill=Response_firstLine, scatter=Tumor_Sample_Barcode))+
  geom_col(position = 'dodge')+
  xlab('')+
  ylab('nb LINE-1 events')+
  theme_bw()+
  theme(text=element_text(size = 20), legend.title = element_blank())

compare_means(n ~ Response_firstLine, nbMut, method = "wilcox.test")
cor.test(nbMut$n, nbMut$nbTD2_New150pb_1.0.2, method = 'spearman')

##############################
# 1. General summaries
##############################

annotColors <- list(qualResponse_firstLine=c('PD' = '#FC8D62', 'SD' = '#8DA0CB', 'PR' = '#66C2A5', 'CR' = '#A6D854'),
                    response = c("responder"= "#B3B3B3", "non_responder" = "#666666"),
                    sex = c("male" = "#B2DF8A", "female" = "#6A3D9A"),
                    smoking_status = c("smoker" = "#084081", "former_smoker" = "#4EB3D3", "never_smoker" = "#CCEBC5"))

# --- ALL VARIANTS oncoplot
oncoplot(maf = all, top = 40, draw_titv = T, clinicalFeatures = c('qualResponse_firstLine', 'response', 'sex', 'smoking_status'),
         annotationColor = annotColors, sortByAnnotation = F, annotationOrder = c('NE', 'PD', 'SD','PR','CR'),
         legendFontSize = 1.2, annotationFontSize = 1.2, bgCol = "white", borderCol = "#CCCCCC", titleText = "\n ONCOPLOT all variants")

# write.table(all_samples[all_samples$Hugo_Symbol%in%m6A.genes,], "20241016_m6A-genes_allVariants_MIBC.txt", sep = "\t", dec = ",", quote = F)

# --- ALL VARIANTS, filter out VAF < 0.075 & VAF > 0.45
maf_filter <- subsetMaf(maf=all,query="t_AF >= 0.075 & t_AF <= 0.45")

oncoplot(maf = maf_filter, top = 40, draw_titv = T, clinicalFeatures = c('qualResponse_firstLine', 'response', 'sex', 'smoking_status'),
         annotationColor = annotColors, sortByAnnotation = F, annotationOrder = c('NE', 'PD', 'SD','PR','CR'), legendFontSize = 1.2,
         annotationFontSize = 1.2, bgCol = "white", borderCol = "#CCCCCC", titleText = "ONCOPLOT filtering out variants with VAF < 0.075 & VAF > 0.45")
dev.off()

# --- ALL VARIANTS, filter out VAF < 0.075 & VAF > 0.45 & DP > 10
maf_filter <- subsetMaf(maf=all,query="t_AF >= 0.075 & t_AF <= 0.45 & t_depth > 10")

# --- ALL VARIANTS, filter out VAF < 0.075 & VAF > 0.45 & DP > 10 & variant attributes
maf_filter <- subsetMaf(maf=all,query="t_AF >= 0.075 & t_AF <= 0.45 & t_depth > 10 & IMPACT != 'LOW' & IMPACT != 'MODIFIER' & Consequence != 'synonymous_variant' & Variant_Classification != 'Silent'")




#################################
# 2. Mutated genes between groups
#################################


## RESPONDERS

resp_mutated_genes <- responders_samples %>% group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>% count()
colnames(resp_mutated_genes)[3] <- paste('nbVariants')
nbSamples_gene <- resp_mutated_genes %>% group_by(Hugo_Symbol) %>% count()
nbVariants_gene <- resp_mutated_genes %>% group_by(Hugo_Symbol) %>% summarise(nbVariants = sum(nbVariants))

resp_mutated_genes <- merge(nbSamples_gene, nbVariants_gene, by="Hugo_Symbol")
write.table(resp_mutated_genes, 'Results/allMutatedGenes_responders.tsv', col.names = T, row.names = F, quote=F, sep='\t')

## NON RESPONDERS

nonResp_mutated_genes <- nonResponders_samples %>% group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>% count()
colnames(nonResp_mutated_genes)[3] <- paste('nbVariants')
nbSamples_gene <- nonResp_mutated_genes %>% group_by(Hugo_Symbol) %>% count()
nbVariants_gene <- nonResp_mutated_genes %>% group_by(Hugo_Symbol) %>% summarise(nbVariants = sum(nbVariants))

nonResp_mutated_genes <- merge(nbSamples_gene, nbVariants_gene, by="Hugo_Symbol")
write.table(nonResp_mutated_genes, 'Results/allMutatedGenes_nonResponders.tsv', col.names = T, row.names = F, quote=F, sep='\t')

## COMPARISSON RESPONDERS VS NON RESPONDERS

multipleResp <- resp_mutated_genes[resp_mutated_genes$n > 1,]
multipleNonResp <- nonResp_mutated_genes[nonResp_mutated_genes$n > 1,]


justResp <- multipleResp[!(multipleResp$Hugo_Symbol %in% multipleNonResp$Hugo_Symbol),]
justResp_strict <- multipleResp[!(multipleResp$Hugo_Symbol %in% nonResp_mutated_genes$Hugo_Symbol),]
write.table(justResp_strict, 'Results/strictlyExclusivelyMutatedGenes_responders.tsv', sep='\t', quote=F, col.names=T, row.names=F)

justNonResp <- multipleNonResp[!(multipleNonResp$Hugo_Symbol %in% multipleResp$Hugo_Symbol),]
justNonResp_strict <- multipleNonResp[!(multipleNonResp$Hugo_Symbol %in% resp_mutated_genes$Hugo_Symbol),]
write.table(justNonResp_strict, 'Results/strictlyExclusivelyMutatedGenes_nonResponders.tsv', sep='\t', quote=F, col.names=T, row.names=F)


#############################
# 3. Variants between groups
#############################

## RESPONDERS

resp_variants <- responders_samples %>% group_by(variantID_noSample, Tumor_Sample_Barcode) %>% count()
colnames(resp_variants)[3] <- paste('nbVariants')
nbSamples_variants <- resp_variants %>% group_by(variantID_noSample) %>% count()
nbVariants_variant <- resp_variants %>% group_by(variantID_noSample) %>% summarise(nbVariants = sum(nbVariants))

resp_variants <- merge(nbSamples_variants, nbVariants_variant, by="variantID_noSample")
write.table(resp_variants, 'Results/allVariants_responders.tsv', col.names = T, row.names = F, quote=F, sep='\t')

## NON RESPONDERS

nonResp_variants <- nonResponders_samples %>% group_by(variantID_noSample, Tumor_Sample_Barcode) %>% count()
colnames(nonResp_variants)[3] <- paste('nbVariants')
nbSamples_variants <- nonResp_variants %>% group_by(variantID_noSample) %>% count()
nbVariants_variant <- nonResp_variants %>% group_by(variantID_noSample) %>% summarise(nbVariants = sum(nbVariants))

nonResp_variants <- merge(nbSamples_variants, nbVariants_variant, by="variantID_noSample")
write.table(nonResp_variants, 'Results/allVariants_nonResponders.tsv', col.names = T, row.names = F, quote=F, sep='\t')

# COMPARISON RESP VS NONRESP

multipleResp <- resp_variants[resp_variants$n > 1,]
multipleNonResp <- nonResp_variants[nonResp_variants$n > 1,]


varjustResp <- multipleResp[!(multipleResp$variantID_noSample %in% multipleNonResp$variantID_noSample),]
varjustResp_strict <- multipleResp[!(multipleResp$variantID_noSample %in% nonResp_variants$variantID_noSample),]

varjustResp_strict_annot <- merge(varjustResp_strict, all_samples, by=c("variantID_noSample"))
varjustRespSamples <- varjustResp_strict_annot %>% group_by(variantID_noSample) %>% summarize(samples=toString(Tumor_Sample_Barcode), 
                                                                                              qualResp=toString(qualResponse_firstLine))
varjustResp_strict_annot <- varjustResp_strict_annot[!duplicated(varjustResp_strict_annot$variantID_noSample),]
varjustResp_strict_annot <- merge(varjustResp_strict_annot, varjustRespSamples, by=c("variantID_noSample"))
varjustResp_strict_annot <- varjustResp_strict_annot[,c(1:3,5,13,14,18,19,39,40,41,42,51:105,121,151,152)]

write.table(varjustResp_strict_annot, 'Results/strictlyExclusivelyMutatedVariants_responders.tsv', sep='\t', quote=F, col.names=T, row.names=F)

mafjustResp_strict <- all@data[all@data$variantID_noSample %in% varjustResp_strict$variantID_noSample,]
mafjustResp_strict1 <- read.maf(mafjustResp_strict, vc_nonSyn = vc_nonsyn)

tnm <- trinucleotideMatrix(maf = mafjustResp_strict1, prefix = 'chr', add = T, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
sign <- estimateSignatures(mat = tnm, nTry = 5)
plotCophenetic(res = sign)
sig = extractSignatures(mat = tnm, n = 2) #adjust n depending on the last plot (you have to choose the one befoore the drope)

maftools::plotSignatures(nmfRes = sig)





varjustNonResp <- multipleNonResp[!(multipleNonResp$variantID_noSample %in% multipleResp$variantID_noSample),]
varjustNonResp_strict <- multipleNonResp[!(multipleNonResp$variantID_noSample %in% resp_variants$variantID_noSample),]

varjustNonResp_strict_annot <- merge(varjustNonResp_strict, all_samples, by=c("variantID_noSample"))
varjustNonRespSamples <- varjustNonResp_strict_annot %>% group_by(variantID_noSample) %>% summarize(samples=toString(Tumor_Sample_Barcode), 
                                                                                                    qualResp=toString(qualResponse_firstLine))
varjustNonResp_strict_annot <- varjustNonResp_strict_annot[!duplicated(varjustNonResp_strict_annot$variantID_noSample),]
varjustNonResp_strict_annot <- merge(varjustNonResp_strict_annot, varjustNonRespSamples, by=c("variantID_noSample"))
varjustNonResp_strict_annot <- varjustNonResp_strict_annot[,c(1:3,5,13,14,18,19,39,40,41,42,51:105,121,151,152)]

write.table(varjustNonResp_strict_annot, 'Results/strictlyExclusivelyMutatedVariants_nonResponders.tsv', sep='\t', quote=F, col.names=T, row.names=F)

mafjustNonResp_strict <- all@data[all@data$variantID_noSample %in% varjustNonResp_strict$variantID_noSample,]
mafjustNonResp_strict1 <- read.maf(mafjustNonResp_strict, vc_nonSyn = vc_nonsyn)

tnm <- trinucleotideMatrix(maf = mafjustNonResp_strict1, prefix = 'chr', add = T, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
sign <- estimateSignatures(mat = tnm, nTry = 3)
plotCophenetic(res = sign)
sig = extractSignatures(mat = tnm, n = 2) #adjust n depending on the last plot (you have to choose the one befoore the drope)

maftools::plotSignatures(nmfRes = sig)


## Some co-plots

coBarplot(m1 = responders, m2 = nonResponders, m1Name = 'Responders', m2Name = 'nonResponders', genes = c('KRAS', 'BCR', 'LUC7L', 'FOXF1', 'ETS2'))




###########################################
##########################################
########################################
## LIKELY PATHOGENIC VARIANTS
######################################
####################################

pathogenic <- all_samples[((grepl('deleterious', all_samples$SIFT) == TRUE & grepl('damaging', all_samples$PolyPhen) == TRUE)
                           | all_samples$IMPACT == 'HIGH') & (all_samples$gnomAD_AF < 0.01 | is.na(all_samples$gnomAD_AF)),]


write.table(pathogenic, 'Results/potentiallyPathogenic/all_variants.maf', quote=F, col.names=T, row.names=F, sep='\t')



pall <- read.maf('Results/potentiallyPathogenic/all_variants.maf', vc_nonSyn = vc_nonsyn, clinicalData = clinical)


presponders <- subsetMaf(pall, clinQuery = "Response_firstLine == 'Responder'")
presponders_samples <- presponders@data

pnonResponders <- subsetMaf(pall, clinQuery = "Response_firstLine == 'NonResponder'")
pnonResponders_samples <- pnonResponders@data


#################################
# 1. General Summaries
#################################

oncoplot(maf = pall, top = 40, draw_titv = T, clinicalFeatures = c('qualResponse_firstLine', 'nbTD2_New150pb_1.0.2'),
         annotationColor = annotColors, sortByAnnotation = T, annotationOrder = c('PD', 'SD', 'PR', 'CR'), legendFontSize = 0.8,
         annotationFontSize = 0.8)

#################################
# 2. Mutated genes between groups
#################################


## RESPONDERS

presp_mutated_genes <- presponders_samples %>% group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>% count()
colnames(presp_mutated_genes)[3] <- paste('nbVariants')
nbSamples_gene <- presp_mutated_genes %>% group_by(Hugo_Symbol) %>% count()
nbVariants_gene <- presp_mutated_genes %>% group_by(Hugo_Symbol) %>% summarise(nbVariants = sum(nbVariants))

presp_mutated_genes <- merge(nbSamples_gene, nbVariants_gene, by="Hugo_Symbol")
write.table(presp_mutated_genes, 'Results/potentiallyPathogenic/allMutatedGenes_responders.tsv', col.names = T, row.names = F, quote=F, sep='\t')


## NON RESPONDERS

pnonResp_mutated_genes <- pnonResponders_samples %>% group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>% count()
colnames(pnonResp_mutated_genes)[3] <- paste('nbVariants')
nbSamples_gene <- pnonResp_mutated_genes %>% group_by(Hugo_Symbol) %>% count()
nbVariants_gene <- pnonResp_mutated_genes %>% group_by(Hugo_Symbol) %>% summarise(nbVariants = sum(nbVariants))

pnonResp_mutated_genes <- merge(nbSamples_gene, nbVariants_gene, by="Hugo_Symbol")
write.table(pnonResp_mutated_genes, 'Results/potentiallyPathogenic/allMutatedGenes_nonResponders.tsv', col.names = T, row.names = F, quote=F, sep='\t')


## COMPARISSON RESPONDERS VS NON RESPONDERS

pmultipleResp <- presp_mutated_genes[presp_mutated_genes$n > 1,]
pmultipleNonResp <- pnonResp_mutated_genes[pnonResp_mutated_genes$n > 1,]


pjustResp <- pmultipleResp[!(pmultipleResp$Hugo_Symbol %in% pmultipleNonResp$Hugo_Symbol),]
pjustResp_strict <- pmultipleResp[!(pmultipleResp$Hugo_Symbol %in% pnonResp_mutated_genes$Hugo_Symbol),]
write.table(pjustResp_strict, 'Results/potentiallyPathogenic/strictlyExclusivelyMutatedGenes_responders.tsv', sep='\t', quote=F, col.names=T, row.names=F)


pjustNonResp <- pmultipleNonResp[!(pmultipleNonResp$Hugo_Symbol %in% pmultipleResp$Hugo_Symbol),]
pjustNonResp_strict <- pmultipleNonResp[!(pmultipleNonResp$Hugo_Symbol %in% presp_mutated_genes$Hugo_Symbol),]
write.table(pjustNonResp_strict, 'Results/potentiallyPathogenic/strictlyExclusivelyMutatedGenes_nonResponders.tsv', sep='\t', quote=F, col.names=T, row.names=F)

penrichr_results <- enrichr(pjustResp_strict$Hugo_Symbol, c('MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures','MSigDB_Computational','BioPlanet_2019', 'KEGG_2019_Human', 'WikiPathways_2019_Human', 'TRANSFAC_and_JASPAR_PWMs', 'ChEA_2016', 'GO_Molecular_Function_2018', 'GO_Biological_Process_2018'))

write.xlsx(penrichr_results$MSigDB_Oncogenic_Signatures, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'MSigDB_Hallmark_2020', row.names = F)
write.xlsx(penrichr_results$MSigDB_Oncogenic_Signatures, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'MSigDB_Oncogenic_Signatures',append = T, row.names = F)
write.xlsx(penrichr_results$MSigDB_Computational, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'MSigDB_Computational',append = T, row.names = F)
write.xlsx(penrichr_results$BioPlanet_2019, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'BioPlanet_2019',append = T, row.names = F)
write.xlsx(penrichr_results$KEGG_2019_Human, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'KEGG_2019_Human', append = T, row.names = F)
write.xlsx(penrichr_results$WikiPathways_2019_Human, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'WikiPathways_2019_Human', append = T, row.names = F)
write.xlsx(penrichr_results$TRANSFAC_and_JASPAR_PWMs, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'TRANSFAC_and_JASPAR_PWMs', append = T, row.names = F)
write.xlsx(penrichr_results$ChEA_2016, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'Chea_2016', append = T, row.names = F)
write.xlsx(penrichr_results$GO_Molecular_Function_2018, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'GO_Molecular_Function_2018', append = T, row.names = F)
write.xlsx(penrichr_results$GO_Biological_Process_2018, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_responders_strict.xlsx', sheetName = 'GO_Biological_Process_2018', append = T, row.names = F)


penrichr_results <- enrichr(pjustNonResp_strict$Hugo_Symbol, c('MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures','MSigDB_Computational','BioPlanet_2019', 'KEGG_2019_Human', 'WikiPathways_2019_Human', 'TRANSFAC_and_JASPAR_PWMs', 'ChEA_2016', 'GO_Molecular_Function_2018', 'GO_Biological_Process_2018'))

write.xlsx(penrichr_results$MSigDB_Oncogenic_Signatures, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'MSigDB_Hallmark_2020', row.names = F)
write.xlsx(penrichr_results$MSigDB_Oncogenic_Signatures, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'MSigDB_Oncogenic_Signatures',append = T, row.names = F)
write.xlsx(penrichr_results$MSigDB_Computational, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'MSigDB_Computational',append = T, row.names = F)
write.xlsx(penrichr_results$BioPlanet_2019, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'BioPlanet_2019',append = T, row.names = F)
write.xlsx(penrichr_results$KEGG_2019_Human, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'KEGG_2019_Human', append = T, row.names = F)
write.xlsx(penrichr_results$WikiPathways_2019_Human, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'WikiPathways_2019_Human', append = T, row.names = F)
write.xlsx(penrichr_results$TRANSFAC_and_JASPAR_PWMs, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'TRANSFAC_and_JASPAR_PWMs', append = T, row.names = F)
write.xlsx(penrichr_results$ChEA_2016, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'Chea_2016', append = T, row.names = F)
write.xlsx(penrichr_results$GO_Molecular_Function_2018, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'GO_Molecular_Function_2018', append = T, row.names = F)
write.xlsx(penrichr_results$GO_Biological_Process_2018, 'Results/potentiallyPathogenic/enrich_strictlyExclusivelyMutated_genes_nonResponders_strict.xlsx', sheetName = 'GO_Biological_Process_2018', append = T, row.names = F)


#############################
# 3. Variants between groups
#############################

# RESPONDERS

presp_variants <- presponders_samples %>% group_by(variantID_noSample, Tumor_Sample_Barcode) %>% count()
colnames(presp_variants)[3] <- paste('nbVariants')
nbSamples_variants <- presp_variants %>% group_by(variantID_noSample) %>% count()
nbVariants_variant <- presp_variants %>% group_by(variantID_noSample) %>% summarise(nbVariants = sum(nbVariants))

presp_variants <- merge(nbSamples_variants, nbVariants_variant, by="variantID_noSample")
write.table(presp_variants, 'Results/potentiallyPathogenic/allVariants_responders.tsv', col.names = T, row.names = F, quote=F, sep='\t')


# NONRESPONDERS

pnonResp_variants <- pnonResponders_samples %>% group_by(variantID_noSample, Tumor_Sample_Barcode) %>% count()
colnames(pnonResp_variants)[3] <- paste('nbVariants')
nbSamples_variants <- pnonResp_variants %>% group_by(variantID_noSample) %>% count()
nbVariants_variant <- pnonResp_variants %>% group_by(variantID_noSample) %>% summarise(nbVariants = sum(nbVariants))

pnonResp_variants <- merge(nbSamples_variants, nbVariants_variant, by="variantID_noSample")
write.table(pnonResp_variants, 'Results/potentiallyPathogenic/allVariants_nonResponders.tsv', col.names = T, row.names = F, quote=F, sep='\t')


# COMPARISON RESP VS NONRESP
pmultipleResp <- presp_variants[presp_variants$n > 1,]
pmultipleNonResp <- pnonResp_variants[pnonResp_variants$n > 1,]


pvarjustResp_strict <- pmultipleResp[!(pmultipleResp$variantID_noSample %in% pnonResp_variants$variantID_noSample),]
write.table(pvarjustResp_strict, 'Results/potentiallyPathogenic/strictlyExclusivelyMutatedVariants_responders.tsv', sep='\t', quote=F, col.names=T, row.names=F)


pvarjustNonResp_strict <- pmultipleNonResp[!(pmultipleNonResp$variantID_noSample %in% presp_variants$variantID_noSample),]
write.table(pvarjustNonResp_strict, 'Results/potentiallyPathogenic/strictlyExclusivelyMutatedVariants_nonResponders.tsv', sep='\t', quote=F, col.names=T, row.names=F)
