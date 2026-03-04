################################################################################################################################################################
################################################################################################################################################################ 
### here comparing the GWAS HITS :


library("cowplot")
library("dplyr")

library("gmodels")
library("gplots")

library("ggplot2")
library("gridExtra")

library("stringr")

library("VennDiagram") 
library("Vennerable")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# list.files()

# ASTIGMATISM = "astigmatism.STRONGEST.associations.txt"   
# HYPEROPIA = "hyperopia.STRONGEST.associations.txt"     
# MYOPIA = "myopia.STRONGEST.associations.txt"        
# MYOPIA_UNDER10 = "myopia.under10.STRONGEST.associations.txt"

ASTIGMATISM = "astigmatism.STRONGEST.associations.txt.sep.txt.processed.txt"
HYPEROPIA = "hyperopia.STRONGEST.associations.txt.sep.txt.processed.txt"
MYOPIA = "myopia.STRONGEST.associations.txt.sep.txt.processed.txt"
MYOPIA_UNDER10 = "myopia.under10.STRONGEST.associations.txt.sep.txt.processed.txt"

#######################################################################

astigmatism <- read.delim(ASTIGMATISM, header=T, sep="\t", stringsAsFactors=F)

hyperopia <-   read.delim(HYPEROPIA, header=T, sep="\t", stringsAsFactors=F)

myopia <-      read.delim(MYOPIA, header=T, sep="\t", stringsAsFactors=F)

myopia.under10 <-  read.delim(MYOPIA_UNDER10, header=T, sep="\t", stringsAsFactors=F)

################################################################################################################################################################
################################################################################################################################################################ 

dim(astigmatism)
dim(hyperopia)
dim(myopia)
dim(myopia.under10)

# colnames(myopia.under10)
# [1] "cytoband"        "assay.name"      "position"        "alleles"        
# [5] "src"             "dose.b"          "pvalue"          "OR"             
# [9] "X95..CI"         "gene.context"    "gene_upstream"   "p1"             
# [13] "intra_gene"      "p2"              "gene_downstream"

# dim(astigmatism)
# [1] 126  15
# dim(hyperopia)
# [1] 191  15
# dim(myopia)
# [1] 601  15
# dim(myopia.under10)
# [1] 32 15

head(astigmatism)
head(hyperopia)
head(myopia)
head(myopia.under10)

################################################################################################################################################################
################################################################################################################################################################ 

# here to do a very quick intersections of the SNP ID :

#  grep(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE,
#               fixed = FALSE, useBytes = FALSE, invert = FALSE)

#      1. ‘fixed = TRUE’: use exact matching.
#      2. ‘perl = TRUE’: use Perl-style regular expressions.
#      3. ‘fixed = FALSE, perl = FALSE’: use POSIX 1003.2 extended
#          regular expressions (the default)

################################################################################################################################################################
################################################################################################################################################################ 
# to separate the LAST COLUMN into GENE and NEIGHBOURING GENES in order to be able to do the comparisons : we did it in a PERL SCRIPT ...

# here to INTERSECT THE SNPs : we RE-LABEL the extragenic SNP

################################################################################################################################################################
################################################################################################################################################################ 

astigmatism$inter_gene <- ifelse( ( astigmatism$gene_upstream != ""  | astigmatism$gene_downstream != "" ) , 
                                    paste(astigmatism$gene_upstream, astigmatism$gene_downstream, sep="::"), 
                                   "") 

# write.table(astigmatism, file=paste(ASTIGMATISM, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

astigmatism$gene_neighbourhood <-  ifelse( ( astigmatism$inter_gene != ""), 
                                             astigmatism$inter_gene, 
                                             astigmatism$intra_gene  ) 

write.table(astigmatism, file=paste(ASTIGMATISM, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

################################################################################################################################################################
################################################################################################################################################################ 

hyperopia$inter_gene <- ifelse( ( hyperopia$gene_upstream != ""  | hyperopia$gene_downstream != "" ) , 
                                  paste(hyperopia$gene_upstream, hyperopia$gene_downstream, sep="::"), 
                                  "") 

# write.table(hyperopia, file=paste(HYPEROPIA, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

hyperopia$gene_neighbourhood <-  ifelse( ( hyperopia$inter_gene != ""), 
                                           hyperopia$inter_gene, 
                                           hyperopia$intra_gene  ) 

write.table(hyperopia, file=paste(HYPEROPIA, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

################################################################################################################################################################
################################################################################################################################################################ 

myopia$inter_gene <- ifelse( ( myopia$gene_upstream != ""  | myopia$gene_downstream != "" ) , 
                               paste(myopia$gene_upstream, myopia$gene_downstream, sep="::"), 
                               "") 


# write.table(myopia, file=paste(MYOPIA, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

myopia$gene_neighbourhood <-  ifelse( ( myopia$inter_gene != ""), 
                                        myopia$inter_gene, 
                                        myopia$intra_gene  ) 

write.table(myopia, file=paste(MYOPIA, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

################################################################################################################################################################
################################################################################################################################################################ 

myopia.under10$inter_gene <- ifelse( ( myopia.under10$gene_upstream != ""  | myopia.under10$gene_downstream != "" ) , 
                                       paste(myopia.under10$gene_upstream, myopia.under10$gene_downstream, sep="::"), 
                                        "") 


# write.table(myopia.under10, file=paste(MYOPIA_UNDER10, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

myopia.under10$gene_neighbourhood <-  ifelse( ( myopia.under10$inter_gene != ""), 
                                                myopia.under10$inter_gene, 
                                                myopia.under10$intra_gene  ) 

write.table(myopia.under10, file=paste(MYOPIA_UNDER10, "with.gene.info", sep="."), quote=F, sep="\t", row.names = FALSE)

################################################################################################################################################################
################################################################################################################################################################ 
################################################################################################################################################################

# astigmatism 
# hyperopia 
# myopia 
# myopia.under10 

library("VennDiagram") 
library("Vennerable")
library("gplots")

################################################################### an example to compute the elements that are in the intersections :

# oneName <- function() paste(sample(LETTERS,5,replace=TRUE),collapse="")
# geneNames <- replicate(1000, oneName())

# GroupA <- sample(geneNames, 400, replace=FALSE)
# GroupB <- sample(geneNames, 750, replace=FALSE)
# GroupC <- sample(geneNames, 250, replace=FALSE)
# GroupD <- sample(geneNames, 300, replace=FALSE)

# input  <-list(GroupA,GroupB,GroupC,GroupD)

# tmp <- venn(input)
# isect <- attr(tmp, "intersections")
# isect <- attr(tmp, "intersection")
# isect

################################################################################################################################################################ 
################################################################################################################################################################
################################################################################################################################################################ 
################################################################################################################################################################
################################################################################################################################################################ intersect SNP POSITION

sets <- list( sets_astigmatism  = astigmatism$position,
              sets_hyperopia  = hyperopia$position, 
              # sets_myopia.under10  = myopia.under10$position,
              sets_myopia = myopia$position) 

sets_venn <- venn.diagram(sets, filename=NULL)

png("figure.intersecting.SNP.POSITIONS.between.astigmatism.myopia.hyperopia.png" )
grid.newpage()
grid.draw(sets_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets_overlap <- calculate.overlap(sets)

################################################################################################################################################################

sets_venn_gplots <- venn(sets)
attr(sets_venn_gplots, "intersections")

png("figure.intersecting.SNP.POSITIONS.between.astigmatism.myopia.hyperopia.gplots.png" )
venn(sets)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including myopia_under_10

sets_m <- list(  sets_myopia.under10  = myopia.under10$position,
                 sets_myopia = myopia$position) 

sets_m_venn <- venn.diagram(sets_m, filename=NULL)

png("figure.intersecting.SNP.POSITIONS.between.myopia.adult.children.png" )
grid.newpage()
grid.draw(sets_m_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets_m_overlap <- calculate.overlap(sets_m)

################################################################################################################################################################

sets_m_venn_gplots <- venn(sets_m)
attr(sets_m_venn_gplots, "intersections")

png("figure.intersecting.SNP.POSITIONS.between.myopia.adult.children.gplots.png" )
venn(sets_m)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including all 4 diagnoses

sets_all <- list( sets_astigmatism  = astigmatism$position,
                  sets_hyperopia  = hyperopia$position, 
                  sets_myopia.under10  = myopia.under10$position,
                  sets_myopia = myopia$position) 

sets_all_venn <- venn.diagram(sets_all, filename=NULL)

png("figure.intersecting.SNP.POSITIONS.between.astigmatism.myopia.hyperopia.and.myopia.kids.png", width = 1024, height = 1024)
grid.newpage()
grid.draw(sets_all_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets_all_overlap <- calculate.overlap(sets_all)

################################################################################################################################################################

sets_all_venn_gplots <- venn(sets_all)
attr(sets_all_venn_gplots, "intersections")

png("figure.intersecting.SNP.POSITIONS.between.astigmatism.myopia.hyperopia.and.myopia.kids.gplots.png", width = 1024, height = 1024)
venn(sets_all)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ intersect SNP ID
################################################################################################################################################################
################################################################################################################################################################

sets2 <- list( sets_astigmatism  = astigmatism$assay.name,
               sets_hyperopia  = hyperopia$assay.name, 
               # sets_myopia.under10  = myopia.under10$assay.name,
               sets_myopia = myopia$assay.name) 
 
sets2_venn <- venn.diagram(sets2, filename=NULL)

png("figure.intersecting.SNP.ID.between.astigmatism.myopia.hyperopia.png" )
grid.newpage()
grid.draw(sets2_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets2_overlap <- calculate.overlap(sets2)

################################################################################################################################################################

sets2_venn_gplots <- venn(sets2)
attr(sets2_venn_gplots, "intersections")

png("figure.intersecting.SNP.ID.between.astigmatism.myopia.hyperopia.gplots.png" )
venn(sets2)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including myopia_under_10

sets2_m <- list(  sets_myopia.under10  = myopia.under10$assay.name,
                  sets_myopia = myopia$assay.name) 

sets2_m_venn <- venn.diagram(sets2_m, filename=NULL)

png("figure.intersecting.SNP.ID.between.myopia.adult.children.png" )
grid.newpage()
grid.draw(sets2_m_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets2_m_overlap <- calculate.overlap(sets2_m)

################################################################################################################################################################

sets2_m_venn_gplots <- venn(sets2_m)
attr(sets2_m_venn_gplots, "intersections")

png("figure.intersecting.SNP.ID.between.myopia.adult.children.gplots.png" )
venn(sets2_m)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including all 4 diagnoses

sets2_all <- list( sets_astigmatism  = astigmatism$assay.name,
                  sets_hyperopia  = hyperopia$assay.name, 
                  sets_myopia.under10  = myopia.under10$assay.name,
                  sets_myopia = myopia$assay.name) 

sets2_all_venn <- venn.diagram(sets2_all, filename=NULL)

png("figure.intersecting.SNP.ID.between.astigmatism.myopia.hyperopia.and.myopia.kids.png", width = 1024, height = 1024)
grid.newpage()
grid.draw(sets2_all_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets2_all_overlap <- calculate.overlap(sets2_all)

################################################################################################################################################################

sets2_all_venn_gplots <- venn(sets2_all)
attr(sets2_all_venn_gplots, "intersections")

png("figure.intersecting.SNP.ID.between.astigmatism.myopia.hyperopia.and.myopia.kids.gplots.png", width = 1024, height = 1024)
venn(sets2_all)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ INTERSECT GENE names
################################################################################################################################################################ gene NEIGHBOURHOODS
################################################################################################################################################################
################## in order to intersect the HITS based on GENE NEIGHBOURHOODS that were computed above ...

sets3 <- list( sets_astigmatism  = astigmatism$gene_neighbourhood,
               sets_hyperopia  = hyperopia$gene_neighbourhood, 
               # sets_myopia.under10  = myopia.under10$gene_neighbourhood,
               sets_myopia = myopia$gene_neighbourhood) 
 
sets3_venn <- venn.diagram(sets3, filename=NULL)

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.astigmatism.myopia.hyperopia.png" )
grid.newpage()
grid.draw(sets3_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets3_overlap <- calculate.overlap(sets3)

################################################################################################################################################################

sets3_venn_gplots <- venn(sets3)
attr(sets3_venn_gplots, "intersections")

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.astigmatism.myopia.hyperopia.gplots.png" )
venn(sets3)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including myopia_under_10

sets3_m <- list(  sets_myopia.under10  = myopia.under10$gene_neighbourhood,
                  sets_myopia = myopia$gene_neighbourhood) 

sets3_m_venn <- venn.diagram(sets3_m, filename=NULL)

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.myopia.adult.children.png" )
grid.newpage()
grid.draw(sets3_m_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets3_m_overlap <- calculate.overlap(sets3_m)

################################################################################################################################################################

sets3_m_venn_gplots <- venn(sets3_m)
attr(sets3_m_venn_gplots, "intersections")

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.myopia.adult.children.gplots.png" )
venn(sets3_m)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
#################### including all 4 diagnoses

sets3_all <- list( sets_astigmatism  = astigmatism$gene_neighbourhood,
                  sets_hyperopia  = hyperopia$gene_neighbourhood, 
                  sets_myopia.under10  = myopia.under10$gene_neighbourhood,
                  sets_myopia = myopia$gene_neighbourhood) 

sets3_all_venn <- venn.diagram(sets3_all, filename=NULL)

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.astigmatism.myopia.hyperopia.and.myopia.kids.png", width = 1024, height = 1024)
grid.newpage()
grid.draw(sets3_all_venn)
dev.off()

# In VennDiagram package, it has a function called "calculate.overlap".
sets3_all_overlap <- calculate.overlap(sets3_all)

################################################################################################################################################################

sets3_all_venn_gplots <- venn(sets3_all)
attr(sets3_all_venn_gplots, "intersections")

png("figure.intersecting.GENE.NEIGHBOURHOOD.between.astigmatism.myopia.hyperopia.and.myopia.kids.gplots.png", width = 1024, height = 1024)
venn(sets3_all)
dev.off()

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ to print 
################################################################################################################################################################ GENE NEIGHBOURHOODS
################################################################################################################################################################ in each category

# str(sets3_all)
#List of 4
# $ sets_astigmatism   : chr [1:126] "PRSS56" "RBFOX1" "LAMA2" "TSPAN10" ...
# $ sets_hyperopia     : chr [1:191] "LAMA2" "GOLGA8B::GJD2" "KCNQ5" "RBFOX1" ...
# $ sets_myopia.under10: chr [1:32] "LAMA2" "GOLGA8B::GJD2" "RBFOX1" "TOX::CA8" ...
# $ sets_myopia        : chr [1:601] "RBFOX1" "GOLGA8B::GJD2" "KCNQ5" "LAMA2" ...

# length(sets3_all$sets_astigmatism)
#[1] 126
# length(sets3_all$sets_hyperopia)
#[1] 191
# length(sets3_all$sets_myopia.under10)
#[1] 32
# length(sets3_all$sets_myopia)
#[1] 601

# sets3_all_venn_gplots <- venn(sets3_all)
# attr(sets3_all_venn_gplots, "intersections")

# here just to print the elements that are in each category :

# $sets_myopia                                                    : chr [1:457] "PZP" "RGR::CCSER2" "METAP1D" "TFAP2B" ...
# $sets_hyperopia                                                 : chr [1:68] "CAMKMT::SIX3" "APOE" "ZKSCAN8::ZSCAN9" "ICA1" ...
# $sets_astigmatism                                               : chr [1:51] "NEGR1::" "MEF2C" "PPM1A::C14orf39" "PER2::TRAF3IP1" ...
# $sets_myopia.under10                                            : chr [1:12] "::ADGRL3" "NRG3::GHITM" "PCBP3" "ARHGAP15" ...
# $sets_hyperopia:sets_myopia                                     : chr [1:68] "GNB3" "BICC1" "RASGRF1" "SLC39A8" ...
# $sets_astigmatism:sets_myopia                                   : chr [1:22] "TSPAN10" "::CADM2" "NT5DC1,COL10A1" "GABRA6::GABRA1" ...
# $sets_astigmatism:sets_hyperopia                                : chr [1:3] "TFAP2B::PKHD1" "RGR" "CD55::CR2"
# $sets_myopia.under10:sets_myopia                                : chr [1:2] "PRDM13::MCHR2" "AKAP6"
# $sets_hyperopia:sets_myopia.under10                             : chr "NPLOC4"
# $sets_astigmatism:sets_myopia.under10                           : chr [1:2] "SLC25A12" "HLA-F::HLA-G"
# $sets_astigmatism:sets_hyperopia:sets_myopia                    : chr [1:34] "ZIC2::PCCA" "BMP3" "SOX4::PRL" "GRM7::LMCD1" ...
# $sets_hyperopia:sets_myopia.under10:sets_myopia                 : chr "BMP4::CDKN3"
# $sets_astigmatism:sets_hyperopia:sets_myopia.under10            : chr "A2M::PZP"
# $sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia: chr [1:13] "PRSS56" "RBFOX1" "LAMA2" "KCNQ5" .

sets3_all_venn_gplots_attr <- attr(sets3_all_venn_gplots, "intersections")

### PRINTING ALL THESE CATEGORIES : 

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism", 
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_hyperopia",
file="fig.sets3_all_venn_gplots_attr.sets_hyperopia.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F)


write.table(sets3_all_venn_gplots_attr$"sets_myopia.under10", 
file="fig.sets3_all_venn_gplots_attr.sets_myopia.under10.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_myopia", 
file="fig.sets3_all_venn_gplots_attr.sets_myopia.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F)

###

write.table(sets3_all_venn_gplots_attr$"sets_hyperopia:sets_myopia", 
file="fig.sets3_all_venn_gplots_attr.sets_hyperopia:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_myopia", 
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia", 
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_myopia.under10:sets_myopia", 
file="fig.sets3_all_venn_gplots_attr.sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_hyperopia:sets_myopia.under10", 
file="fig.sets3_all_venn_gplots_attr.sets_hyperopia:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_myopia.under10",      
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F)

###

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia",                 
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_hyperopia:sets_myopia.under10:sets_myopia",                
file="fig.sets3_all_venn_gplots_attr.sets_hyperopia:sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia.under10",           
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F)

###

write.table(sets3_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia",
file="fig.sets3_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ to print the SNP IDs
################################################################################################################################################################ in each category
################################################################################################################################################################

str(sets_all)

length(sets_all$sets_astigmatism)
# [1] 126
length(sets_all$sets_hyperopia)
# [1] 191
length(sets_all$sets_myopia.under10)
# [1] 32
length(sets_all$sets_myopia)
# [1] 601

sets_all_venn_gplots <- venn(sets_all)
attr(sets_all_venn_gplots, "intersections")

# here just to print the elements that are in each category :

# $sets_myopia                                                    
# $sets_hyperopia                                                 
# $sets_astigmatism                                               
# $sets_myopia.under10                                            
# $sets_hyperopia:sets_myopia                                     
# $sets_astigmatism:sets_myopia                                   
# $sets_astigmatism:sets_hyperopia                                
# $sets_myopia.under10:sets_myopia                                
# $sets_hyperopia:sets_myopia.under10                             
# $sets_astigmatism:sets_myopia.under10                           
# $sets_astigmatism:sets_hyperopia:sets_myopia                    
# $sets_hyperopia:sets_myopia.under10:sets_myopia                 
# $sets_astigmatism:sets_hyperopia:sets_myopia.under10            
# $sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia

sets_all_venn_gplots_attr <- attr(sets_all_venn_gplots, "intersections")

### PRINTING ALL THESE CATEGORIES : 

write.table(sets_all_venn_gplots_attr$"sets_astigmatism", 
file="fig.sets_all_venn_gplots_attr.sets_astigmatism.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_hyperopia",
file="fig.sets_all_venn_gplots_attr.sets_hyperopia.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)


write.table(sets_all_venn_gplots_attr$"sets_myopia.under10", 
file="fig.sets_all_venn_gplots_attr.sets_myopia.under10.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_myopia", 
file="fig.sets_all_venn_gplots_attr.sets_myopia.ONLY.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

###

write.table(sets_all_venn_gplots_attr$"sets_hyperopia:sets_myopia", 
file="fig.sets_all_venn_gplots_attr.sets_hyperopia:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_myopia", 
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia", 
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_myopia.under10:sets_myopia", 
file="fig.sets_all_venn_gplots_attr.sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_hyperopia:sets_myopia.under10", 
file="fig.sets_all_venn_gplots_attr.sets_hyperopia:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_myopia.under10",      
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

###

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia",                 
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_hyperopia:sets_myopia.under10:sets_myopia",                
file="fig.sets_all_venn_gplots_attr.sets_hyperopia:sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia.under10",           
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.under10.txt",  quote = FALSE, sep = "\t", row.names=F, col.names=F)

###

write.table(sets_all_venn_gplots_attr$"sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia",
file="fig.sets_all_venn_gplots_attr.sets_astigmatism:sets_hyperopia:sets_myopia.under10:sets_myopia.txt",  quote = FALSE, sep = "\t", row.names=F,  col.names=F)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
