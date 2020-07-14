# differential expression

library(dplyr)
library(edgeR)
library(limma)

#folder that contains the output of rsem  isoform level expression estimates

isoRsemDir <- ".../TCGA/Liver_LIHC/01.rsem/02.isoforms_rsem_liver_tcga/"

#folder that contains the output of rsem  gene level expression estimates

geneRsemDir <- ".../TCGA/Liver_LIHC/01.rsem/01.genes_rsem_liver_tcga/"

isoFiles <- list.files(isoRsemDir)
geneFiles <- list.files(geneRsemDir)

# get the 249 tcga sample ids
samples249 <- read.csv("data/tcga_sample_ids_249.csv",stringsAsFactors=F)[,1]
clin <- read.csv("data/afp_clin.csv",stringsAsFactors = F)
clin <- filter(clin,bcr_patient_barcode %in% samples249)
clin$afp_group <- ifelse(test = clin$fetoprotein_outcome_value > 20,yes = "high", no = "low")
afp_low_ids <- clin$bcr_patient_barcode[clin$afp_group == "low"]
afp_high_ids <- clin$bcr_patient_barcode[clin$afp_group == "high"]

fclin <- read.delim("data/TCGA_LIHC_liver_clinical_06282017.txt", 
                    stringsAsFactors=FALSE,header=T) %>%
  filter(bcr_patient_barcode %in% clin$bcr_patient_barcode)

# rsem  differential expression

isoformMat <- read.delim(file = paste(isoRsemDir,samples249[1],"-01A.rsem.isoforms.results",sep=''),
                       h=T,stringsAsFactors=F)[,c("gene_id","transcript_id","expected_count")]
geneMat <- read.delim(file = paste(geneRsemDir,samples249[1],"-01A.rsem.genes.results",sep=''),
                        h=T,stringsAsFactors=F)[,c("gene_id","expected_count")]



colnames(isoformMat) <- c("gene_id","transcript_id",samples249[1])
colnames(geneMat) <- c("gene_id",samples249[1])

for(name in samples249[-1]){ # for each isoform file except the first:
  # read the next sample
  #nextSamp <- read.delim(file = paste(isoRsemDir,name,"-01A.rsem.isoforms.results",sep=''),h=T,
                       #  stringsAsFactors=F)[,c("gene_id","transcript_id","expected_count")]
  nextSampG <- read.delim(file = paste(geneRsemDir,name,"-01A.rsem.genes.results",sep=''),h=T,
                          stringsAsFactors = F)[,c("gene_id","expected_count")]
  
 #nextSampG <- read.delim(file = paste(geneRsemDir,name,"-01A.rsem.genes.results",sep=''),h=T,
  #                        stringsAsFactors = F)[,c("gene_id","TPM")]
#colnames(nextSamp)[3:4] <- c(name,name) # rename the TPM column with the sample ID

  #colnames(nextSamp)[3] <- name

  colnames(nextSampG)[2] <- name

  
  isoformMat <- left_join(isoformMat,nextSamp,by=c("transcript_id","gene_id"))
  geneMat <- left_join(geneMat,nextSampG,by="gene_id")
}


#############################
# limma differential expression section
################################


# limma design matrix

# sample conditionA conditionB
# tcgaxx 1          0
# tcgayy 0          1

# etc


designMat <- mutate(clin,low=ifelse(afp_group == "low",1,0)) %>%
  mutate(high=ifelse(low == 1,0,1)) %>%
  select(bcr_patient_barcode,low,high)

designMat <- left_join(designMat,fclin,by="bcr_patient_barcode") %>%
  select(bcr_patient_barcode,low,high,gender,vascular_invasion)

row.names(designMat) <- designMat[,1]
designMat <- designMat[,-1]

designMatPerm <- designMat %>%
  mutate(low=sample(low,replace=F))%>%
  mutate(high=ifelse(low==1,0,1))

row.names(designMatPerm) <- row.names(designMat)

row.names(geneMat) <- geneMat[,1]
geneMat <- select(geneMat,-gene_id)
# order columsn to be same as design mat rows
geneMat <- select(geneMat,row.names(designMat))

groups <- as.factor(designMat$high)
gender <- as.factor(designMat$gender)
vinv <- as.factor(ifelse(designMat$vascular_invasion %in% c("Micro","Macro"),1,0))
#mvi status factor: tbd

design <- model.matrix(~ groups+gender+vinv)

#DGE

y <- DGEList(counts = geneMat)

A <- rowSums(y$counts)
isExpr <- A > 50
rm(A)

y <- y[isExpr, , keep.lib.size = FALSE]

y <- calcNormFactors(y)

AFP <- clin$afp_group

plotMDS(y, labels=AFP, col=ifelse(AFP=="low","blue","red"),gene.selection="common",
        prior.count = 5)



v <- voom(y, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

summary(decideTests(fit,lfc = 1))


de <- topTable(fit, coef="groups1", n = Inf, sort = "p", p = .05)
de_gender <- topTable(fit,coef="genderMALE",n = Inf, p =.05)

# genes which sig for afp group

afp_sig_genes <- de %>% mutate(ensg = row.names(de)) %>%
  filter(logFC>=1 | logFC <= -1)

gen_sig_genes <- de_gender %>% mutate(ensg = row.names(de_gender)) %>%
  filter(logFC>=1 | logFC <= -1)

#write.csv(x = afp_sig_genes,file = "outputdata/difexp_genes_byafp.csv")
#write.csv(x = gen_sig_genes,file = "outputdata/difexp_genes_bygender.csv")


sum(de$adj.P.Val>.05)
sum(de$adj.P.Val<.05)


afp_ensg = "AFP"
afp_fit_position <- ifelse(rownames(fit)==afp_ensg,"AFP Gene","Other Genes")
sig_gene <- row.names(de)[de$logFC > 1 | de $ logFC < -1]

sig_fit_position <- ifelse(rownames(fit) %in% sig_gene,"fold change > 2","fold change < 2")

plotMD(fit,column=2,status = sig_fit_position,hl.col=c("red","black"),
       main = "DGE by AFP status")

abline(h=c(-1,0,1),col=c("red","darkgray","blue"))

afp_low <- filter(clin,afp_group == "low")$bcr_patient_barcode
afp_high <- filter(clin,afp_group == "high")$bcr_patient_barcode


res <- de
head(res)

# volcano plot
with(res, plot(logFC, -log10(adj.P.Val), pch=20, main="DGE by AFP Status"))

with(subset(res, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="red"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, adj.P.Val<.05 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs=afp_ensg, pch=30, cex=.8,offset = 0.8))
abline(h=-log10(0.05),v=c(-1,1), lty=3,col="blue")

