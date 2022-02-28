library(foreach)
library(dplyr)
data.X<- read.table('/storage/home/x/xxw5315/GEUVADIS_chr20_processed.traw', header=TRUE)
data.Y<- read.table('/storage/home/x/xxw5315/GEUVADIS_normalized_expression_chr20.txt', header = TRUE)
data.X<- as.data.frame(data.X)
rownames(data.Y)<- data.Y$gene_id

###NMF
expr<- as.data.frame(t(data.Y[,-(1:4)]))
nmf.expr<- as.matrix(expr-min(expr))
nmf<- nmf(nmf.expr, 50, method = 'brunet', seed=1)
exprs.nmf<- fitted(nmf)
exprs.nmf<- t(exprs.nmf)

#VQ
set.seed(100)
expr<- as.data.frame((data.Y[,-(1:4)]))
vq<- kmeans(expr, 10, nstart=1)
View(vq$centers)
exprs.kmeans<- vq$centers[vq$cluster,]

#PCA
expr<- as.data.frame(t(data.Y[,-(1:4)]))
pca<- prcomp(expr, center = FALSE, scale=FALSE)
loading<- pca$rotation
exprs.pca<- as.data.frame(pca$x[,1:80] %*% t(loading[,1:80]))

require(doMC)
nrow(data.Y[,-(1:4)])
registerDoMC(cores=6)
eqtl<- foreach(i = 1:(nrow(data.Y[,-(1:4)]))) %dopar%{
  genotype<- data.X %>% filter(CHR==data.Y$chromosome[i] & POS >= (data.Y$start[i]-1e5) & POS <=(data.Y$end[i]+1e5))
  exprss<- as.data.frame(exprs.pca[,i])
  colnames(exprss)<- 'expression'
  exprss$sample<- colnames(data.Y[,-(1:4)])
  exprss$sample<- paste0(exprss$sample, '_', exprss$sample)
  lmmodel= foreach(j = 1:nrow(genotype))%dopar%{
    geno<- as.data.frame(t(genotype[j,-(1:6)]))
    colnames(geno)<- 'geno'
    geno$sample<- row.names(geno)
    rownames(geno)<- NULL
    data<- merge(exprss, geno, by = 'sample')
    
    model<- summary(lm(expression~geno, data=data))
    
    if(nrow(model$coefficients)>1){
      
      est<- data.frame(gene=data.Y$gene_id[i], snp=genotype$SNP[j], pvalue=model$coefficients[2,4])
      return(est)}else{ est<- data.frame(gene=data.Y$gene_id[i], snp=genotype$SNP[j], pvalue=NA)
      return(est)}
    
  }
  lmmodel$pvalue.adj<- p.adjust(lmmodel$pvalue, method = 'BH')
  return(lmmodel)
}
save(eqtl, file='/storage/home/x/xxw5315/results/pca.results.RData')
