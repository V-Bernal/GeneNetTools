#===========================
# - Title/Project: GeneNetTools: Tests for Gaussian graphical models with shrinkage
# - Author            : Victor Bernal, Venustiano Soancatl-Aguilar, Jonas Bulthuis, Victor Guryev, Peter Horvatovich, Marco Grzegorczyk.
# - Date created      : 22 JUN 2022
#------------------
# Description:
# This scripts contains the analysis for GeneNetTools: Tests for Gaussian graphical models with shrinkage
#------------------
# Variables:
# - data = data matrix with p columns and n rows
# - p = number of variables (e.g. genes)
# - n = sample size
# - lambda = shrinkage value
# - number= Montecarlo iterations
# - rep = Times to repeat the simulation (for fixed p, n)
# - etaA = proportion of TP
# - alpha = significance level
#-----------------
# Notes: 
# - the data is a n x p matrix of data 
# - the p columns are the random variables 
# - the n rows are the samples
#-----------------
# References
# [1] Schafer,J. and Strimmer,K. A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Stat. Appl. Genet. Mol. Biol.(2005a), 4, 1175-1189.
# [2] Bernal, Victor, et al. "Exact hypothesis testing for shrinkage-based Gaussian graphical models." Bioinformatics 35.23 (2019): 5011-5017.
#===========================

sessionInfo()

#============================
# Shrunk tests
#============================
rm(list=ls())

library(GeneNetTools)
library(GeneNet)
library(ggplot2)

#============================
# simulate a network and data
#============================
p = 100; n = 30;
true.pcor <- GeneNet::ggm.simulate.pcor(num.nodes = p, etaA = 0.2)
dat.sim <- GeneNet::ggm.simulate.data(sample.size = n, true.pcor)

#============================
# GeneNet [1] and GeneNet-tools retrieve same correlations
#============================
estimated.cor.genenet <- cor.shrink(dat.sim) # vector
estimated.cor <- corr.shrunk(dat.sim) # vector
data.frame(estimated.cor,(sm2vec(estimated.cor.genenet)))

#============================
# GeneNet [1] and GeneNeTools retrieve same partial correlations
#============================
estimated.pcor.genenet <- ggm.estimate.pcor(dat.sim) # matrix
estimated.pcor <- GGM.shrunk(dat.sim) # vector

data.frame(estimated.pcor,(sm2vec(estimated.pcor.genenet)))

#==========================
# GeneNet [1] and GeneNeTools are compatible
#============================
# same attributes
attributes(estimated.pcor)
attributes(estimated.pcor.genenet)

# same optimal shrinkage in both packages
GeneNetTools::optimal.shrinkage( x = dat.sim )
corpcor::estimate.lambda(x = dat.sim)

#============================
# p values: test H_0: pcor = 0
#============================

# p-values using the shrunk probability density [2] 
pval.s <- pval.shrunk(x = estimated.pcor)

# p-values using montecarlos [2]
pval.monte <- pval.montecarlo(r = estimated.pcor,
                              number = 100)

# p-values using the new implementation Equation 6
pval.ttest <- ttest.shrunk(x = estimated.pcor)

#===========================
# The new t-test has the same results but 10 times faster
#===========================
# Same result 
data.frame(pval.s , pval.ttest)

# time t -test
start_time.new <- Sys.time()
    ttest.shrunk(x = estimated.pcor)
end_time.new <- Sys.time()
end_time.new - start_time.new

# time integral of prob density
start_time.old <- Sys.time()
  pval.shrunk(x = estimated.pcor)
end_time.old <- Sys.time()
end_time.old - start_time.old

# difference
(as.numeric(end_time.old - start_time.old)/
    as.numeric(end_time.new - start_time.new))

#============================
# Confidence intervals partial correlations Equation 9
#===========================
ci.pcor <- confint.GGM(x = estimated.pcor, alpha = 0.05)
head(ci.pcor)

#============================
# Compare 2 partial correlations Equation 10
#===========================
estimated.pcor1 <- GGM.shrunk(dat.sim, lambda = 0.3)
estimated.pcor2 <- GGM.shrunk(dat.sim, lambda = 0.6)

z <- compare.GGM(x1 = estimated.pcor1 ,
                 x2 = estimated.pcor2 )
head(z)

plot( abs(z[,'z-score']), pch=20, ylim=c(0,4) , cex=0.5)
abline(h=1.96)

#============================
library(manipulate)

mani = function(l2 = NULL){
  estimated.pcor1 <- GGM.shrunk(dat.sim, lambda = 0.1)
  estimated.pcor2 <- GGM.shrunk(dat.sim, lambda = l2)
  z <- compare.GGM(x1 = estimated.pcor1 ,x2 = estimated.pcor2 )
  plot( (z[,'z-score']), ylim=c(-4,4), pch=20, cex=0.75, col = rgb(l2,0,0.5,0.75))
}

manipulate( mani(l2) , l2 = slider(min = 0.1,max =  0.95))


#============================
# Results
#============================

#============================
# Z is normal Equation 10
#============================
pp = 100; nn = 30;

true.pcor <- GeneNet::ggm.simulate.pcor(num.nodes = pp, etaA = 0.00)

d1 <- GeneNet::ggm.simulate.data(sample.size = nn, true.pcor)
d2 <- GeneNet::ggm.simulate.data(sample.size = nn, true.pcor)

pcor1 <- GGM.shrunk( d1 )
pcor2 <- GGM.shrunk( d2 )

z <- compare.GGM(pcor1, pcor2)

qqnorm(y = z[,'z-score'])
qqline(y = z[,'z-score'], col = 2, lwd = 3)


#=============================
# ROC and PR curves
#=============================
rm(list=ls())

library(GeneNetTools)
library(GeneNet)
library(ggplot2)

set.seed(1)

#==========
# Source functions for ROC and PR curves
simple_roc<-function(labels,scores){
  labels<-labels[order(scores,decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels),FPR=cumsum(!labels)/sum(!labels),labels)
}

simple_PR<-function(labels,scores){
  labels<-labels[order(scores,decreasing=TRUE)]
  TPR=cumsum(labels)/sum(labels)
  FPR=cumsum(!labels)/sum(!labels)
  PPV=cumsum(labels)/sum(labels)/((cumsum(labels)/sum(labels))+(cumsum(!labels)/sum(!labels)))
  data.frame(TPR, PPV, labels)
}

simple_auc<-function(TPR,FPR){
  #inputs already sorted, best scores first
  dFPR<-c(diff(FPR),0)
  dTPR<-c(diff(TPR),0)
  sum(TPR*dFPR)+sum(dTPR*dFPR)/2
}

#===========================
# Create two different networks of p nodes. 
# Simulate data with varying sample size
pp <- 100; nn <- seq(from = 30, to = 100, by = 10);
true.pcor <- GeneNet::ggm.simulate.pcor(num.nodes = pp, etaA = 0.01)
true.pcor2 <- GeneNet::ggm.simulate.pcor(num.nodes = pp, etaA = 0.03)

# ROC (or PR)
ZZ <- sapply(X = nn,
             FUN = function(X){

               dat1 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor)
               pcor1 <- GGM.shrunk(dat1)

               while(attributes(pcor1)$lambda == 1){
                 dat1 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor)
                 pcor1 <- GGM.shrunk(dat1)
                 }

               sapply(X = nn,
                      FUN = function(X){

                        dat2 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor2)
                        pcor2 <- GGM.shrunk(dat2)

                        while(attributes(pcor2)$lambda == 1){
                          dat2 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor2)
                          pcor2 <- GGM.shrunk(dat2)
                        }

                        z <- compare.GGM(pcor1, pcor2)
                        z <- abs(z[,'z-score'])
                        
                        #==========
                        # ROC
                        #==========
                        #roc <- simple_roc(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                        #                  scores = z)
                        #out <- simple_auc(TPR = roc$TPR,
                        #                  FPR = roc$FPR)
                        #==========
                        # PR (Un-comment for PR)
                        #==========
                         roc <- simple_PR(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                                         scores = z)
                         out <- simple_auc(TPR = roc$PPV,
                                          FPR = roc$TPR)
                        #==========
                        
                        return( out )

                        }
               )
             }
)

image(x = 1:ncol(ZZ), y = 1:nrow(ZZ), z = t(ZZ) , zlim = c(0,1),
      col = 'white', axes = FALSE, main = 'AUPRC', ylab = 'n2', xlab = 'n1')
axis(1, 1:ncol(ZZ), paste(nn), las = 1)# Update main title 
axis(2, 1:nrow(ZZ), paste(nn), las = 1)
for (x in 1:ncol(ZZ))
  for (y in 1:nrow(ZZ))
    text(x, y, round(ZZ[y,x], 2),cex = 0.75)

#====================
# Compare against DiffNet.FDR
#===================
#
# Install DiffNet.FDR
# Step 1. Install the devtools package. Invoke R and then type
#install.packages("devtools")

# Step 2. Load the devtools package.
#library("devtools")

# Step 3. Install the DiffNetFDR package from GitHub.
#install_github("Zhangxf-ccnu/DiffNetFDR", subdir="pkg")

#====================
# ROC GeneNeTools minus DiffNet.FDR
#===================
library(DiffNetFDR)
library(GeneNet)
library(GeneNetTools)

set.seed(1)

pp <- 100; nn <- seq(from = 30, to = 100, by = 10);
true.pcor <- GeneNet::ggm.simulate.pcor(num.nodes = pp, etaA = 0.03)
true.pcor2 <- GeneNet::ggm.simulate.pcor(num.nodes = pp, etaA = 0.01)

ZZ <- sapply(X = nn,
             FUN = function(X){
               
               dat1 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor) 
               pcor1 <- GGM.shrunk(dat1)
               
               while(attributes(pcor1)$lambda == 1){
                 dat1 <- GeneNet::ggm.simulate.data(sample.size = X, true.pcor) 
                 pcor1 <- GGM.shrunk(dat1)
               }
               
               sapply(X = nn,
                      FUN = function(X){
                        
                        dat2 <- GeneNet::ggm.simulate.data(sample.size = X, 
                                                           true.pcor2) 
                        pcor2 <- GGM.shrunk(dat2)                             
                        
                        while(attributes(pcor2)$lambda == 1){
                          dat2 <- GeneNet::ggm.simulate.data(sample.size = X, 
                                                             true.pcor2) 
                          pcor2 <- GGM.shrunk(dat2) 
                        }
                        
                        # W DiffNetFDR
                        dat3 <- rbind(dat1, dat2)
                        
                        w_stat = DiffNet.FDR(dat3 ,
                                             c(rep('1', nrow(dat1)),rep( '2', nrow(dat2))), 
                                             alpha = 0.05, test.type = "pcor")
                        w_stat = abs(sm2vec(w_stat$W))
                        
                        # Z shrunk
                        z <- compare.GGM(pcor1, pcor2)
                        adj.pval <- abs(z[,'z-score'])
                        
                        #--------------------------
                        # ROC
                        rocS <- simple_roc(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                                           scores = adj.pval)
                        rocW <- simple_roc(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                                           scores = w_stat)
                        outS <- simple_auc(TPR = rocS$TPR,
                                           FPR = rocS$FPR)
                        outW <- simple_auc(TPR = rocW$TPR,
                                           FPR = rocW$FPR)
                        #--------------------------
                        
                        #--------------------------
                        # Un-comment Precision Recall 
                        # rocS <- simple_PR(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                        #                   scores = adj.pval)
                        # rocW <- simple_PR(labels = c(sm2vec(true.pcor2) != sm2vec(true.pcor)),
                        #                 scores = w_stat)
                        # outS <- simple_auc(TPR = rocS$PPV,
                        #                    FPR = rocS$TPR)
                        # outW <- simple_auc(TPR = rocW$PPV,
                        #                    FPR = rocW$TPR)
                        #--------------------------
                        
                        return( (outS - outW)/outW )
                      }
               )
             }
)

image(x = 1:ncol(ZZ), y = 1:nrow(ZZ), z = t(ZZ) , zlim = c(0,1),
      col = 'white', axes = FALSE, main = 'Improved AUROC (%)', ylab = 'n2', xlab = 'n1')
axis(1, 1:ncol(ZZ), paste(nn), las = 1)
axis(2, 1:nrow(ZZ), paste(nn), las = 1)
for (x in 1:ncol(ZZ))
  for (y in 1:nrow(ZZ))
    text(x, y, round(ZZ[y,x], 2),cex = 0.75)
#paste(100*round(ZZ[y,x], 2), '%')
#legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
#       c("0-0.4","0.5","0.6-1"), fill = cm.colors(3), xpd = NA, cex =0.5)

#============================
# Real data: E. coli
#============================
rm(list=ls())

library(GeneNetTools)
library(GeneNet)
library(ggplot2)

data(ecoli)
colnames(ecoli)

# partial correlation pcor
estimated.pcor.ecoli <- GGM.shrunk(ecoli)

# p-values
pval.s.ecoli <- pval.shrunk(x = estimated.pcor.ecoli)
pval.tts.ecoli <- ttest.shrunk(x = estimated.pcor.ecoli)

# time t -test
start_time.new <- Sys.time()
  ttest.shrunk(x = estimated.pcor.ecoli)
end_time.new <- Sys.time()
end_time.new - start_time.new

# time integral of prob density
start_time.old <- Sys.time()
  pval.shrunk(x = estimated.pcor.ecoli)
end_time.old <- Sys.time()
end_time.old - start_time.old

# difference 
(as.numeric(end_time.old - start_time.old)/
    as.numeric(end_time.new - start_time.new))


#--------------------

# Both find same edges at 0.05
sum(pval.tts.ecoli < 0.05)
sum(pval.s.ecoli < 0.05)
sum(which(pval.tts.ecoli < 0.05) %in% which(pval.s.ecoli < 0.05))

plot(pval.s.ecoli, pval.tts.ecoli)
abline(0,1, col=2, lwd =3)

#--------------------
# Bland Altman E. coli
diff.pval <- ( pval.tts.ecoli - pval.s.ecoli)
ave.pval <- 0.5*(pval.tts.ecoli + pval.s.ecoli)
df = data.frame('diffpval' = diff.pval, 'avepval' = ave.pval)
colnames(df) <- c('diffpval','avepval')

ggplot(data=df, aes(y = diffpval, x = avepval)) +
  geom_point() +
  labs(title='', x='ave p-values', y = 'diff p-values')+
  geom_abline(slope = 0, intercept = 0, color='black', linetype='dashed', alpha=.5) +
  theme(panel.border = element_rect(fill = "transparent", color = "black"), panel.background = element_rect(fill = "white"))

#--------------------
# Confidence intervals
ci.ecoli <- confint.GGM(x = estimated.pcor.ecoli , alpha = 0.05)
ci.ecoli

sum(abs(ci.ecoli[,1])>0.1)
sum(abs(ci.ecoli[,2])>0.1)

ci.ecoli <- ci.ecoli[order(ci.ecoli[,1], decreasing = T),]
ci.ecoli

# CI forest plot
ci.ecoli.df = data.frame(ci.ecoli[ci.ecoli[,2]>0, ])

# un-comment for the complete list
ci.ecoli.df = ci.ecoli.df[1:20,]

ggplot(data=ci.ecoli.df, aes(y= 1:nrow(ci.ecoli.df),
                             x=ci.ecoli.df[,1], xmin=ci.ecoli.df[,2],
                             xmax=ci.ecoli.df[,3])) +
  geom_point() +
  geom_errorbarh(height=.4) +
  labs(title='', x='scaled pcors (95% CI)', y = '') +
  geom_vline(xintercept=c(0,.1,.3), color='black', linetype='dashed', alpha=.5) +
  scale_y_continuous(name = "", breaks=1:nrow(ci.ecoli.df), labels=row.names(ci.ecoli.df))+
  theme(panel.border = element_rect(fill = "transparent", color = "black"), panel.background = element_rect(fill = "white"))


# PCA across probsets
pc = prcomp(x = ecoli, scale = T)
barplot(pc$sdev^2/sum(pc$sdev^2))
biplot(x = pc,cex=0.5)

#========================
# Real data 2: M. musculus (Botomly)
#========================
rm(list=ls())

#install.packages("BiocManager")
#BiocManager::install('biomaRt' )
#BiocManager::install('Biobase' )
#BiocManager::install('limma' )

Packages <- c("GeneNet",
              "ggplot2",
              "Biobase",
              "igraph",
              "stats4",
              "devEMF",
              "reshape",
              "STRINGdb",
              "limma")


library(Biobase)
library(limma)
library(GeneNetTools)
library(GeneNet)
library(ggplot2)

cohen_criteria = 0.1
pval.cutoff = 0.05

# load data
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)

# Preprocessing
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)

# filter low mean, low median per strain
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)

# first 10 vs last 11
pdata_bot$strain

# Differential expression at p- value threshold with limma
mod = model.matrix( ~ pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)

# adjust p value
limma_output = topTable(ebayes_limma,
                        number = dim(edata)[1],
                        adjust.method="BH",
                        sort="none")
names(limma_output)
limma_pvals_adj = limma_output$adj.P.Val
limma_pvals_adj[1:10]

hist(limma_pvals_adj, col = 2)
hist(p.adjust(ebayes_limma$p.value[, 2], 'BH'), add=T, col='blue')
plot(limma_pvals_adj, p.adjust(ebayes_limma$p.value[, 2], 'BH'))

sum(limma_pvals_adj < pval.cutoff)

genes = as.integer(limma_pvals_adj < pval.cutoff)
names(genes) = rownames(edata)
sum(is.na(genes))
not_na = !is.na(genes)
genes = genes[not_na]
head(genes)
sum(genes)


p = sum(genes)
n = ncol(edata)
names = names(genes)[which(limma_pvals_adj < pval.cutoff)]

##  Upper Quantile normalization
# Normalisation is dividing the values by the column-based upper quartile statistic.
# The matrix needs to be transposed first (with columns as genes) to get the arithmetic

# Un-comment for quantile normalization
#norm_edata = normalize.quantiles(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ]))
#norm_edata = t(normalize.quantiles(t(as.matrix(norm_edata))))
boxplot(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ]))
data.quantileExpressed <- apply(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ])
                                , 2, function(x){quantile(x, 0.75)});

data.norm <- t(t(as.matrix(edata[ which(limma_pvals_adj < pval.cutoff) , ])) / data.quantileExpressed);
boxplot(as.matrix(data.norm))

## Strain effects
plot(data.norm[1, ],col=as.numeric(pdata_bot$strain)) # level shift?
svd1 = svd(data.norm - rowMeans(data.norm))
plot(svd1$v[, 1],svd1$v[, 2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata_bot$strain))

## Split data
data1 = data.norm[ , pdata_bot$strain == 'C57BL/6J' ]
data2 = data.norm[ , pdata_bot$strain == 'DBA/2J' ]

## Remove zero variance probe
sum(apply(data1,1,sd)==0)
sum(apply(data2,1,sd)==0)

idcc = which(apply(data1,1,sd)==0)
data1 = data1[-idcc, ]
data2 = data2[-idcc, ]

# PCA across probesets
temp.dat= cbind(data1,data2)
colnames(temp.dat) = pdata_bot$strain
rownames(temp.dat) = NULL
pc = prcomp(x=  t(temp.dat), scale = T)
barplot(pc$sdev^2/sum(pc$sdev^2))
biplot(x = pc,cex=0.5)

# Heatmap
heatmap(x = temp.dat, col = cm.colors(256),
        Rowv = T, Colv = NA, scale = 'row',
        labCol = pdata_bot$strain,
        ColSideColors = c('black', 'grey')[1+(pdata_bot$strain=='C57BL/6J')])
#legend(x= 2 , y = 0, legend=c("min", "ave", "max"), cex = 1,fill=cm.colors(3))

# Clustering
d = dist(x = t(temp.dat))
hc = hclust(d)
plot(hc, cex=0.5, hang = -1)
rect.hclust(hc, k = 2, border = 2:5)

#------------------------------
# Network analysis
#------------------------------
estimated.pcor.1 <- GGM.shrunk(t(data1))
estimated.pcor.2 <- GGM.shrunk(t(data2))

ci.1 <- confint.GGM(x = estimated.pcor.1 , alpha = 0.05)
ci.2 <- confint.GGM(x = estimated.pcor.2 , alpha = 0.05)

# compare
z <- compare.GGM(x1 = estimated.pcor.1 ,
                 x2 = estimated.pcor.2  )
cut = 1.96
id = ( z[ , 'z-score'] > cut) + 2*( z[ , 'z-score'] < -cut)
table(id)

#---------------
# BiomaRt
#---------------
library(biomaRt)

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
mouse_gene_ids = names
all_new_gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                      filters = 'ensembl_gene_id', values = mouse_gene_ids, mart = ensembl)
all_new_gene

# Un-comment to check the ordering
#data.frame(all_new_gene,all_new_gene$ensembl_gene_id )

#---------------------
# Bland Altman plot
#---------------------
# Matrix index of the top pcors
id0 = order(abs(z[ , 'z-score']), decreasing = T)[1:9]
id0

idxs = matrix(data = 1:(221*221), nrow = 221, ncol = 221)
idxs2 = matrix(data = idxs %in% id0,
               nrow = 221, ncol = 221)
idxs3 = which(idxs2, arr.ind = T)

# 9 top names
res = data.frame( 'probe1'= all_new_gene$external_gene_name[idxs3[,1]],
                  'probe2' = all_new_gene$external_gene_name[idxs3[,2]] )
# plot
df = data.frame('diff' = (estimated.pcor.1 - estimated.pcor.2),
                'ave' = 0.5*(estimated.pcor.1 + estimated.pcor.2))
colnames(df) <- c('diff' , 'ave' )

plot( x = df$ave[id00], y = df$diff[id00],
      pch = 20, cex= 0.5 ,
      cex.lab = 1,  cex.axis = 1,
      ylim = c(-0.1,0.1),
      col = c('black'),
      xlab = 'ave pcors B6 D2',
      ylab = 'diff pcors B6 D2', las = 2)
abline(a = 0, b = 0)
text(x = df$ave[id0] ,
     y = df$diff[id0]*1.1 ,
     labels = paste(res[,1],res[,2],sep = '-'),
     cex = 0.75, adj = 1)

#---------------------
# Our edges include 4 previously reported genes
# 'Alad', 'Glo1', 'Gabra2', 'Cox7a2l'
#---------------------
id00 = which(abs( z[ , 'z-score']) > cut)
idxs = matrix(data = 1:(221*221), nrow = 221, ncol = 221)
idxs2 = matrix(data = idxs %in% id00,
               nrow = 221, ncol = 221)
idxs3 = which(idxs2, arr.ind = T)

res = data.frame( 'probe1'= all_new_gene$external_gene_name[ idxs3[,1] ],
                  'probe2' = all_new_gene$external_gene_name[ idxs3[,2] ] )


c('Alad', 'Glo1', 'Gabra2', 'Cox7a2l') %in% res[,1]
c('Alad', 'Glo1', 'Gabra2', 'Cox7a2l') %in% res[,2]

#-----------------------------
# B6 only
cut = 1.96
sum((z[ , 'z-score'])>cut )/length(z[ , 'z-score'])
which((z[ , 'z-score'])>cut )
idxs = matrix(data = 1:(221*221), nrow = 221, ncol = 221)
idxs2 = matrix(data = idxs %in% which((z[ , 'z-score']) > cut ), nrow = 221,
               ncol = 221)
idxs3 = which(idxs2, arr.ind = T)

res = data.frame( 'probe1'= all_new_gene$external_gene_name[idxs3[,1]],
                  'probe2' = all_new_gene$external_gene_name[idxs3[,2]] )

c('ALAD', 'GLO1', 'GABRA2', 'COX7A2L') %in% res
#write(x = unique(c(res[,1],res[,2])), file = 'MMUSB6only.txt')
#-----------------------------

#-----------------------------
# D2 only
sum((z[ , 'z-score'])< -cut )/length(z[ , 'z-score'])
which((z[ , 'z-score'])< -cut )
idxs = matrix(data = 1:(221*221), nrow = 221, ncol = 221)
idxs2 = matrix(data = idxs %in% which((z[ , 'z-score'])< -cut ), nrow = 221, ncol = 221)
idxs3 = which(idxs2, arr.ind = T)

res = data.frame( 'probe1'= all_new_gene$external_gene_name[idxs3[,1]],
                  'probe2' = all_new_gene$external_gene_name[idxs3[,2]] )

c('ALAD', 'GLO1', 'GABRA2', 'COX7A2L') %in% res
#write(x = unique(c(res[,1],res[,2])), file = 'MMUSD2only.txt')
#-----------------------------

#------------------
# Graph
#------------------
library(igraph)

pcor1 = abs(estimated.pcor.1)/(1-(attr(x = estimated.pcor.1, 'lambda')) )
pcor2 = abs(estimated.pcor.2)/(1-(attr(x = estimated.pcor.2, 'lambda')) )

summary(pcor1)
summary(pcor2)

cut = 1.3*1.96
id = ( z[ , 'z-score'] > cut) + 2*( z[ , 'z-score'] < -cut)
table(id)

PCOR = vec2sm(id)
diag(PCOR) = 1

g.pcor = graph_from_adjacency_matrix(
  PCOR ,
  mode = c("undirected"),
  weighted = NULL,
  diag = FALSE,
  add.colnames = NA,
  add.rownames = NA
)

l <- layout_nicely(g.pcor,dim = 2)

plot(g.pcor,
     edge.width = 2,
     edge.color = c(rgb(0.5, 0, 0,0.5), rgb(0,0,0.5,0.5))[ id] ,
     vertex.color = "black", vertex.size = 2,
     vertex.frame.color = "black",
     vertex.label.color = "black",
     vertex.label.cex = 0.45,
     vertex.label.dist = 0.5,
     vertex.label = all_new_gene$external_gene_name,
     layout = l)
legend('bottomright', c("only C57BL", "only 6J-DBA/2J"), pch=20,
       col= c(rgb(0.5, 0, 0,0.5), rgb(0,0,0.5,0.5)),
       pt.bg = c(rgb(0.5, 0, 0,0.5), rgb(0,0,0.5,0.5)),
       pt.cex=0.75, cex = 0.5, bty = "n", ncol = 1)

#floor(which(abs(z$z.score)>1.96) / 223)
#(which(abs(z$z.score)>1.96) %% 223) +1

#Row = I / LineLenght; // Integer division
#Column = I % LineLenght; // Reminder of the division of I by LineLenght

