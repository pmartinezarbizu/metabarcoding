
# welcome back ...
###################################################################
#import community table with blast hits
#working directory
#setwd('/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio/vsearchR/')
setwd("C:/Users/pmartinez/Documents/projects/mangan/2018/mangan18")

#load libraries
library(vegan)

#install_github("pmartinezarbizu/RFtools")
library(RFtools)

library(pairwiseAdonis)

#install from github dada2pp
#library(devtools)
#install_github("pmartinezarbizu/dada2pp")
library(dada2pp)

#library(ggridges)
#library(ggplot2)
#library(iNEXT)
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)




#read community table
mf <- read.table('Taxon_table.V1V2.0.97.txt',header=TRUE)
# NOTE 5 line s were without blast result, removed manually

cont <- c('beetles','crickets','moths','monocots','slime_nets','centipedes',
'cryptomonads','roaches','chytrids','green_plants','spiders','springtails','conifers',
'eudicots','bugs','flies','haptophytes','scorpions','black_corals','fungi','primates',
'birds','flowering_plants','high_GC_Gram+','booklice','brown_algae','pelagophytes')              

notbemet <- c( 'eukaryotes','dinoflagellates','cercozoans','ascomycetes',
'apicomplexans','mesozoans','ciliates','tunicates', 'choanoflagellates','green_algae',
'golden_algae','diatoms','animals','oomycetes','basidiomycetes','N/A','arrow_worms',
'glomeromycetes','jellyfishes','arthropods','rotifers','mat_anemones',
'deep-sea_limpets','sea_spiders','micrognathozoans','sea_urchins',
'starfish','protozoa','bacteria','chitons','tunicata',
'sea_anemones','brittle_stars','bony_fishes','horsehair_worms')

# get statistics and reduce file to target group only
st.mf <- stats.tt(mf,mf[,9:45],mf$Group,cont,notbemet)

#write stats to file
write.table(st.mf$stats,file='target.stats.csv',sep=';',row.names=TRUE,col.names=NA)

#keep and rename table as meio
meio <- st.mf$target.table

#check length distribution
pdf(file='distlength.pdf')
plot(table(meio$length),ylab='number of  OTUs',xlab='read length')
dev.off()

#assign isopoda to crustaceans
meio$Group[meio$Group == 'isopods'] <- 'crustaceans'

#create a color table
meiocol <- pal2table(c(as.character(meio$Group),'Branchiopoda','Calanoida','Cyclopoida',
						'Harpacticoida','Isopoda','Misophrioida','Ostracoda',
						'Siphonostomatoida','Tanaidacea','Tantulocarida'),pal='c25bro2')

#Apply possible filters
# 1. remove OTUs shorter than 300bp
meio <- meio[meio$length>=250,]

# 2. remove station with very low number of reads
# first see number of reads per station
#apply(meio[,9:45],FUN=sum,MARGIN=2)

# station X24 has only 821 reads. Remove that station
meio <- meio[,-which(names(meio) == 'X24')]
						
#keep only groups with more than 0 OTUs
t1 <- table(meio$Group) > 0
or <- order(table(meio$Group)[t1],decreasing=TRUE)
tt <- table(meio$Group)[t1][or]
#write table of OTU counts
write.table(data.frame(n_otus=tt,file='OTUsbyGroup.csv',sep=';', row.names=TRUE,col.names=NA)

# create a vector of working areas			  
area <- c("WA1","WA1","WA1","PRZ", "PRZ", "WA1", "PRZ", "PRZ", "PRZ", "PRZ", "WA1", "WA1", "WA1",    
 "WA1", "WA1", "WA1", "WA1", "WA1", "WA1", "WA1", "WA1","WA1", "WA1", "WA1", "WA1", "WA1",    
 "WA1", "WA1", "IRA", "IRA", "IRA", "IRA", "IRA", "IRA", "WA1", "WA1")

#order the vector
or.wa <- order(area)
# order columns of meio according to working area
meio <- cbind(meio[,1:8],meio[,c(9:44)[or.wa]])
#rename colums
colnames(meio)[9:44] <- paste(area[or.wa],'_',names(meio[9:44]),sep='')

pdf(file='totalreads_taxon.pdf')
par(oma = c(4, 0, 0, 0),xpd=NA)
barplot(tt,col=match2table(names(tt),meiocol,'col'),las=2)
dev.off()

bp <-barp.table(meio[,9:44],meio$Group,sum)	 
pdf('barplots_meiofauna.pdf',height=10,width=12)
plot(bp,'vsearchsum',pal2table=meiocol)
dev.off()

## Community analysis
##MEIO##############################

#create color palette for plots
colarea <- pal2table(area,pal='bro2',alpha=200)

dataset <- meio
dat <- t(dataset[,c(9:44)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area[or.wa],dat)

#########BRAY-CURTIS
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')

ad <- adonis(dat[-1] ~ area[or.wa])

ad
Call:
adonis(formula = dat[-1] ~ area[or.wa]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
area[or.wa]  2    0.5106 0.25530  1.9171 0.10409  0.001 ***
Residuals   33    4.3947 0.13317         0.89591           
Total       35    4.9054                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 

pairwise.adonis(dat[-1],area[or.wa])
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 1.571746 0.13582615   0.008      0.024   .
2 IRA vs WA1 1.865395 0.06246007   0.001      0.003   *
3 PRZ vs WA1 2.101890 0.06982586   0.001      0.003   *
> 
be_bray <- permutest(betadisper(vegdist(dat[-1]),area[or.wa] ))

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)  
Groups     2 0.009875 0.0049377 4.605    999  0.016 *
Residuals 33 0.035384 0.0010723                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 
#if permutest(betadisper) is significant then the dispersion is different between groups

> permutest(betadisper(vegdist(dat[-1]),area[or.wa] ),pairwise=TRUE)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)  
Groups     2 0.009875 0.0049377 4.605    999  0.027 *
Residuals 33 0.035384 0.0010723                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
           IRA        PRZ   WA1
IRA            0.55300000 0.001
PRZ 0.53586042            0.132
WA1 0.00025317 0.12891865      
> 


#match colors and symbols
col <- match2table(area[or.wa],colarea,'col')
pch <- match2table(area[or.wa],colarea,'pch')
pch2 <- match2table(area[or.wa],colarea,'pch2')


pdf(file='mds_meiofauna_brayc_vsearch.pdf')
plot(mds$points,pch=pch, bg=col,cex=2,xlab='',ylab='',main='Meiofauna, Bray-Curtis (vsearch)')
legend(0.45,0.35,legend=unique(area),pch=colarea$pch2,col=as.vector(colarea$col),cex=1.3)
legend(0.45,0.35,legend=unique(area),pch=colarea$pch,cex=1.3)
dev.off()

################PRESENCE ABSENCE
#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa,distance='euclidean')
pdf(file='mds_meiofauna_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=col,cex=2,xlab='',ylab='',main='Meiofauna, p/a (vsearch)')
legend(-15,-4,legend=unique(area),pch=colarea$pch2,col=as.vector(colarea$col),cex=1.3)
legend(-15,-4,legend=unique(area),pch=colarea$pch,cex=1.3)
dev.off()

ad_pa <- adonis(datpa ~ area[or.wa])


Call:
adonis(formula = datpa ~ area[or.wa]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
area[or.wa]  2    0.5001 0.25007  2.0288 0.10949  0.001 ***
Residuals   33    4.0675 0.12326         0.89051           
Total       35    4.5677                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pairwise.adonis(datpa,area[or.wa])

> pairwise.adonis(datpa,area[or.wa])
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 1.547532 0.13401410   0.022      0.066    
2 IRA vs WA1 2.087992 0.06939619   0.001      0.003   *
3 PRZ vs WA1 2.158122 0.07156022   0.001      0.003   *

> pairwise.adonis(datpa , area[or.wa],p.adjust.m='holm')
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 1.547532 0.13401410   0.016      0.016   .
2 IRA vs WA1 2.087992 0.06939619   0.001      0.003   *
3 PRZ vs WA1 2.158122 0.07156022   0.001      0.003   *


> pairwise.adonis(datpa , area[or.wa],p.adjust.m='BH')
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 1.547532 0.13401410   0.020     0.0200   .
2 IRA vs WA1 2.087992 0.06939619   0.001     0.0015   *
3 PRZ vs WA1 2.158122 0.07156022   0.001     0.0015   *

permutest(betadisper(vegdist(datpa,method='euclidean'),area[or.wa] ),pairwise=TRUE)

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
Groups     2 16.933  8.4667 4.5144    999  0.018 *
Residuals 33 61.891  1.8755                       
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
          IRA       PRZ   WA1
IRA           0.2140000 0.006
PRZ 0.2453486           0.221
WA1 0.0048939 0.2212497      
> 

#calculate Diversity evenness rarefaction
dat_d <- t(meio[9:44])

dd <- divers.tt(dat_d,area[or.wa],pal=colarea)

pdf(file='meiofauna_divers.pdf')
plot(dd)
dev.off()

sac <- specaccum(dat_d,method='rarefaction')
sac_ira <- specaccum(dat_d[area[or.wa]=='IRA',],method='rarefaction')
sac_prz <- specaccum(dat_d[area[or.wa]=='PRZ',],method='rarefaction')
sac_wa1 <- specaccum(dat_d[area[or.wa]=='WA1',],method='rarefaction')

pdf(file='meiofauna_sac_sites.pdf')
plot(sac,lwd=3,ci.type='bar',ci=2,main='Meiofauna sites-rarefaction curve')
plot(sac_wa1,col=match2table('WA1',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
plot(sac_ira,col=match2table('IRA',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
plot(sac_prz,col=match2table('PRZ',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
legend(28,700,legend=c('total',as.character(colarea$class)),pch='-',col=c('black',as.vector(colarea$col)),lwd=3)
dev.off()

pdf(file='meiofauna_rarefaction_indiviuals.pdf')
plot(sac$individuals,sac$richness,type='l',lwd=3,ylim=c(500,max(sac$richness)),main='Meiofauna rarefaction',xlab='Number of reads',ylab='Number of OTUs')
lines(sac_wa1$individuals,sac_wa1$richness,type='l',col=match2table('WA1',colarea,'col'),lwd=3)
lines(sac_ira$individuals,sac_ira$richness,type='l',col=match2table('IRA',colarea,'col'),lwd=3)
lines(sac_prz$individuals,sac_prz$richness,type='l',col=match2table('PRZ',colarea,'col'),lwd=3)
legend(2300000,1000,legend=c('total',as.character(colarea$class)),pch='-',col=c('black',as.vector(colarea$col)),lwd=3)
dev.off()



##########
######### Rarefied community of 5000 reads
perm <- 20
drar <- c()
for(i in 1:perm){
drar <- rbind(drar, rrarefy(dat_d,5000))
 }
dat <- decostand(log(drar+1),method='hellinger')
dat <- data.frame(rep(area[or.wa],perm),dat)

#BRAY-CURTIS
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')

ad <- adonis(dat[-1] ~ rep(area[or.wa],perm))
> ad

Call:
adonis(formula = dat[-1] ~ rep(area[or.wa], perm)) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
rep(area[or.wa], perm)   2    10.278  5.1392  28.846 0.07447  0.001 ***
Residuals              717   127.738  0.1782         0.92553           
Total                  719   138.016                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 25.56739 0.09700512   0.001      0.003   *
2 IRA vs WA1 26.23115 0.04202154   0.001      0.003   *
3 PRZ vs WA1 32.78671 0.05197749   0.001      0.003   *
>  

#match colors and symbols
colrar <- match2table(rep(area[or.wa],perm),colarea,'col')
pchrar <- match2table(rep(area[or.wa],perm),colarea,'pch')
pch2rar <- match2table(rep(area[or.wa],perm),colarea,'pch2')


pdf(file='mds_meiofauna_brayc_vsearch_rare.pdf')
plot(mds$points,pch=pchrar, bg=colrar,cex=2,xlab='',ylab='',main='Meiofauna, Bray-Curtis (vsearch) rarefied(5000,x20)')
legend(0.8,0.4,legend=unique(area),pch=colarea$pch2,col=as.vector(colarea$col),cex=1.3)
legend(0.8,0.4,legend=unique(area),pch=colarea$pch,cex=1.3)
dev.off()

##Rarefied P/A

datparar <- decostand(drar,method='pa')
mds <- metaMDS(datparar,distance='euclidean')
pdf(file='mds_meiofauna_pa_vsearch_rar.pdf')
plot(mds$points,pch=pchrar, bg=colrar,cex=2,xlab='',ylab='',main='Meiofauna, p/a (vsearch) rarefied(5000,x20)')
legend(-10.5,7.5,legend=unique(area),pch=colarea$pch2,col=as.vector(colarea$col),cex=1.3)
legend(-10.5,7.5,legend=unique(area),pch=colarea$pch,cex=1.3)
dev.off()

adparar <- adonis(datparar ~ rep(area[or.wa],perm))

Call:
adonis(formula = datparar ~ rep(area[or.wa], perm)) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                        Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
rep(area[or.wa], perm)   2     9.859  4.9296  28.137 0.07277  0.001 ***
Residuals              717   125.615  0.1752         0.92723           
Total                  719   135.475                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 

pairwise.adonis(datparar,rep(area[or.wa],perm),sim.method='euclidean')
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 18.05613 0.07051630   0.001      0.003   *
2 IRA vs WA1 21.12628 0.03412272   0.001      0.003   *
3 PRZ vs WA1 20.33454 0.03288599   0.001      0.003   *
> 


#veen diagrams

venn1 <- aggregate(dat_d,list(area[or.wa]),sum)
vennpa <- decostand(venn1[,-1],'pa')
rownames(vennpa) <- venn1[,1]
vc <-  vennCounts(t(vennpa))
pdf(file='vennd_meiofauna.pdf')
vennDiagram(vc,circle.col=match2table(colnames(vc),colarea,'col'),main='Shared OTUs between sites')
dev.off()



#############################
####just the great groups
gg <- t(bp$table)
dat <- decostand(log(gg+1),method='hellinger')
dat <- data.frame(area[or.wa],dat)

#########BRAY-CURTIS
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')

ad <- adonis(dat[-1] ~ area[or.wa])

Call:
adonis(formula = dat[-1] ~ area[or.wa]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
area[or.wa]  2   0.02653 0.013263 0.87888 0.05057  0.612
Residuals   33   0.49800 0.015091         0.94943       
Total       35   0.52453                  1.00000       
> 
> pairwise.adonis(dat[-1],area[or.wa])
       pairs    F.Model          R2 p.value p.adjusted sig
1 IRA vs PRZ 0.09401154 0.009313596   0.996      1.000    
2 IRA vs WA1 1.48571164 0.050387512   0.138      0.414    
3 PRZ vs WA1 0.60020378 0.020985997   0.810      1.000    
#############################################

################
####NEMATODES AND CRUSTACEANS

#reduce dataset
necru <- meio[meio$Group %in% c('nematodes','crustaceans'),]
dat_nc <- t(necru[9:44])
#veen diagrams with nema&crust only

venn_nc <- aggregate(dat_nc,list(area[or.wa]),sum)
vennpa_nc <- decostand(venn_nc[,-1],'pa')
rownames(vennpa_nc) <- venn_nc[,1]
vc <-  vennCounts(t(vennpa_nc))

pdf(file='vennd_necrust.pdf')
vennDiagram(vc,circle.col=match2table(colnames(vc),colarea,'col'),main='Shared OTUs between sites: Nematodes and Crustaceans')
dev.off()

##Diversity
ddnc <- divers.tt(dat_nc,area[or.wa],pal=colarea)

#write results
write.table(ddnc,file='diversity.metabarcoding.csv',sep=';',row.names=TRUE,col.names=NA)


pdf(file='necrust_n2s.pdf')
plot(ddnc$N,ddnc$S,pch=ddnc$pch,bg=as.vector(ddnc$col),cex=1.5,xlab='Number of reads',ylab='Number of OTUs',main='Nematoda-Crustacea')
dev.off()



pdf(file='necrus_divers.pdf')
plot(ddnc)
dev.off()

sac <- specaccum(dat_nc,method='rarefaction')
sac_ira <- specaccum(dat_nc[area[or.wa]=='IRA',],method='rarefaction')
sac_prz <- specaccum(dat_nc[area[or.wa]=='PRZ',],method='rarefaction')
sac_wa1 <- specaccum(dat_nc[area[or.wa]=='WA1',],method='rarefaction')

pdf(file='necrus_sac_sites.pdf')
plot(sac,lwd=3,ci.type='bar',ci=2,main='Nematodes-Crustaceans sites-rarefaction curve')
plot(sac_wa1,col=match2table('WA1',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
plot(sac_ira,col=match2table('IRA',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
plot(sac_prz,col=match2table('PRZ',colarea,'col'),add=TRUE,lwd=3,ci.type='bar',ci=2)
legend(28,500,legend=c('total',as.character(colarea$class)),pch='-',col=c('black',as.vector(colarea$col)),lwd=3)
dev.off()

pdf(file='necrus_rarefaction_indiviuals.pdf')
plot(sac$individuals,sac$richness,type='l',lwd=3,ylim=c(500,max(sac$richness)),main='Nematodes-Crustaceans rarefaction',xlab='Number of reads',ylab='Number of OTUs')
lines(sac_wa1$individuals,sac_wa1$richness,type='l',col=match2table('WA1',colarea,'col'),lwd=3)
lines(sac_ira$individuals,sac_ira$richness,type='l',col=match2table('IRA',colarea,'col'),lwd=3)
lines(sac_prz$individuals,sac_prz$richness,type='l',col=match2table('PRZ',colarea,'col'),lwd=3)
legend(600000,800,legend=c('total',as.character(colarea$class)),pch='-',col=c('black',as.vector(colarea$col)),lwd=3)
dev.off()

#species turnover
library(betapart)
benc <- betapart.core(dat_nc)
bencmulti <- beta.multi(benc)
benc_sample <- beta.sample(benc,site=10,sample=4000)

dist.nc <- benc_sample$sampled.values

plot(density(dist.nc$beta.SOR), xlim=c(0,0.8), ylim=c(0, 19), xlab='Beta diversity', main='', lwd=3)
lines(density(dist.nc$beta.SNE), lty=1, lwd=2)
lines(density(dist.nc$beta.SIM), lty=2, lwd=2)
#divide by group and plot together.
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00224.x
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12029

#community nema crust
dat <- decostand(log(dat_nc+1),method='hellinger')
dat <- data.frame(area[or.wa],dat)

#########BRAY-CURTIS
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')

ad <- adonis(dat[-1] ~ area[or.wa])
all:
adonis(formula = dat[-1] ~ area[or.wa]) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
area[or.wa]  2    0.4383 0.21915  1.7957 0.09815  0.001 ***
Residuals   33    4.0273 0.12204         0.90185           
Total       35    4.4656                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 
> pairwise.adonis(dat[-1],area[or.wa])
       pairs  F.Model         R2 p.value p.adjusted sig
1 IRA vs PRZ 1.427708 0.12493386   0.022      0.066    
2 IRA vs WA1 1.632814 0.05510155   0.003      0.009   *
3 PRZ vs WA1 2.101998 0.06982918   0.001      0.003   *
> 
> be_bray <- permutest(betadisper(vegdist(dat[-1]),area[or.wa] ))
> be_bray

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.005420 0.0027099 1.7331    999  0.182
Residuals 33 0.051598 0.0015636                     

pdf(file='mds_nema_crust_brayc_vsearch.pdf')
plot(mds$points,pch=pch, bg=col,cex=2,xlab='',ylab='',main='Nematoda-Crustacea, Bray-Curtis (vsearch)')
legend(0.45,0.35,legend=unique(area),pch=colarea$pch2,col=as.vector(colarea$col),cex=1.3)
legend(0.45,0.35,legend=unique(area),pch=colarea$pch,cex=1.3)
dev.off()







#Multivariate structure test
###
library(RFtools)
dataset <- meio
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
colnames(dat) <- paste('v',1:ncol(dat),sep='')
mv.test <- MVSF.test(area,dat)

> mv.test
      oob.err  null.err P.null   q05.null q95.null smooth.err      p.smooth
IRA 0.4285714 0.8002288  0.999 0.42857143      1.0          0 6.661722e-132
PRA 0.0000000 0.2759426  0.999 0.08333333      0.5          0           NaN
> 


pdf('mvsf.vsearch.pdf',height=4,width=10)
plot(mv.test)
dev.off()




##multipatt
library(indicspecies) 

 area <- as.factor(c(rep('IRA',7),rep('PRA',12)))
 dataset <- meio
 dat <- data.frame(t(dataset[,c(9:27)]))
 colnames(dat) <- paste(1:ncol(dat),'_',meio$taxon,sep='')


 mp <- multipatt(dat,area)

> summary(mp,indvalcomp=TRUE)

 Multilevel pattern analysis
 ---------------------------

 Association function: IndVal.g
 Significance level (alpha): 0.05

 Total number of species: 3127
 Selected number of species: 98 
 Number of species associated to 1 group: 98 

 List of species associated to each combination: 

 Group IRA  #sps.  89 
                                                   A      B  stat p.value   
3057_Sosane_wahrbergi                         1.0000 1.0000 1.000   0.005 **
33_Cerviniella                                0.9881 1.0000 0.994   0.005 **
1603_Anguillosyllis_capensis                  0.9690 1.0000 0.984   0.005 **
310_Ameriidae_sp                              0.9505 1.0000 0.975   0.010 **
2872_Torwolia_sp._VTDes024                    0.9332 1.0000 0.966   0.005 **
2769_Bathyeurystomina_sp._MC124               0.8958 1.0000 0.946   0.010 **
474_Cyclopinella_sp._DZMB477                  0.8689 1.0000 0.932   0.010 **
1829_Bathyarca_glomerula                      1.0000 0.8571 0.926   0.005 **
2775_Crenopharynx_sp._MC123                   0.9999 0.8571 0.926   0.005 **
1109_Aricidea_fragilis                        0.9998 0.8571 0.926   0.015 * 
2791_Argestidae_sp._DZMB629                   0.9121 0.8571 0.884   0.045 * 
1446_Halalaimus_sp._BC056                     0.8926 0.8571 0.875   0.005 **
2864_Cervonema_sp._28L20C15                   0.8877 0.8571 0.872   0.005 **
928_Caulleriella_parva                        0.8818 0.8571 0.869   0.015 * 
520_Cyclopinella_Cyclopinid                   0.8797 0.8571 0.868   0.020 * 
528_Cerviniella                               0.8679 0.8571 0.863   0.010 **
100_Cyclopinella_Cyclopinid                   0.7350 1.0000 0.857   0.010 **
421_Acantholaimus_sp._K41                     1.0000 0.7143 0.845   0.005 **
721_Parabradya_dilatata                       1.0000 0.7143 0.845   0.005 **
1349_Cerviniella                              1.0000 0.7143 0.845   0.005 **
1633_Xylora_bathyalis                         1.0000 0.7143 0.845   0.005 **
2532_Emplectonema_gracile                     1.0000 0.7143 0.845   0.010 **
413_Xylora_bathyalis                          0.9998 0.7143 0.845   0.010 **
228_Archinotodelphys_nsp                      0.8219 0.8571 0.839   0.050 * 
535_Theristus_sp._TN2consensus                0.9803 0.7143 0.837   0.005 **
1119_Bathyeurystomina_sp._MC124               0.6925 1.0000 0.832   0.020 * 
643_Parabradya_dilatata                       0.6876 1.0000 0.829   0.040 * 
2541_Parabradya_dilatata                      0.9445 0.7143 0.821   0.015 * 
2986_Acantholaimus_sp._K41                    0.9375 0.7143 0.818   0.020 * 
64_Anguillosyllis_capensis                    0.9335 0.7143 0.817   0.015 * 
1481_Stradorhynchus_sp._WW2004                0.7731 0.8571 0.814   0.020 * 
1534_Cerviniella                              0.9236 0.7143 0.812   0.030 * 
1177_Acantholaimus_sp._K41                    0.9231 0.7143 0.812   0.010 **
290_Stenocopiinae_sp._DZMB560                 0.9068 0.7143 0.805   0.050 * 
212_Acantholaimus_sp._K41                     0.8802 0.7143 0.793   0.010 **
2814_Rhabdocoma_sp._TCR125                    0.8733 0.7143 0.790   0.035 * 
1327_Cyclopinella                             0.8727 0.7143 0.790   0.025 * 
1522_Parabradya_dilatata                      0.8522 0.7143 0.780   0.040 * 
1121_Parabradya_dilatata                      0.8372 0.7143 0.773   0.020 * 
90_Nerillidium_sp._KW2005                     0.8352 0.7143 0.772   0.025 * 
1300_Halalaimus_sp._51L26C15                  0.8327 0.7143 0.771   0.040 * 
292_Troglochaetus_beranecki                   1.0000 0.5714 0.756   0.015 * 
833_Aphelochaeta_marioni                      1.0000 0.5714 0.756   0.015 * 
1220_Ectinosomatidae_sp_2                     1.0000 0.5714 0.756   0.025 * 
1248_Zosimeidae_sp                            1.0000 0.5714 0.756   0.015 * 
1261_Nematoda_environmental_sample            1.0000 0.5714 0.756   0.010 **
1430_Acantholaimus_sp._M137                   1.0000 0.5714 0.756   0.015 * 
1450_Thoracostomopsidae_sp._BCA14             1.0000 0.5714 0.756   0.015 * 
2756_Halalaimus_sp._MC121                     1.0000 0.5714 0.756   0.015 * 
2774_Halalaimus_sp._89L30C15                  1.0000 0.5714 0.756   0.015 * 
2996_Astomonema_sp._NCM2006                   1.0000 0.5714 0.756   0.020 * 
3111_Speleiothona_sp                          1.0000 0.5714 0.756   0.005 **
2788_Nuculana_pernula                         0.9999 0.5714 0.756   0.025 * 
1230_Gadila_sp._NHMUK_20170076                0.9991 0.5714 0.756   0.015 * 
1532_Cerviniella                              0.9903 0.5714 0.752   0.015 * 
1399_Cyclopinella_sp._DZMB477                 0.9668 0.5714 0.743   0.020 * 
1485_Argestigens_sp._GreenlandRJH2007         0.9496 0.5714 0.737   0.035 * 
1451_Acantholaimus_sp._K4                     0.7595 0.7143 0.737   0.035 * 
544_Anticoma_sp._TCR44                        0.9488 0.5714 0.736   0.010 **
1283_Anguillosyllis_capensis                  0.9437 0.5714 0.734   0.025 * 
703_Cyclopinella_Cyclopinid                   0.9421 0.5714 0.734   0.035 * 
407_Acantholaimus_sp._K41                     0.9369 0.5714 0.732   0.025 * 
665_Acantholaimus_sp._K4                      0.9300 0.5714 0.729   0.035 * 
3024_Acantholaimus_sp._M137                   0.9050 0.5714 0.719   0.025 * 
2773_Eurycletodes_laticauda                   0.9005 0.5714 0.717   0.020 * 
2974_Acantholaimus_sp._K4                     0.8961 0.5714 0.716   0.030 * 
1216_Trigonostomum_penicillatum               0.8896 0.5714 0.713   0.045 * 
162_Acantholaimus_sp._M137                    0.8842 0.5714 0.711   0.025 * 
73_Axinella_verrucosa                         1.0000 0.4286 0.655   0.045 * 
204_Quasitetrastemma_stimpsoni                1.0000 0.4286 0.655   0.040 * 
282_Tharyx_sp._THS2012                        1.0000 0.4286 0.655   0.040 * 
718_Stradorhynchus_sp._WW2004                 1.0000 0.4286 0.655   0.030 * 
731_Vieitezia_luzmurubeae                     1.0000 0.4286 0.655   0.040 * 
739_Acantholaimus_sp._K41                     1.0000 0.4286 0.655   0.040 * 
741_Cerviniella                               1.0000 0.4286 0.655   0.035 * 
1260_Metacyclopina_sp._DZMB717                1.0000 0.4286 0.655   0.035 * 
1316_Anguillosyllis_capensis                  1.0000 0.4286 0.655   0.030 * 
1338_Acantholaimus_sp._K41                    1.0000 0.4286 0.655   0.030 * 
1419_Lumbrineris_inflata                      1.0000 0.4286 0.655   0.045 * 
1604_Macrostylis_sp._ML08                     1.0000 0.4286 0.655   0.030 * 
2661_Halalaimus_sp._BC056                     1.0000 0.4286 0.655   0.030 * 
3016_Baicalellia_canadensis                   1.0000 0.4286 0.655   0.035 * 
475_Crenopharynx_sp._MC123                    0.9997 0.4286 0.655   0.045 * 
910_Gyratrix_hermaphroditus                   0.9995 0.4286 0.654   0.045 * 
2523_Stradorhynchus_sp._WW2004                0.7325 0.5714 0.647   0.040 * 
161_Acantholaimus_sp._K41                     0.9730 0.4286 0.646   0.030 * 
2897_Halalaimus_sp._51L26C15                  0.9626 0.4286 0.642   0.035 * 
399_Cerviniella                               0.9536 0.4286 0.639   0.040 * 
2945_Schminkepinellidae_Einslepinella_mediana 0.9075 0.4286 0.624   0.045 * 

 Group PRA  #sps.  9 
                                        A      B  stat p.value   
1191_Amphicorina_mobilis           0.9998 1.0000 1.000   0.005 **
403_Micrura_ignea                  0.9998 0.9167 0.957   0.040 * 
2605_Theristus_sp._1268            0.9891 0.9167 0.952   0.010 **
1269_Perarella_schneideri          1.0000 0.7500 0.866   0.030 * 
2204_Argestidae_sp._DZMB629        0.9873 0.7500 0.860   0.045 * 
1270_Cytheromorpha_acupunctata     1.0000 0.5833 0.764   0.045 * 
2976_Kotoracythere_inconspicua     1.0000 0.5833 0.764   0.035 * 
2880_Barantolla_lepte              0.9997 0.5833 0.764   0.050 * 
2991_Metacyatholaimus_sp._90L31C15 0.9437 0.5833 0.742   0.050 * 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> 










rfnames <- c()
for(i in 1:4571){rfnames<-c(rfnames, paste(sample(c(letters,LETTERS),6,replace=TRUE),collapse=''))}

colnames(meionb) <- c('area',rfnames)

rf <- randomForest(area~.,data=meionb)

datlogc <- decostand(log(dat[-1])+1,method='total')
mds <- metaMDS(datlogc[-1])
plot(mds$points,pch=21, bg = as.factor(area))


datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1])
plot(mds$points,pch=21, bg = as.factor(area))
