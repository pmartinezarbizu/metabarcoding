
# welcome back ...
###################################################################
#import community table with blast hits
#working directory
setwd('/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio/dadaR/')
load('dada2.RData')

#load libraries
library(dada2)
library(vegan)
library(RFtools)
library(indicspecies) 

#read community table
TTdada <- read.table('Taxon_table.V1V2.dada2.txt',header=FALSE)
TTnames <- c('OTU','pident','qcovs','eval','length','accn','Group','Species',sample.names)
colnames(TTdada) <- TTnames

mf <- TTdada

pdf(file='pident.qcovs.dada2.pdf')
plot(mf$pident,mf$qcovs,xlim=c(70,100),ylim=c(20,100),xlab='percent identity',ylab='query coverage',pch=21,bg='orange')
abline(v=97,h=90)
dev.off()

## define contamination
colnames(mf) <- c('OTU','pident','qcovs','eval','qlength','accn','Group','taxon',sample.names)
unique(mf$Group)
> unique(mf$Group)
 [1] segmented_worms   crustaceans       ribbon_worms      nematodes        
 [5] isopods           cercozoans        hemichordates     hydrozoans       
 [9] bivalves          bryozoans         flatworms         gastrotrichs     
[13] gastropods        loriciferans      animals           protozoa         
[17] spoonworms        eukaryotes        ascomycetes       dinoflagellates  
[21] mesozoans         sea_cucumbers     tunicates         tusk_shells      
[25] golden_algae      solenogasters     jellyfishes       sponges          
[29] centipedes        peanut_worms      crickets          moths            
[33] ciliates          monocots          apicomplexans     tardigrades      
[37] chytrids          green_plants      oomycetes         slime_nets       
[41] mud_dragons       rotifers          diatoms           beetles          
[45] green_algae       chitons           fungi             brachiopods      
[49] springtails       eudicots          basidiomycetes    choanoflagellates
[53] haptophytes       flies             tunicata          soft_corals      
56 Levels: animals apicomplexans ascomycetes basidiomycetes ... tusk_shells

   
cont <- c('beetles','crickets','moths','monocots','slime_nets','centipedes',
'cryptomonads','roaches','chytrids','green_plants','spiders','springtails','conifers',
'eudicots','bugs','flies','haptophytes','scorpions','black_corals','fungi')              

nrow(mf)
#[1] 3615
 
sum(apply(mf[,-c(1:8)],MARGIN=1,FUN=sum))
#[1] 4843916
 
nrow(mf)
#[1] 3615
 
sum(apply(mf[,-c(1:8)],MARGIN=1,FUN=sum))
#[1] 4843916

length(unique(mf$Group))
#[1] 56
  
length(cont)
#[1] 20
 
sum(apply(mf[mf$Group %in% cont,-c(1:8)],MARGIN=1,FUN=sum)) 
#[1] 1543

####### mf2 is the datset without contaminants
mf2 <- mf[!(mf$Group %in% cont),]
nrow(mf2)
[1] 4500

unique(mf2$Group)
 [1] segmented_worms   crustaceans       ribbon_worms      nematodes        
 [5] isopods           cercozoans        hemichordates     hydrozoans       
 [9] bivalves          bryozoans         flatworms         gastrotrichs     
[13] gastropods        loriciferans      animals           protozoa         
[17] spoonworms        eukaryotes        ascomycetes       dinoflagellates  
[21] mesozoans         sea_cucumbers     tunicates         tusk_shells      
[25] golden_algae      solenogasters     jellyfishes       sponges          
[29] peanut_worms      ciliates          apicomplexans     tardigrades      
[33] oomycetes         mud_dragons       rotifers          diatoms          
[37] green_algae       chitons           brachiopods       basidiomycetes   
[41] choanoflagellates tunicata          soft_corals      
56 Levels: animals apicomplexans ascomycetes basidiomycetes ... tusk_shells
 

length(unique(mf2$Group))
#[1] 43
 
sum(apply(mf2[-c(1:8)],MARGIN=1,FUN=sum))
#[1]  4842373

#define not benthic meiofauna
notbemet <- c( 'eukaryotes','dinoflagellates','cercozoans','ascomycetes',
'apicomplexans','mesozoans','ciliates','tunicates', 'choanoflagellates','green_algae',
'cryptomonads','golden_algae','diatoms','animals','oomycetes',
'basidiomycetes','N/A','haptophytes','arrow_worms','glomeromycetes','fungi','jellyfishes',
'arthropods','rotifers','mat_anemones','deep-sea_limpets','sea_spiders',
'micrognathozoans','sea_urchins','starfish','protozoa','bacteria','chitons','tunicata')

#assign isopoda to crustacea
mf2$Group[mf2$Group == 'isopods'] <- 'crustaceans' 


unique(mf3$Group)
 [1] segmented_worms crustaceans     ribbon_worms    nematodes      
 [5] hemichordates   hydrozoans      bivalves        bryozoans      
 [9] flatworms       gastrotrichs    gastropods      loriciferans   
[13] spoonworms      sea_cucumbers   tusk_shells     solenogasters  
[17] sponges         peanut_worms    tardigrades     mud_dragons    
[21] brachiopods     soft_corals  
56 Levels: animals apicomplexans ascomycetes basidiomycetes ... tusk_shells

nrow(mf3)
#[1] 2500
length(unique(mf3$Group))
#[1] 22
sum(apply(mf3[-c(1:8)],MARGIN=1,FUN=sum))
#[1] 4628206

#check length distribution
plot(table(mf3$qlength))
#remove OTUs shorter than 300bp
mf4 <- mf3[mf3$qlength>=300,]

nrow(mf4)
#2310

t1 <- table(mf4$Group) > 0
 
table(mf4$Group)[t1]

       bivalves       bryozoans     crustaceans       flatworms      gastropods 
             28               9             769             134              12 
   gastrotrichs   hemichordates      hydrozoans    loriciferans     mud_dragons 
              8               9              20              10               4 
      nematodes    peanut_worms    ribbon_worms   sea_cucumbers segmented_worms 
            914               1              50               1             322 
    soft_corals   solenogasters         sponges      spoonworms     tardigrades 
              1               1               6               2               8 
    tusk_shells 
              1 


mf4$Group <- droplevels(mf4$Group)
 
levels(mf4$Group)
  [1] "bivalves"        "bryozoans"       "crustaceans"     "flatworms"      
 [5] "gastropods"      "gastrotrichs"    "hemichordates"   "hydrozoans"     
 [9] "loriciferans"    "mud_dragons"     "nematodes"       "peanut_worms"   
[13] "ribbon_worms"    "sea_cucumbers"   "segmented_worms" "soft_corals"    
[17] "solenogasters"   "sponges"         "spoonworms"      "tardigrades"    
[21] "tusk_shells"    
> 


library(ggridges)
library(ggplot2)

pdf(file='ridge_normal.dada.pdf')
#ridges plots make probability density function with length coverage
  ggplot(mf4, aes(x = qlength, y = Group, fill = Group)) +
  geom_density_ridges() +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()

#calculate some grade statistics
grade2 <- (mf4$pident + mf4$qcovs + (mf4$qlength/max(mf4$qlength)*100))/3 
grade_gen <- (mf4$pident*0.5) + (mf4$qcovs*0.25) + ((mf4$qlength/max(mf4$qlength)*100)*0.25)

taxa <- unique(mf4$Group)
grp <- c()
pid <- c()
for(tax in taxa){
prob.ca <- mf4$pident[mf4$Group %in% tax]

        # first calculate lower 5% quantile from ecdf
        ed <- ecdf(prob.ca)
        # delete beta values
        beta1 = NULL
        beta2 = NULL
        # try optimization of empirical beta distribution of the probability of correct
        # assigment
        try(beta1 <- MASS::fitdistr((prob.ca-0.00001)/100, "beta", start = list(shape1 = 1, shape2 = 1),lower=c(0,1)),silent = TRUE)
        # now try to generate 5000 random numbers with the parameters of fitdistr
        try(beta2 <- rbeta(5000, beta1$estimate[1], beta1$estimate[2]), silent = TRUE)
        if(!is.null(beta2)){ed <- ecdf(beta2)}
        # write to list
prob <- runif(500, 1e-04, 1)
vy <- quantile(ed, prob)

grp <- c(grp,rep(tax,500))
pid <- c(pid,vy*100)
}

ridges <- data.frame(Group=grp,pident=pid)

pdf(file='ridges_beta.dada.pdf')
#ridges plots 
  ggplot(ridges, aes(x = pident, y = Group, fill = Group)) +
  geom_density_ridges() +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  geom_density_ridges(scale = 3) +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()

#beta qlength '''IN PROGRESSS'''
###########################################################
for(tax in taxa){
prob.ca <- mf4$qlength[mf4$Group %in% tax]

        # first calculate lower 5% quantile from ecdf
        ed <- ecdf(prob.ca)
        # delete beta values
        beta1 = NULL
        beta2 = NULL
        # try optimization of empirical beta distribution of the probability of correct
        # assigment
        try(beta1 <- MASS::fitdistr(prob.ca, "beta", start = list(shape1 = 1, shape2 = 1),lower=c(0,1)),silent = TRUE)
        # now try to generate 5000 random numbers with the parameters of fitdistr
        try(beta2 <- rbeta(5000, beta1$estimate[1], beta1$estimate[2]), silent = TRUE)
        if(!is.null(beta2)){ed <- ecdf(beta2)}
        # write to list
prob <- runif(500, 1e-04, 1)
vy <- quantile(ed, prob)

grp <- c(grp,rep(tax,500))
pid <- c(pid,vy)
}

ridges <- data.frame(Group=grp,qlength=pid)

pdf(file='ridges_beta_qlength.dada.pdf')
#ridges plots 
  ggplot(ridges, aes(x = qlength, y = Group, fill = Group)) +
  geom_density_ridges() +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  geom_density_ridges(scale = 3) +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()
##############################################################

##color palette from Kevin Wright https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c("#6A3D9A",# purple
	"#E31A1C", # red
	"green4",
	"dodgerblue2", 
	"#FF7F00", # orange
	"black",
	"gold1",
	"skyblue2",
	"#FB9A99", # lt pink
	"palegreen2",
	"#CAB2D6", # lt purple
	"#FDBF6F", # lt orange
	"gray70",
	"khaki2",
	"maroon",
	"orchid1",
	"deeppink1",
	"blue1",
	"steelblue4",
	"darkturquoise",
	"green1",
	"yellow4",
	"yellow3",
	"darkorange4",
	"brown")



##new fuction barp.sample	
###########################################################
barp.sample <- function(x, col, by , fun,...) {
res <- matrix(ncol=1, nrow=length(unique(by)))

for(i in colnames(col)){
res <- cbind(res,aggregate(x[,i] ~ by, x, fun)[2])
}

rnames <- aggregate(x[,i] ~ by, x, fun)[1]
res <- as.matrix(res[,2:ncol(res)]) 
colnames(res) <- names(col)
rownames(res) <- rnames[,1]
res <- res[order(apply(res,1,sum),decreasing=TRUE),]
return(res)
}   
###########################################################

#match colors
#create linking table
colc25 <- data.frame(unique(mf4$Group),c25[1:length(unique(mf4$Group))])

#plots for meiofauna total
pdf('barplots_meiofauna.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(mf4,mf4[,9:27],mf4$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,  col=col, main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

mf4pa <- cbind(mf4[,1:8],decostand(mf4[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(mf4pa,mf4pa[,9:27],mf4pa$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=2,legend=colc25[1:8,1],bty='n',pch=15,col=as.vector(colc25[1:8,2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=2,legend=colc25[1:8,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

legend(-12,-0.4,ncol=2,legend=colc25[9:16,1],bty='n',pch=15,col=as.vector(colc25[9:16,2]),pt.cex=2.5, x.intersp=1 )
legend(-12,-0.4,ncol=2,legend=colc25[9:16,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

legend(8,-0.4,ncol=2,legend=colc25[17:21,1],bty='n',pch=15,col=as.vector(colc25[17:21,2]),pt.cex=2.5, x.intersp=1 )
legend(8,-0.4,ncol=2,legend=colc25[17:21,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )
##########
dev.off()

meio <- mf4


##crustaceans
unique(mf4[mf4$Group=='crustaceans',8])
  

Isopoda <-c('Storthyngurella_menziesi','Macrostylis_sp._ML08','Janirella_sp._JW2004','Eugerdella_sp._VTDes136','Eugerda_sp._JW2004',
'Eurycope_sp._MB_I31'  ) 
 
Tanaidacea <- c('Pseudotanaid')

Harpacticoida <- c('Microarthridion_littorale','Tisbe_furcata','Diarthrodes_sp._2_New_CaledoniaRJH2007','Thalestridae_sp._DZMB531',
 'Eurycletodes_laticauda','Amphiascoides_atopus','Ameira_scotti','Laophontodes_sp._DZMB547','Paramphiascella_fulvofasciata',
'Idyanthidae_sp','Cerviniopsis','Tisbe_furcata','Phyllothalestris_sp._New_CaledoniaRJH2007','Mesochra_rapiens','Parastenhelia_sp._New_CaledoniaRJH2007',
'Stenhelia_sp._GreenlandRJH2007','Microarthridion_fallax','Hase','Leptomesochra_sp','Stratiopontotes_meditrraneus','Dactylopusia_sp._New_CaledoniaRJH2007',
 'Cletodidae_sp._DZMB535','Volkmannia_attennata','Tetragonicepsidae_sp','Idyanthidae_sp._DZMB628','Ameriidae_sp',
'Mesocletodes_sp._DZMB330','Argestigens_sp._GreenlandRJH2007','Xylora_bathyalis','Delavalia_palustris','Ectinosomatidae_sp',
 'Mesochra_sp', 'Pseudotachidius_bipartitus','Itunella_muelleri', 'Parameiropsis_sp','Aegisthidae','Stenocopiinae_sp._DZMB560',
 'Parabradya_dilatata','Cerviniella_hyunsui','Paranannopus_sp','Argestidae_sp._DZMB629','Sarsameira_sp._DZMB468','Cerviniella',
 'Normanellidae_sp','Bradya_sp._GreenlandRJH2004','Ectinosomatidae_sp_2','Laophontodes_sp._DZMB547' )

Cyclopoida <- c('Euryte', 'Oncaeidae_Oncaea_sp3','Speleoithona_bermudensis','Oncaeidae_Oncaea_sp1','Oromiina_fortunata',
'Cyclopicina', 'Cyclopinella', 'Cyclopoida','Cyclopinodes_nsp','Archinotodelphys_sp._New_CaledoniaRJH2001',
'Schminkepinellidae_Einslepinella_mediana','Speleiothona_sp','Cyclopinodes_sp._DZMB060','Euryte_sp._DZMB480','Euryte_sp._New_CaledoniaRJH2004',
'Cerviniopsis_longicaudata','Oithona_sp','Cyclopina', 'Metacyclopina_sp._DZMB717', 'Poecilostomatoida_Erebonasteridae_Centobnaster',
 'Speleoithona_sp','Archinotodelphys_nsp','Cyclopinella_Cyclopinid','Cyclopinella_sp._DZMB477','Paracyclopina_nana')     

Siphonostomatoida <- c('Aphotopontius_mammillatus','Rhizorhina_soyoae','Hatschekia_pagrosomi','Kroyeria_longicauda')

Misophrioida <- c('Expansophria_sp', 'Speleophriopsis_sp','Misophriopsis_sp')

Ostracoda <- c( 'Hirsutocythere_hanaii', 'Cytheromorpha_acupunctata', 'Pontocypris_mytiloides','Bythoceratina_hanejiensis',
 'Polycope_sp._OC2001','Kotoracythere_inconspicua','Cobanocythere_japonica','Uncinocythere_occidentalis','Parapolycope_spiralis',
'Parapolycope_uncata','Spinileberis_quadriaculeata','Cytheropteron_subuchioi','Ostracoda_Conchoecia_sp_91.3','Neomonoceratina_microreticulata',
'Neomonoceratina_microreticulata','Bairdiocopina_sp._New_Zealand'  )

Tantulocarida <- c('Microdajus_tchesunovi' )

Calanoida <- c( 'Mimocalanus_nudus','Acartia_bifilosa','Mimocalanus_crassus', 'Monacilla_typica', 'Paraeuchaeta_norvegica',
'Gaetanus_tenuispinus' )


 	
# divide dataset only crustaceans
crust <- mf4[mf4$Group=='crustaceans',]

# reassign Groups
grcrust <- as.vector(crust$Group)
grcrust[crust$taxon %in% Isopoda] <- 'Isopoda'
grcrust[crust$taxon %in% Calanoida] <- 'Calanoida'
grcrust[crust$taxon %in% Harpacticoida] <- 'Harpacticoida'
grcrust[crust$taxon %in% Siphonostomatoida] <- 'Siphonostomatoida'
grcrust[crust$taxon %in% Misophrioida] <- 'Misophrioida'
grcrust[crust$taxon %in% Cyclopoida] <- 'Cyclopoida'
grcrust[crust$taxon %in% Tantulocarida] <- 'Tantulocarida'
grcrust[crust$taxon %in% Ostracoda] <- 'Ostracoda'
grcrust[crust$taxon %in% Tanaidacea] <- 'Tanaidacea'
crust$Group <- grcrust

#check correct assigment
 table(crust$Group)

        Calanoida        Cyclopoida     Harpacticoida           Isopoda 
               49                82               548                14 
     Misophrioida         Ostracoda Siphonostomatoida        Tanaidacea 
                3                59                 8                 5 
    Tantulocarida 
                1 


> 

#match colors
#create linking table
c25r <- c('cadetblue2','darkolivegreen2','darkred','orange','bisque3','cyan4','chocolate','darkorchid1','brown1')
colc25 <- data.frame(unique(crust$Group),c25r)

pdf(file='barplot_crust.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(crust,crust[,9:27],crust$Group,sum)
col <- c25r[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

crustpa <- cbind(crust[,1:8],decostand(crust[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(crustpa,crustpa[,9:27],crustpa$Group,sum)
col <- c25r[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=5,legend=colc25[1:9,1],bty='n',pch=15,col=as.vector(colc25[1:9,2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=5,legend=colc25[1:9,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

dev.off()

# divide dataset only Harpacticoids Cyclopoids and and Nematodes
hacyne <- mf4[mf4$Group %in% c('crustaceans','nematodes'),]
table(hacyne$Group)

      bivalves       bryozoans     crustaceans       flatworms      gastropods 
              0               0             769               0               0 
   gastrotrichs   hemichordates      hydrozoans    loriciferans     mud_dragons 
              0               0               0               0               0 
      nematodes    peanut_worms    ribbon_worms   sea_cucumbers segmented_worms 
            914               0               0               0               0 
    soft_corals   solenogasters         sponges      spoonworms     tardigrades 
              0               0               0               0               0 
    tusk_shells 

grhacyne <- as.vector(hacyne$Group)             
grhacyne[hacyne$taxon %in% Harpacticoida] <- 'Harpacticoida'
grhacyne[hacyne$taxon %in% Cyclopoida] <- 'Cyclopoida'
grhacyne[grhacyne == 'nematodes'] <- 'Nematoda'
hacyne$Group <- grhacyne

#exclude 'crustaceans'
hacyne <- hacyne[hacyne$Group %in% c('Harpacticoida','Cyclopoida','Nematoda'),]
table(hacyne$Group)

   Cyclopoida Harpacticoida      Nematoda 
           82           548           914 
###

#match colors
#create linking table
c25r <- c(c25[4],'darkolivegreen2','orange')
colc25 <- data.frame(unique(hacyne$Group),c25r)

pdf(file='barplot_Har_Cyc_Nema.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(hacyne,hacyne[,9:27],hacyne$Group,sum)
col <- c25r[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

hacynepa <- cbind(hacyne[,1:8],decostand(hacyne[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(hacynepa,hacynepa[,9:27],hacynepa$Group,sum)
col <- c25r[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=3,legend=colc25[1:3,1],bty='n',pch=15,col=as.vector(colc25[1:3,2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=3,legend=colc25[1:3,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

dev.off()

# divide dataset only Worms
worms <- mf4[mf4$Group %in% c('segmented_worms', 'peanut_worms','ribbon_worms','spoonworms','hemichordates','flatworms'),]

#create linking table
colc25 <- data.frame(unique(mf4$Group),c25[1:length(unique(mf4$Group))])

#plots
pdf('barplots_worms.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(worms,worms[,9:27],worms$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,  col=col, main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

wormspa <- cbind(worms[,1:8],decostand(worms[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(wormspa,wormspa[,9:27],wormspa$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=3,legend=colc25[c(1,3,5,9,13,18),1],bty='n',pch=15,col=as.vector(colc25[c(1,3,5,9,13,18),2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=3,legend=colc25[c(1,3,5,9,13,18),1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )
##########
dev.off()


molluscs <- mf4[mf4$Group %in% c('bivalves','gastropods','solenogasters','tusk_shells'),]

#create linking table
colc25 <- data.frame(unique(mf4$Group),c25[1:length(unique(mf4$Group))])

#plots
pdf('barplots_molluscs.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(molluscs,molluscs[,9:27],molluscs$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,  col=col, main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

molluscsspa <- cbind(worms[,1:8],decostand(molluscs[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(molluscspa,molluscspa[,9:27],molluscspa$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=4,legend=colc25[c(7,11,15,16),1],bty='n',pch=15,col=as.vector(colc25[c(7,11,15,16),2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=4,legend=colc25[c(7,11,15,16),1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )
##########
dev.off()

# we have datasets meio, crust and hacyne, worms, molluscs

#start community analysis meio
area <- as.factor(c(rep('IRA',7),rep('PRA',12)))
colmds <- c(rgb(245,90,50,max=255,alpha=120),rgb(45,225,200,max=255,alpha=120))

#par(mfrow=c(3,2))
#pdf(file='mds_dada2.pdf')

 

##MEIO##############################
dataset <- meio
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')

pdf(file='mds_meiofauna_brayc_dada2.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Meiofauna, Bray-Curtis (dada2)')
legend(-0.4,0.7,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.4,0.7,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_meiofauna_pa_dada2.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Meiofauna, p/a (dada2)')
legend(-10,6,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-10,6,legend=unique(area),pch=1,cex=1.3)
dev.off()

## CRUST ###########################
dataset <- crust
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_crust_brayc_dada2.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Crustacea, Bray-Curtis (dada2)')
legend(-0.6,0.4,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.6,0.4,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_crust_pa_dada2.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Crustacea, p/a (dada2)')
legend(-5,5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-5,5,legend=unique(area),pch=1,cex=1.3)
dev.off()

## HACYNE ###########################
dataset <- hacyne
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_hacyne_brayc_dada2.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Nematoda-Copepoda, Bray-Curtis (dada2)')
legend(-0.4,0.6,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.4,0.6,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_hacyne_pa_dada2.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Nematoda-Copepoda, p/a (dada2)')
legend(-10,6,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-10,6,legend=unique(area),pch=1,cex=1.3)
dev.off()

##
#dev.off()

## WORMS ###########################
dataset <- worms
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_worms_brayc_dada2.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Worms, Bray-Curtis (dada2)')
legend(-0.8,0.5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.8,0.5,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_worms_pa_dada2.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Worms, p/a (dada2)')
legend(-5.5,3,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-5.5,3,legend=unique(area),pch=1,cex=1.3)
dev.off()


## MOLLUSCS ###########################
dataset <- molluscs
dat <- t(dataset[,c(9:16,18:24,26:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area[1:17],dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_molluscs_brayc_dada2.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Molluscs, Bray-Curtis (dada2)')
legend(-0.7,0.65,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.7,0.65,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_molluscs_pa_dada2.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='molluscs, p/a (dada2)')
legend(-2,1.25,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-2,1.25,legend=unique(area),pch=1,cex=1.3)
dev.off()


##adonis permanova table

datasets <- list(meio,crust,hacyne, worms, molluscs)
names <- c('meio','crust','hacyne', 'worms', 'molluscs')
ad.F <- c()
ad.p <- c()
bd.F <- c()
bd.p <- c()
adpa.F <- c()
adpa.p <- c()
bdpa.F <- c()
bdpa.p <- c()


for (i in 1:length(datasets)){

dat <- t(datasets[[i]][,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat2 <- dat[!(apply(dat,sum,MARGIN=1)==0),]
ad <- adonis(dat2 ~ area[!(apply(dat,sum,MARGIN=1)==0) ])
ad.F <- c(ad.F,ad$aov.tab[1,4])
ad.p <- c(ad.p,ad$aov.tab[1,6])
bd <- permutest(betadisper(vegdist(dat2),area[!(apply(dat,sum,MARGIN=1)==0) ]))
bd.F <- c(bd.F,bd[[1]][1,4])
bd.p <- c(bd.p,bd[[1]][1,6])

#now pa
datpa <- decostand(dat2,'pa')
ad <- adonis(datpa ~ area[!(apply(dat,sum,MARGIN=1)==0) ],distance='euclidean')
adpa.F <- c(adpa.F,ad$aov.tab[1,4])
adpa.p <- c(adpa.p,ad$aov.tab[1,6])
bd <- permutest(betadisper(vegdist(datpa,'euclidean'),area[!(apply(dat,sum,MARGIN=1)==0) ]))
bdpa.F <- c(bdpa.F,bd[[1]][1,4])
bdpa.p <- c(bdpa.p,bd[[1]][1,6])
}

res <- data.frame(names,ad.F = round(ad.F,2),ad.p= ad.p,bd.F=round(bd.F,2),bd.p,adpa.F=round(adpa.F,2),adpa.p,bdpa.F=round(bdpa.F,2),bdpa.p)
    names ad.F  ad.p  bd.F  bd.p adpa.F adpa.p bdpa.F bdpa.p
1     meio 1.81 0.001 11.33 0.003   1.84  0.001   2.29  0.147
2    crust 1.82 0.002  7.52 0.013   1.91  0.001   1.44  0.242
3   hacyne 1.83 0.001 12.62 0.002   1.81  0.001   1.61  0.227
4    worms 1.71 0.002  7.65 0.018   1.78  0.002   0.00  0.989
5 molluscs 1.48 0.061  1.79 0.210   1.42  0.100   2.04  0.175
> 

#mvstest
library(RFtools)
dataset <- meio
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
colnames(dat) <- paste('v',1:ncol(dat),sep='')
mv.test <- MVSF.test(dat,area)

> print(mv.test)
      oob.err  null.err P.null  q05.null q95.null smooth.err      p.smooth
IRA 0.4285714 0.8684399  0.999 0.5714286      1.0          0 6.661722e-132
PRA 0.0000000 0.2352352  0.999 0.0750000      0.5          0           NaN
> 
pdf('mvsf.dada2.pdf',height=4,width=10)
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

 Total number of species: 2310
 Selected number of species: 39 
 Number of species associated to 1 group: 39 

 List of species associated to each combination: 

 Group IRA  #sps.  35 
                                                 A      B  stat p.value   
82_Cerviniella                              1.0000 0.8571 0.926   0.005 **
13_Cyclopinella_sp._DZMB477                 0.8941 0.8571 0.875   0.005 **
124_Bathyeurystomina_sp._MC124              0.7149 1.0000 0.846   0.025 * 
866_Janirella_sp._JW2004                    1.0000 0.7143 0.845   0.010 **
35_Acantholaimus_sp._K41                    0.8325 0.8571 0.845   0.030 * 
190_Halalaimus_sp._MC1212                   0.6747 1.0000 0.821   0.040 * 
5_Vieitezia_luzmurubeae                     0.9409 0.7143 0.820   0.025 * 
128_Cerviniella                             1.0000 0.5714 0.756   0.015 * 
533_Anticoma_sp._TCR44                      1.0000 0.5714 0.756   0.010 **
1365_Achromadora_cf_terricola_JH2004        1.0000 0.5714 0.756   0.020 * 
81_Bathyeurystomina_sp._MC124               0.9849 0.5714 0.750   0.010 **
49_Cyclopinella_Cyclopinid                  0.9608 0.5714 0.741   0.010 **
39_Anguillosyllis_capensis                  0.9516 0.5714 0.737   0.045 * 
60_Parabradya_dilatata                      0.9465 0.5714 0.735   0.025 * 
421_Neochromadora_sp._RF6R09                0.9066 0.5714 0.720   0.040 * 
347_Trigonostomum_penicillatum              0.8967 0.5714 0.716   0.025 * 
54_Vieitezia_luzmurubeae                    1.0000 0.4286 0.655   0.020 * 
188_Archinotodelphys_nsp                    1.0000 0.4286 0.655   0.050 * 
199_Halalaimus_sp._MC1212                   1.0000 0.4286 0.655   0.035 * 
480_Normanellidae_sp                        1.0000 0.4286 0.655   0.020 * 
524_Metacyclopina_sp._DZMB717               1.0000 0.4286 0.655   0.040 * 
534_Acantholaimus_sp._K84                   1.0000 0.4286 0.655   0.035 * 
683_Acantholaimus_sp._K84                   1.0000 0.4286 0.655   0.020 * 
686_Gadila_sp._NHMUK_20170076               1.0000 0.4286 0.655   0.020 * 
703_Parabradya_dilatata                     1.0000 0.4286 0.655   0.020 * 
761_Halalaimus_sp._89L30C15                 1.0000 0.4286 0.655   0.020 * 
784_Stenocopiinae_sp._DZMB560               1.0000 0.4286 0.655   0.020 * 
818_Parabradya_dilatata                     1.0000 0.4286 0.655   0.035 * 
888_Acantholaimus_sp._K4                    1.0000 0.4286 0.655   0.040 * 
1038_Cerviniella                            1.0000 0.4286 0.655   0.020 * 
1329_Diarthrodes_sp._2_New_CaledoniaRJH2007 1.0000 0.4286 0.655   0.050 * 
974_Halalaimus_sp._51L26C15                 0.9574 0.4286 0.641   0.040 * 
653_Argonemertes_australiensis              0.9471 0.4286 0.637   0.050 * 
716_Neochromadora_BHMM2005                  0.9353 0.4286 0.633   0.035 * 
1239_Draconema_japonicum                    0.9297 0.4286 0.631   0.035 * 

 Group PRA  #sps.  4 
                                     A      B  stat p.value  
26_Theristus_sp._1268           0.9862 0.7500 0.860    0.02 *
926_Spinileberis_quadriaculeata 0.9545 0.6667 0.798    0.02 *
74_Theristus_sp._TN2consensus   1.0000 0.5833 0.764    0.05 *
6_Laonice_norgensis             0.9999 0.5833 0.764    0.05 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> 



















datlogc <- decostand(log(dat[-1])+1,method='total')
mds <- metaMDS(datlogc[-1])
plot(mds$points,pch=21, bg = as.factor(area))


datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1])
plot(mds$points,pch=21, bg = as.factor(area))


####
# naive bayes model
library(e1071)
library(vegan)
setwd("/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio")
dat <- t(TTdada[,c(9:27)])
area <- as.factor(c(rep('IRA',7),rep('PRA',12)))
dat <- data.frame(area,dat)
mv <- MVSF.test(dat$area,dat[,2:ncol(dat)])
print(mv)
> print(mv)
      oob.err  null.err P.null  q05.null  q95.null smooth.err      p.smooth
IRA 0.4285714 0.8624339  0.999 0.4285714 1.0000000          0 6.661722e-132
PRA 0.0000000 0.2348182  0.999 0.0000000 0.5833333          0           NaN

pdf(file='mvsf.pdf')
plot(mv)
dev.off()


nb <- naiveBayes(area ~.,data=dat)


meio <- seqtab.nochim
colnames(meio) <- 1:ncol(meio)
meiologchord <- decostand(log(meio+1),'total')
colnames(meiologchord) <- 1:ncol(meio)

meionb <- data.frame(area,meiologchord)
colnames(meionb) <- c('area',1:ncol(meio))
nb <- naiveBayes(area ~.,data=meionb)

naiv = predict(nb, meiologchord)

unl = unlist(nb$tables)
unlmat =matrix(unl,nrow=4571,ncol=4,byrow=T)
unlmat = data.frame(unlmat)
colnames(unlmat) <- c('IRA','PRA','IRA_sd','PRA_sd')
unlmat = cbind(OTU=colnames(meiologchord),unlmat)
write.csv(unlmat,file='naive_bayes.csv')

plot(unlmat[,2:3],pch=16)
plot(unlmat[,2:3],pch=21,bg=rgb(200,100,150,max=255,alpha=150),cex=1.7,xlim=c(0,0.005),ylim=c(0,0.005))
abline(0,1,lwd=2,col='red')

##
meiopa <- decostand(meio,'pa')
colnames(meiologchord) <- 1:ncol(meio)

meiopa <- data.frame(area,meiopa)
colnames(meiopa) <- c('area',1:ncol(meio))
nbpa <- naiveBayes(area ~.,data=meiopa)

test <- sample(c(0,1),4571)

rfnames <- c()
for(i in 1:4571){rfnames<-c(rfnames, paste(sample(c(letters,LETTERS),6,replace=TRUE),collapse=''))}


naiv = predict(nbpa, decostand(meiologchord,'pa')

###
rfnames <- c()
for(i in 1:4571){rfnames<-c(rfnames, paste(sample(c(letters,LETTERS),6,replace=TRUE),collapse=''))}

colnames(meionb) <- c('area',rfnames)

rf <- randomForest(area~.,data=meionb)


