
# welcome back ...
###################################################################
#import community table with blast hits
#working directory
setwd('/home/pmartinez/projects/metabarcoding/vsearch/mangan16meio/vsearchR/')

#load libraries
library(vegan)
library(RFtools)

#read community table
TTvsearch <- read.table('Taxon_table.V1V2.0.97.txt',header=TRUE)
# NOTE 5 line s were without blast result, removed manually

mf <- TTvsearch

pdf(file='pident.qcovs.vsearch.pdf')
plot(mf$pident,mf$qcovs,xlim=c(70,100),ylim=c(20,100),xlab='percent identity',ylab='query coverage',pch=21,bg='orange')
abline(v=97,h=90)
dev.off()

## define contamination
colnames(mf) <- c('accn', 'Group', 'OTU','taxon','pident','qcovs','eval','qlength',"IRA43c1",   "IRA43c9",   "IRA44c5" ,  "IRA54c12" ,
 "IRA54c5",   "IRA56c1" ,  "IRA56c2" ,  "PRA100c12", "PRA102c11" ,"PRA103c4", 
 "PRA81c4" ,  "PRA84c12" , "PRA85c8" , "PRA86c8"  , "PRA87c8" , "PRA89c8" , 
 "PRA90c11" , "PRA94c3" ,  "PRA97c6")
unique(mf$Group)
> unique(mf$Group)
 [1] segmented_worms   crustaceans       nematodes         cercozoans       
 [5] diatoms           flatworms         ribbon_worms      bivalves         
 [9] goblet_worms      dinoflagellates   tunicates         protozoa         
[13] slime_nets        loriciferans      eukaryotes        hemichordates    
[17] centipedes        springtails       beetles           sponges          
[21] hydrozoans        chytrids          mud_dragons       bryozoans        
[25] tardigrades       chitons           rotifers          animals          
[29] isopods           ciliates          apicomplexans     brachiopods      
[33] mesozoans         gastrotrichs      glomeromycetes    green_algae      
[37] moths             cryptomonads      golden_algae      oomycetes        
[41] gastropods        choanoflagellates monocots          primates         
[45] fungi             basidiomycetes    ascomycetes       birds            
[49] sea_cucumbers     spoonworms        tusk_shells       haptophytes      
[53] sea_anemones      Annelida          peanut_worms      arrow_worms      
[57] phoronid_worms    black_corals      flies             flowering_plants 
[61] tunicata          high_GC_Gram+     solenogasters     soft_corals      
[65] jellyfishes       crickets          sea_urchins      
67 Levels: animals Annelida apicomplexans arrow_worms ... tusk_shells

   
cont <- c('beetles','crickets','moths','monocots','slime_nets','centipedes',
'cryptomonads','roaches','chytrids','green_plants','spiders','springtails','conifers',
'eudicots','bugs','flies','haptophytes','scorpions','black_corals','fungi','primates','birds','flowering_plants','high_GC_Gram+')              

 nrow(mf)
#[1] 4520
 
 sum(apply(mf[,-c(1:8)],MARGIN=1,FUN=sum))
#[1] 6338485

 nrow(mf)
#[1] 4520
  
 sum(apply(mf[,-c(1:8)],MARGIN=1,FUN=sum))
#[1] 6338485

 length(unique(mf$Group))
#[1] 67
  
 length(cont)
#[1] 24
 
 sum(apply(mf[mf$Group %in% cont,-c(1:8)],MARGIN=1,FUN=sum)) 
#[1] 1249


####### mf2 is the datset without contaminants
mf2 <- mf[!(mf$Group %in% cont),]
nrow(mf2)
[1] 4447


unique(mf2$Group)
  [1] segmented_worms   crustaceans       nematodes         cercozoans       
 [5] diatoms           flatworms         ribbon_worms      bivalves         
 [9] goblet_worms      dinoflagellates   tunicates         protozoa         
[13] loriciferans      eukaryotes        hemichordates     sponges          
[17] hydrozoans        mud_dragons       bryozoans         tardigrades      
[21] chitons           rotifers          animals           isopods          
[25] ciliates          apicomplexans     brachiopods       mesozoans        
[29] gastrotrichs      glomeromycetes    green_algae       golden_algae     
[33] oomycetes         gastropods        choanoflagellates basidiomycetes   
[37] ascomycetes       sea_cucumbers     spoonworms        tusk_shells      
[41] sea_anemones      Annelida          peanut_worms      arrow_worms      
[45] phoronid_worms    tunicata          solenogasters     soft_corals      
[49] jellyfishes       sea_urchins      
67 Levels: animals Annelida apicomplexans arrow_worms ... tusk_shells

 

length(unique(mf2$Group))
#[1] 50
 
sum(apply(mf2[-c(1:8)],MARGIN=1,FUN=sum))
#[1]  6337236


#define not benthic meiofauna
notbemet <- c( 'eukaryotes','dinoflagellates','cercozoans','ascomycetes',
'apicomplexans','mesozoans','ciliates','tunicates', 'choanoflagellates','green_algae',
'cryptomonads','golden_algae','diatoms','animals','oomycetes',
'basidiomycetes','N/A','haptophytes','arrow_worms','glomeromycetes','fungi','jellyfishes',
'arthropods','rotifers','mat_anemones','deep-sea_limpets','sea_spiders',
'micrognathozoans','sea_urchins','starfish','protozoa','bacteria','chitons','tunicata','sea_anemones')

#assign isopoda to crustacea
mf2$Group[mf2$Group == 'isopods'] <- 'crustaceans' 
#assign Annelida to segmented worms
mf2$Group[mf2$Group == 'Annelida'] <- 'segmented_worms' 


mf3 <- mf2[!(mf2$Group %in% notbemet),]
unique(mf3$Group)

 [1] segmented_worms crustaceans     nematodes       flatworms      
 [5] ribbon_worms    bivalves        goblet_worms    loriciferans   
 [9] hemichordates   sponges         hydrozoans      mud_dragons    
[13] bryozoans       tardigrades     brachiopods     gastrotrichs   
[17] gastropods      sea_cucumbers   spoonworms      tusk_shells    
[21] peanut_worms    phoronid_worms  solenogasters   soft_corals    
67 Levels: animals Annelida apicomplexans arrow_worms ... tusk_shells
 
 nrow(mf3)
#[1] 3213
 length(unique(mf3$Group))
#[1] 24
 sum(apply(mf3[-c(1:8)],MARGIN=1,FUN=sum))
# 6123963
> 

#check length distribution
plot(table(mf3$qlength))
#remove OTUs shorter than 300bp
mf4 <- mf3[mf3$qlength>=300,]

nrow(mf4)
#3127


t1 <- table(mf4$Group) > 0
 
table(mf4$Group)[t1]

       bivalves     brachiopods       bryozoans     crustaceans       flatworms 
             42              29              13             986             151 
     gastropods    gastrotrichs    goblet_worms   hemichordates      hydrozoans 
             20              11               1               9              29 
   loriciferans     mud_dragons       nematodes    peanut_worms  phoronid_worms 
              8               3             990               4               1 
   ribbon_worms   sea_cucumbers segmented_worms     soft_corals   solenogasters 
            122               1             683               2               1 
        sponges      spoonworms     tardigrades     tusk_shells 
             12               2               6               1 


mf4$Group <- droplevels(mf4$Group)
 
levels(mf4$Group)

 [1] "bivalves"        "brachiopods"     "bryozoans"       "crustaceans"    
 [5] "flatworms"       "gastropods"      "gastrotrichs"    "goblet_worms"   
 [9] "hemichordates"   "hydrozoans"      "loriciferans"    "mud_dragons"    
[13] "nematodes"       "peanut_worms"    "phoronid_worms"  "ribbon_worms"   
[17] "sea_cucumbers"   "segmented_worms" "soft_corals"     "solenogasters"  
[21] "sponges"         "spoonworms"      "tardigrades"     "tusk_shells"    
> 
  



library(ggridges)
library(ggplot2)

pdf(file='ridge_normal.vsearch.pdf')
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

pdf(file='ridges_beta.vsearch.pdf')
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

pdf(file='ridges_beta_qlength.vsearch.pdf')
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
load('colc25.RData')
levels(colc25[,1]) <- c(levels(colc25[,1]),'phoronid_worms','goblet_worms','brachiopods')
levels(colc25[,2]) <- c(levels(colc25[,2]),'orange','grey30','tan3')

colc25 <- rbind(colc25,c('phoronid_worms','orange'))
colc25 <- rbind(colc25,c('goblet_worms','grey30'))
colc25 <- rbind(colc25,c('brachiopods','tan3'))
 
#plots for meiofauna total
pdf('barplots_meiofauna.dada.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(mf4,mf4[,9:27],mf4$Group,sum)
col <- as.vector( colc25[,2][match(rownames(bp),colc25[,1])])
barplot(bp,las=2,  col=col, main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

mf4pa <- cbind(mf4[,1:8],decostand(mf4[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(mf4pa,mf4pa[,9:27],mf4pa$Group,sum)
col <- c25[match(rownames(bp),colc25[,1])]
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=2,legend=colc25[1:10,1],bty='n',pch=15,col=as.vector(colc25[1:10,2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=2,legend=colc25[1:10,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

legend(-12,-0.4,ncol=2,legend=colc25[11:20,1],bty='n',pch=15,col=as.vector(colc25[11:20,2]),pt.cex=2.5, x.intersp=1 )
legend(-12,-0.4,ncol=2,legend=colc25[11:20,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

legend(8,-0.4,ncol=1,legend=colc25[21:24,1],bty='n',pch=15,col=as.vector(colc25[21:24,2]),pt.cex=2.5, x.intersp=1 )
legend(8,-0.4,ncol=1,legend=colc25[21:24,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )
##########
dev.off()

meio <- mf4


##crustaceans
unique(mf4[mf4$Group=='crustaceans',4])
 
Isopoda <-c('Mirabilicoxa_sp._VTDes146','Torwolia_sp._VTDes024','Munneurycope_murrayi', 'Syneurycope_sp._MR2008','Chelator_rugosus','Eurycope_sp._MB_I31','Janirella_sp._JW2004', 'Eugerdella_sp._VTDes011','Macrostylis_sp._ML08','Storthyngurella_menziesi','Ilyarachna_triangulata' ) 
 
Tanaidacea <- c('Pseudotanaid')

Harpacticoida <- c('Mesochra_sp','Tisbe_furcata','Thalestridae_sp._DZMB531','Ectinosomatidae_sp',
 'Sewellia_tropica','Aegisthus','Ameridae_sp','Harpacticus_nipponicus','Lourinia_armata','Malacopsyllus_sp',
'Canthocamptus_staphylinoides','Stenocaris_sp._DZMB019','Superornatiremiidae_sp','Stratiopontotes','Cerviniopsis',
'Euryte_sp._DZMB480','Tachidius_discipes','Cletodidae_sp._DZMB535','sentiropsis_sp',
'Attheyella_crassa','Parastenocarididae_sp._DZMB669','Cerviniopsis_longicaudata',
'Aegisthidae_sp._DZMB331','Canthocamptus_kitaurensis','Idyanthidae_sp._DZMB628','Cerviniella_sp._DZMB428',
'Stenhelia_sp._GreenlandRJH2007','Parameiropsis_sp._DZMB582','Idyanthidae_sp','Tetragonicepsidae_sp',
'Paramphiascella_fulvofasciata','Diarthrodes_sp._2_New_CaledoniaRJH2007','Volkmannia_attennata', 'Bryocamptus_pygmaeus',
'Eurycletodes_laticauda', 'Amphiascoides_atopus','Enhydrosoma_gariene','Zosimeidae_sp','Mesocletodes_sp._DZMB330',
 'Phyllothalestris_sp._New_CaledoniaRJH2007','Hase', 'Parastenhelia_sp._New_CaledoniaRJH2007','Ameira_scotti',
 'Danielsseniidae_sp','Paranannopus_sp','Xylora_bathyalis','Mesochra_rapiens','Dactylopusia_sp._New_CaledoniaRJH2007',
 'Sarsameira_sp._DZMB468','Parameiropsis_sp','Stratiopontotes_meditrraneus','Cerviniella_hyunsui','Itunella_muelleri',
 'Argestigens_sp._GreenlandRJH2007','Argestidae_sp._DZMB62','Stenocopiinae_sp._DZMB560','Microarthridion_fallax',  
 'Pseudotachidius_bipartitus','Ameriidae_sp', 'Aegisthidae', 'Cerviniella','Ectinosomatidae_sp_2', 'Bradya_sp._GreenlandRJH2004',
 'Delavalia_palustris', 'Parabradya_dilatata','Argestidae_sp._DZMB629','Normanellidae_sp')

Cyclopoida <- c('Cyclopina','Oromiina_fortunata','Oithonidae_sp._DZMB624','Speleoithona_sp','Oncaeidae_Oncaea_sp3',
 'Cyclopinella_sp','Cyclopinodes_sp._DZMB060','Metacyclopina_sp._DZMB717','Cyclopicina','Euryte_sp._New_CaledoniaRJH2004',
'Cyclopinella', 'Cyclopinidae_Cyclopinoides','Schminkepinellidae_Einslepinella_mediana', 'Speleiothona_sp','Cyclopinella_sp._DZMB477', 'Paracyclopina_nana','Cyclopinodes_nsp', 'Pterinopsyllus_sp',
'Archinotodelphys_sp._New_CaledoniaRJH2001','Poecilostomatoida_Erebonasteridae_Centobnaster',
'Archinotodelphys_nsp','Oncaeidae_Oncaea_sp1','Speleoithona_bermudensis','Cyclopinella_Cyclopinid')     

Siphonostomatoida <- c('Ecbathyrion_prolixicauda','Ceuthoecetes_sp._RJH2006','Rhizorhina_soyoae','Aphotopontius_mammillatus',
 'Hatschekia_pagrosomi' )

Misophrioida <- c('Speleophriopsis_sp','Palpophria_sp','Expansophria_sp','Misophriopsis_sp')

Ostracoda <- c('Ostracoda_Conchoecia_sp_91.3','Bythoceratina_hanejiensis','Cytheropteron_subuchioi','Kroyeria_longicauda',
 'Spinileberis_quadriaculeata','Cytheromorpha_acupunctata','Parapolycope_uncata', 'Pontocypris_mytiloides','Polycope_sp._OC2001','Cobanocythere_japonica','Kotoracythere_inconspicua','Parapolycope_spiralis',
'Uncinocythere_occidentalis', 'Nenesidea_oligodentata' )

Tantulocarida <- c( 'Arcticotantulus_pertzovi','Microdajus_tchesunovi' )

Calanoida <- c('Eucalanus_elongatus','Gaetanus_environmental_sample','Scolecithrix_danae','Euchaeta_acuta','Aetideopsis_minor',
'Mesocalanus_tenuicornis','Temora_longicornis','Calanoida_sp','Acartia_bifilosa','Neocalanus_cristatus','Paraeuchaeta_norvegica',
 'Mimocalanus_nudus','Gaetanus_tenuispinus', 'Ctenocalanus_vanus', 'Pseudocalanus_elongatus','Mimocalanus_crassus')

Branchiopoda <- c('Lepidurus_packardi')

 	
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
grcrust[crust$taxon %in% Branchiopoda] <- 'Branchiopoda'
crust$Group <- grcrust

#check correct assigment
 table(crust$Group)

     Branchiopoda         Calanoida        Cyclopoida     Harpacticoida 
                1                58               127               717 
          Isopoda      Misophrioida         Ostracoda Siphonostomatoida 
               20                 9                41                10 
       Tanaidacea     Tantulocarida 
                1                 2 
> 



> 

#match colors
#create linking table
#match colors
#create linking table
load('colc25cru.RData')
levels(colc25cru[,1]) <- c(levels(colc25cru[,1]),'Branchiopoda')
levels(colc25cru[,2]) <- c(levels(colc25cru[,2]),'yellow3')

colc25cru <- rbind(colc25cru,c('Branchiopoda','yellow3'))
 

pdf(file='barplot_crust.vsearch.pdf',width=8,height=9)
par(mfrow=c(2,2), oma = c(7, 0, 0, 0),xpd=NA)

bp <- barp.sample(crust,crust[,9:27],crust$Group,sum)
col <- as.vector(colc25cru[,2][match(rownames(bp),colc25cru[,1])])
barplot(bp,las=2,col=col,main='absolute number of reads per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of reads per Taxon')

crustpa <- cbind(crust[,1:8],decostand(crust[,9:27],method='pa',MARGIN=2) )
bp <- barp.sample(crustpa,crustpa[,9:27],crustpa$Group,sum)
col <- as.vector(colc25cru[,2][match(rownames(bp),colc25cru[,1])])
barplot(bp,las=2,col=col,main='absolute number of OTUS per Taxon')
barplot(decostand(bp,method='total',MARGIN=2),las=2,col=col,main='relative number of OTUS per Taxon')

legend(-34,-0.4,ncol=5,legend=colc25cru[1:10,1],bty='n',pch=15,col=as.vector(colc25cru[1:10,2]),pt.cex=2.5, x.intersp=1 )
legend(-34,-0.4,ncol=5,legend=colc25cru[1:10,1],bty='n',pch=0,pt.cex=2.5, x.intersp=1 )

dev.off()

# divide dataset only Harpacticoids Cyclopoids and and Nematodes
hacyne <- mf4[mf4$Group %in% c('crustaceans','nematodes'),]
table(hacyne$Group)

       bivalves     brachiopods       bryozoans     crustaceans       flatworms 
              0               0               0             986               0 
     gastropods    gastrotrichs    goblet_worms   hemichordates      hydrozoans 
              0               0               0               0               0 
   loriciferans     mud_dragons       nematodes    peanut_worms  phoronid_worms 
              0               0             990               0               0 
   ribbon_worms   sea_cucumbers segmented_worms     soft_corals   solenogasters 
              0               0               0               0               0 
        sponges      spoonworms     tardigrades     tusk_shells 
              0               0               0               0 
> 


grhacyne <- as.vector(hacyne$Group)             
grhacyne[hacyne$taxon %in% Harpacticoida] <- 'Harpacticoida'
grhacyne[hacyne$taxon %in% Cyclopoida] <- 'Cyclopoida'
grhacyne[grhacyne == 'nematodes'] <- 'Nematoda'
hacyne$Group <- grhacyne

#exclude 'crustaceans'
hacyne <- hacyne[hacyne$Group %in% c('Harpacticoida','Cyclopoida','Nematoda'),]
table(hacyne$Group)

   Cyclopoida Harpacticoida      Nematoda 
          127           717           990 

###

#match colors
#create linking table
c25r <- c(c25[4],'darkolivegreen2','orange')
colc25 <- data.frame(c('Nematoda','Cyclopoida','Harpacticoida'),c25r)

pdf(file='barplot_Har_Cyc_Nema.vsearch.pdf',width=8,height=9)
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
load('colc25.RData')
levels(colc25[,1]) <- c(levels(colc25[,1]),'phoronid_worms','goblet_worms','brachiopods')
levels(colc25[,2]) <- c(levels(colc25[,2]),'orange','grey30','tan3')

colc25 <- rbind(colc25,c('phoronid_worms','orange'))
colc25 <- rbind(colc25,c('goblet_worms','grey30'))
colc25 <- rbind(colc25,c('brachiopods','tan3'))
 

#plots
pdf('barplots_worms.vsearch.pdf',width=8,height=9)
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
#colc25 <- data.frame(unique(mf4$Group),c25[1:length(unique(mf4$Group))])

#plots
pdf('barplots_molluscs.vsearch.pdf',width=8,height=9)
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

pdf(file='mds_meiofauna_brayc_vsearch.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Meiofauna, Bray-Curtis (vsearch)')
legend(-0.3,-0.4,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.3,-0.4,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_meiofauna_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Meiofauna, p/a (vsearch)')
legend(-10,-7,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-10,-7,legend=unique(area),pch=1,cex=1.3)
dev.off()

## CRUST ###########################
dataset <- crust
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_crust_brayc_vsearch.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Crustacea, Bray-Curtis (vsearch)')
legend(-0.6,-0.15,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.6,-0.15,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_crust_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Crustacea, p/a (vsearch)')
legend(-12,-5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-12,-5,legend=unique(area),pch=1,cex=1.3)
dev.off()

## HACYNE ###########################
dataset <- hacyne
dat <- t(dataset[,c(9:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area,dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_hacyne_brayc_vsearch.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Nematoda-Copepoda, Bray-Curtis (vsearch)')
legend(-0.4,0.5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.4,0.5,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_hacyne_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Nematoda-Copepoda, p/a (vsearch)')
legend(-11,-5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-11,-5,legend=unique(area),pch=1,cex=1.3)
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
pdf(file='mds_worms_brayc_vsearch.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Worms, Bray-Curtis (vsearch)')
legend(-0.6,0.5,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.6,0.5,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_worms_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Worms, p/a (vsearch)')
legend(-9.5,4,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-9.5,4,legend=unique(area),pch=1,cex=1.3)
dev.off()


## MOLLUSCS ###########################
dataset <- molluscs
dat <- t(dataset[,c(9:16,18:24,26:27)])
dat <- decostand(log(dat+1),method='hellinger')
dat <- data.frame(area[1:17],dat)

# bray curtis
#mds <- metaMDS(dat[-1])
mds <- metaMDS(dat[-1],distance='euclidean')
pdf(file='mds_molluscs_brayc_vsearch.pdf')
plot(mds$points,pch=21, bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='Molluscs, Bray-Curtis (vsearch)')
legend(-0.47,0.47,legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(-0.47,0.47,legend=unique(area),pch=1,cex=1.3)
dev.off()

#pa
datpa <- decostand(dat[-1],method='pa')
mds <- metaMDS(datpa[-1],distance='euclidean')
pdf(file='mds_molluscs_pa_vsearch.pdf')
plot(mds$points,pch=21,  bg=colmds[as.factor(area)],cex=2,xlab='',ylab='',main='molluscs, p/a (vsearch)')
legend(2.7,-1, legend=unique(area),pch=16,col=colmds,cex=1.3)
legend(2.7,-1, legend=unique(area),pch=1,cex=1.3)
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

#
