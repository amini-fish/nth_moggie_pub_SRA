library(dartRverse)
library(devtools)

# set working directory:
setwd("C:/Users/Samue/Desktop/Research/Mogurnda")

#gl <- gl.read.dart(filename="Report_DMog23-7998_4_moreOrders_SNP_2.csv", covfilename="mogcov14_PJU.csv")

#save(gl, file = "Northern_Mogurnda.RData") #Updated as at 26/08/2023

## Lets load in our data - be sure to upload this R file to github so others can use subsetted data for reproducibility where needed

gl <- get(load("C:/Users/Samue/Desktop/Research/Mogurnda/Northern_Mogurnda.RData")); gl

## Just keep our northern individuals of interest 

index <- gl@other$ind.metrics$Northern_ALL !="delete"

gl <- gl[ index, ]; gl

## Quickly visualise how shit our data may be

#gl.smearplot(gl)

################################################################################

## Starting no SNPs is 41384 SNPS from 438 inds 

## Start filtering our SNP data - remember the wise words of Renee Catullo 

gl.report.callrate(gl, method = "loc") # lets set a threshold of 0.90 

gl <- gl.filter.callrate(gl, method = "loc", t = 0.9, v = 5) #7244 Loci left 

gl.report.reproducibility(gl)

gl <- dartR.base::gl.filter.reproducibility(gl, t = 0.99); gl # 6432 loci left 

gl.report.secondaries(gl)

## For PCOA datasets only - phylo + FD needs secondaries

gl <- gl.filter.secondaries(gl); gl #3635 loci left 

gl.report.rdepth(gl)

gl <- gl.filter.rdepth(gl, lower = 15, upper = 200) #2781 loci left

gl.report.maf(gl)

# MAF removes a bucnh of loci....no bueno

#gl <- gl.filter.maf(gl, t = 0.1); gl

gl.report.callrate(gl, method = "ind")

gl <- gl.filter.callrate(gl, method = "ind", threshold=0.85, v=3)

gl <- gl.filter.monomorphs(gl, v=3) # 2207 is final no loci 

## Save the clean data - new starting point unlocked LOL 

#save(gl, file = "mog_nth_filtered.RData") #Updated as at 26/08/2023

# and load it up for future use 

clean_gl <- get(load("C:/Users/Samue/Desktop/Research/Mogurnda/mog_nth_filtered.RData")); clean_gl

##################################################################################

#Assign individuals to putative populations
pop(clean_gl) <- clean_gl@other$ind.metrics$species
table(pop(clean_gl))

#Get the individual het values and chuck them in a csv...then load it back in 

het <- gl.report.heterozygosity(clean_gl, method = "ind")
write.csv(het, file = "24.3.23_Northern_heterozygosity.csv")

#het.ind <- read.csv(file = "24.3.23_Northern_heterozygosity.csv", stringsAsFactors = T)

str(het)
ggplot(data = het, aes(x = Ho)) +
  geom_histogram(bins = 200)

het.ind

## Heterozygosity by population to get rough ideas of Ho and FIS

het.pop <- gl.report.heterozygosity(clean_gl, method = "pop")

write.csv(het.pop, file = "24.3.23_nth_pop_het.csv")
het.pop <- read.csv(file = "24.3.23_nth_pop_het.csv")
View(het.pop)

## Key takeways for me are the very high Fis for all pops with > 10 inds - especially the Kimberley mogs 

###################################################################################

#Run a preliminary PCoA
pc <- gl.pcoa(clean_gl, nfactors = 5)

pcoa.1 <- gl.pcoa.plot(glPca = pc, clean_gl, 
                       pop.labels = "pop",
                       label.size = 1.1,
                       pt.size = 3)


gl.pcoa.plot(glPca = pc, gl, interactive = T)

pcoa.1 <- gl.pcoa.plot(glPca = pc, gl, 
                       pop.labels = "pop",
                       label.size = 1.1,
                       pt.size = 3,
                       y = 2,
                       z = 3
                       )
pcoa.3 <- gl.pcoa.plot(glPca = pc, gl, 
                       pop.labels = "pop",
                       label.size = 1.1,
                       pt.size = 3,
                       y = 3, interactive = T
)



#Based on the PCoA clustering we can seperate the northern range of Mogurnda (i.e. Mogurnda Mogrunda, and Oligolepis into three main groups)
#North West Mogs - Regent, Mitchell, Fotzroy, Frank and Glenelg
#North East - Gulf_East, Gulf_West, Mitchell, Arnhem, Cape_York 
#Northern - Snowdrop, Yali, mog_type, Drysdale_admixed, Bindoola, Bindoola_admixed, Drys_olig_type, mogurnda, Victoria, NT_Victoria admixed


NthxNthWest <- gl.drop.pop(clean_gl, pop.list = c("Gulf_East", "Gulf_West", "Mithcell", "Arnhem", "Cape_York", "larapintae", "thermophila", "clivicola"))

NthxNthWest <- gl.filter.monomorphs(NthxNthWest, v = 5) #1928 loci 264 inds

NthWest.pc <- gl.pcoa(NthxNthWest, nfactors = 5)

Nthwest.plot <- gl.pcoa.plot(glPca = NthWest.pc, NthxNthWest, 
                             pop.labels = "pop",
                             label.size = 1.1,
                             pt.size = 3, 
                             ellipse = F)

#I'm fairly confident I can remove the North Western group from the analysis now so we can focus on the central species although I'm including the Drys and Olig type based on placement in the PCoA (See above)

North <- gl.drop.pop(NthxNthWest, 
                     pop.list = c("Regent", "Mitchell", "Fitzroy", "Frank", "Glenelg", "Drysdale_Admixed", "Bindoola_admixed", "Drys_olig_type", "Snowdrop_admixed"))

save(North, file = "Nth_mog_strict.RData")

North <- get(load("C:/Users/Samue/Desktop/Research/Mogurnda/Nth_mog_strict.RData")); North

rivers <- North@other$ind.metrics$river

North.pc <- gl.pcoa(North, nfactors = 5)

pop(North) <- North$other$ind.metrics$Location

gl.pcoa.plot(North.pc, North, interactive = F)

pop(North) <- North$other$ind.metrics$species

plot2 <- gl.pcoa.plot(North.pc, North, interactive = F)


#Quick FD's 
D <- gl.fixed.diff(North, 
                   test = T, 
                   mono.rm = T, 
                   pb = T,
                   reps = 1000)

D[2]$fd #raw fixed differences
D[3]$pcfd #percentage fixed differences
D[4]$nobs #number of individuals used in each comparison
D[6]$expfpos #expected number of false positives

D$pval

North.merged <- gl.merge.pop(North, old = "Drys_olig_type", new = "Drysdale")
North.merged <- gl.merge.pop(North.merged, old = "mog_type", new = "mogurnda")

D.2 <- gl.fixed.diff(North.merged, 
                   test = T, 
                   mono.rm = T, 
                   pb = T,
                   reps = 1000)

Mog_Vic <- gl.keep.pop(gl,
                       pop.list = c("mogurnda", "Victoria", "NT_Victoria_admixed"))

popNames(Mog_Vic) <- gsub(" ","_",popNames(Mog_Vic))
indNames(Mog_Vic) <- gsub(" ","_",indNames(Mog_Vic))

Mog_Vic_str <- Mog_Vic
indNames(Mog_Vic_str) <- as.character(1:length(indNames(Mog_Vic)))


require(parallel)

detectCores()
cl <- makeCluster(6)

Mog_Vic_str <- gl.filter.callrate(Mog_Vic_str, 
                                  method = "loc", 
                                  threshold = 1)

Mog_Vic_str

STR <- gl.run.structure(Mog_Vic_str, 
                        k.range = 3:8,
                        num.k.rep = 3, 
                        burnin = 5000, 
                        numreps = 5000,
                        noadmix = F,
                        exec = "C:/Users/samue/Downloads/structure_windows_console/console/structure.exe", 
                        plot.out = F, 
                        plot_theme = theme_bw())


ev <- gl.evanno(STR)
ev

qmat <- gl.plot.structure(STR, K= 4, 
                          ind_name = F, 
                          plot_theme = theme_bw(), 
                          clumpak = T,
                          iter_clumpp = 10000)


str(qmat)
qmat$ord <- paste(indNames(Mog_Vic))


sr.k4 <- data.frame(qmat); sr

K4 <- sr.k4[,c(2:5,7,9)]; K4

library(tidyverse)

K4 <- K4[, c(6,5,1,2,3,4)]
K4


K4 <- K4 %>% pivot_longer(cols = 3:6, 
               names_to = "clus", 
               values_to = "value")

K4 <- data.frame(K4)

K4
  
K4.plot <- ggplot(K4, aes(x = ord, y = value, fill = clus)) +
           geom_bar(stat = "identity", show.legend = F) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Dark2") +
  facet_grid("K4" ~ X1.orig.pop, scales = "free", space = "free")
  
  
K4.plot

qmat.7 <- gl.plot.structure(STR, K= 7, 
                          ind_name = F, 
                          plot_theme = theme_bw(), 
                          clumpak = T,
                          iter_clumpp = 10000)

qmat.7$ord <- paste(indNames(Mog_Vic))

sr.k7 <- data.frame(qmat.7); sr.k7

ncol(sr.k7)

K7 <- sr.k7[,c(2:8,10,12)]; K7


K7 <- K7[, c(9,8,1:7)]
K7


K7 <- K7 %>% pivot_longer(cols = 3:9, 
                          names_to = "cluster", 
                          values_to = "value")

K7 <- data.frame(K7)

K7

K7.plot <- ggplot(K7, aes(x = ord, y = value, fill = cluster)) +
  geom_bar(stat = "identity", show.legend = F) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Paired") +
  facet_grid("K7" ~ X1.orig.pop, scales = "free", space = "free") +
  theme(strip.text.x.top = element_text(angle = 90, size = 15))

K7.plot

pc.m.v <- gl.pcoa(Mog_Vic,
                  nfactors = 5)

gl.pcoa.plot(pc.m.v, Mog_Vic, 
             interactive = T)

#It appreas the NT_Victoria moggies are hybrids with "mogurnda" specifically the population including inds MOMoyle.1 (Moyle River) MOMH1627.1 (Reynolds River)


mogurnda <- gl.keep.pop(Mog_Vic, pop.list = "mogurnda")

pc.mog <- gl.pcoa(mogurnda, nfactors = 5)

gl.pcoa.plot(pc.mog, mogurnda, interactive = T)


qmat.6 <- gl.plot.structure(STR, 
                            K = 6, 
                            ind_name = T, 
                            plot_theme = theme_bw(), 
                            clumpak = T,
                            iter_clumpp = 10000)

qmat.6$names <- paste(indNames(Mog_Vic))

sr.k6 <- data.frame(qmat.6); sr.k6

ncol(sr.k6)

K6 <- sr.k6[,c(2:7,9,11)]; K6

library(tidyverse)

K6 <- K6[, c(8,7,1:6)]
K6


K6 <- K6 %>% pivot_longer(cols = 3:8, 
                          names_to = "cluster", 
                          values_to = "value")

K6 <- data.frame(K6)

K6

library(viridis)

K6.plot <- ggplot(K6, aes(x = names, y = value, fill = cluster)) +
  geom_bar(stat = "identity", show.legend = F) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_brewer(palette = "Paired") +
  facet_grid("K6" ~ X1.orig.pop, scales = "free", space = "free") +
  theme(strip.text.x.top = element_text(angle = 90, size = 15))

K6.plot


Gulf <- gl.keep.pop()


North.pc <- gl.pcoa(North, nfactors = 5)
North.plot <- gl.pcoa.plot(glPca = North.pc, North, 
                           pop.labels = "pop",
                           label.size = 1.1,
                           pt.size = 3, 
                           ellipse = F)


North.plot <- gl.pcoa.plot(glPca = North.pc, North, 
                           pop.labels = "pop",
                           label.size = 1.1,
                           pt.size = 3,
                           y = 3)

#Interactive
North.plot <- gl.pcoa.plot(glPca = North.pc, North, 
                           pop.labels = "pop",
                           label.size = 1.1,
                           pt.size = 3, 
                           interactive = T,
                           ellipse = F)


North <- gl.drop.ind(North, ind.list = c("MODrysdale.1", "MODrysdale.2", 
                                "MODrysdaleheadw.1", "MODrysdaleheadw.2"))

Central <- gl.drop.pop(North, 
                       pop.list = c("Gulf_West", "Arnhem", "Cape_York", "Gulf_East", "clivicola"))

Central.pc <- gl.pcoa(Central, nfactors = 5)
Central.plot <- gl.pcoa.plot(glPca = Central.pc, Central, 
                           pop.labels = Central@other$ind.metrics$river,
                           label.size = 1.1,
                           pt.size = 3, 
                           ellipse = F)

Central.interactive <- gl.pcoa.plot(glPca = Central.pc, Central, 
                             pop.labels = "pop",
                             label.size = 1.1,
                             pt.size = 3,
                             interactive = T,
                             ellipse = F)
library(dartR)
t1 <- testset.gl
t1 <- gl.keep.pop(testset.gl,pop.list = popNames(t1)[1:10])
res <- gl.report.pa(t1,loc_names = TRUE)
#getting loc names from fixed differences
fd_loc_names <- lapply(res$names_loci,"[","fd")
#getting trimmed sequences
fd_seq <- lapply(fd_loc_names, function(x){
  as.character(t1$other$loc.metrics[which(locNames(t1) %in% unname(unlist(x))),"TrimmedSequence"])
})


library(dartR.captive)



               
