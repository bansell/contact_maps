# Contact maps

#Author: Brendan Ansell
#email: ansell.b@wehi.edu.au

library(tidyverse); 
library(plotly)
library(Rpdb)
#library(rgl)
#library(gplots);
#library(ggplot2)
#library(data.table)
#devtools::install_github('ropensci/plotly') #https://community.plot.ly/t/ggplotly-giving-error-unused-argument-environment/3538/2

#TBC: recode angst. distance as variable

setwd("~/Dropbox/protein_modelling/contact_maps/bundle/")

rm(list=ls()); gc()

#user input:
ang_thresh <- 10

set.seed(1234)
nonSyn <- sample(125,10,replace=FALSE) #vector of mutated/substituted AAs, by position

#hit_pdb <- "2WLB";
hit_pdb <- "6VOB"; 
#hit_pdb <- "4O75";


########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ##########

#AAs_20 <- readRDS('AminoAcids_abv.Rds')

hit_chain <- 'A'

system(paste0("curl -s https://files.rcsb.org/view/", hit_pdb, ".pdb > ",hit_pdb, ".pdb"))

pdb_import <- Rpdb::read.pdb(paste0(hit_pdb,".pdb")); 
names(pdb_import$atom)

head(pdb_import$atoms)
pdb_import$atoms %>% filter(elename == "CA") %>% filter(chainid == hit_chain) %>% select(x1,x2,x3) %>% head()

meta <- pdb_import$atoms %>% filter(elename == "CA") %>% filter(chainid == hit_chain) %>% select(resid,resname) %>% distinct()
xyz  <- pdb_import$atoms %>% filter(elename=="CA") %>% select(x1,x2,x3)

min(meta$resid)

# #remove MODEL line from file:
# system("cat model1.pdb | grep -v MODEL | tr -s ' ' ',' > mod_model1.pdb")
# 
# query <- "mod_model1.pdb" 
# 
# query_pdb <- read.csv(query,header=FALSE); 
# head(query_pdb)
# 
# names(query_pdb) <- c('recname','eleid','elename','resname','chainid','x1','x2','x3','occ','temp')
# 
# xyz <- query_pdb %>% filter(elename=="CA") %>% select(x1,x2,x3)


xyz[,1] <- as.numeric(as.character(xyz[,1]))
xyz[,2] <- as.numeric(as.character(xyz[,2]))
xyz[,3] <- as.numeric(as.character(xyz[,3]))

as.matrix(dist(xyz[,1], upper=TRUE,diag=TRUE)) %>% as_tibble(rownames = 'Var1') %>% gather(key=Var2,value,-Var1) 

mat1 <- as.matrix(dist(xyz[,1], upper=TRUE,diag=TRUE)) %>% as_tibble(rownames = 'Var1') %>% gather(key=Var2,value,-Var1) 
mat2 <- as.matrix(dist(xyz[,2], upper=TRUE,diag=TRUE)) %>% as_tibble(rownames = 'Var1') %>% gather(key=Var2,value,-Var1) 
mat3 <- as.matrix(dist(xyz[,3], upper=TRUE,diag=TRUE)) %>% as_tibble(rownames = 'Var1') %>% gather(key=Var2,value,-Var1) 

#bound <- cbind(mat1, mat2, mat3)
bound <- bind_cols(mat1, mat2, mat3) %>% mutate_all(as.numeric)
names(bound) <- c('x1','x2','xdist','y1','y2','ydist','z1','z2','zdist')

#To limit needless computation, consider only AAs within a n Å cubic space:

bound_filtered <- bound %>% 
  filter(xdist< ang_thresh & ydist< ang_thresh & zdist< ang_thresh) 

#Euclidean distance from xyz coords: 
#https://stackoverflow.com/questions/39671579/compute-euclidean-distance-matrix-from-x-y-z-coordinates

#index of AA pairs of interest:
ind <- bound_filtered %>% select(x1,x2)

Angstrom = apply(ind, 1, function(z){
  sqrt(sum((xyz[z[1],] - xyz[z[2], ])^2))
})
distances <- cbind(data.frame(ind), Angstrom)


#exclude data within 5 AAs from each other in the primary peptide sequence:
distances_filt <- distances %>% mutate(AAdist = x2-x1) %>% filter(abs(AAdist)>5) 

#mirror df and append to create square DF for plotting:
distancesMirror <- distances_filt %>% rename(x2t = x1, x1t=x2) %>% rename(x2=x2t, x1=x1t) %>% select(2,1,3)
#select fields 1:3 to conform to shape of distancesMirror
comb <- bind_rows(distances_filt %>% select(1:3), distancesMirror) 


#Add n to all x and y axes to represent the AA position from which Chain A was crystallized
comb_recal <- comb %>% mutate(x1=x1 + min(meta$resid), x2= x2 + min(meta$resid) ) %>% 
  left_join(meta,by=c('x1'='resid')) %>% rename(res_x1 = resname) %>% 
  left_join(meta,by=c('x2'='resid')) %>% rename(res_x2 = resname) 

comb_recal %>% filter(Angstrom < ang_thresh) %>% ggplot(aes(x=x1,y=x2)) + geom_point()

comb_recal %>% filter(Angstrom < ang_thresh) %>% select(1,2) %>% summary() #125 

dat_for_plot <- comb_recal %>% as_tibble() %>%  filter(Angstrom < ang_thresh) %>% 
  mutate(SNPx1 = ifelse(x1 %in% nonSyn,'y','n'),
         SNPx2 = ifelse(x2 %in% nonSyn,'y','n')) %>% 
  distinct()

myPlot <- dat_for_plot %>% 
  ggplot(aes(x=x1,y=x2)) + geom_point(col= 'dodger blue') +
  geom_abline() +
  geom_point(data = . %>% filter(SNPx1=='y'|SNPx2=='y' ), col='#F8766d',
             aes(text=paste0(x1,'. ',res_x1,' meets ',x2,'. ',res_x2))) 

tidyExt::default_GG_col(2)

ggplotly(myPlot, tooltip = 'text')

sessionInfo()
