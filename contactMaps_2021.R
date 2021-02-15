# Contact maps

#Author: Brendan Ansell
#email: ansell.b@wehi.edu.au

library(tidyverse); 
#library(plotly) #required
#library(Rpdb)   #required


setwd("")

#User input:
ang_thresh <- 10   #Maximum 3D distance in Å between AA residues in the plot. Default = 10; Increasing this number makes a denser plot.
primary_dist <- 10 #Distance between residues in primary AA sequence. Default = 10 (essentially controls the density around x-y line)


# for testing:
source('generate_random_EBA_nonSyn.R') #for testing with EBA175
subst_tbl  <- read_csv('EBA175RII_new_subst_table.csv')  #two column csv file containing the position (RES_POS) and substituted AA residue (RES_NAME)

hit_chain <- 'A' #pdb chain to render. default should be A.

#EITHER get pdb files from local .pdb file e.g
local_file <- "EBA175RII_new.pdb"
pdb_import <- Rpdb::read.pdb(local_file) # if pdb file is directly supplied, read from .pdb file

#OR hit the PDB API (un-comment next 4 lines). NB this will fail if also using test file EBA175RII_new_subst_table.csv

#hit_pdb <- "2WLB"; #hit_pdb <- "6VOB";  #hit_pdb <- "4O75";
#Download pdb file from rcsb:
#system(paste0("curl -s https://files.rcsb.org/view/", hit_pdb, ".pdb > ",hit_pdb, ".pdb"))
#pdb_import <- Rpdb::read.pdb(paste0(hit_pdb,".pdb")); 


#TBD: detect whether local_file, or a PDB code is present in input arguments. run ifelse statement.


############ ############ ############ ############ ############ ############ ############ ############ ############
############ ############ ############ ############ ############ ############ ############ ############ ############


names(pdb_import$atom)

head(pdb_import$atoms)
summary(pdb_import$atoms$resid); min(pdb_import$atoms$resid)

pdb_import$atoms %>% filter(elename == "CA") %>% filter(chainid == hit_chain) %>% select(x1,x2,x3) %>% head()



meta <- pdb_import$atoms %>% filter(elename == "CA") %>% filter(chainid == hit_chain) %>% select(resid,resname) %>% distinct()
xyz  <- pdb_import$atoms %>% filter(elename=="CA") %>% select(x1,x2,x3)

min(meta$resid)

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
distances_filt <- distances %>% mutate(AAdist = x2-x1) %>% filter(abs(AAdist) > primary_dist) 

#mirror df and append to create square DF for plotting:
distancesMirror <- distances_filt %>% rename(x2t = x1, x1t=x2) %>% rename(x2=x2t, x1=x1t) %>% select(2,1,3)
#select fields 1:3 to conform to shape of distancesMirror
comb <- bind_rows(distances_filt %>% select(1:3), distancesMirror) 


#Add n to all x and y axes to represent the AA position from which Chain A was crystallized. Then annotate with AA IDs. 
comb_recal <- comb %>%
  mutate(x1=x1 + min(meta$resid), x2= x2 + min(meta$resid) ) %>% 
  left_join(meta,by=c('x1'='resid')) %>% rename(res_x1 = resname) %>% 
  left_join(meta,by=c('x2'='resid')) %>% rename(res_x2 = resname) 


if(any(ls() =='subst_tbl')){
  comb_recal <- comb_recal %>% 
    left_join(subst_tbl,by=c('x1'='RES_POS')) %>% rename(x1_RES_mut =RES_NAME) %>% 
    left_join(subst_tbl,by=c('x2'='RES_POS')) %>% rename(x2_RES_mut =RES_NAME ) %>% 
    rowwise() %>% 
    mutate(res_x1 = ifelse(!is.na(x1_RES_mut), paste0(res_x1,'>',x1_RES_mut),res_x1),
           res_x2 = ifelse(!is.na(x2_RES_mut), paste0(res_x2,'>',x2_RES_mut),res_x2))
                              
  
  myPlot <- comb_recal %>%  filter(Angstrom < ang_thresh) %>% 
    mutate(SNPx1 = ifelse(x1 %in% nonSyn,'y','n'),
           SNPx2 = ifelse(x2 %in% nonSyn,'y','n')) %>% 
    ggplot(aes(x=x1,y=x2)) + geom_point(col= 'dodger blue') +
    geom_abline() +
    geom_point(data = . %>% filter(SNPx1=='y'|SNPx2=='y' ), col='#F8766d',
               aes(text=paste0(x1,'. ',res_x1,' meets ',x2,'. ',res_x2))) 
  
}else{
  comb_recal <- comb_recal
  
  myPlot <- comb_recal %>%  filter(Angstrom < ang_thresh) %>% 
    mutate(SNPx1 = ifelse(x1 %in% nonSyn,'y','n'),
           SNPx2 = ifelse(x2 %in% nonSyn,'y','n')) %>% 
    ggplot(aes(x=x1,y=x2)) + geom_point(col= 'dodger blue') +
    geom_abline() +
    geom_point(data = . %>% filter(SNPx1=='y'|SNPx2=='y' ), col='#F8766d',
               aes(text=paste0(x1,'. ',res_x1,' meets ',x2,'. ',res_x2))) 
  
}
  

ggplotly(myPlot, tooltip = 'text')
