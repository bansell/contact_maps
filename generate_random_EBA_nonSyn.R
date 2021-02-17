
#slave of contactMaps~.R scripts

#Randomly generate substituted AAs for testing with EBA175:
AAs_20 <- readRDS('AminoAcids_abv.Rds')

#nonSyn <- sample(125,10,replace=FALSE) #vector of mutated/substituted AAs, by position

set.seed(1234)
nonSyn <- sample(152:745, 30, replace=FALSE) #vector of mutated/substituted AAs, by position
set.seed(1234)
sample(AAs_20,length(nonSyn), replace=TRUE)

subst_tbl <- tibble(RES_POS = nonSyn, RES_NAME = sample(AAs_20,length(nonSyn), replace=TRUE))

write_csv(subst_tbl,'EBA175RII_new_subst_table.csv')