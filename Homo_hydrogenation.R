### HOMOGENEOUS PAIR CORRELATIONS ###

library(data.table)
library(mltools)
library(stringi)

data <- read.csv('data/cleaned_datasets/homo_hydrogenation.csv', fill = TRUE)
data$X <- NULL
data <- data.table(data)
data <- data[is.na(Reaction_T) == FALSE]
data <- data[is.na(Reaction_Time_hrs) == FALSE]
data <- data[, .(Product_Yield_PCT_Area_UV, PRODUCT_STRUCTURE, Solvent_1_Name, Reaction_T, Reactant_1_SMILES, reactant_2_SMILES, catalyst, ligand, Reagent_1_Short_Hand)]
data <- data[Reagent_1_Short_Hand == "", Reagent_1_Short_Hand := 'None']
data <- data[ligand == "R-JosiphosSL-J009-1", ligand := "R-Josiphos SL-J009-1"]

for (i in 1:nrow(data)) {
  halide <- toString(data[i, Reactant_1_SMILES])
  catalyst <- toString(data[i, catalyst])
  lig <- toString(data[i, ligand])
  y <- paste(catalyst, lig, sep='+')
  data[i, pair := halide]
  data[i, cat_lig := y]
}

data$cat_lig <- as.factor(data$cat_lig)
data$pair <- as.factor(data$pair)

d <- copy(data)


zscore <- function(dt, product){
  rows_remove <- list()
  p_dt <- dt[PRODUCT_STRUCTURE == product]
  if (nrow(p_dt) > 1){
    mean <- mean(p_dt[, Product_Yield_PCT_Area_UV])
    std <- sd(p_dt[, Product_Yield_PCT_Area_UV])
    dt[PRODUCT_STRUCTURE == product, z_score := ((Product_Yield_PCT_Area_UV - mean)/(std + 0.01))]
  }
  else {
    rows_remove <- append(rows_remove, product)
  }
  dt <- dt[PRODUCT_STRUCTURE %in% setdiff(dt[,PRODUCT_STRUCTURE], rows_remove)]
  
  
  return (dt)
  
}

un_pairs <- unique(d[, pair])

for (p in un_pairs) {
  if (p == "c12c(ncnc1N1CCN(CC1)C(OC(C)(C)C)=O)C=C(CC2)c1c(ccc2c1cnn2[C@@H]1OCCCC1)C" | p == "C1C(C(=C1c1ccccc1)C)=O" | p == "c1cccc(c1)C(\\C(F)(F)F)=C/C(O)=O" | p == "c1nccnc1C(=NCc1ccc(cc1)OC)C" | p == "n1(nc(c(c1)C=1CCCN1)C)C" | p == "C(C(F)(F)F)(C(O)=O)=C" | p == "c1(cnc2c(c1C)c(cn2C)C1=CCN(C(C1)(C)C)C(OC(C)(C)C)=O)[N+]([O-])=O" | p == "c12c(C(NCC1=C)=O)cccc2" | p == "C1C=C(CC=C1O[Si](C)(C)C(C)(C)C)C" | p == "N1(C(=N[C@](C(S1(=O)=O)=C)(C)c1scc(n1)NC(c1ccc(cn1)OC(F)F)=O)NC(OC(C)(C)C)=O)C" | p == "c1(cc(c(nn1)N)C(c1ccccc1)=C)Cl") {
    d[pair == p, rxn_type := 'alkene']
  }
  if (p == "c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)OC)=O)C)C(=O)[C@@H]1COCC1" | p == "c1cc(c(c(c1I)C(C)=O)Cl)F" | p == "O1[C@H]([C@H]2[C@@H]([C@H]1C(=O)c1cc(c(cc1)F)F)OC(O2)(C)C)n1cc(c2c1ncnc2Cl)F" | p == "c1(c(ccc(c1)F)OC)C(C)=O" | p == "c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)C)=O)Cl)C([C@@H]1CCOC1)=O" | p == "c1(cc(c2c(c1Cl)C(N(CC2)Cc1c(cc(nc1OCc1ccccc1)C)C)=O)Cl)C(=O)C1CCN(CC1)C(=O)OC(C)(C)C") {
    d[pair == p, rxn_type := 'CO_reduction']
  
  }
}

for (prdt in unique(d[, PRODUCT_STRUCTURE])) {
  d <- zscore(d, prdt)
}

random_forest_files <- function(d){
  d$Reactant_1_SMILES <- NULL
  d$reactant_2_SMILES <- NULL
  d$catalyst <- NULL
  d$ligand <- NULL
  d$rxn_type <- NULL
  d$PRODUCT_STRUCTURE <- NULL
  
  d <- one_hot(d)
  
  for (i in unique(d[, rxn_type])) {
    d_i <- d[rxn_type == i]
    d_i <- one_hot(d_i)
    d_i <- d_i[, rxn_type := NULL]
    d_i <- d_i[, .SD, .SDcols = apply(d_i, 2, var) != 0] # No constant value columns.
    n <- paste('one_hot_rxn_type_homo_hydrog', i, sep='_')
    n <- paste(n, 'txt', sep='.')
    write.table(d_i, n)
  }
  return(0)
}


rf_files <- random_forest_files(d)

# correlations: - Corrleations between solvents and certain pairs/reagents, but this is stripped out in the rxn class saving.
correlations <- function(d){
  d$PRODUCT_STRUCTURE <- NULL
  d$Reactant_1_SMILES <- NULL
  d$reactant_2_SMILES <- NULL
  d$catalyst <- NULL
  d$ligand <- NULL
  d$rxn_type <- NULL
  
  d <- one_hot(d)
  
  d <- d[, .SD, .SDcols = colSums(d) > 0] # No all zero columns.
  coln <- colnames(d)
  coln <- coln[coln != 'Product_Yield_PCT_Area_UV']
  
  rows_to_remove <- list()
  for (i in coln){
    for (j in coln){
      if (i != j){
        c <- cor(d[, c(..i)], d[, c(..j)])
        if (abs(c) > 0.85) {
          print(paste(i, j, 'have correlation of', c))
          rows_to_remove <- append(rows_to_remove, i)
          rows_to_remove <- append(rows_to_remove, j)
        }
      }
    }
  }
  return(rows_to_remove)
}


###### ALKENE ######
dalkene <- d[rxn_type == "alkene"]

dalkene.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=dalkene) # only reagents important
tukey <- TukeyHSD(dalkene.aov.factor)
tr <- tukey$Reagent_1_Short_Hand
tr <- as.data.frame(tr)
tr$factor <- rownames(tr)
tr_sig <- tr[tr[, "p adj"] < 0.05,]
tr_sig <- data.table(tr_sig)

for (i in unique(tr_sig[,factor])) {
  
  if (grepl("P1-t-Bu-tris", i, fixed=T) ==T) {
    f1 <- "P1-t-Bu-tris"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P1-t-Bu", i, fixed=T) == T) {
    f1 <- "P1-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P2-Et", i, fixed=T) == T) {
    f1 <- "P2-Et"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P2-t-Bu", i, fixed=T) == T) {
    f1 <- "P2-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P4-t-Bu", i, fixed=T) == T) {
    f1 <- "P2-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else {
    f1 <- strsplit(i, "\\-")[[1]][1]
    f2 <- strsplit(i, "\\-")[[1]][2]
  }
  
  tr_sig[factor == i, factor_1 := f1]
  tr_sig[factor == i, factor_2 := f2]
  
}

r_un_factors <- unique(union(tr_sig[, factor_1], tr_sig[ ,factor_2]))

r_un_factors_dt <- data.table(factor = as.character(), count = as.numeric())

for (i in r_un_factors) {
  n <- nrow(tr_sig[factor_1 == i]) + nrow(tr_sig[factor_2 == i])
  r_un_factors_dt <- rbind(r_un_factors_dt, list(i, n))
}

r_un_factors_dt <- r_un_factors_dt[order(count, decreasing=T)]

homo_alkene_top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- dalkene[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  homo_alkene_top_15_outlier_r <- rbind(homo_alkene_top_15_outlier_r, list(i, avg))
  
}

write.table(homo_alkene_top_15_outlier_r, 'top_15_outlier_reagents_alkene_homo.txt')

###### CO REDUCTIONS ######
dco <- d[rxn_type == "CO_reduction"]

dco.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=dco) # both important
tukey <- TukeyHSD(dco.aov.factor)
tcat <- tukey$cat_lig
tcat <- as.data.frame(tcat)
tcat$factor <- rownames(tcat)
tcat_sig <- tcat[tcat[, "p adj"] < 0.05,]
tcat_sig <- data.table(tcat_sig)
tcat <- data.table(tcat)

# making factor 1 and factor 2

t <- copy(tcat_sig)
for (i in unique(t[,factor])) {
  
  
  if (grepl("C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.[Fe].c1ccc(cc1)P(c1ccccc1)c1ccccc1C1CCCC1[C@H](C)P(C1[C@H]2CC[C@@H](C1)C2)C1[C@@H]2CC[C@H](C1)C2-", i, fixed=T) ==T) {
    f1 <- "C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.[Fe].c1ccc(cc1)P(c1ccccc1)c1ccccc1C1CCCC1[C@H](C)P(C1[C@H]2CC[C@@H](C1)C2)C1[C@@H]2CC[C@H](C1)C2"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[B-](F)(F)(F)F.[Rh+].C1=CCCC=CCC1.C1=CCCC=CCC1+P(C1C([C@H](P(C(C)(C)C)C(C)(C)C)C)CCC1)(c1ccccc1)c1ccccc1.C1CCCC1.[Fe]-", i, fixed=T) ==T) {
    f1 <- "[B-](F)(F)(F)F.[Rh+].C1=CCCC=CCC1.C1=CCCC=CCC1+P(C1C([C@H](P(C(C)(C)C)C(C)(C)C)C)CCC1)(c1ccccc1)c1ccccc1.C1CCCC1.[Fe]"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C1CC=CCCC=C1.C1CC=CCCC=C1.Cl[Ir].Cl[Ir]+c1(c2c(P(c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)ccc3c2OCO3)c(P(c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)ccc2c1OCO2-", i, fixed=T) ==T) {
    f1 <- "C1CC=CCCC=C1.C1CC=CCCC=C1.Cl[Ir].Cl[Ir]+c1(c2c(P(c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)ccc3c2OCO3)c(P(c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)ccc2c1OCO2"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R)-(+)-Cl-MeO-BIPHEP-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R)-(+)-Cl-MeO-BIPHEP"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R)-(+)-MeO-BIPHEP, SL-A101-1-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R)-(+)-MeO-BIPHEP, SL-A101-1"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R)-(+)-XylBINAP-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R)-(+)-XylBINAP"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R)-DM-SEGPHOS-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R)-DM-SEGPHOS"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+R-Josiphos SL-J009-1-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+R-Josiphos SL-J009-1"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R)-Xylyl-P-Phos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R)-Xylyl-P-Phos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R,R)-DIPAMP-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R,R)-DIPAMP"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R,R)-Me-BPE-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R,R)-Me-BPE"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(R,R)-Me-DuPhos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(R,R)-Me-DuPhos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(Ra,S)-DTB-Bn-SIPHOX-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(Ra,S)-DTB-Bn-SIPHOX"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(Ra,S)-Ph-Bn-SIPHOX-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(Ra,S)-Ph-Bn-SIPHOX"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(S)-Phanephos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(S)-Phanephos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(S)-Xylyl-P-Phos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(S)-Xylyl-P-Phos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+(S,S,R,R)-TangPhos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+(S,S,R,R)-TangPhos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+[Rh(cod)2]BF4 /(S)-Phanephos-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+[Rh(cod)2]BF4 /(S)-Phanephos"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+CTH-(R)-BINAM-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+CTH-(R)-BINAM"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Rh(cod)2]BF4+Mandyphos SL-M002-1-", i, fixed=T) ==T) {
    f1 <- "[Rh(cod)2]BF4+Mandyphos SL-M002-1"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("-[B-](F)(F)(F)F.[Rh+].C1=CCCC=CCC1.C1=CCCC=CCC1+P(C1C([C@H](P(C(C)(C)C)C(C)(C)C)C)CCC1)(c1ccccc1)c1ccccc1.C1CCCC1.[Fe]", i, fixed=T) ==T) {
    f2 <- "[B-](F)(F)(F)F.[Rh+].C1=CCCC=CCC1.C1=CCCC=CCC1+P(C1C([C@H](P(C(C)(C)C)C(C)(C)C)C)CCC1)(c1ccccc1)c1ccccc1.C1CCCC1.[Fe]"
    f1 <- substr(i, 1, nchar(i) - nchar(f2) - 1)
  }
  
  else if (grepl("C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.CC(P(C(C)(C)C)[C@H](C1CCCC1P(=O)c1ccccc1)C)(C)C.[Fe]-", i, fixed=T) ==T) {
    f1 <- "C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.CC(P(C(C)(C)C)[C@H](C1CCCC1P(=O)c1ccccc1)C)(C)C.[Fe]"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("-C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.[Fe].c1ccc(cc1)P(c1ccccc1)c1ccccc1C1CCCC1[C@H](C)P(C1[C@H]2CC[C@@H](C1)C2)C1[C@@H]2CC[C@H](C1)C2", i, fixed=T) ==T) {
    f2 <- "C(S(=O)(=O)[O-])(F)(F)F.[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1+C1CCCC1.[Fe].c1ccc(cc1)P(c1ccccc1)c1ccccc1C1CCCC1[C@H](C)P(C1[C@H]2CC[C@@H](C1)C2)C1[C@@H]2CC[C@H](C1)C2"
    f1 <- substr(i, 1, nchar(i) - nchar(f2) - 1)
  }
  
  else if (grepl("C1-358+-", i, fixed=T) ==T) {
    f1 <- "C1-358+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C1-360+-", i, fixed=T) ==T) {
    f1 <- "C1-360+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("-C1CC=CCCC=C1.C1CC=CCCC=C1.Cl[Ir].Cl[Ir]+c1(c2c(P(c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)ccc3c2OCO3)c(P(c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)ccc2c1OCO2", i, fixed=T) ==T) {
    f2 <- "C1CC=CCCC=C1.C1CC=CCCC=C1.Cl[Ir].Cl[Ir]+c1(c2c(P(c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)c3cc(c(c(c3)C(C)(C)C)OC)C(C)(C)C)ccc3c2OCO3)c(P(c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)c2cc(c(c(c2)C(C)(C)C)OC)C(C)(C)C)ccc2c1OCO2"
    f1 <- substr(i, 1, nchar(i) - nchar(f2) - 1)
  }
  
  else if (grepl("-[Rh(cod)2]BF4+R-Josiphos SL-J009-1", i, fixed=T) ==T) {
    f2 <- "[Rh(cod)2]BF4+R-Josiphos SL-J009-1"
    f1 <- substr(i, 1, nchar(i) - nchar(f2) - 1)
  }
  
  else if (grepl("Ru-721+-", i, fixed=T) ==T) {
    f1 <- "Ru-721+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C1(C(P(c2ccccc2)c2ccccc2)CCC1)C1=N[C@@H](CO1)C(C)C.P(c1ccccc1)(c1ccccc1)c1ccccc1.[Ru](Cl)Cl.C1CCCC1.[Fe]+-", i, fixed=T) ==T) {
    f1 <- "C1(C(P(c2ccccc2)c2ccccc2)CCC1)C1=N[C@@H](CO1)C(C)C.P(c1ccccc1)(c1ccccc1)c1ccccc1.[Ru](Cl)Cl.C1CCCC1.[Fe]+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else {
    f1 <- strsplit(i, "\\-")[[1]][1]
    f2 <- strsplit(i, "\\-")[[1]][2]
  }
  
  t[factor == i, factor_1 := f1]
  t[factor == i, factor_2 := f2]
  
}

un_factors <- unique(union(t[, factor_1], t[ ,factor_2]))

un_factors_dt <- data.table(factor = as.character(), count = as.numeric())

for (i in un_factors) {
  n <- nrow(t[factor_1 == i]) + nrow(t[factor_2 == i])
  un_factors_dt <- rbind(un_factors_dt, list(i, n))
}

un_factors_dt <- un_factors_dt[order(count, decreasing=T)]

homo_co_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in un_factors_dt[, factor][1:15]) {
  zscores <- dco[cat_lig == i, z_score]
  avg <- mean(zscores)
  homo_co_outlier_cats <- rbind(homo_co_outlier_cats, list(i, avg))
}

tr <- tukey$Reagent_1_Short_Hand
tr <- as.data.frame(tr)
tr$factor <- rownames(tr)
tr_sig <- tr[tr[, "p adj"] < 0.05,]
tr_sig <- data.table(tr_sig)
#### NO SIGNIFICANT REAGENTS ####
for (i in unique(tr_sig[,factor])) {
  
  if (grepl("P1-t-Bu-tris", i, fixed=T) ==T) {
    f1 <- "P1-t-Bu-tris"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P1-t-Bu", i, fixed=T) == T) {
    f1 <- "P1-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P2-Et", i, fixed=T) == T) {
    f1 <- "P2-Et"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P2-t-Bu", i, fixed=T) == T) {
    f1 <- "P2-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("P4-t-Bu", i, fixed=T) == T) {
    f1 <- "P2-t-Bu"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else {
    f1 <- strsplit(i, "\\-")[[1]][1]
    f2 <- strsplit(i, "\\-")[[1]][2]
  }
  
  tr_sig[factor == i, factor_1 := f1]
  tr_sig[factor == i, factor_2 := f2]
  
}

r_un_factors <- unique(union(tr_sig[, factor_1], tr_sig[ ,factor_2]))

r_un_factors_dt <- data.table(factor = as.character(), count = as.numeric())

for (i in r_un_factors) {
  n <- nrow(tr_sig[factor_1 == i]) + nrow(tr_sig[factor_2 == i])
  r_un_factors_dt <- rbind(r_un_factors_dt, list(i, n))
}

r_un_factors_dt <- r_un_factors_dt[order(count, decreasing=T)]

homo_alkene_top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- dalkene[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  homo_alkene_top_15_outlier_r <- rbind(homo_alkene_top_15_outlier_r, list(i, avg))
  
}

write.table(homo_alkene_top_15_outlier_r, 'top_15_outlier_reagents_CO_homo.txt')
write.table(homo_co_outlier_cats, 'top_15_outlier_reagents_CO_homo.txt')




