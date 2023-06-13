### HETEROGENEOUS HYDROGENATION ANALYSIS ###

library(data.table)
library(mltools)
library(stringi)

data <- read.csv('../data/cleaned_datasets/hetero_hydrogenation.csv', fill = TRUE) # Your path here.
data$X <- NULL
data <- data.table(data)
data <- data[is.na(Reaction_T) == FALSE]
data <- data[is.na(Reaction_Time_hrs) == FALSE]
data <- data[, .(Product_Yield_PCT_Area_UV, PRODUCT_STRUCTURE, Solvent_1_Name, Reaction_T, Reactant_1_SMILES, reactant_2_SMILES, catalyst, ligand, Reagent_1_Short_Hand)]

data <- data[Reagent_1_Short_Hand == "", Reagent_1_Short_Hand := 'None']
data <- data[catalyst == "RaNi 4100", catalyst := "[Ni]"]
data <- data[catalyst == "RaNi 4200", catalyst := "[Ni]"]

# data <- data[Year == 2018]
# data <- data[Product_Yield_PCT_Area_UV > 0]

for (i in 1:nrow(data)) {
  halide <- toString(data[i, Reactant_1_SMILES])
  nuc <- toString(data[i, reactant_2_SMILES])
  catalyst <- toString(data[i, catalyst])
  lig <- toString(data[i, ligand])
  x <- paste(halide, nuc, sep='+')
  y <- paste(catalyst, lig, sep='+')
  data[i, pair := x]
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

un_pairs = unique(d[, pair])

for (p in un_pairs) {
  if (p == "c12c(ncnc1N1CCN(CC1)C(OC(C)(C)C)=O)C=C(CC2)c1c(ccc2c1cnn2[C@@H]1OCCCC1)C+NA" | p == "N1(C(=N[C@](C(S1(=O)=O)=C)(C)c1scc(n1)NC(c1ccc(cn1)OC(F)F)=O)NC(OC(C)(C)C)=O)C+NA" | p == "C=1(C(CCCC1NC(C)=O)=O)C+NA") {
    d[pair == p, rxn_type := 'alkene']
  }
  if (p == "c1c(cc2c(c1)C(=C([C@H](O2)c1ccc(cc1)OCc1ccccc1)c1ccc(cc1)OC)C)OC+NA" | p == "C1(CC(C1)OCc1ccccc1)N1CCOCC1+NA" | p == "N1(CC(CC1)(C(=O)OCC)C(=O)OCC)Cc1ccccc1+NA" | p == "n1c(nc2c(c1c1cnc(nc1)N)CCN2[C@]1(CCN(C1)C(=O)OCc1ccccc1)C)N1CCOCC1+NA" | p == "c1cc(c2c(c1OCc1ccccc1)c(n(cc2)C(OC(C)(C)C)=O)=O)F+NA" | p == "C1C2(CC1C2)N(Cc1ccccc1)Cc1ccccc1+NA" | p == "[nH]1c(nc2c(c1=O)C[C@H](CN2)CCc1sc(cc1C)C(=O)N[C@H](C(OCc1ccccc1)=O)CCc1nnnn1Cc1ccccc1)N+NA" | p == "[nH]1c(nc2c(c1=O)c(c[nH]2)CCc1ccc(cc1)C(=O)N[C@H](C(=O)OCc1ccccc1)CCc1[nH]nnn1)N+NA" | p == "n1c(nc2c(c1c1cnc(nc1)N)CCN2[C@]1(CCN(C1)Cc1ccccc1)C)N1CCOCC1+NA") {
    d[pair == p, rxn_type := 'deprotection']
  }
  if (p == "c1c(cncc1O)C+NA" | p == "c1(cc(ccc1)CO)CC(O)=O+NA" | p == "n1(nc(cc1)c1nccc(c1)C)C+NA" | p == "c1c(ccc(c1)C[C@H](C(OC)=O)NC(OC(C)(C)C)=O)O+NA" | p == "c12nc([nH]c2cccc1C(=O)O)C+NA") {
    d[pair == p, rxn_type := 'dearomatization']
  }
  if (p == "C1([C@@H]([C@@H]([C@H](O1)COC(=O)c1ccccc1)OC(c1ccccc1)=O)OC(c1ccccc1)=O)N=[N+]=[N-]+NA" | p == "c1(sc2c(c1)c(nn2C(C)OCC)Br)C(=O)OCC+NA" | p == "c1(cc(c2c(c1)c([nH]cc2)=O)[N+](=O)[O-])c1c(nnn1C)C+NA" | p == "c1cc2n(c(c1)C)n(c(c2C#N)=O)Cc1ccccc1+NA") {
    d[pair == p, rxn_type := 'other_FG']
  }
  if (p == "c12c(ncc(c1)C#Cc1sc(cc1)C(=O)OCC)nc([nH]c2=O)NC(C(C)(C)C)=O+NA") {
    d[pair == p, rxn_type := 'alkyne']
  }
  
}

for (prdt in unique(d[, PRODUCT_STRUCTURE])) {
  d <- zscore(d, prdt)
}

#### ANOVA ON INDIVIDUAL RXN CLASSES ####

# ALKENE
dalkene <- d[rxn_type == "alkene"]

dalkene.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=dalkene) # only cat-lig important
tukey <- TukeyHSD(dalkene.aov.factor)
tcat <- tukey$cat_lig
tcat <- as.data.frame(tcat)
tcat$factor <- rownames(tcat)
tcat_sig <- tcat[tcat[, "p adj"] < 0.05,]
tcat_sig <- data.table(tcat_sig)

t <- copy(tcat_sig)
for (i in unique(t[,factor])) {
  if (grepl("C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", i, fixed=T) ==T) {
    f1 <- "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]", i, fixed=T) == T) {
    f1 <- "C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+", i, fixed=T) == T) {
    f1 <- "c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+", i, fixed=T) == T) {
    f1 <- "c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O-[Cu]I+c1(C(=NC#N)N)ccccn1", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
    f2 <- "[Cu]I+c1(C(=NC#N)N)ccccn1"
  }
  
  else if (grepl("[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]-[Cu]I+c1(cc(ccn1)OC)C(N)=N", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]"
    f2 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
  }
  
  else if (grepl("[Cu]I+c1(cc(ccn1)OC)C(N)=N-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
  }
  
  else if (grepl("[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
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

hetero_alkene_top_15_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in un_factors_dt[, factor]) {
  zscores <- dalkene[cat_lig == i, z_score]
  avg <- mean(zscores)
  hetero_alkene_top_15_outlier_cats <- rbind(hetero_alkene_top_15_outlier_cats, list(i, avg))
}

write.table(hetero_alkene_top_15_outlier_cats, 'top_15_outlier_catalysts_alkene_hetero.txt')

# DEPROTECTION
ddeprot <- d[rxn_type == "deprotection"]


ddeprot.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=ddeprot) # both important
tukey <- TukeyHSD(ddeprot.aov.factor)
tcat <- tukey$cat_lig
tcat <- as.data.frame(tcat)
tcat$factor <- rownames(tcat)
tcat_sig <- tcat[tcat[, "p adj"] < 0.05,]
tcat_sig <- data.table(tcat_sig)

t <- copy(tcat_sig)
for (i in unique(t[,factor])) {
  if (grepl("C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", i, fixed=T) ==T) {
    f1 <- "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]", i, fixed=T) == T) {
    f1 <- "C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+", i, fixed=T) == T) {
    f1 <- "c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+", i, fixed=T) == T) {
    f1 <- "c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O-[Cu]I+c1(C(=NC#N)N)ccccn1", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
    f2 <- "[Cu]I+c1(C(=NC#N)N)ccccn1"
  }
  
  else if (grepl("[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]-[Cu]I+c1(cc(ccn1)OC)C(N)=N", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]"
    f2 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
  }
  
  else if (grepl("[Cu]I+c1(cc(ccn1)OC)C(N)=N-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
  }
  
  else if (grepl("[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
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

hetero_deprot_top_15_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in un_factors_dt[, factor][1:15]) {
  zscores <- ddeprot[cat_lig == i, z_score]
  avg <- mean(zscores)
  hetero_deprot_top_15_outlier_cats <- rbind(hetero_deprot_top_15_outlier_cats, list(i, avg))
}


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

hetereo_deprot_top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- ddeprot[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  hetereo_deprot_top_15_outlier_r <- rbind(hetereo_deprot_top_15_outlier_r, list(i, avg))
}

write.table(hetereo_deprot_top_15_outlier_r, 'top_15_outlier_catalysts_deprotection_hetero.txt')
write.table(hetero_deprot_top_15_outlier_cats, 'top_15_outlier_catalysts_deprotection_hetero.txt')

####### DEAROMATIZATION ########

ddearo <- d[rxn_type == "dearomatization"]


ddearo.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=ddearo) # both important
tukey <- TukeyHSD(ddearo.aov.factor)
tcat <- tukey$cat_lig
tcat <- as.data.frame(tcat)
tcat$factor <- rownames(tcat)
tcat_sig <- tcat[tcat[, "p adj"] < 0.05,]
tcat_sig <- data.table(tcat_sig)

t <- copy(tcat_sig)
for (i in unique(t[,factor])) {
  if (grepl("C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", i, fixed=T) ==T) {
    f1 <- "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]", i, fixed=T) == T) {
    f1 <- "C(COC)OC.[Ni](Cl)Cl+[N+]1(=CN(c2c(C(C)C)cccc2C(C)C)CC1)c1c(C(C)C)cccc1C(C)C.[Cl-]"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+", i, fixed=T) == T) {
    f1 <- "c1(c(c(ccc1OC)OC)P(C1CCCCC1)C1CCCCC1)c1c(cc(cc1C(C)C)C(C)C)C(C)C.c1(ccccc1c1ccccc1N)[Pd+].[O-]S(C)(=O)=O+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+", i, fixed=T) == T) {
    f1 <- "c1c(c(ccc1)[Pd+])c1ccccc1N.[O-]S(C)(=O)=O.c1c(c(c(cc1C(C)C)C(C)C)c1c(ccc(c1P(C(C)(C)C)C(C)(C)C)OC)OC)C(C)C+"
    f2 <- substr(i, nchar(f1) + 2, nchar(i))
  }
  
  else if (grepl("[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O-[Cu]I+c1(C(=NC#N)N)ccccn1", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
    f2 <- "[Cu]I+c1(C(=NC#N)N)ccccn1"
  }
  
  else if (grepl("[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]-[Cu]I+c1(cc(ccn1)OC)C(N)=N", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(ccc1cccn2)cccn3.[Cu+2].[N+](=O)([O-])[O-].[N+](=O)([O-])[O-]"
    f2 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
  }
  
  else if (grepl("[Cu]I+c1(cc(ccn1)OC)C(N)=N-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c1(cc(ccn1)OC)C(N)=N"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
  }
  
  else if (grepl("[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC-[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O", i, fixed=T) == T) {
    f1 <- "[Cu]I+c12c3c(c(ccn3)OC)ccc1c(ccn2)OC"
    f2 <- "[Cu]I+c1(c(ccc(c1)OC)OC)NC(Nc1c(cc(cc1)[N+](=O)[O-])OC)=O"
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

hetero_dearo_top_15_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in un_factors_dt[, factor][1:15]) {
  zscores <- ddearo[cat_lig == i, z_score]
  avg <- mean(zscores)
  hetero_dearo_top_15_outlier_cats <- rbind(hetero_dearo_top_15_outlier_cats, list(i, avg))
}

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

hetereo_dearo_top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- ddearo[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  hetereo_dearo_top_15_outlier_r <- rbind(hetereo_dearo_top_15_outlier_r, list(i, avg))
}

write.table(hetereo_dearo_top_15_outlier_r, 'top_15_outlier_reagents_deoaromatization_hetero.txt')
write.table(hetero_dearo_top_15_outlier_cats, 'top_15_outlier_catalysts_deoaromatization_hetero.txt')


