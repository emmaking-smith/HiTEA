### ULLMANN ANOVA-TUKEY ANALYSIS ###

library(data.table)
library(mltools)
library(stringi)

data <- read.csv('../data/cleaned_datasets/ullmann.csv', fill = TRUE)
data$X <- NULL
data <- data.table(data)
data <- data[is.na(Reaction_T) == FALSE]
data <- data[is.na(Reaction_Time_hrs) == FALSE]
data <- data[, .(Product_Yield_PCT_Area_UV, PRODUCT_STRUCTURE, Solvent_1_Name, Reaction_T, halide, nuc, catalyst, ligand, Reagent_1_Short_Hand)]

# data <- data[Year == 2016]
# data <- data[Product_Yield_PCT_Area_UV > 0]

for (i in 1:nrow(data)) {
  halide <- toString(data[i, halide])
  nuc <- toString(data[i, nuc])
  catalyst <- toString(data[i, catalyst])
  lig <- toString(data[i, ligand])
  x <- paste(halide, nuc, sep='+')
  y <- paste(catalyst, lig, sep='+')
  data[i, pair := x]
  data[i, cat_lig := y]
}

data$cat_lig <- as.factor(data$cat_lig)
data$pair <- as.factor(data$pair)

# Calculating the z scores.

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

d <- copy(data)

un_pairs <- unique(d[, pair])
for (p in un_pairs) {
  if (p == "Clc1ccc(I)cc1+CCOC(=O)c1c[nH]nc1O" | p == "CCOc1cc(Cl)c(C#N)cc1I+c1cn[nH]c1" | p == "Cc1sc2nc(NCc3ccccc3)cnc2c1I+Clc1ncnc2[nH]ccc12"){
    d[pair == p, rxn_type := 'i+aro_N']
  }
  else if (p == "Nc1ncncc1I+CC(C)(C)OC(=O)N1CCC(CO)CC1" | p == "CC1=C(c2ccc(OC3CCCCO3)cc2)C(c2ccc(I)cc2)Oc2ccc(OC3CCCCO3)cc21+C[C@@H](CO)N1CC[C@@H](C)C1") {
    d[pair == p, rxn_type := 'i+pri_alc']
  }
  else if (p == "COc1ccc(C2=C(C)c3ccc(OC)cc3O[C@@H]2c2ccc(I)cc2)cc1+CCCN1CC(O)C1" | p == "CCOc1ccccc1I+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1") {
    d[pair == p, rxn_type := 'i+sec_alc']
  }
  else if (p == "CCOC(=O)C(F)Br+OB(O)c1ccc(-c2ccccc2)cc1" | p == "CCOC(=O)C(F)Br+OB(O)c1ccc(-c2cccnc2)cc1") {
    d[pair == p, rxn_type := 'br+boronic_acid']
  }
  else if (p == "Cc1cc(Br)nn1C+C=CCO" | p == "Nc1ncncc1Br+CC(C)(C)OC(=O)N1CCC(CO)CC1") {
    d[pair == p, rxn_type := 'br+pri_alc']
  }
  
}

for (prdt in unique(d[, PRODUCT_STRUCTURE])) {
  d <- zscore(d, prdt)
}

### ANOVA WITH TUKEY ###

d.aov.factor <- aov(z_score ~ cat_lig + Reagent_1_Short_Hand, data=d) 
tukey <- TukeyHSD(d.aov.factor)
tcat <- tukey$cat_lig
tcat <- as.data.frame(tcat)
tcat$factor <- rownames(tcat)
tcat_sig <- tcat[tcat[, "p adj"] < 0.05,]
tcat_sig <- data.table(tcat_sig)
tcat <- data.table(tcat)

# making factor 1 and factor 2

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

ull_top_15_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in un_factors_dt[, factor][1:15]) {
  zscores <- d[cat_lig == i, z_score]
  avg <- mean(zscores)
  ull_top_15_outlier_cats <- rbind(ull_top_15_outlier_cats, list(i, avg))
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

ull_top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- d[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  ull_top_15_outlier_r <- rbind(ull_top_15_outlier_r, list(i, avg))
}

write.table(ull_top_15_outlier_r, 'top_15_outlier_reagents_ullmann.txt')
write.table(ull_top_15_outlier_cats, 'top_15_outlier_catalysts_ullmann.txt')


