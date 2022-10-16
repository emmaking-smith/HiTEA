### BUCHWALD ANALYSIS ###
library(data.table)
library(mltools)
library(stringi)
data <- read.csv('data/cleaned_datasets/buchwald.csv', fill = TRUE)
data$X <- NULL
data <- data.table(data)
data <- data[is.na(Reaction_T) == FALSE]
data <- data[is.na(Reaction_Time_hrs) == FALSE]
data <- data[, .(Product_Yield_PCT_Area_UV, PRODUCT_STRUCTURE, Solvent_1_Name, Reaction_T, halide, nuc, catalyst, ligand, Reagent_1_Short_Hand)]
data[Reagent_1_Short_Hand == "KOPnt", Reagent_1_Short_Hand := "KOtPn"]
data[Reagent_1_Short_Hand == "LiOBut", Reagent_1_Short_Hand := "LiOBut 1M in Hexanes"]

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

data[cat_lig == "c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+", cat_lig := "BrettPhos Palladacycle+"]
data[cat_lig == "C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+", cat_lig := "Xantphos Pd G4+"]

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

correlations <- function(d){
  d$PRODUCT_STRUCTURE <- NULL
  d$halide <- NULL
  d$nuc <- NULL
  d$catalyst <- NULL
  d$ligand <- NULL
  
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

correlations_of_d <- correlations(d)
# remove reagent_1_Short_Hand_Special in onehot.

un_pairs <- unique(d[, pair])
for (p in un_pairs) {
  if (p == "Brc1ccc2ncccc2c1+O=C1NCc2c(OCCCCN3CCNCC3)cccc21" | p == "Brc1ccccn1+Cn1cc(Nc2nc(O[C@H]3C[C@H](N(C)C(=O)OC(C)(C)C)C3)c3ccn(COCC[Si](C)(C)C)c3n2)cn1" | p == "CCc1ccn2ccnc2c1Br+CC(C)(C)OC(=O)N1CCCNCC1"){
    d[pair == p, rxn_type := 'br+sec_amine']
  }
  else if (p == "CN(C(=O)OC(C)(C)C)[C@H]1C[C@H](Oc2nc(Cl)nc3[nH]cc(-c4ccccn4)c23)C1+Cn1cc(N)cn1" | p == "Fc1c(Cl)nc(N2CCOCC2)nc1Cl+CC(C)(C)OC(=O)N1CC[C@](N)(CO)C1" | p == "Cc1cc(C)c(CN2CCc3ccc(OS(=O)(=O)C(F)(F)F)c(Cl)c3C2=O)c(OCc2ccccc2)n1+CCN.Cl" | p == "Cn1ccnc1Cl+Nc1cccnc1" | p == "Nc1ccc(-c2cc(Cl)nc(N3CCOCC3)n2)cn1+CN(C(=O)OC(C)(C)C)[C@H]1C[C@H](N)C1" | p == "CC(C)(C)OC(=O)N1C[C@H](COc2nc(Cl)nc3[nH]ccc23)[C@@H](C(F)(F)F)C1+Nc1cnoc1") {
    d[pair == p, rxn_type := 'cl+pri_amine']
  }
  else if (p == "Cn1cc(Nc2nc(Cl)nc3[nH]cnc23)cn1+CC(C)(C)OC(=O)N1C[C@@H]2CNC[C@@H]2C1" | p == "Cc1ncc(Cl)nc1Cl+C1COCCN1" | p == "Cc1cc(C)c(CN2CCc3ccc(OS(=O)(=O)C(F)(F)F)c(Cl)c3C2=O)c(OCc2ccccc2)n1+CCNC1CCOCC1"){
    d[pair == p, rxn_type := 'cl+sec_amine']
  }
  else if (p == "Brc1ccc(N2CCOCC2)nc1+CC(C)(C)OC(=O)N1CC[C@@H](N)C1" | p == "CC(C)(C)OC(=O)N1Cc2c(Br)nn(C(=O)OC(C)(C)C)c2C1(C)C+Cn1cc(N)cn1" | p == "Cn1ccnc1Br+Nc1cccnc1" | p == "CC(C)(C)OC(=O)N1CC[C@@H](N)C1"){
    d[pair == p, rxn_type := 'br+pri_amine']
  }
  else if (p == "Brc1cnsc1+CC(C)(C)OC(N)=O" | p == "Brc1ccc2ncccc2c1+CC(C)(C)OC(=O)NNC(=O)OC(C)(C)C"){
    d[pair == p, rxn_type := 'br+amide']
  }
  else if (p == "O=C1OCCn2c1ccc(Br)c2=O+Cc1c[nH]cn1" | p == "Brc1ccc2ncccc2c1+CCOC(=O)c1c[nH]nc1O"){
    d[pair == p, rxn_type := 'br+aro_N']
  }
  else if (p == "Cn1cc(Nc2nc(Cl)c3c(Cl)c[nH]c3n2)cn1+CC(C)(C)OC(=O)N1CC[C@@H](CO)C1"){
    d[pair == p, rxn_type := 'cl+pri_alc']
  }
  else if (p == "Cc1sc2nc(NCc3ccccc3)cnc2c1I+Cc1ncnc2[nH]ccc12"){
    d[pair == p, rxn_type := 'i+aro_N']
  }
  else if (p == "Cn1cc(I)cn1+CCOC(=O)n1nc(N)c2c1C(C)(C)N(C(=O)OC(C)(C)C)C2" | p == "Ic1cccnc1+Cn1ccnc1N"){
    d[pair == p, rxn_type := 'i+pri_amine']
  }
  else if (p == "Clc1ccc(I)cc1+Cc1nnc2n1-c1ccccc1NC(=O)C2"){
    d[pair == p, rxn_type := 'i+amide']
  }
  else if (p == "Cc1cc(Br)n(C)n1+CC(C)(C)OC(=O)CO" | p == "Nc1ncncc1Br+CC(C)(C)OC(=O)N1CCC(CO)CC1"){
    d[pair == p, rxn_type := 'br+pri_alc']
  }
  else if (p == "Nc1ncncc1I+CC(C)(C)OC(=O)N1CCC(CO)CC1" | p == "CC1=C(c2ccc(OC3CCCCO3)cc2)C(c2ccc(I)cc2)Oc2ccc(OC3CCCCO3)cc21+C[C@@H](CO)N1CC[C@@H](C)C1"){
    d[pair == p, rxn_type := 'i+pri_alc']
  }
  else if (p == "COc1ncccc1-c1cnc(N)c(O[C@H](C)c2cc(F)cnc2Cl)c1+c1c[nH]nn1"){
    d[pair == p, rxn_type := 'cl+aro_N']
  }
  else if (p == "CCOc1ccccc1I+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1"){
    d[pair == p, rxn_type := 'i+sec_alc']
  }
  else if (p == "CCOc1ccccc1Br+CC(C)(C)OC(=O)N1CCC[C@@H](O)C1"){
    d[pair == p, rxn_type := 'br+sec_alc']
  }
}

##### RANDOM FOREST FILES ######
random_forest_files <- function(d){
  # Remove superfulous rows.
  d$halide <- NULL
  d$nuc <- NULL
  d$catalyst <- NULL
  d$ligand <- NULL
  d$PRODUCT_STRUCTURE <- NULL
  d$Product_Yield_PCT_Area_UV <- NULL
  
  # One hot encoding & removal of duplicates.
  d <- one_hot(d)
  d$Reagent_1_Short_Hand_Special <- NULL
  d[, Reagent_1_Short_Hand_KOtPn := Reagent_1_Short_Hand_KOtPn + Reagent_1_Short_Hand_KOPnt]
  d$Reagent_1_Short_Hand_KOPnt <- NULL
  d$'cat_lig_c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+' <- d$'cat_lig_c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+' + d$'cat_lig_BrettPhos Palladacycle+'
  d$'cat_lig_BrettPhos Palladacycle+' <- NULL
  d$Reagent_1_Short_Hand_LiOBut <- d$Reagent_1_Short_Hand_LiOBut + d$'Reagent_1_Short_Hand_LiOBut 1M in Hexanes'
  d$'Reagent_1_Short_Hand_LiOBut 1M in Hexanes' <- NULL
  d$'cat_lig_C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+' <- d$'cat_lig_C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+' + d$'cat_lig_Xantphos Pd G4+'
  d$'cat_lig_Xantphos Pd G4+' <- NULL
  
  # Writing to text files
  for (i in unique(d[, rxn_type])) {
    d_i <- d[rxn_type == i]
    print(i)
    print(nrow(d_i))
    d_i <- one_hot(d_i)
    d_i <- d_i[, rxn_type := NULL]
    d_i <- d_i[, .SD, .SDcols = apply(d_i, 2, var) != 0] # No constant value columns.
    n <- paste('one_hot_rxn_type', i, sep='_')
    n <- paste(n, 'txt', sep='.')
    write.table(d_i, n)
  }
  
  # For any overall dataset comparisons, save out entire dataset.
  write.table(d, 'one_hot_RF_zscores.txt')
  
  return(0)
}

rf_files <- random_forest_files(d)



### ANOVA WITH TUKEY ###

d$rxn_type <- NULL
d$Product_Yield_PCT_Area_UV <- NULL
setnames(d, 'cat_lig_c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+', 'BrettPhos_Pd_G1')

for (i in colnames(d)) {
  d$i <- as.factor(d$i)
}

# Removing identical rows.
d <- d[Reagent_1_Short_Hand != 'Special']
d <- d[Reagent_1_Short_Hand == 'KOPnt', Reagent_1_Short_Hand := 'KOtPn']
d <- d[cat_lig == 'c1(cc(cc(c1c1c(ccc(c1P(C1CCCCC1)C1CCCCC1)OC)OC)C(C)C)C(C)C)C(C)C.c1(c(cccc1)[Pd]Cl)CCN+', cat_lig := 'BrettPhos Palladacycle+']
d <- d[Reagent_1_Short_Hand == 'LiOBut 1M in Hexanes', Reagent_1_Short_Hand := 'LiOBut']
d <- d[cat_lig == 'C1(C)(C)c2c(c(P(c3ccccc3)c3ccccc3)ccc2)Oc2c1cccc2P(c1ccccc1)c1ccccc1.c1ccc(c(c1)[Pd+])c1ccccc1NC.[O-]S(C)(=O)=O+', cat_lig := 'Xantphos Pd G4+']

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

top_15_outlier_cats <- data.table(factor = as.character(), avg_zscore = as.numeric())
for (i in un_factors_dt[, factor][1:15]) {
  zscores <- d[cat_lig == i, z_score]
  avg <- mean(zscores)
  top_15_outlier_cats <- rbind(top_15_outlier_cats, list(i, avg))
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

top_15_outlier_r <- data.table(factor = as.character(), avg_zscore = as.numeric())

for (i in r_un_factors_dt[, factor][1:15]) {
  zscores <- d[Reagent_1_Short_Hand == i, z_score]
  avg <- mean(zscores)
  top_15_outlier_r <- rbind(top_15_outlier_r, list(i, avg))
}


#### THE TOP OUTLIER CATALYSTS AND REGENTS DATATABLES #####
write.table(top_15_outlier_r, 'top_15_outlier_reagents_buchwald.txt')
write.table(top_15_outlier_cats, 'top_15_outlier_catalysts_buchwald.txt')