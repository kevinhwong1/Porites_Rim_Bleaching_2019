

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")

library(ropls)


#### PCA

data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
head(foods)
foodMN <- as.matrix(foods[, colnames(foods) != "Country"])
rownames(foodMN) <- foods[, "Country"]
head(foodMN)
foo.pca <- opls(foodMN)

#### PLS with a single response

data(cornell) ## see Tenenhaus, 1998
head(cornell)
cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                    cornell[, "y"])

## Complementary graphics

plot(cornell.pls, typeVc = c("outlier", "predict-train", "xy-score", "xy-weight"))

#### PLS with multiple (quantitative) responses

data(lowarp) ## see Eriksson et al. (2001); presence of NAs
head(lowarp)
lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                   as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                      grepl("^st", colnames(lowarp))]))

#### PLS-DA

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
plot(sacurine.plsda, parPaletteVc = c("green4", "magenta"))

#### OPLS-DA

sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)

#### Application to an ExpressionSet

sacSet <- Biobase::ExpressionSet(assayData = t(dataMatrix), 
                                 phenoData = new("AnnotatedDataFrame", 
                                                 data = sampleMetadata), 
                                 featureData = new("AnnotatedDataFrame", 
                                                   data = variableMetadata),
                                 experimentData = new("MIAME", 
                                                      title = "sacurine"))

sacPlsda <- opls(sacSet, "gender")
sacSet <- getEset(sacPlsda)
head(Biobase::pData(sacSet))
head(Biobase::fData(sacSet))

detach(sacurine)

______________

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")