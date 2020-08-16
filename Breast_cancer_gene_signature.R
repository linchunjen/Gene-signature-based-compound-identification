library(GEOquery)
library(ggplot2)
library(tidyr)
library(dplyr)
library(limma)
library(factoextra)
library(cluster)
library(ggfortify)
library(caret)
library(ROCR)
library(e1071)


#=======================================
# Download and Create Data Tables
#=======================================

# Download GSE data with GEO query library
gse20194 <- getGEO("GSE20194", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Create data tables
gse.pdata <- tbl_df(pData(gse20194$GSE20194_series_matrix.txt.gz))
gse.dat <- tbl_df(exprs(gse20194$GSE20194_series_matrix.txt.gz))
gse.ano <- tbl_df(fData(gse20194$GSE20194_series_matrix.txt.gz))
gse.pdata.m <- rename(gse.pdata, age = "age:ch1", bmngrd = "bmngrd:ch1", er_status = "er_status:ch1",
                      er = "er:ch1", her2.fish = "her2 fish:ch1", her2.ihc = "her2 ihc:ch1", 
                      her2_status = "her2 status:ch1", histology = "histology:ch1", 
                      nbefore = "nbefore:ch1", pcr_vs_rd = "pcr_vs_rd:ch1", pr_status = "pr_status:ch1", 
                      race = "race:ch1", tbefore = "tbefore:ch1", tissue = "tissue:ch1", 
                      treatment.code = "treatment code:ch1", treatments.comments = "treatments comments:ch1")

#Modify class types
gse.pdata.m.f <- mutate_if(gse.pdata.m, sapply(gse.pdata.m, is.character), as.factor)
gse.pdata.m.f.t <- gse.pdata.m.f
gse.pdata.m.f.t$age <- as.numeric(as.character(gse.pdata.m.f.t$age)) 
gse.pdata.m.f.t$er <- as.numeric(as.character(gse.pdata.m.f.t$er)) 
gse.pdata.m.f.t$her2.fish <- as.numeric(as.character(gse.pdata.m.f.t$her2.fish))


#=========================================
# Exploratory Data Analysis
#=========================================

## Age distribution vs races
ggplot(gse.pdata.m.f.t, aes(x = race, y = age)) + geom_boxplot() + labs(x = "Race", y = "AGE", 
                                                                        title = "Race vs AGE")
## Age distribution vs treatment response
ggplot(gse.pdata.m.f.t, aes(x = pcr_vs_rd, y = age)) + geom_boxplot() + labs(x = "pCR vs RD", y = "AGE", 
                                                                             title = "pCR vs RD vs AGE")

## Age distribution vs pathological characteristics
ggplot(gse.pdata.m.f.t, aes(x = bmngrd, y = age)) + geom_boxplot() + labs(x = "BMN grade", y = "AGE", 
                                                                          title = "BMN grade vs AGE")
ggplot(gse.pdata.m.f.t, aes(x = nbefore, y = age)) + geom_boxplot() + labs(x = "N-status before", y = "AGE", 
                                                                           title = "N stage before vs AGE")
ggplot(gse.pdata.m.f.t, aes(x = tbefore, y = age)) + geom_boxplot() + labs(x = "T stage", y = "AGE", 
                                                                           title = "T stage vs AGE")

## Age distribution vs molecular characteristics
ggplot(gse.pdata.m.f.t, aes(x = er_status, y = age)) + geom_boxplot() + labs(x = "ER status", y = "AGE", 
                                                                             title = "ER status vs AGE")
ggplot(gse.pdata.m.f.t, aes(x = pr_status, y = age)) + geom_boxplot() + labs(x = "PR status", y = "AGE", 
                                                                             title = "PR status vs AGE")
ggplot(gse.pdata.m.f.t, aes(x = her2_status, y = age)) + geom_boxplot() + labs(x = "HER2 status", y = "AGE", 
                                                                               title = "HER2 status vs AGE")

#================================================
# Array Data Normalization and Data Annotation
#================================================

# Array data normalization
plotDensities(gse.dat, main = "Before normalization") 
gse.dat.n <- tbl_df(normalizeBetweenArrays(gse.dat, method = "quantile"))
plotDensities(gse.dat.n, main = "After normalization")

# Data annotation
gse.dat.n.t <- t(gse.dat.n)
colnames(gse.dat.n.t) <- gse.ano$ID
gse.dat.n.t <- tbl_df(gse.dat.n.t)
gse.dat.n.pdata <- bind_cols(gse.pdata.m.f.t, gse.dat.n.t)

#================================================
# Clustering and PCA
#================================================

## K-mean cluster
fviz_nbclust(gse.dat.n.pdata[68:22350], kmeans, method = "silhouette") + 
  labs(subtitle = "Silhouette method") #find optimal number of cluster for K-mean
gdnp.km <- kmeans(gse.dat.n.pdata[68:22350], 2, nstart = 10)
fviz_cluster(gdnp.km, data = gse.dat.n.pdata[68:22350], geom = "point", ellipse.type = "norm", 
             palette = "jco", repel = TRUE, ggtheme = theme_minimal()) + labs(subtitle = "K-mean cluster")

# Clustering Large Applications
fviz_nbclust(gse.dat.n.pdata[68:22350], cluster::clara, method = "silhouette") + 
  labs(subtitle = "Silhouette method")
gdnp.clara <- clara(gse.dat.n.pdata[68:22350], 6)
fviz_cluster(gdnp.clara, data = gse.dat.n.pdata[68:22350], geom = "point", ellipse.type = "norm", 
             palette = "jco", repel = TRUE, ggtheme = theme_minimal()) + labs(subtitle = "Clustering Large Applications")

#PCA analysis
gdnp.pca <- prcomp(gse.dat.n.pdata[68:22350])
autoplot(gdnp.pca, data = gse.dat.n.pdata, colour = "description", 
         shape = "description", frame = TRUE, frame.type = 'norm')


#====================================================================
# Data visualization based on treatment and ER/PR/HER2 status
#====================================================================

#Subest data
gse.dat.n.pdata.train <- filter(gse.dat.n.pdata, description =="MAQC_Distribution_Status: MAQC_T -- Training")
gse.dat.n.pdata.test <- filter(gse.dat.n.pdata, description =="MAQC_Distribution_Status: MAQC_V -- Validation")
gdnp.trim <- bind_rows(gse.dat.n.pdata.train, gse.dat.n.pdata.test) #remove not used data

#Plot ER/PR/HER2 status base on the selected group
gdnp.trim.pca <- prcomp(gdnp.trim[68:22350])
autoplot(gdnp.trim.pca, data = gdnp.trim, colour = "pcr_vs_rd", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "pCR vs RD")
autoplot(gdnp.trim.pca, data = gdnp.trim, colour = "er_status", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "ER Status")
autoplot(gdnp.trim.pca, data = gdnp.trim, colour = "pr_status", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "PR Status")
autoplot(gdnp.trim.pca, data = gdnp.trim, colour = "her2_status", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "HER2 Status")

#===================================================================
# Machine learning classification (based on Whole transcriptome)
#===================================================================

## k-Nearest Neighbor Classification (based on ER status)
ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
gdnp.trim.knnFit <- train(as.data.frame(gse.dat.n.pdata.train[68:ncol(gse.dat.n.pdata.train)]), 
                          gse.dat.n.pdata.train$er_status, 
                          method = "knn", preProcess = c("center", "scale"), tuneLength = 10, 
                          trControl = ctrl)

gdnp.trim.knn.train <- predict(gdnp.trim.knnFit,
                               newdata = as.data.frame(gse.dat.n.pdata.train[68:ncol(gse.dat.n.pdata.train)]),
                               type = "prob")

gdnp.trim.knn.train.rocr <- prediction(gdnp.trim.knn.train[,2], gse.dat.n.pdata.train$er_status)
auc.1 <- performance(gdnp.trim.knn.train.rocr, "auc")@y.values
gdnp.trim.knn.train.rocr.tf <- performance(gdnp.trim.knn.train.rocr, "tpr", "fpr")
plot(gdnp.trim.knn.train.rocr.tf, main = "Training Data Set_K-NN")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.1),2)), bty = "n")

gdnp.trim.knn.pred <- predict(gdnp.trim.knnFit,
                              newdata = as.data.frame(gse.dat.n.pdata.test[68:ncol(gse.dat.n.pdata.test)]),
                              type = "prob")

gdnp.trim.knn.pred.rocr <- prediction(gdnp.trim.knn.pred[,2], gse.dat.n.pdata.test$er_status)
auc.2 <- performance(gdnp.trim.knn.pred.rocr, "auc")@y.values
gdnp.trim.knn.pred.rocr.tf <- performance(gdnp.trim.knn.pred.rocr, "tpr", "fpr")
plot(gdnp.trim.knn.pred.rocr.tf, main = "Validation Data Set_K-NN")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.2),2)), bty = "n")

## Support Vector Machine Classification (based on ER status)
gdnp.trim.svmFit <- svm(as.data.frame(gse.dat.n.pdata.train[68:ncol(gse.dat.n.pdata.train)]), 
                        gse.dat.n.pdata.train$er_status, probability=TRUE)

gdnp.trim.svm.train <- predict(gdnp.trim.svmFit, newdata = as.data.frame(gse.dat.n.pdata.train[68:ncol(gse.dat.n.pdata.train)]), probability=TRUE)


gdnp.trim.svm.train.rocr <- prediction(attr(gdnp.trim.svm.train, "prob")[,1], gse.dat.n.pdata.train$er_status)

auc.3 <- performance(gdnp.trim.svm.train.rocr, "auc")@y.values
gdnp.trim.svm.train.rocr.tf <- performance(gdnp.trim.svm.train.rocr, "tpr", "fpr")
plot(gdnp.trim.svm.train.rocr.tf, main = "Training Data Set_SVM")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.3),2)), bty = "n")

gdnp.trim.svm.pred <- predict(gdnp.trim.svmFit,
                              newdata = as.data.frame(gse.dat.n.pdata.test[68:ncol(gse.dat.n.pdata.test)]), 
                              probability=TRUE)

gdnp.trim.svm.pred.rocr <- prediction(attr(gdnp.trim.svm.pred, "prob")[,1], gse.dat.n.pdata.test$pcr_vs_rd)
auc.4 <- performance(gdnp.trim.svm.pred.rocr, "auc")@y.values
gdnp.trim.svm.pred.rocr.tf <- performance(gdnp.trim.svm.pred.rocr, "tpr", "fpr")
plot(gdnp.trim.svm.pred.rocr.tf, main = "Validation Data Set_SVM")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.4),2)), bty = "n")

                               
#==============================================================================
# Differential gene expression analysis (limma) for generate pCR signature
#==============================================================================

Group <- factor(gdnp.trim$pcr_vs_rd, levels=c("RD","pCR"))
design <- model.matrix(~Group)
colnames(design) <- c("WT","MUvsWT")
pCRvsRD.fit <- lmFit(as.data.frame(t(gdnp.trim[68:22350])), design=design)
pCRvsRD.fit.eb <- eBayes(pCRvsRD.fit)
pCRvsRD.de <- topTable(pCRvsRD.fit.eb, coef=2, p.value = 1e-04, number = Inf) #Select adj.p-value < 0.01 as pCR signature

pCRvsRD.sig <- as.vector(rownames(pCRvsRD.de))
gdnp.trim.pCRvsRD.sig <- select(gdnp.trim, pCRvsRD.sig)
gdnp.trim.pCRvsRD.sig.p <- bind_cols(gdnp.trim[1:67], gdnp.trim.pCRvsRD.sig)
gdnp.trim.pCRvsRD.sig.pca <- prcomp(gdnp.trim.pCRvsRD.sig.p[68:ncol(gdnp.trim.pCRvsRD.sig.p)])
autoplot(gdnp.trim.pCRvsRD.sig.pca, data = gdnp.trim.pCRvsRD.sig.p, colour = "pcr_vs_rd", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "pCR vs RD based on pCR signature")
autoplot(gdnp.trim.pCRvsRD.sig.pca, data = gdnp.trim.pCRvsRD.sig.p, colour = "er_status", 
         shape = "pcr_vs_rd", frame = TRUE, frame.type = 'norm') + labs(title = "ER status based on pCR signature")


#==================================================================================
# Examine whether gene signature can improve the prediction of treatment response. 
#==================================================================================

# Subset data
pCRvsRD.sig.p.train <- filter(gdnp.trim.pCRvsRD.sig.p, description =="MAQC_Distribution_Status: MAQC_T -- Training")
pCRvsRD.sig.p.test <- filter(gdnp.trim.pCRvsRD.sig.p, description =="MAQC_Distribution_Status: MAQC_V -- Validation")

#K-Nearest Neighbor Classification on pCR signature data (based on ER status)
ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
pCRvsRD.sig.knnFit <- train(as.data.frame(pCRvsRD.sig.p.train[68:ncol(pCRvsRD.sig.p.train)]), 
                            pCRvsRD.sig.p.train$er_status, 
                            method = "knn", preProcess = c("center", "scale"), tuneLength = 10, 
                            trControl = ctrl)
pCRvsRD.sig.knn.pred <- predict(pCRvsRD.sig.knnFit,
                                newdata = as.data.frame(pCRvsRD.sig.p.test[68:ncol(pCRvsRD.sig.p.test)]),
                                type = "prob")
pCRvsRD.sig.knn.pred.rocr <- prediction(pCRvsRD.sig.knn.pred[,2], pCRvsRD.sig.p.test$er_status)
auc.5 <- performance(pCRvsRD.sig.knn.pred.rocr, "auc")@y.values
pCRvsRD.sig.knn.pred.rocr.tf <- performance(pCRvsRD.sig.knn.pred.rocr, "tpr", "fpr")
plot(pCRvsRD.sig.knn.pred.rocr.tf, main = "Validation Data Set_K-NN_pCR signature")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.5),2)), bty = "n")


## Support Vector Machine Classification on pCR signature data  (based on ER status)
pCRvsRD.sig.svmFit <- svm(as.data.frame(pCRvsRD.sig.p.train[68:ncol(pCRvsRD.sig.p.train)]), 
                          pCRvsRD.sig.p.train$pcr_vs_rd, probability=TRUE)
pCRvsRD.sig.svm.pred <- predict(pCRvsRD.sig.svmFit,
                                newdata = as.data.frame(pCRvsRD.sig.p.test[68:ncol(pCRvsRD.sig.p.test)]), 
                                probability=TRUE)
pCRvsRD.sig.svm.pred.rocr <- prediction(attr(pCRvsRD.sig.svm.pred, "prob")[,1], pCRvsRD.sig.p.test$er_status)
auc.6 <- performance(pCRvsRD.sig.svm.pred.rocr, "auc")@y.values
pCRvsRD.sig.svm.pred.rocr.tf <- performance(pCRvsRD.sig.svm.pred.rocr, "tpr", "fpr")
plot(pCRvsRD.sig.svm.pred.rocr.tf, main = "Validation Data Set_SVM_pCR signature")
abline(a = c(0,0), b = c(1,1), lty = 2)
legend(0.6, 0.2, legend = paste("AUC =", round(unlist(auc.6),2)), bty = "n")

