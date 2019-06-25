#Load the data
data = read.table("3algos_2116seq.txt", sep = "\t", header = TRUE)
names(data)
# create value labels 
nrs = subset(data, category == "NRS") #Non-riboSNitches
rs = subset(data, category == "RS") # RiboSNitches

###ROC analysis
if(!require(pROC)) install.packages("pROC")
library(pROC)
#https://rpubs.com/Wangzf/pROC
##plotting ROCs

rocobj1 <- plot.roc(data$category, data$SNPfold3,
                    main="ROC curves on benchmark data",
                    col="#0000FF")
rocobj2 <- lines.roc(data$category, data$remuRNA, 
                     col="#FF0000")
rocobj3 <- lines.roc(data$category, data$RNAsnp1, 
                     col="#006400")

legend("top", legend=c("                                                      ",
                           paste("SNPfold  (AUC:",round(as.double(auc(rocobj1)),digits = 3),")"),
                           paste("remuRNA  (AUC:",round(as.double(auc(rocobj2)),digits=3),")"),
                           paste("RNAsnp  (AUC:",round(as.double(auc(rocobj3)),digits=3),")")
                           ),
col=c("#FFFFFF", "#0000FF", "#FF0000","#006400")
,lwd=2, box.col = "black")


