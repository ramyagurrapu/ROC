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

#for top 5%
five_percent_top_SNPfold3 = subset(data, data$SNPfold3 > quantile(data$SNPfold3, 0.95))
five_percent_bottom_SNPfold3 = subset(data, data$SNPfold3 < quantile(data$SNPfold3, 0.05))
five_percent_SNPfold3 = rbind(five_percent_top_SNPfold3, five_percent_bottom_SNPfold3)
str(five_percent_SNPfold3)

five_percent_top_remuRNA = subset(data, data$remuRNA > quantile(data$remuRNA, 0.95))
five_percent_bottom_remuRNA = subset(data, data$remuRNA < quantile(data$remuRNA, 0.05))
five_percent_remuRNA = rbind(five_percent_top_remuRNA, five_percent_bottom_remuRNA)
str(five_percent_remuRNA)

five_percent_top_RNAsnp1 = subset(data, data$RNAsnp1 > quantile(data$RNAsnp1, 0.95))
five_percent_bottom_RNAsnp1 = subset(data, data$RNAsnp1 < quantile(data$RNAsnp1, 0.05))
five_percent_RNAsnp1 = rbind(five_percent_top_RNAsnp1, five_percent_bottom_RNAsnp1)
str(five_percent_RNAsnp1)
#ROC for 5% tails
rocobj_1 <- plot.roc(five_percent_SNPfold3$category, five_percent_SNPfold3$SNPfold3,
                     main="ROC curves on benchmark data - 5% tails",
                     col="#0000FF")
rocobj_2 <- lines.roc(five_percent_remuRNA$category, five_percent_remuRNA$remuRNA, 
                      col="#FF0000")
rocobj_3 <- lines.roc(five_percent_RNAsnp1$category, five_percent_RNAsnp1$RNAsnp1, 
                      col="#006400")

legend("topleft", legend=c("                                                  ",
                           paste("SNPfold  AUC(",round(as.double(auc(rocobj_1)),digits = 3),")"),
                           paste("remuRNA  AUC(",round(as.double(auc(rocobj_2)),digits=3),")"),
                           paste("RNAsnp  AUC(",round(as.double(auc(rocobj_3)),digits=3),")")
),
col=c("#FFFFFF", "#0000FF", "#FF0000","#006400")
,lwd=2, box.col = "black")

###subcategory plots##

asr = data[grep("_ASR", data$old_name),]
str(asr)
sr = data[grep("_SR", data$old_name),]
vr = data[grep("_VR", data$old_name),]
pr = data[grep("_PR", data$old_name),]

for_asr = rbind(asr, nrs)
rocobj_asr_1 <- plot.roc(for_asr$category, for_asr$SNPfold3,
                         main="ROC curves for three algorithms - Asymmetric set",
                         col="#0000FF")
rocobj_asr_2 <- lines.roc(for_asr$category, for_asr$remuRNA,
                          col="#008000")
rocobj_asr_3 <- lines.roc(for_asr$category, for_asr$RNAsnp1,
                          col="#FF0000")
legend("topleft", legend=c("                                                  ",
                           paste("SNPfold  (AUC:",round(as.double(auc(rocobj_asr_1)),digits=3),")"),
                           paste("remuRNA  (AUC:",round(as.double(auc(rocobj_asr_2)),digits=3),")"),
                           paste("RNAsnp  (AUC:",round(as.double(auc(rocobj_asr_3)),digits=3),")")
),
col=c("#FFFFFF","#0000FF", "#008000","#FF0000")
,lwd=2, box.col = "black")

#####################
for_sr = rbind(sr, nrs)
rocobj_sr1 <- plot.roc(for_sr$category, for_sr$SNPfold3,
                       main="ROC curves for three algorithms - Symmetric set",
                       col="#0000FF")
rocobj_sr2 <- lines.roc(for_sr$category, for_sr$remuRNA,
                        col="#008000")
rocobj_sr3 <- lines.roc(for_sr$category, for_sr$RNAsnp1,
                        col="#FF0000")
legend("topleft", legend=c("                                                  ",
                           paste("SNPfold  (AUC:",round(as.double(auc(rocobj_sr1)),digits=3),")"),
                           paste("remuRNA  (AUC:",round(as.double(auc(rocobj_sr3)),digits=3),")"), 
                           paste("RNAsnp  (AUC:",round(as.double(auc(rocobj_sr2)),digits=3),")")
),
col=c("#FFFFFF","#0000FF", "#008000","#FF0000")
,lwd=2, box.col = "black")


#####################
for_vr = rbind(vr, nrs)
rocobj_vr1 <- plot.roc(for_vr$category, for_vr$SNPfold3,
                       main="ROC curves for three algorithms - Validated set",
                       col="#0000FF")
rocobj_vr2 <- lines.roc(for_vr$category, for_vr$remuRNA,
                        col="#008000")
rocobj_vr3 <- lines.roc(for_vr$category, for_vr$RNAsnp1,
                        col="#FF0000")
legend("topleft", legend=c("                                                  ",
                           paste("SNPfold  (AUC:",round(as.double(auc(rocobj_vr1)),digits=3),")"),
                           paste("remuRNA  (AUC:",round(as.double(auc(rocobj_vr2)),digits=3),")"),
                           paste("RNAsnp  (AUC:",round(as.double(auc(rocobj_vr3)),digits=3),")")
),
col=c("#FFFFFF", "#0000FF", "#008000","#FF0000")
,lwd=2, box.col = "black")
##############################
for_pr = rbind(pr, nrs)
rocobj_pr1 <- plot.roc(for_pr$category, for_pr$SNPfold3,
                       main="ROC curves for three algorithms - Probed set",
                       col="#0000FF")
rocobj_pr2 <- lines.roc(for_pr$category, for_pr$remuRNA,
                        col="#008000")
rocobj_pr3 <- lines.roc(for_pr$category, for_pr$RNAsnp1,
                        col="#FF0000")
legend("topleft", legend=c("                                                  ",
                           paste("SNPfold  (AUC:",round(as.double(auc(rocobj_pr1)),digits=3),")"),
                           paste("remuRNA  (AUC:",round(as.double(auc(rocobj_pr2)),digits=3),")"),
                           paste("RNAsnp  (AUC:",round(as.double(auc(rocobj_pr3)),digits=3),")")
),
col=c("#FFFFFF", "#0000FF", "#008000","#FF0000")
,lwd=2, box.col = "black")

