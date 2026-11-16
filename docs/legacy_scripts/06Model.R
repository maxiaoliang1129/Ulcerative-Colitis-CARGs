

rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\06Model\\01model'
setwd(path)

#data
rt<-read.table('diffGeneExp_p0.05.txt',header = T,sep='\t',row.names=1)
hub<-read.table('venn.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Type')

rowsid<-rownames(rt)
rowsid<-toupper(rowsid)
rownames(rt)<-rowsid

rt<-rt[hub$V1,]
rt<-t(rt)
rt<-data.frame(rt)

rt$Sample<-rownames(rt)

rt<-merge(group,rt)
rownames(rt)<-rt[,1]
rt<-rt[,-1]
colnames(rt)[1]<-'Group'
table(rt$Group)
#control      UC 
#455     865 

#训练最终模型
rt1<-rt
dim(rt1)

#处理分组
rt1[which(rt1$Group == 'UC'),]$Group<-1
rt1[which(rt1$Group == 'control'),]$Group<-0
rt1$Group<-as.numeric(rt1$Group)

### ========== 训练集：全模型（保留原逻辑，但改列名避免覆盖） ==========
f_full <- glm(Group~., data = rt1, family = binomial())

library(ROCR)
prob_full <- predict(f_full, newdata = rt1, type = "response")
pred_full <- prediction(prob_full, rt1$Group)
performance(pred_full, "auc")@y.values[[1]]   # 0.8217

rt1$prob_full <- prob_full        # 改名为 prob_full，不占用 prob

library(timeROC)
library(pROC)

roc1_full <- roc(rt1$Group, rt1$prob_full)
roc2<- roc(rt1$Group, rt1$CCL11)
roc3<- roc(rt1$Group, rt1$CLDN1)
roc4<- roc(rt1$Group, rt1$POSTN)
roc5<- roc(rt1$Group, rt1$THBS2)
roc6<- roc(rt1$Group, rt1$IL1RN)
roc7<- roc(rt1$Group, rt1$IL1B)
roc8<- roc(rt1$Group, rt1$FAP)
roc9<- roc(rt1$Group, rt1$SOCS1)
roc10<- roc(rt1$Group, rt1$MMP12)

#ALL
pdf(file="ALL.pdf", width =5, height =5)
plot(roc2, col="orange")
plot(roc3, add=TRUE, col="yellow")
plot(roc4, add=TRUE, col="green")
plot(roc5, add=TRUE, col="blue")
plot(roc6, add=TRUE, col="purple")
plot(roc7, add=TRUE, col="pink")
plot(roc8, add=TRUE, col="brown")
plot(roc9, add=TRUE, col="#BD6263")
plot(roc10, add=TRUE, col="#8EA325")

legend('bottomright',
       c(paste0('CCL11: ', round(roc2$auc, 3)),
         paste0('CLDN1: ', round(roc3$auc, 3)),
         paste0('POSTN: ', round(roc4$auc, 3)),
         paste0('THBS2: ', round(roc5$auc, 3)),
         paste0('IL1RN: ', round(roc6$auc, 3)),
         paste0('IL1B: ', round(roc7$auc, 3)),
         paste0('FAP: ', round(roc8$auc, 3)),
         paste0('SOCS1: ', round(roc9$auc, 3)),
         paste0('MMP12: ', round(roc10$auc, 3))
       ),
       col = c("orange","yellow",'green','blue','purple','pink','brown',
        "#BD6263","#8EA325"),
       lwd = 2, bty = 'n')
dev.off()

### ========== 最终模型：CCL11 + CLDN1（prob 改为 prob_final） ==========
f <- glm(Group~CCL11+CLDN1, data = rt1, family = binomial())
summary(f)

# ---------- 导出模型系数（β、SE、OR、95% CI、P） ----------
coef_tab <- data.frame(
  Term      = names(coef(f)),
  Beta      = as.numeric(coef(f)),
  SE        = as.numeric(sqrt(diag(vcov(f)))),
  OR        = exp(as.numeric(coef(f))),
  OR_low    = exp(as.numeric(coef(f)) - 1.96 * sqrt(diag(vcov(f)))),
  OR_high   = exp(as.numeric(coef(f)) + 1.96 * sqrt(diag(vcov(f)))),
  P         = summary(f)$coefficients[, 4]
)
print(coef_tab)
'''
                   Term       Beta         SE         OR     OR_low    OR_high            P
(Intercept) (Intercept) -3.9003923 0.30972588 0.02023397 0.01102649 0.03713002 2.305992e-36
CCL11             CCL11  0.4179810 0.03846414 1.51889187 1.40859299 1.63782762 1.659677e-27
CLDN1             CLDN1  0.4093874 0.05071363 1.50589499 1.36340970 1.66327094 6.885467e-16
'''
write.table(coef_tab, "model_coefficients.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)



prob_final <- predict(f, newdata = rt1, type = "response")
pred_final <- prediction(prob_final, rt1$Group)
performance(pred_final, "auc")@y.values[[1]]   # 0.777275

rt1$prob_final <- prob_final      # 关键：不再用 prob，避免与 prob_full 冲突

roc1<- roc(rt1$Group, rt1$prob_final)   # 现在明确是 CCL11+CLDN1 的预测值
roc2<- roc(rt1$Group, rt1$CCL11)
roc3<- roc(rt1$Group, rt1$CLDN1)

### ---------- 新增：训练集 95% CI（DeLong 法） ----------
ci_model <- ci.auc(roc1, method = "delong")
ci_CCL11 <- ci.auc(roc2, method = "delong")
ci_CLDN1 <- ci.auc(roc3, method = "delong")

cat("\n===== 训练集 AUC 及 95% CI（DeLong） =====\n")
cat(sprintf("CCL11+CLDN1: AUC = %.3f (95%% CI %.3f–%.3f)\n",
            ci_model[2], ci_model[1], ci_model[3]))
#CCL11+CLDN1: AUC = 0.777 (95% CI 0.751–0.803)
cat(sprintf("CCL11      : AUC = %.3f (95%% CI %.3f–%.3f)\n",
            ci_CCL11[2], ci_CCL11[1], ci_CCL11[3]))
#CCL11      : AUC = 0.741 (95% CI 0.713–0.769)
cat(sprintf("CLDN1      : AUC = %.3f (95%% CI %.3f–%.3f)\n",
            ci_CLDN1[2], ci_CLDN1[1], ci_CLDN1[3]))
#CLDN1      : AUC = 0.723 (95% CI 0.695–0.751)

pdf(file="model_CI.pdf", width =4.5, height =4.5)
plot(roc1,col="red", lwd=2, legacy.axes=TRUE)
plot(roc2,add=TRUE,  col="green", lwd=2)
plot(roc3, add=TRUE, col="blue",  lwd=2)
legend('bottomright',
       c(paste0('Training set: ', round(roc1$auc, 3),
                ' (95% CI ', round(ci_model[1], 3), '–', round(ci_model[3], 3), ')'),
         paste0('CCL11: ', round(roc2$auc, 3),
                ' (95% CI ', round(ci_CCL11[1], 3), '–', round(ci_CCL11[3], 3), ')'),
         paste0('CLDN1: ', round(roc3$auc, 3),
                ' (95% CI ', round(ci_CLDN1[1], 3), '–', round(ci_CLDN1[3], 3), ')')
       ),
       col = c("red", "green", "blue"),
       lwd = 2, bty = 'n', cex = 0.85)
dev.off()



### ========== 验证集：GSE206285 + GSE87466 + GSE128682 ==========
library(pROC)

val_files <- c(
  GSE206285 = "GSE206285_data.txt",
  GSE87466  = "GSE87466_data.txt",
  GSE128682 = "GSE128682_data.txt"
)

val_list  <- list()
ci_list   <- list()

for (ds in names(val_files)) {
  cat("\n========================================\n")
  cat("处理验证集：", ds, "\n")
  
  val <- read.table(val_files[ds], header = TRUE, sep = "\t", row.names = 1)
  
  # --- 检查必需列 ---
  need_cols <- c("Group", "CCL11", "CLDN1")
  miss <- setdiff(need_cols, colnames(val))
  if (length(miss) > 0) {
    stop(ds, " 缺少列：", paste(miss, collapse = ", "),
         "。请检查基因名大小写（应为 CCL11、CLDN1）和 Group 列是否存在")
  }
  
  # --- 统一 Group 编码为 0/1（UC=1, control=0） ---
  if (is.character(val$Group) || is.factor(val$Group)) {
    val$Group_num <- ifelse(val$Group %in% c("UC", "uc"), 1, 0)
  } else {
    val$Group_num <- as.numeric(val$Group)
  }
  cat("样本数：", nrow(val),
      "  UC =", sum(val$Group_num == 1),
      "  control =", sum(val$Group_num == 0), "\n")
  
  # --- 用训练好的 f（CCL11+CLDN1）预测 ---
  pred <- predict(f, newdata = val, type = "response")
  
  roc_val   <- roc(val$Group_num, pred,          quiet = TRUE)
  roc_CCL11 <- roc(val$Group_num, val$CCL11,     quiet = TRUE)
  roc_CLDN1 <- roc(val$Group_num, val$CLDN1,     quiet = TRUE)
  
  ci_val    <- ci.auc(roc_val,    method = "delong")
  ci_CCL11v <- ci.auc(roc_CCL11,  method = "delong")
  ci_CLDN1v <- ci.auc(roc_CLDN1,  method = "delong")
  
  cat(sprintf("CCL11+CLDN1: AUC = %.3f (95%% CI %.3f–%.3f)\n",
              ci_val[2], ci_val[1], ci_val[3]))
  cat(sprintf("CCL11      : AUC = %.3f (95%% CI %.3f–%.3f)\n",
              ci_CCL11v[2], ci_CCL11v[1], ci_CCL11v[3]))
  cat(sprintf("CLDN1      : AUC = %.3f (95%% CI %.3f–%.3f)\n",
              ci_CLDN1v[2], ci_CLDN1v[1], ci_CLDN1v[3]))
  
  # --- 每个验证集单独出 ROC 图 ---
  pdf(file = paste0(ds, "_model_ROC.pdf"), width = 4.5, height = 4.5)
  plot(roc_val, col = "red", lwd = 2, legacy.axes = TRUE,
       main = paste0("Validation Set ROC (", ds, ")"))
  plot(roc_CCL11, add = TRUE, col = "green", lwd = 2)
  plot(roc_CLDN1, add = TRUE, col = "blue",  lwd = 2)
  legend("bottomright",
         legend = c(
           paste0("Model: ", round(ci_val[2], 3),
                  " (95% CI ", round(ci_val[1], 3), "–", round(ci_val[3], 3), ")"),
           paste0("CCL11: ", round(ci_CCL11v[2], 3),
                  " (95% CI ", round(ci_CCL11v[1], 3), "–", round(ci_CCL11v[3], 3), ")"),
           paste0("CLDN1: ", round(ci_CLDN1v[2], 3),
                  " (95% CI ", round(ci_CLDN1v[1], 3), "–", round(ci_CLDN1v[3], 3), ")")
         ),
         col = c("red", "green", "blue"),
         lwd = 2, bty = "n", cex = 0.8)
  dev.off()
  
  # --- 收集结果 ---
  val_list[[ds]] <- list(roc_val = roc_val, roc_CCL11 = roc_CCL11, roc_CLDN1 = roc_CLDN1,
                         ci_val = ci_val, ci_CCL11 = ci_CCL11v, ci_CLDN1 = ci_CLDN1v)
  ci_list[[ds]] <- data.frame(
    Dataset = ds,
    Model   = c("CCL11+CLDN1", "CCL11", "CLDN1"),
    AUC     = c(ci_val[2],    ci_CCL11v[2], ci_CLDN1v[2]),
    CI_low  = c(ci_val[1],    ci_CCL11v[1], ci_CLDN1v[1]),
    CI_high = c(ci_val[3],    ci_CCL11v[3], ci_CLDN1v[3])
  )
}
'''
========================================
处理验证集： GSE206285 
样本数： 568   UC = 550   control = 18 
CCL11+CLDN1: AUC = 0.996 (95% CI 0.991–1.000)
CCL11      : AUC = 0.990 (95% CI 0.975–1.000)
CLDN1      : AUC = 0.980 (95% CI 0.965–0.995)

========================================
处理验证集： GSE87466 
样本数： 108   UC = 87   control = 21 
CCL11+CLDN1: AUC = 0.985 (95% CI 0.966–1.000)
CCL11      : AUC = 0.939 (95% CI 0.894–0.985)
CLDN1      : AUC = 0.936 (95% CI 0.890–0.982)

========================================
处理验证集： GSE128682 
样本数： 30   UC = 14   control = 16 
CCL11+CLDN1: AUC = 0.906 (95% CI 0.800–1.000)
CCL11      : AUC = 0.746 (95% CI 0.558–0.933)
CLDN1      : AUC = 0.920 (95% CI 0.794–1.000)
'''
### ========== 三个验证集合并 ROC 图（每个验证集一条模型 ROC） ==========
pdf("All_validation_model_ROC.pdf", width = 5, height = 5)
cols <- c("red", "blue", "darkgreen")
plot(val_list[[1]]$roc_val, col = cols[1], lwd = 2, legacy.axes = TRUE,
     main = "Validation cohorts (CCL11+CLDN1 model)")
for (i in 2:length(val_list)) {
  plot(val_list[[i]]$roc_val, col = cols[i], lwd = 2, add = TRUE)
}
legend("bottomright",
       legend = sapply(names(val_list), function(ds) {
         ci <- val_list[[ds]]$ci_val
         paste0(ds, ": ", round(ci[2], 3),
                " (95% CI ", round(ci[1], 3), "–", round(ci[3], 3), ")")
       }),
       col = cols[1:length(val_list)], lwd = 2, bty = "n", cex = 0.85)
dev.off()

### ========== 汇总表：训练集 + 三个验证集 ==========
ci_tab_all <- rbind(
  data.frame(
    Dataset = rep("Training (GSE193677)", 3),
    Model   = c("CCL11+CLDN1", "CCL11", "CLDN1"),
    AUC     = c(ci_model[2], ci_CCL11[2], ci_CLDN1[2]),
    CI_low  = c(ci_model[1], ci_CCL11[1], ci_CLDN1[1]),
    CI_high = c(ci_model[3], ci_CCL11[3], ci_CLDN1[3])
  ),
  do.call(rbind, ci_list)
)
print(ci_tab_all)
'''
                         Dataset       Model       AUC    CI_low   CI_high
1           Training (GSE193677) CCL11+CLDN1 0.7772750 0.7513978 0.8031522
2           Training (GSE193677)       CCL11 0.7408982 0.7128463 0.7689500
3           Training (GSE193677)       CLDN1 0.7229651 0.6952806 0.7506496
GSE206285.1            GSE206285 CCL11+CLDN1 0.9962626 0.9906899 1.0000000
GSE206285.2            GSE206285       CCL11 0.9900000 0.9746460 1.0000000
GSE206285.3            GSE206285       CLDN1 0.9795455 0.9645282 0.9945628
GSE87466.1              GSE87466 CCL11+CLDN1 0.9846743 0.9664964 1.0000000
GSE87466.2              GSE87466       CCL11 0.9392447 0.8936235 0.9848658
GSE87466.3              GSE87466       CLDN1 0.9359606 0.8903232 0.9815980
GSE128682.1            GSE128682 CCL11+CLDN1 0.9062500 0.7995273 1.0000000
GSE128682.2            GSE128682       CCL11 0.7455357 0.5582055 0.9328659
GSE128682.3            GSE128682       CLDN1 0.9196429 0.7943188 1.0000000
'''
write.table(ci_tab_all, "AUC_with_95CI_all_cohorts.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


#其他基于的AUC
# 在 model_95%CI.R 的 roc4–roc10 之后加：
cat("POSTN  AUC =", round(roc4$auc, 3), "\n")
#POSTN  AUC = 0.641 
cat("THBS2  AUC =", round(roc5$auc, 3), "\n")
#THBS2  AUC = 0.672 
cat("IL1RN  AUC =", round(roc6$auc, 3), "\n")
#IL1RN  AUC = 0.646 
cat("IL1B   AUC =", round(roc7$auc, 3), "\n")
#IL1B   AUC = 0.634 
cat("FAP    AUC =", round(roc8$auc, 3), "\n")
#FAP    AUC = 0.65 
cat("SOCS1  AUC =", round(roc9$auc, 3), "\n")
#SOCS1  AUC = 0.682
cat("MMP12  AUC =", round(roc10$auc, 3), "\n")
#MMP12  AUC = 0.601 

roc10<- roc(rt1$Group, rt1$MMP12)

# ---------- 导出 9 个基因的单基因 AUC（含 95% CI） ----------
gene_auc <- data.frame(
  Gene   = c("CCL11","CLDN1","POSTN","THBS2","IL1RN","IL1B","FAP","SOCS1","MMP12"),
  AUC    = c(roc2$auc, roc3$auc, roc4$auc, roc5$auc, roc6$auc,
             roc7$auc, roc8$auc, roc9$auc, roc10$auc),
  CI_low = NA, CI_high = NA
)

# 用 DeLong 法补 95% CI
for (i in 1:9) {
  r <- list(roc2, roc3, roc4, roc5, roc6, roc7, roc8, roc9, roc10)[[i]]
  ci <- ci.auc(r, method = "delong")
  gene_auc$AUC[i]     <- ci[2]
  gene_auc$CI_low[i]  <- ci[1]
  gene_auc$CI_high[i] <- ci[3]
}

print(gene_auc)
'''
   Gene       AUC    CI_low   CI_high
1 CCL11 0.7408982 0.7128463 0.7689500
2 CLDN1 0.7229651 0.6952806 0.7506496
3 POSTN 0.6412628 0.6116600 0.6708655
4 THBS2 0.6716763 0.6426786 0.7006740
5 IL1RN 0.6455313 0.6162616 0.6748011
6  IL1B 0.6337547 0.6038070 0.6637024
7   FAP 0.6499956 0.6203922 0.6795989
8 SOCS1 0.6821025 0.6530706 0.7111345
9 MMP12 0.6010811 0.5706750 0.6314873
'''
write.table(gene_auc, "single_gene_AUC.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


# ==========================================================
# 诉求 5：预设 cutoff 下的 Sensitivity / Specificity + PR 曲线
# ==========================================================

# ---------- 5.1 训练集约登指数点作为 cutoff ----------
cut_train <- coords(roc1, "best", best.method = "youden",
                    ret = "threshold", transpose = FALSE)
cut_value <- as.numeric(cut_train$threshold[1])
cat("Training Youden cutoff =", round(cut_value, 4), "\n")
#Training Youden cutoff = 0.6396

# ---------- 5.2 各队列在 cutoff 下的敏感度/特异度 ----------
perf_list <- list()

# 训练集
cm_train <- coords(roc1, cut_value, ret = c("sensitivity","specificity"),
                   transpose = FALSE)
perf_list[["Training"]] <- data.frame(
  Cohort = "Training (GSE193677)", Cutoff = cut_value,
  Sensitivity = cm_train$sensitivity, Specificity = cm_train$specificity
)

# 三个验证集
for (ds in names(val_list)) {
  r <- val_list[[ds]]$roc_val
  cm <- coords(r, cut_value, ret = c("sensitivity","specificity"),
               transpose = FALSE)
  perf_list[[ds]] <- data.frame(
    Cohort = ds, Cutoff = cut_value,
    Sensitivity = cm$sensitivity, Specificity = cm$specificity
  )
}

perf_tab <- do.call(rbind, perf_list)
print(perf_tab)
'''
                        Cohort    Cutoff Sensitivity Specificity
Training  Training (GSE193677) 0.6396434   0.7560694   0.6659341
GSE206285            GSE206285 0.6396434   0.9909091   0.8888889
GSE87466              GSE87466 0.6396434   1.0000000   0.0000000
GSE128682            GSE128682 0.6396434   1.0000000   0.2500000
'''
write.table(perf_tab, "sensitivity_specificity_at_youden.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

# ---------- 5.2b 各队列内部的 Youden cutoff（补充：平台特异常数） ----------
perf_internal <- list()
for (nm in c("Training", names(val_list))) {
  r <- if (nm == "Training") roc1 else val_list[[nm]]$roc_val
  cm <- coords(r, "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity"),
               transpose = FALSE)
  perf_internal[[nm]] <- data.frame(
    Cohort = if (nm == "Training") "Training (GSE193677)" else nm,
    Internal_cutoff = as.numeric(cm$threshold[1]),
    Sensitivity = cm$sensitivity,
    Specificity = cm$specificity
  )
}
perf_internal_tab <- do.call(rbind, perf_internal)
print(perf_internal_tab)
'''
                        Cohort Internal_cutoff Sensitivity Specificity
Training  Training (GSE193677)       0.6396434   0.7560694   0.6659341
GSE206285            GSE206285       0.7635756   0.9545455   1.0000000
GSE87466              GSE87466       0.8699931   0.9425287   0.9523810
GSE128682            GSE128682       0.8368034   0.8571429   0.8750000
'''
write.table(perf_internal_tab, "sensitivity_specificity_internal_youden.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


# ---------- 5.3 PR 曲线（GSE206285，类别不平衡 550:18） ----------
library(PRROC)

val206285 <- read.table("GSE206285_data.txt", header = TRUE, sep = "\t", row.names = 1)
pred206285 <- predict(f, newdata = val206285, type = "response")

# --- 关键修复：用与主循环完全相同的编码逻辑 ---
cat("\n===== PR 曲线诊断 =====\n")
cat("Group 原始类型:", class(val206285$Group), "\n")
#Group 原始类型: integer 
cat("Group 原始取值:\n"); print(table(val206285$Group))
#Group 原始取值:
#0   1 
#18 550 

if (is.character(val206285$Group) || is.factor(val206285$Group)) {
  y206285 <- ifelse(val206285$Group %in% c("UC", "uc"), 1, 0)
} else {
  y206285 <- as.numeric(val206285$Group)
}

cat("编码后分布:\n"); print(table(y206285))
#编码后分布:
#y206285
#0   1 
#18 550 

cat("pred NA 数:", sum(is.na(pred206285)), "\n")
#pred NA 数: 0 
cat("UC 组 pred 范围:", range(pred206285[y206285 == 1], na.rm = TRUE), "\n")
#UC 组 pred 范围: 0.4838315 0.9679872 
cat("control 组 pred 范围:", range(pred206285[y206285 == 0], na.rm = TRUE), "\n")
#control 组 pred 范围: 0.3212253 0.7629161 
cat("是否完美分离:",
    min(pred206285[y206285 == 1], na.rm = TRUE) >
      max(pred206285[y206285 == 0], na.rm = TRUE), "\n")
#是否完美分离: FALSE

# --- 计算 PR 曲线 ---
pr206285 <- pr.curve(scores.class0 = pred206285[y206285 == 1],
                     scores.class1 = pred206285[y206285 == 0],
                     curve = TRUE)

auprc_integral <- pr206285$auc.integral
auprc_dg       <- pr206285$auc.davis.goadrich
baseline       <- mean(y206285 == 1)

cat("\nAUPRC (integral)       =", round(auprc_integral, 4), "\n")
#AUPRC (integral)       = 0.9999 
cat("AUPRC (Davis-Goadrich) =", round(auprc_dg, 4), "\n")
#AUPRC (Davis-Goadrich) = 0.9999 
cat("Baseline precision     =", round(baseline, 4), "\n")
#Baseline precision     = 0.9683 

# --- 最终报告值 ---
auprc_final  <- ifelse(is.na(auprc_dg), auprc_integral, auprc_dg)
auprc_method <- ifelse(is.na(auprc_dg), "integral", "Davis-Goadrich")

if (is.na(auprc_final)) {
  pr_df <- data.frame(recall    = pr206285$curve[, 1],
                      precision = pr206285$curve[, 2])
  pr_df <- pr_df[order(pr_df$recall), ]
  auprc_final <- sum(diff(pr_df$recall) *
                       (head(pr_df$precision, -1) + tail(pr_df$precision, -1)) / 2,
                     na.rm = TRUE)
  auprc_method <- "manual trapezoid"
}

cat("\n最终报告 AUPRC =", round(auprc_final, 4),
    "（方法:", auprc_method, "）\n")
#最终报告 AUPRC = 0.9999 （方法: Davis-Goadrich ）
cat("相对 baseline 提升 =", round(auprc_final - baseline, 4), "\n")
#相对 baseline 提升 = 0.0316 

# --- 出图 ---
pdf("PR_curve_GSE206285.pdf", width = 4.5, height = 4.5)
plot(pr206285,
     main = paste0("PR curve (GSE206285)\nAUPRC = ", round(auprc_final, 3),
                   " (baseline = ", round(baseline, 3), ")"),
     col = "red", lwd = 2)
abline(h = baseline, lty = 2, col = "grey40")
legend("bottomleft",
       legend = c(paste0("AUPRC = ", round(auprc_final, 3)),
                  paste0("Baseline = ", round(baseline, 3))),
       col = c("red", "grey40"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
dev.off()

write.table(data.frame(Cohort      = "GSE206285",
                       AUPRC       = auprc_final,
                       AUPRC_method = auprc_method,
                       Baseline    = baseline,
                       Improvement = auprc_final - baseline),
            "AUPRC_GSE206285.txt", row.names = FALSE, sep = "\t", quote = FALSE)




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\06Model\\02nom'
setwd(path)

#data
rt<-read.table('diffGeneExp_p0.05.txt',header = T,sep='\t',row.names=1)
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Type')

rowsid<-rownames(rt)
rowsid<-toupper(rowsid)
rownames(rt)<-rowsid

rt<-rt[hub$V1,]
rt<-t(rt)
rt<-data.frame(rt)

rt$Sample<-rownames(rt)

rt<-merge(group,rt)
rownames(rt)<-rt[,1]
rt<-rt[,-1]
colnames(rt)[1]<-'Group'
table(rt$Group)
#control  UC
#455     865 

rt1<-rt
dim(rt1)  #[1] 1320   3

#plot
library(Hmisc)
library(rms)

Mydata<-rt1

#设置二元分类变量
Mydata[which(Mydata$Group == 'UC'),'Group']<-1
Mydata[which(Mydata$Group == 'control'),'Group']<-0
Mydata$Group<-as.factor(Mydata$Group)
table(Mydata$Group)

#清理无限值
anyNA(Mydata) # 是否有NA
sapply(Mydata, function(x) sum(!is.finite(x))) # 是否有Inf或-Inf

#绘制列线图
library(rms)
Mydata_complete <- Mydata[complete.cases(Mydata[, c('Group', 'CCL11', 'CLDN1')]), ]
dd <- datadist(Mydata_complete)
options(datadist='dd')

##CCL11+CLDN1
f_lrm <- lrm(Group ~ CCL11 + CLDN1, data = Mydata_complete, x=TRUE, y=TRUE)
summary(f_lrm)

#plot1
#以后先观察刻度信息在绘图
nomogram(f_lrm, fun=function(x)1/(1+exp(-x)), fun.at=c(0.1, 0.9), abbrev=F, conf.int=F, verbose=T)

nom <- nomogram(f_lrm, fun=function(x)1/(1+exp(-x)),
                fun.at=c(0.1, 0.9),
                funlabel="Prob of UC",
                conf.int=FALSE,
                abbrev=FALSE)
plot(nom)

#save
pdf(file="nom.pdf", width =7.5, height =5)
plot(nom)
dev.off()



rm(list=ls())

# ==========================================================
# 路径设置
# ==========================================================
path <- 'D:\\LZ-25234\\UC\\06Model\\07nom2'
setwd(path)

model_dir  <- 'D:\\LZ-25234\\UC\\06Model\\01model'

val_files <- c(
  GSE206285 = 'D:\\LZ-25234\\UC\\06Model\\01model\\GSE206285_data.txt',
  GSE87466  = 'D:\\LZ-25234\\UC\\06Model\\01model\\GSE87466_data.txt',
  GSE128682 = 'D:\\LZ-25234\\UC\\06Model\\01model\\GSE128682_data.txt'
)

library(data.table)
library(rms)
library(rmda)
library(pROC)

# ==========================================================
# 1. 重建训练模型（与 dca.R / model_95%CI.R 一致）
# ==========================================================
rt    <- read.table(file.path(model_dir, 'diffGeneExp_p0.05.txt'),
                    header=TRUE, sep='\t', row.names=1)
hub   <- read.table(file.path(model_dir, 'venn.txt'), header=FALSE, sep='\t')
group <- read.table(file.path(model_dir, 'group.txt'), header=TRUE, sep='\t')
colnames(group) <- c('Sample','Type')

rownames(rt) <- toupper(rownames(rt))
rt <- rt[hub$V1, , drop=FALSE]
rt <- t(rt); rt <- data.frame(rt); rt$Sample <- rownames(rt)
rt <- merge(group, rt, by='Sample')
rownames(rt) <- rt[,1]; rt <- rt[,-1]; colnames(rt)[1] <- 'Group'

rt1 <- rt
rt1[which(rt1$Group == 'UC'),]$Group <- 1
rt1[which(rt1$Group == 'control'),]$Group <- 0
rt1$Group <- as.numeric(rt1$Group)

f <- glm(Group ~ CCL11 + CLDN1, data = rt1, family = binomial())
cat('Training model:\n'); print(coef(f))

# ==========================================================
# 工具函数：稳健 PDF 绘图（保证设备关闭）
# ==========================================================
safe_plot_pdf <- function(file, width, height, expr) {
  if (dev.cur() > 1) dev.off()            # 清掉可能残留的设备
  pdf(file, width = width, height = height)
  tryCatch({
    expr
  }, error = function(e) {
    cat('  ⚠ 绘图失败：', conditionMessage(e), '\n')
  }, finally = {
    if (dev.cur() > 1) dev.off()          # 关键：无论成败都关闭
  })
}

# ==========================================================
# 2. 循环处理三个验证集
# ==========================================================
cal_stats_all <- list()
dca_stats_all <- list()

for (ds in names(val_files)) {
  
  cat('\n========================================\n')
  cat('处理验证集：', ds, '\n')
  cat('========================================\n')
  
  res <- tryCatch({
    
    val <- read.table(val_files[ds], header=TRUE, sep='\t', row.names=1)
    
    need <- c('Group','CCL11','CLDN1')
    miss <- setdiff(need, colnames(val))
    if (length(miss) > 0) stop('缺少列：', paste(miss, collapse=', '))
    
    # --- Group 编码 ---
    if (is.character(val$Group) || is.factor(val$Group)) {
      val$y <- ifelse(val$Group %in% c('UC','uc'), 1, 0)
    } else {
      val$y <- as.numeric(val$Group)
    }
    val$prob <- predict(f, newdata = val, type = 'response')
    
    roc_val <- roc(val$y, val$prob, quiet = TRUE)
    ci_auc  <- ci.auc(roc_val, method = 'delong')
    
    cat('样本数：', nrow(val),
        '  UC =', sum(val$y == 1),
        '  control =', sum(val$y == 0), '\n')
    cat('预测概率范围：', round(range(val$prob), 4), '\n')
    cat(sprintf('原始 AUC = %.3f (95%% CI %.3f–%.3f)\n',
                ci_auc[2], ci_auc[1], ci_auc[3]))
    
    # ======================================================
    # 2.1 Platt 重新校准
    # ======================================================
    platt <- glm(y ~ prob, data = val, family = binomial())
    val$prob_recal <- predict(platt, newdata = val, type = 'response')
    cat('Platt 系数：intercept =', round(coef(platt)[1], 3),
        ' slope =', round(coef(platt)[2], 3), '\n')
    
    # ======================================================
    # 2.2 校准曲线（两面板：原始 | 重新校准）→ PDF
    # ======================================================
    dd <- datadist(val); options(datadist = 'dd')
    B_use <- ifelse(nrow(val) < 50, 200, 1000)
    
    fit_orig  <- lrm(y ~ prob,       data = val, x = TRUE, y = TRUE)
    fit_recal <- lrm(y ~ prob_recal, data = val, x = TRUE, y = TRUE)
    
    cal_orig  <- calibrate(fit_orig,  method = 'boot', B = B_use)
    cal_recal <- calibrate(fit_recal, method = 'boot', B = B_use)
    
    safe_plot_pdf(paste0('calibration_', ds, '_orig_vs_recal.pdf'),
                  width = 8, height = 4.2, {
                    par(mfrow = c(1, 2))
                    
                    plot(cal_orig,
                         xlim = c(0,1), ylim = c(0,1),
                         xlab = 'Predicted Probability',
                         ylab = 'Observed Probability',
                         main = paste0('Original (', ds, ')'),
                         subtitles = FALSE, conf.int = TRUE)
                    abline(0, 1, lty = 2, lwd = 2, col = '#224444')
                    
                    plot(cal_recal,
                         xlim = c(0,1), ylim = c(0,1),
                         xlab = 'Predicted Probability',
                         ylab = 'Observed Probability',
                         main = paste0('Recalibrated (', ds, ')'),
                         subtitles = FALSE, conf.int = TRUE)
                    abline(0, 1, lty = 2, lwd = 2, col = '#224444')
                    
                    par(mfrow = c(1, 1))
                  })
    
    # --- 校准统计量 ---
    s_orig  <- val.prob(val$prob,       val$y, pl = FALSE)
    s_recal <- val.prob(val$prob_recal, val$y, pl = FALSE)
    
    cal_tab <- data.frame(
      Cohort        = rep(ds, 2),
      Calibration   = c('Original', 'Recalibrated'),
      n             = nrow(val),
      n_UC          = sum(val$y == 1),
      n_ctrl        = sum(val$y == 0),
      AUC           = round(ci_auc[2], 4),
      AUC_CI_low    = round(ci_auc[1], 4),
      AUC_CI_high   = round(ci_auc[3], 4),
      Brier         = round(c(s_orig['Brier'],     s_recal['Brier']), 4),
      Slope         = round(c(s_orig['Slope'],     s_recal['Slope']), 4),
      Intercept     = round(c(s_orig['Intercept'], s_recal['Intercept']), 4)
    )
    print(cal_tab)
    
    # ======================================================
    # 2.3 DCA：与训练集 dca.R 完全一致的 rmda 写法 → PDF
    #     每个队列分别输出 original 和 recalibrated 两张图
    # ======================================================
    dca_df <- data.frame(y = val$y,
                         prob = val$prob,
                         prob_recal = val$prob_recal)
    
    dca_orig <- decision_curve(y ~ prob,
                               family = binomial(link='logit'),
                               data = dca_df,
                               thresholds = seq(0, 1, by = 0.01),
                               bootstraps = 500,
                               confidence.intervals = 0.95)
    
    dca_recal <- decision_curve(y ~ prob_recal,
                                family = binomial(link='logit'),
                                data = dca_df,
                                thresholds = seq(0, 1, by = 0.01),
                                bootstraps = 500,
                                confidence.intervals = 0.95)
    
    cat('DCA derived.data 行数：Original =', nrow(dca_orig$derived.data),
        '  Recal =', nrow(dca_recal$derived.data), '\n')
    
    # --- 原始模型 DCA 图 ---
    safe_plot_pdf(paste0('DCA_', ds, '_original.pdf'),
                  width = 5, height = 5, {
                    plot_decision_curve(dca_orig,
                                        curve.names = paste0('Original (', ds, ')'),
                                        confidence.intervals = FALSE,
                                        standardize = FALSE,
                                        legend.position = 'right')
                  })
    
    # --- 重新校准后 DCA 图 ---
    safe_plot_pdf(paste0('DCA_', ds, '_recalibrated.pdf'),
                  width = 5, height = 5, {
                    plot_decision_curve(dca_recal,
                                        curve.names = paste0('Recalibrated (', ds, ')'),
                                        confidence.intervals = FALSE,
                                        standardize = FALSE,
                                        legend.position = 'right')
                  })
    
    # --- DCA 概括指标 ---
    get_nb <- function(dca_fit) {
      nb <- dca_fit$derived.data
      col_nb  <- intersect(c('netbenefit','net_benefit'), colnames(nb))[1]
      col_all <- intersect(c('netbenefit_all','net_benefit_all'), colnames(nb))[1]
      col_thr <- intersect(c('thresholds','threshold'), colnames(nb))[1]
      if (is.na(col_nb) || is.na(col_all) || is.na(col_thr)) {
        return(list(max_model=NA, max_all=NA, thr_above='none'))
      }
      ok <- is.finite(nb[[col_nb]]) & is.finite(nb[[col_all]])
      if (sum(ok) == 0) return(list(max_model=NA, max_all=NA, thr_above='none'))
      above <- nb[[col_thr]][ok][ nb[[col_nb]][ok] > nb[[col_all]][ok] ]
      thr_above <- if (length(above) > 0) {
        paste0(round(min(above), 2), '-', round(max(above), 2))
      } else 'none'
      list(max_model = round(max(nb[[col_nb]][ok]), 4),
           max_all   = round(max(nb[[col_all]][ok]), 4),
           thr_above = thr_above)
    }
    
    s_o <- get_nb(dca_orig)
    s_r <- get_nb(dca_recal)
    
    dca_tab <- data.frame(
      Cohort      = rep(ds, 2),
      Calibration = c('Original', 'Recalibrated'),
      Max_net_benefit_model = c(s_o$max_model, s_r$max_model),
      Max_net_benefit_all   = c(s_o$max_all,   s_r$max_all),
      Threshold_range_above_all = c(s_o$thr_above, s_r$thr_above)
    )
    print(dca_tab)
    
    list(cal = cal_tab, dca = dca_tab)
    
  }, error = function(e) {
    cat('  ⚠ 队列', ds, '处理失败：', conditionMessage(e), '\n')
    if (dev.cur() > 1) dev.off()
    NULL
  })
  
  if (!is.null(res)) {
    cal_stats_all[[ds]] <- res$cal
    dca_stats_all[[ds]] <- res$dca
  }
}

# ==========================================================
# 3. 汇总输出
# ==========================================================
cat('\n========================================\n')
cat('校准统计量汇总\n')
cat('========================================\n')
if (length(cal_stats_all) > 0) {
  cal_all <- do.call(rbind, cal_stats_all)
  print(cal_all)
  write.table(cal_all, 'calibration_stats_all_validation.txt',
              row.names=FALSE, sep='\t', quote=FALSE)
}

cat('\n========================================\n')
cat('DCA 净收益汇总\n')
cat('========================================\n')
if (length(dca_stats_all) > 0) {
  dca_all <- do.call(rbind, dca_stats_all)
  print(dca_all)
  write.table(dca_all, 'DCA_stats_all_validation.txt',
              row.names=FALSE, sep='\t', quote=FALSE)
}

cat('\n完成。输出文件（PDF）：\n')
cat('  - calibration_<cohort>_orig_vs_recal.pdf（两面板）\n')
cat('  - DCA_<cohort>_original.pdf\n')
cat('  - DCA_<cohort>_recalibrated.pdf\n')
cat('  - calibration_stats_all_validation.txt\n')
cat('  - DCA_stats_all_validation.txt\n')




rm(list=ls())

path<-'D:\\LZ-25234\\UC\\06Model\\03dca'
setwd(path)

#data
rt<-read.table('diffGeneExp_p0.05.txt',header = T,sep='\t',row.names=1)
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Type')

rowsid<-rownames(rt)
rowsid<-toupper(rowsid)
rownames(rt)<-rowsid

rt<-rt[hub$V1,]
rt<-t(rt)
rt<-data.frame(rt)

rt$Sample<-rownames(rt)

rt<-merge(group,rt)
rownames(rt)<-rt[,1]
rt<-rt[,-1]
colnames(rt)[1]<-'Group'
table(rt$Group)
#control      UC 
#    455     864 


#训练最终模型
rt1<-rt
dim(rt1)
#[1] 1320 5

head(rt1)

#处理分组
rt1[which(rt1$Group == 'UC'),]$Group<-1
rt1[which(rt1$Group == 'control'),]$Group<-0
rt1$Group<-as.numeric(rt1$Group)


#-----------------------------------
library(rmda)

merge_rt<-rt1
dcaData<-as.data.frame(lapply(merge_rt, function(x) as.numeric(as.character(x))))

###CCL11+CLDN1
#构架基础DCA模型
baseline.model <- decision_curve(Group~CCL11+CLDN1,
                                 family=binomial(link='logit'),
                                 data = dcaData,
                                 thresholds = seq(0,1, by = 0.01),
                                 bootstraps = 1000,
                                 confidence.intervals= 0.95)
#plot
plot_decision_curve(baseline.model,  
                    curve.names = "model",
                    confidence.intervals=F, #不显示置信区间
                    standardize = FALSE,legend.position = 'right')


pdf("dca.pdf", onefile=F,width =5, height =5)
plot_decision_curve(baseline.model,  
                    curve.names = "model",
                    confidence.intervals=F, #不显示置信区间
                    legend.position = 'right',
                    standardize = FALSE)
dev.off()


#矫正曲线
library(rms)
dcaData$Group<-factor(dcaData$Group)

##
#fit <- glm(Group ~ CCL11+CLDN1,
#           data = dcaData,
#           family = binomial)

library(rms)
dd <- datadist(dcaData)
options(datadist = "dd")

fit_lrm <- lrm(Group ~ CCL11+CLDN1, data = dcaData,
               x = TRUE,y = TRUE)
cal <- calibrate(fit_lrm, method = "boot", B = 1000) # 建议设定自助法次数如100或1000
plot(cal, xlab = "Prediced Probability", ylab = "Observed Probability", lwd = 2)


pdf("cal1.pdf", onefile=F,width =4.5, height =5)
plot(cal,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Prediced Probability",
     ylab = "Observed Probability",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.8,
     #subtitles = FALSE,
     legend = FALSE)
lines(cal[,c("predy","calibrated.corrected")], 
      type = 'l', #连线的类型，可以是"p","b","o"
      lwd = 3, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = "#2166AC") #连线的颜色

lines(cal[,c("predy","calibrated.orig")],type="l",pch=16,lwd=3,col="tomato")

abline(0,1,
       lty = 2, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = "#224444") #对角线的颜色
legend(0.6,0.2,
       c("Apparent","Bias-corrected","Ideal"), 
       lty = c(2,1,1), 
       lwd = c(2,3,3), 
       col = c("black","#2166AC","tomato"), 
       bty = "n")
dev.off()


##校准曲线修改
library(rms)

dcaData$Group <- factor(dcaData$Group)
dd <- datadist(dcaData)
options(datadist = "dd")

fit_lrm <- lrm(Group ~ CCL11 + CLDN1, data = dcaData, x = TRUE, y = TRUE)

# 校准，默认自助法和1000次重复
cal <- calibrate(fit_lrm, method = "boot", B = 100)

pdf("cal1.pdf", width = 4.5, height = 4)
plot(cal,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     cex.lab = 1.2,
     cex.axis = 1,
     cex.main = 1.2,
     cex.sub = 0.8,
     subtitles=FALSE, 
     conf.int=TRUE)
abline(0, 1, lty = 2, lwd = 2, col = "#224444")
dev.off()




rm(list=ls())

#path
path<-'D:\\LZ-25234\\UC\\06Model\\05boxplot'
setwd(path)

#data
rt<-read.table('diffGeneExp_p0.05.txt',header = T,sep='\t',row.names=1)
hub<-read.table('final.model.txt',header = F,sep='\t')
group<-read.table('group.txt',header=T,sep='\t')
colnames(group)<-c('Sample','Group')

rt<-rt[hub$V1,]

rt<-t(rt)
rt<-rt[group$Sample,]
data = data.frame()
rows = na.omit(match(rownames(rt), group$Sample))
table(group$Group)
#control      UC  
#    455      865 
for(i in colnames(rt)){
  a = as.data.frame(cbind(expression=rt[,i],gene=i,Group=as.vector(group[rows,2])))
  a[,1] = as.numeric(a[,1])
  test=wilcox.test(expression ~ Group, data = a)
  pValue = test$p.value
  data=rbind(data,a)
}

#save result
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

#plot
library(ggplot2)
library(ggpubr)

data=read.table("data.txt",sep="\t",header=T,check.names=F)

data$Group = factor(data$Group, levels = c("control", "UC"))


##箱线图绘制
#SDC1
#data1<-data[which(data$gene == 'SDC1'),]
p = ggboxplot(data, x="gene", y="expression", fill = "Group",
            ylab = "Gene expression",
            xlab = "", 
            palette = c("#00A087B2","#E64B35FF"),#"#00A087B2","#E64B35FF"
            #add='jitter',
            #notch = TRUE
            )+
  theme(text = element_text(size = 10))
#p = p+rotate_x_text(60)
p=p+stat_compare_means(aes(group=Group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
p

ggsave(filename = "boxplot.pdf", plot = p, dpi = 500, width = 3, height = 3。5)




rm(list=ls())

# ==========================================================
# 路径设置
# ==========================================================
model_path <- 'D:\\LZ-25234\\UC\\06Model\\06model2'    # 模型目录
setwd(model_path)

# 临床信息 CSV（data.R 导出）路径；如已复制到 model_path，改成本地文件名即可
pheno_path <- 'D:\\LZ-25234\\UC\\00data\\00train\\GSE193677\\GSE193677.csv'

# ==========================================================
# 加载包
# ==========================================================
library(data.table)
library(pROC)
library(limma)

# ==========================================================
# 1. 读入表达谱、hub 基因、分组（沿用 model.R 逻辑）
# ==========================================================
rt    <- read.table('D://LZ-25234//UC//06Model//01model//diffGeneExp_p0.05.txt', header=TRUE, sep='\t', row.names=1)
hub   <- read.table('D://LZ-25234//UC//06Model//01model//venn.txt',              header=FALSE, sep='\t')
group <- read.table('D://LZ-25234//UC//06Model//01model//group.txt',             header=TRUE,  sep='\t')
colnames(group) <- c('Sample','Type')

# 基因名统一大写
rownames(rt) <- toupper(rownames(rt))

# 只保留 hub 基因（与 model.R 一致）
rt <- rt[hub$V1, , drop=FALSE]
rt <- t(rt)
rt <- data.frame(rt)
rt$Sample <- rownames(rt)

# 合并分组
rt1 <- merge(group, rt, by='Sample')
rownames(rt1) <- rt1$Sample
rt1 <- rt1[, -1]
colnames(rt1)[1] <- 'Group'
table(rt1$Group)   # 期望 control 455 / UC 865
#control      UC 
#   455     865 

# ==========================================================
# 2. 构建最终模型（CCL11 + CLDN1）
# ==========================================================
if (!all(c('CCL11','CLDN1') %in% colnames(rt1))) {
  stop("CCL11 或 CLDN1 不在 rt1 中，请检查 venn.txt 是否包含这两个基因")
}

rt1$Group_num <- ifelse(rt1$Group == 'UC', 1, 0)

f <- glm(Group_num ~ CCL11 + CLDN1, data = rt1, family = binomial())
summary(f)

# 模型预测概率
rt1$prob <- predict(f, newdata = rt1, type = 'response')

# 训练集总体 AUC（应约 0.777）
roc_all <- roc(rt1$Group_num, rt1$prob, quiet=TRUE)
cat('Overall training AUC =', round(as.numeric(auc(roc_all)), 4), '\n')
#Overall training AUC = 0.7773 

# ==========================================================
# 3. 合并临床信息（typere:ch1 + patient_id）
# ==========================================================
pheno <- fread(pheno_path)
pheno$patient_id <- sub('.*reGRID_([0-9]+)_Biopsy_.*', '\\1', pheno$title)
stopifnot(sum(!grepl('reGRID_[0-9]+_Biopsy_', pheno$title)) == 0)

clin <- pheno[, c('geo_accession', 'patient_id', 'typere:ch1')]
colnames(clin) <- c('Sample', 'Patient', 'typere')

rt1$Sample <- rownames(rt1)
rt1 <- merge(rt1, clin, by='Sample', all.x=TRUE)
rownames(rt1) <- rt1$Sample
rt1 <- rt1[, -1]                       # 去掉 Sample 列

# 检查匹配情况
cat('\nGroup × typere 交叉表（确认没有未匹配的样本）：\n')
print(table(rt1$Group, rt1$typere, useNA='ifany'))
#            I NonI
# control   0  455
# UC      291  574

# ==========================================================
# 4. 分层子集：I+control 与 NonI+control
# ==========================================================
uc_I    <- rt1[rt1$Group == 'UC' & rt1$typere == 'I',    ]
uc_NonI <- rt1[rt1$Group == 'UC' & rt1$typere == 'NonI', ]
ctrl    <- rt1[rt1$Group == 'control', ]

sub_I    <- rbind(uc_I,    ctrl)
sub_NonI <- rbind(uc_NonI, ctrl)

cat('\n子集样本数：\n')
cat('I    UC 活检:', nrow(uc_I),    '来自患者:', length(unique(uc_I$Patient)),    '\n')
#I    UC 活检: 291 来自患者: 229 
cat('NonI UC 活检:', nrow(uc_NonI), '来自患者:', length(unique(uc_NonI$Patient)), '\n')
#NonI UC 活检: 574 来自患者: 375
cat('Control 活检 :', nrow(ctrl),   '来自患者:', length(unique(ctrl$Patient)),    '\n')
#Control 活检 : 455 来自患者: 241 

# ==========================================================
# 5. 点估计 AUC
# ==========================================================
roc_I    <- roc(sub_I$Group_num,    sub_I$prob,    quiet=TRUE)
roc_NonI <- roc(sub_NonI$Group_num, sub_NonI$prob, quiet=TRUE)

auc_I    <- as.numeric(auc(roc_I))
auc_NonI <- as.numeric(auc(roc_NonI))

cat('\n===== 点估计 AUC =====\n')
cat('AUC (I + control)    =', round(auc_I, 4), '\n')
#AUC (I + control)    = 0.8924 
cat('AUC (NonI + control) =', round(auc_NonI, 4), '\n')
#AUC (NonI + control) = 0.7189
cat('ΔAUC (I - NonI)      =', round(auc_I - auc_NonI, 4), '\n')
#ΔAUC (I - NonI)      = 0.1734

# ==========================================================
# 6. 患者聚类 bootstrap（一次 bootstrap 同时估计 I、NonI、ΔAUC）
# ==========================================================
set.seed(20240101)
B <- 2000

pat_I    <- unique(uc_I$Patient)
pat_NonI <- unique(uc_NonI$Patient)
pat_ctrl <- unique(ctrl$Patient)

boot_res <- replicate(B, {
  ctrl_b <- do.call(rbind, lapply(
    sample(pat_ctrl, length(pat_ctrl), TRUE),
    function(id) ctrl[ctrl$Patient == id, , drop=FALSE]))
  I_b <- do.call(rbind, lapply(
    sample(pat_I, length(pat_I), TRUE),
    function(id) uc_I[uc_I$Patient == id, , drop=FALSE]))
  N_b <- do.call(rbind, lapply(
    sample(pat_NonI, length(pat_NonI), TRUE),
    function(id) uc_NonI[uc_NonI$Patient == id, , drop=FALSE]))
  
  d_I <- rbind(I_b, ctrl_b)
  d_N <- rbind(N_b, ctrl_b)
  
  a_I <- as.numeric(auc(roc(d_I$Group_num, d_I$prob, quiet=TRUE)))
  a_N <- as.numeric(auc(roc(d_N$Group_num, d_N$prob, quiet=TRUE)))
  
  c(auc_I=a_I, auc_NonI=a_N, delta=a_I-a_N)
})

cat('\n===== 患者聚类 bootstrap 95% CI =====\n')
ci_I    <- quantile(boot_res['auc_I',],    c(.025, .975), na.rm=TRUE)
ci_NonI <- quantile(boot_res['auc_NonI',], c(.025, .975), na.rm=TRUE)
ci_delta<- quantile(boot_res['delta',],    c(.025, .975), na.rm=TRUE)

cat('AUC_I    95% CI:', round(ci_I,    4), '\n')
#AUC_I    95% CI: 0.868 0.9142 
cat('AUC_NonI 95% CI:', round(ci_NonI, 4), '\n')
#AUC_NonI 95% CI: 0.6893 0.7465 
cat('ΔAUC     95% CI:', round(ci_delta, 4), '\n')
#ΔAUC     95% CI: 0.1431 0.2044 

# ==========================================================
# 7. 结果表输出
# ==========================================================
result_tab <- data.frame(
  Subset           = c('I + control', 'NonI + control'),
  n_UC_biopsy      = c(nrow(uc_I), nrow(uc_NonI)),
  n_UC_patient     = c(length(pat_I), length(pat_NonI)),
  n_ctrl_biopsy    = c(nrow(ctrl), nrow(ctrl)),
  n_ctrl_patient   = c(length(pat_ctrl), length(pat_ctrl)),
  AUC              = c(auc_I, auc_NonI),
  CI_low           = c(ci_I[1], ci_NonI[1]),
  CI_high          = c(ci_I[2], ci_NonI[2])
)
print(result_tab)
'''
          Subset n_UC_biopsy n_UC_patient n_ctrl_biopsy n_ctrl_patient       AUC    CI_low   CI_high
1    I + control         291          229           455            241 0.8923606 0.8680394 0.9142459
2 NonI + control         574          375           455            241 0.7189302 0.6892942 0.7464591
'''
write.table(result_tab, 'stratified_AUC_I_vs_NonI.txt',
            row.names=FALSE, sep='\t', quote=FALSE)

delta_tab <- data.frame(
  Comparison = 'I - NonI',
  delta_AUC  = auc_I - auc_NonI,
  CI_low     = ci_delta[1],
  CI_high    = ci_delta[2],
  crosses_zero = (ci_delta[1] <= 0 & ci_delta[2] >= 0)
)
print(delta_tab)
'''
     Comparison delta_AUC    CI_low   CI_high crosses_zero
2.5%   I - NonI 0.1734304 0.1430844 0.2043599        FALSE
'''
write.table(delta_tab, 'stratified_deltaAUC.txt',
            row.names=FALSE, sep='\t', quote=FALSE)

# ==========================================================
# 8. 分层 ROC 图
# ==========================================================
pdf('ROC_training_I_vs_NonI.pdf', width=5, height=5)
plot(roc_I, col='red',
     main='Training set stratified by inflammation status')
plot(roc_NonI, col='blue', add=TRUE)
legend('bottomright',
       legend=c(
         paste0('I + control (n=', nrow(sub_I), ', AUC=', round(auc_I,3), ')'),
         paste0('NonI + control (n=', nrow(sub_NonI), ', AUC=', round(auc_NonI,3), ')')
       ),
       col=c('red','blue'), lwd=2, bty='n')
dev.off()

cat('\n完成。输出文件：\n')
cat('  - stratified_AUC_I_vs_NonI.txt\n')
cat('  - stratified_deltaAUC.txt\n')
cat('  - ROC_training_I_vs_NonI.pdf\n')


# ==========================================================
# 9. UC vs CD（利用 GSE193677 中的 CD 样本）
#    预处理严格对齐 data.R：log2 → impute.knn → normalizeBetweenArrays
# ==========================================================
cat('\n========================================\n')
cat('9. UC vs CD\n')
cat('========================================\n')

library(impute)

counts_file <- 'D:/LZ-25234/UC/00data/00train/GSE193677/deduplication_counts.txt'
if (!file.exists(counts_file)) {
  stop('未找到 deduplication_counts.txt，请修改 counts_file 路径')
}

# ---- 9.1 读入原始 counts（已 avereps 合并重复基因，但未做任何数值转换） ----
exp_all <- read.table(counts_file, header=TRUE, sep='\t', row.names=1)
cat('全队列 counts 维度：', dim(exp_all), '\n')
#全队列 counts 维度： 39374 2471 
cat('counts 中 NA 数：', sum(is.na(exp_all)), '\n')
3counts 中 NA 数： 0 

# ---- 9.2 三步预处理，与训练集 data.R 完全一致 ----
# (1) log2(count + 1)
exp_all <- log2(exp_all + 1)

# (2) KNN 填补（训练集用的就是 impute.knn 默认参数 k=10）
mat <- impute.knn(as.matrix(exp_all))
exp_all <- mat$data
cat('KNN 后 NA 数：', sum(is.na(exp_all)), '\n')
#KNN 后 NA 数： 0 

# (3) 分位数归一化（对样本/列）
exp_all <- normalizeBetweenArrays(exp_all)
cat('预处理后范围：', round(range(exp_all), 3), '\n')
#预处理后范围： 0 19.306 

# ---- 9.3 提取 CCL11、CLDN1 ----
need_genes <- c('CCL11','CLDN1')
if (!all(need_genes %in% rownames(exp_all))) {
  stop('全队列 counts 中缺少 CCL11 或 CLDN1')
}

sub_expr <- t(exp_all[need_genes, , drop=FALSE])
sub_expr <- data.frame(sub_expr, check.names = FALSE)
sub_expr$Sample <- rownames(sub_expr)

# ---- 9.4 合并疾病标签（UC / CD） ----
clin_all <- pheno[, c('geo_accession','patient_id','ibd_disease:ch1')]
colnames(clin_all) <- c('Sample','Patient','disease')

cd_uc <- merge(clin_all[clin_all$disease %in% c('UC','CD'), ],
               sub_expr, by='Sample')
cd_uc$y_CD <- ifelse(cd_uc$disease == 'CD', 1, 0)

cat('UC / CD 样本数：', sum(cd_uc$y_CD == 0), '/', sum(cd_uc$y_CD == 1), '\n')
#UC / CD 样本数： 865 / 1151 
cat('UC / CD 患者数：',
    length(unique(cd_uc$Patient[cd_uc$y_CD == 0])), '/',
    length(unique(cd_uc$Patient[cd_uc$y_CD == 1])), '\n')
#UC / CD 患者数： 419 / 495

# ---- 9.5 固定系数 f 预测 ----
cd_uc$prob_fixed <- predict(f, newdata = cd_uc, type = 'response')
roc_CD_fixed <- roc(cd_uc$y_CD, cd_uc$prob_fixed, quiet = TRUE)
ci_CD_fixed  <- ci.auc(roc_CD_fixed, method = 'delong')
cat(sprintf('固定系数：UC vs CD AUC = %.3f (95%% CI %.3f–%.3f)\n',
            ci_CD_fixed[2], ci_CD_fixed[1], ci_CD_fixed[3]))
#固定系数：UC vs CD AUC = 0.568 (95% CI 0.543–0.593)

# ---- 9.6 重新拟合（UC vs CD 判别任务） ----
f_CD <- glm(y_CD ~ CCL11 + CLDN1, data = cd_uc, family = binomial())
cd_uc$prob_refit <- predict(f_CD, newdata = cd_uc, type = 'response')
roc_CD_refit <- roc(cd_uc$y_CD, cd_uc$prob_refit, quiet = TRUE)
ci_CD_refit  <- ci.auc(roc_CD_refit, method = 'delong')
cat(sprintf('重新拟合：UC vs CD AUC = %.3f (95%% CI %.3f–%.3f)\n',
            ci_CD_refit[2], ci_CD_refit[1], ci_CD_refit[3]))
#重新拟合：UC vs CD AUC = 0.569 (95% CI 0.544–0.594)

# ---- 9.7 患者聚类 bootstrap（固定系数版本） ----
set.seed(20240101)
pat_CD <- unique(cd_uc$Patient)
boot_CD <- replicate(2000, {
  p <- sample(pat_CD, length(pat_CD), TRUE)
  d <- do.call(rbind, lapply(p, function(id) cd_uc[cd_uc$Patient == id, , drop=FALSE]))
  as.numeric(auc(roc(d$y_CD, d$prob_fixed, quiet=TRUE)))
})
ci_CD_boot <- quantile(boot_CD, c(.025, .975), na.rm=TRUE)
cat('固定系数 患者聚类 bootstrap 95% CI:',
    round(ci_CD_boot[1], 3), '–', round(ci_CD_boot[2], 3), '\n')
#固定系数 患者聚类 bootstrap 95% CI: 0.542 – 0.595 

# ---- 9.8 出图 ----
pdf('ROC_UC_vs_CD.pdf', width = 5, height = 5)
plot(roc_CD_fixed, col='red', lwd=2, legacy.axes=TRUE,
     main='UC vs CD (CCL11 + CLDN1)')
plot(roc_CD_refit, col='blue', lwd=2, add=TRUE)
legend('bottomright',
       legend=c(
         paste0('Fixed f: AUC=', round(ci_CD_fixed[2],3),
                ' (', round(ci_CD_fixed[1],3), '-', round(ci_CD_fixed[3],3), ')'),
         paste0('Refit: AUC=', round(ci_CD_refit[2],3),
                ' (', round(ci_CD_refit[1],3), '-', round(ci_CD_refit[3],3), ')')
       ),
       col=c('red','blue'), lwd=2, bty='n', cex=0.8)
dev.off()

# ---- 9.9 输出（含 bootstrap CI） ----
cd_tab <- data.frame(
  Analysis     = c('Fixed f', 'Refit'),
  AUC          = c(ci_CD_fixed[2], ci_CD_refit[2]),
  CI_low       = c(ci_CD_fixed[1], ci_CD_refit[1]),
  CI_high      = c(ci_CD_fixed[3], ci_CD_refit[3]),
  Boot_CI_low  = c(ci_CD_boot[1], NA),
  Boot_CI_high = c(ci_CD_boot[2], NA),
  n_UC         = sum(cd_uc$y_CD == 0),
  n_CD         = sum(cd_uc$y_CD == 1)
)
print(cd_tab)
'''
     Analysis       AUC    CI_low   CI_high Boot_CI_low Boot_CI_high n_UC n_CD
2.5%  Fixed f 0.5676471 0.5426387 0.5926555   0.5417083    0.5950409  865 1151
        Refit 0.5687429 0.5436784 0.5938074          NA           NA  865 1151
'''
write.table(cd_tab, 'UC_vs_CD_AUC.txt',
            row.names=FALSE, sep='\t', quote=FALSE)


# ==========================================================
# 10. I vs NonI（UC 内部）
# ==========================================================
cat('\n========================================\n')
cat('10. I vs NonI（UC 内部）\n')
cat('========================================\n')

uc <- rt1[rt1$Group == 'UC' & rt1$typere %in% c('I','NonI'), ]
uc$I_status <- ifelse(uc$typere == 'I', 1, 0)

cat('UC 活检中 I / NonI：', sum(uc$I_status == 1), '/', sum(uc$I_status == 0), '\n')
#UC 活检中 I / NonI： 291 / 574 
cat('患者数：', length(unique(uc$Patient)), '\n')
#患者数： 419 

# ---- 10.1 用现有 UC vs control 模型预测值（prob）区分 I vs NonI ----
roc_I_vs_NonI <- roc(uc$I_status, uc$prob, quiet = TRUE)
ci_I_vs_NonI  <- ci.auc(roc_I_vs_NonI, method = 'delong')
cat(sprintf('现有模型预测值：I vs NonI AUC = %.3f (95%% CI %.3f–%.3f)\n',
            ci_I_vs_NonI[2], ci_I_vs_NonI[1], ci_I_vs_NonI[3]))
#现有模型预测值：I vs NonI AUC = 0.778 (95% CI 0.744–0.812)

# ---- 10.2 UC 内部重新拟合 I vs NonI ----
f_I <- glm(I_status ~ CCL11 + CLDN1, data = uc, family = binomial())
uc$prob_I <- predict(f_I, newdata = uc, type = 'response')
roc_I_refit <- roc(uc$I_status, uc$prob_I, quiet = TRUE)
ci_I_refit  <- ci.auc(roc_I_refit, method = 'delong')
cat(sprintf('UC 内部重新拟合：I vs NonI AUC = %.3f (95%% CI %.3f–%.3f)\n',
            ci_I_refit[2], ci_I_refit[1], ci_I_refit[3]))
#UC 内部重新拟合：I vs NonI AUC = 0.814 (95% CI 0.782–0.845)

# ---- 10.3 患者聚类 bootstrap（现有模型预测值） ----
set.seed(20240101)
pat_uc <- unique(uc$Patient)
boot_IN <- replicate(2000, {
  p <- sample(pat_uc, length(pat_uc), TRUE)
  d <- do.call(rbind, lapply(p, function(id) uc[uc$Patient == id, , drop=FALSE]))
  as.numeric(auc(roc(d$I_status, d$prob, quiet=TRUE)))
})
ci_IN_boot <- quantile(boot_IN, c(.025, .975), na.rm=TRUE)
cat('现有模型预测值 患者聚类 bootstrap 95% CI:',
    round(ci_IN_boot[1], 3), '–', round(ci_IN_boot[2], 3), '\n')
#现有模型预测值 患者聚类 bootstrap 95% CI: 0.742 – 0.811 

# ---- 10.4 出图 ----
pdf('ROC_I_vs_NonI_UC.pdf', width = 5, height = 5)
plot(roc_I_vs_NonI, col='red', lwd=2, legacy.axes=TRUE,
     main='I vs NonI within UC')
plot(roc_I_refit, col='blue', lwd=2, add=TRUE)
legend('bottomright',
       legend=c(
         paste0('Existing model prob: AUC=', round(ci_I_vs_NonI[2],3),
                ' (', round(ci_I_vs_NonI[1],3), '-', round(ci_I_vs_NonI[3],3), ')'),
         paste0('Refit (I~CCL11+CLDN1): AUC=', round(ci_I_refit[2],3),
                ' (', round(ci_I_refit[1],3), '-', round(ci_I_refit[3],3), ')')
       ),
       col=c('red','blue'), lwd=2, bty='n', cex=0.8)
dev.off()

# ---- 10.5 输出 ----
in_tab <- data.frame(
  Analysis = c('Existing model prob', 'Refit within UC'),
  AUC      = c(ci_I_vs_NonI[2], ci_I_refit[2]),
  CI_low   = c(ci_I_vs_NonI[1], ci_I_refit[1]),
  CI_high  = c(ci_I_vs_NonI[3], ci_I_refit[3]),
  n_I      = sum(uc$I_status == 1),
  n_NonI   = sum(uc$I_status == 0),
  n_patient= length(unique(uc$Patient))
)
print(in_tab)
'''
             Analysis       AUC    CI_low  CI_high n_I n_NonI n_patient
1 Existing model prob 0.7778357 0.7437693 0.811902 291    574       419
2     Refit within UC 0.8135889 0.7824497 0.844728 291    574       419
'''
write.table(in_tab, 'I_vs_NonI_UC_AUC.txt',
            row.names=FALSE, sep='\t', quote=FALSE)

cat('\n完成。新增输出文件：\n')
cat('  - UC_vs_CD_AUC.txt\n')
cat('  - ROC_UC_vs_CD.pdf\n')
cat('  - I_vs_NonI_UC_AUC.txt\n')
cat('  - ROC_I_vs_NonI_UC.pdf\n')





