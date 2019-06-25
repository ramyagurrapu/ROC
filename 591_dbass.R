data = read.table("591_dbass.txt", sep = "\t", header = TRUE)
##SNPfold
five_percent_bottom_SNPfold = subset(data, data$SNPfold < quantile(data$SNPfold, 0.05))
str(five_percent_bottom_SNPfold)
five_percent_bottom_SNPfold$SNPfold
five_percent_bottom_SNPfold$UID
plot(density(five_percent_bottom_SNPfold$SNPfold), main = "Top 5% SNPfold results of DBASS splice site mutations (n = 30)"
     ,xlab = "Pearson correlation coefficient values from SNPfold")

##remuRNA
five_percent_top_remuRNA = subset(data, data$remuRNA > quantile(data$remuRNA, 0.95))
str(five_percent_top_remuRNA)
five_percent_top_remuRNA$remuRNA
five_percent_top_remuRNA$UID
plot(density(five_percent_top_remuRNA$remuRNA),main = "Top 5% remuRNA results of DBASS splice site mutations (n = 30)"
     ,xlab = "Relative entropy values from remuRNA")

##RNAsnp
five_percent_top_RNAsnp = subset(data, data$RNAsnp < quantile(data$RNAsnp, 0.05))
str(five_percent_top_RNAsnp)
five_percent_top_RNAsnp$RNAsnp
five_percent_top_RNAsnp$UID
plot(density(five_percent_top_RNAsnp$RNAsnp),main = "Top 5% RNAsnp results of DBASS splice site mutations (n = 30)"
     ,xlab = "Relative entropy values from RNAsnp")


five_percent = rbind(five_percent_bottom_SNPfold,five_percent_top_remuRNA,five_percent_top_RNAsnp)
write.csv(five_percent,file = "591_dbass_fivepercent.csv")

_______________________________________________________________________
## Highlight the 5% tails
##ggplot
y <- data$SNPfold3
cutoff <- quantile(y, probs = 0.05)
hist.y <- density(y) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)
the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  theme(legend.position = "none") +
  annotate(geom = 'text', x = cutoff, y = 0.025, color = 'red', label = '5% tail',hjust = -0.1) +
  labs(title="SNPfold results of DBASS splice site mutations", y="Density", x="Pearson correlation coefficient values from SNPfold")
print(the.plot)
###remuRNA

five_percent_top_remuRNA = subset(data, data$remuRNA > quantile(data$remuRNA, 0.95))
str(five_percent_top_remuRNA)
five_percent_top_remuRNA$remuRNA
plot(density(five_percent_top_remuRNA$remuRNA),main = "Top 5% remuRNA results of DBASS splice site mutations (n = 29)"
     ,xlab = "Relative entropy values from remuRNA")
##ggplot
y <- data$remuRNA
cutoff <- quantile(y, probs = 0.95)
hist.y <- density(y) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)
the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  theme(legend.position = "none") +
  annotate(geom = 'text', x = cutoff, y = 0.025, color = 'red', label = '5% tail',hjust = -0.1) +
  labs(title="remuRNA results of DBASS splice site mutations", y="Density", x="Relative entropy values from remuRNA")
print(the.plot)
###RNAsnp

five_percent_bottom_RNAsnp1 = subset(data, data$RNAsnp1 < quantile(data$RNAsnp1, 0.05))
str(five_percent_bottom_RNAsnp1)
plot(density(five_percent_bottom_RNAsnp1$RNAsnp1), main = "Top 5% RNAsnp results of DBASS splice site mutations (n = 29)"
     ,xlab = "P-values on Euclidean distance scores from RNAsnp")
##ggplot
y <- data$RNAsnp1
cutoff <- quantile(y, probs = 0.05)
hist.y <- density(y) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)
the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  theme(legend.position = "none") +
  annotate(geom = 'text', x = cutoff, y = 0.025, color = 'red', label = '5% tail',hjust = -0.1) +
  labs(title="RNAsnp results of DBASS splice site mutations", y="Density", x="P-values on Euclidean distance scores from RNAsnp")
print(the.plot)
###
str(five_percent_bottom_SNPfold3)
five_percent_bottom_SNPfold3$SNPfold3
str(five_percent_top_remuRNA)
five_percent_top_remuRNA$remuRNA
str(five_percent_bottom_RNAsnp1)
five_percent_bottom_RNAsnp1$RNAsnp1

five_percent = rbind(five_percent_bottom_SNPfold3,five_percent_top_remuRNA,five_percent_bottom_RNAsnp1)
write.csv(five_percent,file = "dbass_fivepercent.csv")

