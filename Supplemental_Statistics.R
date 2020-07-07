
library("phyloseq")
library("parallel")
library("vegan")
library("ggplot2")
library("reshape2")
library("plyr")
library("forcats")
library("dplyr")
library("egg")

# Samples and cateogories for differential abundance tests:
# 
# 	season
# MB-E1	Sum
# MB-E2	Sum
# MB-E3	Sum
# MB-E4	Sum
# MB-E5	Sum
# MB-E6	Sum
# MB-A1b	WS
# MB-A2b	WS
# MB-A3b	WS
# MB-A4b	WS
# MB-A5b	WS
# MB-B2	WS
# MB-B3	WS
# MB-B4	WS
# MB-B5	WS
# MB-B6	WS
# MB-C1	Sum
# MB-C2	Sum
# MB-C3	Sum
# MB-C4	Sum
# MB-C5	Sum
# MB-C6	Sum
# MB-D1	Sum
# MB-D2	Sum
# MB-D3	Sum
# MB-D4	Sum
# MB-D5	Sum
# MB-D6	Sum

# pruning and ordination


wh0 = genefilter_sample(arctic_lakes, filterfun_sample(function(x) x > 2), A=0.1*nsamples(arctic_lakes))
arctic_lakes_pr = prune_taxa(wh0, arctic_lakes)
arctic_lakes_pr
write.csv(otu_table(arctic_lakes_pr), file = "arctic_lakes_pr.csv")
pr_asv <- read.csv("arctic_lakes_pr.csv", header = TRUE, row.names = 1)
d.czm <- cmultRepl(pr_asv,  label=0, method="CZM") # zero replacement
d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) # CLR transformation
a.dist <- dist(d.clr) # Aitchinson Distance Calculation
Aitchinson_PCA <- ordinate(arctic_lakes_pr, 'MDS', distance=a.dist)
aitch_ord1 <- plot_ordination(physeq=arctic_fm_pr, ordination=Aitchinson_PCA, justDF = TRUE)

# Figure
arct_aitch <- ggplot(data=aitch_ord1, aes(y =Axis.2, x = Axis.1)) + geom_point(size = 4, aes(shape=Season, color = Season, fill=Season)) + scale_color_manual(values=nmdscols) + 
scale_fill_manual(values=nmdscols)  +  
scale_shape_manual(values=c(9,7,5,6)) + 
labs(y="Axis-2 [13% Variation]", x = "Axis-1 [16.3% Variation]") +
theme_classic()
arct_aitch

# Differential Abundance
                                                       

write.csv(otu_table(arctic_lakes_pr), file = "lakes_aldex2_pruned.csv")
lakes_asv <- read.csv("lakes_aldex2_pruned.csv", header = TRUE, row.names = 1)
lakes_asv_T <- t(lakes_asv)
write.csv(sample_data(arctic_lakes_pr), file = "lakes_AldMet.csv")
lake_aldex2_meta <- read.csv(file="lakes_AldMet.csv", header = TRUE, row.names=1)
conds <- as.character(lake_aldex2_meta$season)
lake_aldex2 <- aldex.clr(lakes_asv_T, conds = conds, mc.samples=128)
lake_effects <- aldex.effect(lake_aldex2, conds)
lake_stat <- aldex.ttest(lake_aldex2, conds, paired.test= FALSE)
lake_DAT <- data.frame(lake_effects, lake_stat)
aldex.plot(lake_DAT, method = "welch", cutoff=0.05)
lake_DAT = lake_DAT[order(lake_DAT$diff.btw, na.last=NA), ]
alpha = .05
sigtab = lake_DAT[(lake_DAT$we.eBH < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(arctic_lakes_pr)[rownames(sigtab), ], "matrix"))

#Figure
 
 theme_set(theme_bw(base_size=13)) 
#sigtabgen = subset(sigtab1, !is.na(Genus)) this removes things with Genus as NA
sigtabgen = sigtab1
# Phylum order
x = tapply(sigtabgen$Median.CLR.Difference, sigtabgen$Class, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Class = factor(as.character(sigtabgen$Class), levels=names(x))
# Genus order
x = tapply(sigtabgen$Median.CLR.Difference, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
diff_plot <- ggplot(sigtabgen, aes(y=Genus, x=Median.CLR.Difference, color=Class)) +
  theme(text = element_text(face="bold", color="#000000")) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4, shape=9) + 
  geom_errorbarh(aes(xmax = Median.CLR.Difference + 0.5*diff.win, xmin = Median.CLR.Difference - 0.5*diff.win, height = .4)) +
  scale_color_manual(values=pal31) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
diff_plot                                                     

# Regression of qPCR data (https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/)

ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue") +
  labs(size = 2, title = paste("Adjusted R-squared = ",signif(summary(fit)$adj.r.squared, 1),
                     " P-value =",signif(summary(fit)$coef[2,4], 1)))
}


mod3 <- lm(Ordination_Axis1 ~ Log10_Methane)
summary(mod3)

Call:
lm(formula = Ordination_Axis1 ~ Log10_Methane)

Residuals:
    Min      1Q  Median      3Q     Max 
-39.048  -9.996   2.582  15.582  29.278 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)   
(Intercept)     35.331     13.127   2.691  0.01227 * 
Log10_Methane   -9.736      3.474  -2.802  0.00946 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 19.35 on 26 degrees of freedom
Multiple R-squared:  0.232,	Adjusted R-squared:  0.2024 
F-statistic: 7.853 on 1 and 26 DF,  p-value: 0.009457

axis1_regression <- ggplotRegression(mod3)


mod1 <- lm(Log10_PmoA ~ Log10_Methane)
summary(mod1)

Call:
lm(formula = Log10_PmoA ~ Log10_Methane)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.37404 -0.20625 -0.06867  0.09377  1.10010 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)     1.4787     0.3666   4.033 0.000779 ***
Log10_Methane   0.3180     0.1064   2.990 0.007850 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3734 on 18 degrees of freedom
Multiple R-squared:  0.3319,	Adjusted R-squared:  0.2948 
F-statistic: 8.941 on 1 and 18 DF,  p-value: 0.00785
