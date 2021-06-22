
library(dplyr)
library(ggplot2)
library(lm.beta)
library(mediation)
library(tidyverse)
library(broom)
library(memisc)
library(boot)


rm(list = ls())
basedir <- "/your/path/to/ADNI/data"

######===========
## FUNCTIONS
######===========
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

elapsed_year <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  ed$year - sd$year
}

bootCorTest <- function(data, i){
  d <- data[i, ]
  res <- cor.test(d$x, d$y)
  c(stat = res$estimate, p.value = res$p.value)
}
#######################################################################################
## PREPARE DATA
#######################################################################################
# tau PET
# cross sectional
dat_tau <- read.csv(paste0(basedir,"/UCBERKELEYAV1451_05_12_20_PROCESSED.csv")) # data has been SUVr corrected using inferior cerebellar grey matter as reference and log-transformed to approach normal distribution
dat_tau$tau_InfTemp_roi1 <- apply(dat_tau[,c("CTX_LH_INFERIORTEMPORAL_SUVR_INFCEREB", "CTX_RH_INFERIORTEMPORAL_SUVR_INFCEREB")],1,function(x)mean(x))
dat_tau$tau_global_roi2 <- dat_tau$global_volume_weighted_SUVR_INFCEREB
#longitudinal
subject_w_long_tau <- names(table(dat_tau$RID))[table(dat_tau$RID)>1]
dat_tau_long <- dat_tau[dat_tau$RID %in% subject_w_long_tau,]
## select earliest and latest tau PET visit
dat_tau_long   <- dat_tau_long %>% group_by(RID) %>% mutate(
  date_first_tau_visit=min(as.Date(EXAMDATE)),
  first_tau_visit=ifelse(as.Date(EXAMDATE)==date_first_tau_visit, 1, 0),
  years_to_first_tau_visit=as.numeric(elapsed_year(as.Date(EXAMDATE), date_first_tau_visit)),
  months_to_first_tau_visit=as.numeric(elapsed_months(as.Date(EXAMDATE), date_first_tau_visit)),
  last_tau_visit=ifelse(months_to_first_tau_visit==max(months_to_first_tau_visit),1,0)
)
dat_tau_long <- dat_tau_long[dat_tau_long$first_tau_visit==1 | dat_tau_long$last_tau_visit==1,]

dat_tau_long$CTX_LH_INFERIORTEMPORAL_SUVR_INFCEREB_CHANGE <- (dat_tau_long$CTX_LH_INFERIORTEMPORAL_SUVR_INFCEREB[dat_tau_long$last_tau_visit==1] - dat_tau_long$CTX_LH_INFERIORTEMPORAL_SUVR_INFCEREB[dat_tau_long$first_tau_visit==1]) / dat_tau_long$years_to_first_tau_visit[dat_tau_long$last_tau_visit==1] 
dat_tau_long$CTX_RH_INFERIORTEMPORAL_SUVR_INFCEREB_CHANGE <- (dat_tau_long$CTX_RH_INFERIORTEMPORAL_SUVR_INFCEREB[dat_tau_long$last_tau_visit==1] - dat_tau_long$CTX_RH_INFERIORTEMPORAL_SUVR_INFCEREB[dat_tau_long$first_tau_visit==1]) / dat_tau_long$years_to_first_tau_visit[dat_tau_long$last_tau_visit==1] 
dat_tau_long$tau_InfTemp_CHANGE <- rowMeans(dat_tau_long[,c("CTX_LH_INFERIORTEMPORAL_SUVR_INFCEREB_CHANGE", "CTX_RH_INFERIORTEMPORAL_SUVR_INFCEREB_CHANGE")])
dat_tau_long$tau_global_CHANGE <- (dat_tau_long$global_volume_weighted_SUVR_INFCEREB[dat_tau_long$last_tau_visit==1] - dat_tau_long$global_volume_weighted_SUVR_INFCEREB[dat_tau_long$first_tau_visit==1]) / dat_tau_long$years_to_first_tau_visit[dat_tau_long$last_tau_visit==1] 
dat_tau$tau_InfTemp_CHANGE <- dat_tau_long$tau_InfTemp_CHANGE[match(dat_tau$RID, dat_tau_long$RID),]
dat_tau$tau_global_CHANGE <- dat_tau_long$tau_global_CHANGE[match(dat_tau$RID, dat_tau_long$RID),]


# amyloid PET
dat_fbp <- read.csv(paste0(basedir, "/UCBERKELEYAV45_05_12_20_PROCESSED.csv")) # FBP-PET closest to tau PET has been selected
dat_fbp$globalFBP_CL <- (196.9 * dat_fbp$SUMMARYSUVR_WHOLECEREBNORM) - 196.03

dat_fbb <- read.csv(paste0(basedir, "/UCBERKELEYFBB_05_12_20_PROCESSES.csv")) # FBB-PET closest to tau PET has been selected
dat_fbb$globalFBB_CL <- (159.08 * dat_fbb$SUMMARYSUVR_WHOLECEREBNORM) - 151.65

# ADNI sample characteristics
dat_merge <- read.csv(paste0(basedir, "/ADNIMERGE_PROCESSED.csv")) # data closest to tau PET has been selected
dat_merge$APOE <- factor(ifelse(dat_merge$APOE4>0,1,0))

# ADNI cognition
dat_cog <- read.csv(paste0(basedir, "//UWNPSYCHSUM_03_26_20_PROCESSED.csv")) # assessment closest to tau PET has been selected

# KL-VS
dat_klotho <- read.csv("/Klotho_coding.csv") # based on imputed ADNI3 genetic data
dat_klotho$KL_VS <- factor(ifelse(dat_klotho$rs9536314 == "GT" & dat_klotho$rs9527025 == "CG", 1 ,0 ), levels = c(1,0), labels = c("carrier","noncarrier"))

## KL expression from Allen brain atlas
# http://human.brain-map.org
# http://figshare.com/articles/A_FreeSurfer_view_of_the_cortical_transcriptome_generated_from_the_Allen_Human_Brain_Atlas/1439749) 
dat_KL_LH_allen <- read.csv(paste0(basedir, "/Allen_to_DK/lh-KL-ExpressionSummary.csv"))
dat_KL_RH_allen <- read.csv(paste0(basedir, "/Allen_to_DK/rh-KL-ExpressionSummary.csv")) #don't use

# df
df <- data.frame(dat_tau,
                 dat_fbp[match(dat_tau$RID, dat_fbp$RID),],
                 dat_fbb[match(dat_tau$RID, dat_fbb$RID),],
                 dat_merge[match(dat_tau$RID, dat_merge$RID),],
                 dat_cog[match(dat_tau$RID, dat_cog$RID),],
                 dat_klotho[match(dat_tau$RID, dat_klotho$RID),] )

# Amyloid measures
df$abeta_global <- df$globalFBP_CL
df$abeta_global[is.na(df$abeta_global)] <- df$globalFBB_CL[is.na(df$abeta_global)]
df$abeta_status <- df$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF
df$abeta_status[is.na(df$abeta_status)] <- df$SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF[is.na(df$abeta_status)]



######################################################################################################################
## STATISTICS - KL-VS x amyloid interaction on tau PET - CROSS-SECTIONAL- Figure 1a - b
######################################################################################################################
df_tmp <- df[(!is.na(df$KL_VS) & !is.na(df$DX) & !is.na(df$abeta_global)),]; nrow(df_tmp)
#551

LM1 <- lm.beta(lm(tau_InfTemp_roi1 ~ Klotho * abeta_global + Age + Sex + DX + APOE + Edu, data = df_tmp));summary(LM1)
LM2 <- lm.beta(lm(tau_global_roi2 ~ Klotho * abeta_global + Age + Sex + DX + APOE + Edu, data = df_tmp));summary(LM2)


# Figure 1a-b
for_annotation <- paste0("beta=", round(LM1$standardized.coefficients[10],2),
                         ", p=", round(summary(LM1)$coef[10,5],3),
                         ", R2=", round(summary(LM1)$r.squared,2),
                         ", N=", nrow(df_tmp) )
ggplot(df_tmp, aes(x=abeta_global, y=tau_InfTemp_roi1, group=factor(KL_VS), color=factor(KL_VS), shape=factor(APOE) )) + 
  geom_point( size = 2, alpha = 0.7)+
  scale_color_manual(values = c("grey50", "blue"))+
  stat_smooth(method = "lm", se=T)+
  theme_minimal()+
  scale_y_continuous(limits = c(0.7,2.9))+
  scale_x_continuous(limits = c(-30,160))+
  xlab("\namyloid-PET CL")+ylab("tau-PET SUVR\n")+
  annotate("text", x=63,y=0.7, label=for_annotation, size=5, color="black")+
  theme(legend.position = "none", axis.text=element_text(size=16), axis.title = element_text(size = 18),
        axis.line = element_line(size = 0.5), 
        panel.grid.major.x =  element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 0.2), axis.text.x = element_text(angle = 0))



######################################################################################################################
## STATISTICS - KL-VS x amyloid interaction on tau PET - LONGITUDINAL- Figure 1c - d
######################################################################################################################
df_tmp <- df[(!is.na(df$tau_InfTemp_CHANGE) & !is.na(df$KL_VS) & & !is.na(df$DX) & !is.na(df$abeta_global)),]; nrow(df_tmp)
#200

LM1 <- lm.beta(lm(tau_InfTemp_CHANGE ~ KL_VS * abeta_global + Age + Sex + DX + APOE + Edu, data = df_tmp));summary(LM1)
LM2 <- lm.beta(lm(tau_global_CHANGE ~ KL_VS * abeta_global + Age + Sex + DX + APOE + Edu, data = df_tmp));summary(LM2)


# Figure 1c - d
for_annotation <- paste0("beta=", round(LM1$standardized.coefficients[10],2),
                         ", p=", round(summary(LM1)$coef[10,5],3),
                         ", R2=", round(summary(LM1)$r.squared,2),
                         ", N=", nrow(df_tmp) )
ggplot(df_tmp, aes(x=abeta_global, y=tau_InfTemp_CHANGE, group=factor(KL_VS), color=factor(KL_VS), shape=APOE )) + 
  geom_point( size = 2, alpha=0.7)+
  scale_color_manual(values = c("grey40", "blue"))+
  stat_smooth(method = "lm")+
  theme_minimal()+
  scale_y_continuous(limits = c(-0.08,0.))+
  xlab("\namyloid-PET CL")+ylab("tau-PET SUVR\nannual change")+
  annotate("text", x=63,y=-0.08, label=for_annotation, size=5, color="black")+
  theme(legend.position = "none", axis.text=element_text(size=16), axis.title = element_text(size = 18),
        axis.line = element_line(size = 0.5), 
        panel.grid.major.x =  element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 0.2), axis.text.x = element_text(angle = 0),
        plot.margin = unit(c(0.1, 0.15, 0.1, 0.1), "in"))


######################################################################################################################
## STATISTICS - spatial correlation between KL-VS and KL expression - CROSS-SECTIONAL- Figure 2e
######################################################################################################################
df_tmp <- df[(!is.na(df$KL_VS) & !is.na(df$DX) & !is.na(df$abeta_global)),]; nrow(df_tmp)
#551

# only use data from left hemisphere, because KL expression is only complete for this side
regions_of_interest <- names(df_tmp)[intersect(grep("CTX_LH",names(df_tmp)), grep("INFCEREB$",names(df_tmp)))]; length(regions_of_interest)
# 34
covariates <- c("KL_VS", "Age", "Sex", "DX", "APOE", "Edu", "abeta_global")
df_tmp_s <- dplyr::select(df_tmp, c(regions_of_interest, covariates))

# compute KL-VS * amyloid interaction on tau PET for each cortical (lh) Freesurfer region
lm_func <- function(y) lm(y ~ KL_VS * abeta_global + Age + Sex + DX + APOE + Edu, data = df_tmp_s)
fit <- purrr::map_dfr(df_tmp_s[1:length(regions_of_interest)], function(x) lm_func(x) %>% broom::tidy())
fit.out <- fit[fit$term=="KL_VScarrier:abeta_global",]
fit.out$p.value.adj <- p.adjust(fit.out$p.value, method = "fdr")
fit.out$name <- regions_of_interest
knitr::kable(fit.out) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)

# correlation
fit.out$KL_expression <- dat_KL_LH_allen$KL[match(fit.out$name, dat_KL_LH_allen$Label.name)]
LM1 <- cor.test(fit.out$statistic*(-1), fit.out$KL_expression)
tmp_boot <- data.frame(x=fit.out$KL_expression, y=fit.out$statistic*(-1))
b <- boot(tmp_boot, bootCorTest, R = 1000)
ci_boot <- boot.ci(b)


# Figure 2e
for_annotation <- paste0("r=", round(LM1$estimate,2),
                         ", p=",round(LM1$p.value,3),
                         ", CI95=",round(ci_boot$normal[2],2), "-",round(ci_boot$normal[3],2) )
ggplot(fit.out, aes(x=KL_expression, y=statistic*(-1))) + 
  geom_point( size = 2, alpha=0.7, color="blue")+
  stat_smooth(method = "lm", se=T, color="grey40")+
  theme_minimal()+
  xlab("\nKlotho mRNA expression")+ylab("Klotho x amyloid-PET\ninteraction effect")+
  annotate("text", x=3.2,y=0.1, label=for_annotation, size=5, color="black")+
  theme(legend.position = "none", axis.text=element_text(size=16), axis.title = element_text(size = 18),
        axis.line = element_line(size = 0.5), 
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 0.2), axis.text.x = element_text(angle = 0))


######################################################################################################################
## STATISTICS - KL-VS effect on ADNI-MEM - CROSS-SECTIONAL- Figure 3
######################################################################################################################
df_tmp <- df[(!is.na(df$KL_VS) & !is.na(df$adni.mem) & !is.na(df$DX) & df$abeta_status == 1),]; nrow(df_tmp)
#229

LM1 <- lm.beta(lm(adni.mem ~ KL_VS + Age + Sex + DX + APOE + Edu, data = df_tmp));summary(LM1)

# Figure 3
for_annotation <- paste0("beta=", round(LM1$standardized.coefficients[9],2),
                         ", p=", round(summary(LM1)$coef[9,5],3),
                         ", R2=", round(summary(LM1)$r.squared,2),
                         ", N=", nrow(df_tmp) )
ggplot(df_tmp, aes(y=adni.mem, x=KL_VS, group=KL_VS, color=KL_VS, shape=APOE )) + 
  geom_boxplot(notch = T, outlier.colour = NA, fatten=2)+
  geom_jitter(size=2, width = 0.1, alpha=0.7)+
  scale_color_manual(values = c("blue", "grey40"))+
  stat_smooth(method = "lm", se=T)+
  theme_minimal()+
  xlab(bquote("\nKL-VS"^het))+ylab("ADNI-MEM\n")+
  annotate("text", x=1.5,y=-2.5, label=for_annotation, size=4, color="black")+
  theme(legend.position = "none", axis.text=element_text(size=16), axis.title = element_text(size = 18),
        axis.line = element_line(size = 0.5), 
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 0.2), axis.text.x = element_text(angle = 0))

######################################################################################################################
## STATISTICS - mediation - CROSS-SECTIONAL- Figure 4
######################################################################################################################
df_tmp <- df[(!is.na(df$KL_VS) & !is.na(df$adni.mem) & !is.na(df$DX) & df$abeta_status == 1),]; nrow(df_tmp)
#229

med.fit <- lm(tau_global_roi2 ~ KL_VS + Age + Sex + DX + APOE + Edu, data = df_tmp); summary(med.fit)
out.fit <- lm(adni.mem ~ KL_VS + tau_global_roi2 + Age + Sex + DX + APOE + Edu, data = df_tmp); summary(out.fit)
set.seed(323232)
med.out <- mediate(med.fit, out.fit, treat = "KL_VS", mediator = "tau_global_roi2",  boot=T, sims = 10000)
summary(med.out)


