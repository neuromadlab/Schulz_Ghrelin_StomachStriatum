# TUE008 NIMG Study
# Corinna Schulz 

######################################################################
# Read in Ghrelin Blood data from Lab Files
# Make first Plots to inspect acyl and des-acyl data 
# save final csv with all blood data and preprocessed ghrelin values 

######################################################################
# (1) SET UP

# Install and load all required libraries
if (!require("librarian")) install.packages("librarian")
librarian::shelf(readxl, lme4, lmerTest, foreign, MASS, readr,
ggplot2, ggpubr, cowplot, gridExtra, viridis, tidybayes, dplyr, httpgd, 
languageserver, dabestr, stringr)


# Set all themes
theme_set(theme_cowplot(font_size = 12))

# Set Colors 

color_Placebo <- "darkblue"  
color_Ghrelin <- "darkgoldenrod"  

# Set all Paths
setwd(getwd())
path_in <- "./input/" 
path_out <- "./output/" 
sub_prepro <-  "0.preprocessingQC/"

if (file.exists(paste(path_out, sub_prepro, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_prepro, sep = ""))}

######################################################################
# (2) LOAD DATA and PREP DATA

# 2.1. Load Ghrelin Values

# Ghrelin Values ("raw results" and group membership)
d_acylG <- read_excel(path=paste(path_in, "TUE008_NIMG_Acyl.xlsx", sep = ""))
d_desaG <- read_excel(path=paste(path_in, "TUE008_NIMG_Desacyl.xlsx", sep = ""))

d_acylG <- d_acylG %>%
  mutate(`Durchschn. Einzelkonz. (pg/ml)` = na_if(`Durchschn. Einzelkonz. (pg/ml)`, "NoCalc"),
         `Durchschn. Einzelkonz. (pg/ml)` = as.numeric(`Durchschn. Einzelkonz. (pg/ml)`)) %>%
  filter(`Replikations-info` == "1/2") %>%
  mutate(
    is_repeat = str_detect(`Proben-ID`, "W$"),
    base_id = str_remove(`Proben-ID`, "W$")
  )

# Wenn es Wiederholungen gibt, nimm nur diese; sonst Originale
d_acylG <- d_acylG %>%
  filter(!is.na(`Durchschn. Einzelkonz. (pg/ml)`)) %>%
  mutate(
    ID = str_extract(base_id, "P\\d+") %>%
         str_remove("P") %>%
         as.integer() %>%
         sprintf("%06d", .),
    Session = str_extract(base_id, "S\\d") %>%
              str_remove("S") %>%
              as.integer(),
    Timepoint = str_extract(base_id, "T\\d") %>%
                str_remove("T") %>%
                as.integer(),
    F_AG = `Durchschn. Einzelkonz. (pg/ml)`
  ) %>%
  # Group by unique measurements
  group_by(ID, Session, Timepoint) %>%
  # Keep only the repeat (W) if available
  arrange(desc(is_repeat)) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, Session, Timepoint, F_AG,base_id, is_repeat)


# Step 1: count timepoints per ID-session
session_counts <- d_acylG %>%
  count(ID, Session, name = "n_timepoints")

# Step 2: keep only sessions that have 3 timepoints
complete_sessions <- session_counts %>%
  filter(n_timepoints == 3)

# Step 3: count number of such sessions per ID
complete_IDs <- complete_sessions %>%
  count(ID, name = "n_sessions") %>%
  filter(n_sessions == 2)

# Step 4: how many such IDs?
nrow(complete_IDs)


# Clean desacyl (here only W, no noCalc) 
d_desaG <- d_desaG %>%
  filter(`Replikatinfo` == "1/2") %>%
  filter(!is.na(`Durchschn.   Einzelkonz.     pg/ml)`)) %>%
    mutate(
    is_repeat = str_detect(`Proben-ID`, "W$"),
    base_id = str_remove(`Proben-ID`, "W$")
  )

# Wenn es Wiederholungen gibt, nimm nur diese; sonst Originale
d_desaG <- d_desaG %>%
  mutate(
    ID = str_extract(base_id, "P\\d+") %>%
         str_remove("P") %>%
         as.integer() %>%
         sprintf("%06d", .),
    Session = str_extract(base_id, "S\\d") %>%
              str_remove("S") %>%
              as.integer(),
    Timepoint = str_extract(base_id, "T\\d") %>%
                str_remove("T") %>%
                as.integer(),
    F_DG = `Durchschn.   Einzelkonz.     pg/ml)`
  ) %>%
  # Group by unique measurements
  group_by(ID, Session, Timepoint) %>%
  # Keep only the repeat (W) if available
  arrange(desc(is_repeat)) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, Session, Timepoint, F_DG,base_id, is_repeat)

# Step 1: count timepoints per ID-session
session_counts <- d_desaG %>%
  count(ID, Session, name = "n_timepoints")

# Step 2: keep only sessions that have 3 timepoints
complete_sessions <- session_counts %>%
  filter(n_timepoints == 3)

# Step 3: count number of such sessions per ID
complete_IDs <- complete_sessions %>%
  count(ID, name = "n_sessions") %>%
  filter(n_sessions == 2)

# Step 4: how many such IDs?
nrow(complete_IDs)

d_acylG$fID  <- factor(as.numeric(d_acylG$ID))
d_desaG$fID  <- factor(as.numeric(d_desaG$ID))

# Add Condition Info 
# Load Participant data / Conditions 
d_blood <- read.csv(paste(path_in,"blood_preprocessed.csv", sep = ""), header = TRUE)

# Prep Confounds
d_blood$cAge <- d_blood$Age - mean(d_blood$Age)
d_blood$cBMI <- d_blood$BMI - mean(d_blood$BMI)
d_blood$cSex <- d_blood$Sex_female - mean(d_blood$Sex_female)
d_blood$fID  <- factor(as.numeric(d_blood$ID))

# Merge Condition file 
d_ghrelin <- merge(d_acylG, d_blood, by = c("fID","Session"))
d_ghrelin <- merge(d_ghrelin, d_desaG, by = c("fID","Session","Timepoint"))

## Correct ghrelin levels (i.e., log transformation)
d_ghrelin$logF_AG <- log(d_ghrelin$F_AG)
d_ghrelin$logF_DG <- log(d_ghrelin$F_DG)

length(unique(d_ghrelin$ID)) # Check! 

d_ghrelin$fGhrelin <- factor(d_ghrelin$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_ghrelin$cSession <- d_ghrelin$Session - mean(d_ghrelin$Session)

# We need to match right VAS Timepoint to right Blood Timepoint  
d_ghrelin <- d_ghrelin %>%
  mutate(Ghrelin_Timepoint = case_when(
    Timepoint == 0 ~ 1,
    Timepoint == 1 ~ 3,
    Timepoint == 2 ~ 4,
  ))

d_ghrelin$fTimepoint <- factor(d_ghrelin$Ghrelin_Timepoint, labels = c("T0", "T2","T3"))

######################################################################
# (3) Ghrelin Plots 



hist1 <- ggplot(d_ghrelin, aes(x = F_AG)) +
 #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 20, fill = "#76767f") +
  xlab(label = 'Blood levels of acyl ghrelin [pg/ml]') +
  ggtitle("Acyl ghrelin (raw)")

hist1log <- ggplot(d_ghrelin, aes(x = logF_AG)) +
  #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 0.2, fill = "#a3a3c9") +
  xlab(label = 'Blood levels of acyl ghrelin [pg/ml]') +
  ggtitle("Acyl ghrelin (log)")

hist2 <- ggplot(d_ghrelin, aes(x = F_DG)) +
  #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 20, fill = "#76767f") +  
  xlab(label = 'Blood levels of des-acyl ghrelin [pg/ml]') +
  ggtitle("Des-acyl ghrelin (raw)")

hist2log <- ggplot(d_ghrelin, aes(x = logF_DG)) +
  #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 0.2, fill = "#a3a3c9") +
  xlab(label = 'Blood levels of des-acyl ghrelin [log(pg/ml)]') +
  ggtitle("Des-acyl ghrelin (log)")

# Save Grid Plot
G <- grid.arrange(hist1, hist1log, hist2, hist2log, ncol=2, nrow=2, 
                  layout_matrix = rbind(c(1,2), c(3,4)),
                  widths=c(6,6), heights=c(6, 6))

ggsave(paste(path_out, sub_prepro, "0.Ghrelin_histograms.png", sep=""), 
        plot = G,  height = 6, width = 10, units = "in", dpi = 600, bg = "white")

# Test for Normality 
shapiro.test(d_ghrelin$logF_DG)
shapiro.test(d_ghrelin$logF_AG)


# Correlation of Fasting Acyl and Desacyl Ghrelin levels 

p1<-
  ggplot(aes(x = F_AG,y = F_DG, fill = fGhrelin),data = d_ghrelin) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  geom_smooth(aes(group = fGhrelin, color = fGhrelin), method = 'rlm', alpha = 0.1) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Ghrelin,color_Placebo)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Ghrelin,color_Placebo)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Correlation of ghrelin levels") +
  facet_wrap(~Timepoint, scale = "free") +
  xlab(label = 'Fasting levels of acyl ghrelin (raw)') +
  ylab(label = 'Fasting levels of des-acyl ghrelin (raw)')

ggsave(paste(path_out, sub_prepro , "0.QC_Ghrelin_Corr_Raw.png", sep=""), 
              plot = p1,  height =6, width = 15, units = "in", dpi = 600, bg = "white")

p3<-
  ggplot(aes(x = logF_AG,y = logF_DG),data = d_ghrelin) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  geom_smooth(color = 'black', group = 1, method = 'rlm', size = 1.5, alpha = 0.5) +
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Ghrelin,color_Placebo)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Ghrelin,color_Placebo)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Correlation of ghrelin levels") +
    facet_wrap(~Timepoint, scale = "free") +
  xlab(label = 'Fasting levels of acyl ghrelin (log)') +
  ylab(label = 'Fasting levels of des-acyl ghrelin (log)')

ggsave(paste(path_out, sub_prepro , "0.QC_Ghrelin_Corr_Log.png", sep=""), 
              plot = p3,  height =6, width = 15, units = "in", dpi = 600, bg = "white")


contrasts(d_ghrelin$fTimepoint) <-  contr.treatment(levels(d_ghrelin$fTimepoint), base = 1)


lm_F_AG <- lmer(logF_AG ~   cSession + fGhrelin* fTimepoint + cBMI + cSex + cAge + (1 + fGhrelin + fTimepoint |fID ) , data = d_ghrelin)
summary(lm_F_AG)

lm_F_DG <- lm(logF_DG ~   cSession + fGhrelin* fTimepoint + cBMI + cSex + cAge , data = d_ghrelin)
summary(lm_F_DG)

# Residualize the Ghrelin Levels
# i.e., regress out cBMI and cAge and cSex (c = centered)

# Fasting levels
lm_F_AG <- lm(logF_AG ~ cBMI + cAge + cSex   , data = d_ghrelin, na.rm = TRUE)
summary(lm_F_AG)
lm_F_DG <- lm(logF_DG ~ cBMI + cAge + cSex  , data = d_ghrelin, na.rm = TRUE)
summary(lm_F_DG)

length(unique(d_ghrelin$fID))
# Linear Model for Residualization
d_ghrelin$res_logF_AG <- lm_F_AG$residuals
d_ghrelin$res_logF_DG <- lm_F_DG$residuals

# Test for Normality 
shapiro.test(d_ghrelin$res_logF_AG)
shapiro.test(d_ghrelin$res_logF_DG)


hist1res <- ggplot(d_ghrelin, aes(x = res_logF_AG)) +
  #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 0.2, fill = "#bbbbf3") +  
  xlab(label = 'Blood levels of acyl ghrelin [res]') +
  ggtitle("Acyl ghrelin (log + res)")

hist2res <- ggplot(d_ghrelin, aes(x =res_logF_DG)) +
  #facet_wrap(~fPlate) +
  geom_histogram(binwidth = 0.2, fill = "#bbbbf3") +
  xlab(label = 'Blood levels of des-acyl ghrelin [res]') +
  ggtitle("Des-acyl ghrelin (log + res)")

# Save Grid Plot
G <- grid.arrange(hist1, hist1log, hist1res, hist2, hist2log, hist2res, ncol=3, nrow=2, 
                  layout_matrix = rbind(c(1,2,3), c(4,5,6)),
                  widths=c(13,13,13), heights=c(6, 6))

ggsave(paste(path_out, sub_prepro, "0.Ghrelin_histograms_reslog.png", sep=""), 
        plot = G,  height = 6, width = 16, units = "in", dpi = 600, bg = "white")

# Test correlations
cor.test(d_ghrelin$res_logF_AG, d_ghrelin$res_logF_DG)
cor.test(d_ghrelin$logF_AG, d_ghrelin$logF_DG)
cor.test(d_ghrelin$F_AG, d_ghrelin$F_DG)

# Test again but residualized
library(emmeans)
lm_F_AG <- lmer(res_logF_AG ~   cSession + fGhrelin* fTimepoint + cBMI + cSex + cAge  +(1  |fID ), data = d_ghrelin)
summary(lm_F_AG)

emm <- emmeans(lm_F_AG, ~ fTimepoint | fGhrelin)
emmeans::contrast(emm, "trt.vs.ctrl", ref = "T2") 

lm_F_AG <- lmer(res_logF_AG ~   cSession + fGhrelin* fTimepoint * cBMI + cBMI + cSex + cAge + (1 + fGhrelin + fTimepoint |fID ) ,  , data = d_ghrelin)
summary(lm_F_AG)

lm_F_DG <- lm(res_logF_DG ~   cSession + fGhrelin* fTimepoint + cBMI + cSex + cAge , data = d_ghrelin)
summary(lm_F_DG)


## Plot Res Log 
p4<-
  ggplot(aes(x = res_logF_AG,y = res_logF_DG),data = d_ghrelin) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  #geom_smooth(method = 'rlm', alpha = 0.1) +
  geom_smooth(aes(group = fGhrelin, color = fGhrelin), method = 'rlm', size = 1.5, alpha = 0.5) +
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Correlation of ghrelin levels") +
    facet_wrap(~fTimepoint, scale = "free") +
  xlab(label = 'Fasting levels of acyl ghrelin (res log)') +
  ylab(label = 'Fasting levels of des-acyl ghrelin (res log)')

ggsave(paste(path_out, sub_prepro , "0.QC_Ghrelin_Corr_ResLog.png", sep=""), 
              plot = p4,  height =6, width = 12, units = "in", dpi = 600, bg = "white")

##########################################################################################
# Density Plots 

plot_AG1 <-
  ggplot(aes(x = res_logF_AG, group = fGhrelin, fill = fGhrelin),data = d_ghrelin) +
  geom_density(alpha=.5, adjust = 2) +
  #geom_vline(xintercept = 0, color = "#3d3b3b", linewidth = 2, linetype = "dashed") +
  #stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.7) +
  #stat_dotsinterval(side = "bottom", scale = 0.7, slab_linewidth = NA) +  
  scale_fill_manual(guide = guide_legend(title=""),values = c(color_Placebo, color_Ghrelin)) +
  theme(legend.position = "top",text = element_text(face = 'bold',size = 20.0),
        axis.text = element_text(face = 'plain',size = 18.0),axis.text.x = element_text(size = 18.0)) +
 facet_wrap(~fTimepoint, scale = "free") +
xlab(label = 'Fasting acyl ghrelin (res log)') +
  ylab(label = 'Density')

ggsave(paste(path_out, sub_prepro , "0.QC_AG_Timepoints.png", sep=""), 
              plot = plot_AG1,  height =6, width = 10, units = "in", dpi = 600, bg = "white")


plot_DG1 <-
  ggplot(aes(x = res_logF_DG, group = fGhrelin, fill = fGhrelin),data = d_ghrelin) +
  geom_density(alpha=.5, adjust = 2) +
  geom_vline(xintercept = 0, color = "#3d3b3b", linewidth = 2, linetype = "dashed") +
  #stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.7) +
  #stat_dotsinterval(side = "bottom", scale = 0.7, slab_linewidth = NA) +  
  scale_fill_manual(guide = guide_legend(title=""),values = c(color_Placebo, color_Ghrelin)) +
  theme(legend.position = "top",text = element_text(face = 'bold',size = 20.0),
        axis.text = element_text(face = 'plain',size = 18.0),axis.text.x = element_text(size = 18.0)) +
   facet_wrap(~fTimepoint, scale = "free") +
    xlab(label = 'Fasting des-acyl ghrelin (res log)') +
  ylab(label = 'Density')

ggsave(paste(path_out, sub_prepro , "0.QC_DG_Timepoints.png", sep=""), 
              plot = plot_DG1,  height =6, width = 10, units = "in", dpi = 600, bg = "white")


######################################################################
# Clean and Save Ghrelin CSV 

d_ghrelin_clean <- select(d_ghrelin, ID, fID, Session, fTimepoint, F_AG, F_DG, logF_AG, logF_DG, res_logF_AG, res_logF_DG)

# Write final csv file for further analysis

# write.csv(d_ghrelin_clean, paste("./input/","TUE008_ghrelin_summary_preprocessed.csv", sep = ""), row.names=FALSE)

