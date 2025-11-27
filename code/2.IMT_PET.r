###############################################################################
# TUE008 Neuroimaging Study: Instrumental Motivation Task (IMT) during PET-MR 
# Corinna Schulz, April 2024 
# Main Tasks of the Script: 
# IMT analysis (effort lmers, model comparison, sensitivity analyses)
# PET BP Plots + combine PET + Effort 
# IMT + Pulses + matches 
###############################################################################

###############################################################################
# (1) SET UP
###############################################################################

# Install and load all required libraries
if (!require("librarian")) install.packages("librarian")
librarian::shelf(readxl, lme4, lmerTest, foreign, MASS, 
                 sjPlot, stargazer, table1, readr,
                 ggplot2, ggpubr, cowplot, viridis, tidybayes, dplyr, httpgd, 
                 languageserver, dabestr, smplot2, emmeans, tidyr, devtools, flextable, 
                 devtools, formattable, gtsummary,)


# Raincloud library (credit Micah Allen Paper)
# remotes::install_github('jorvlan/raincloudplots')
library(raincloudplots)

#Set all themes
theme_set(theme_cowplot(font_size = 12))

# Set Colors
color_Placebo <- "darkblue"  
color_Ghrelin <- "darkgoldenrod"  

color_Phen <- "plum2"   
color_Nimg_Placebo <- "purple4"  

color_PhenLow <- "aquamarine"   
color_Nimg_PlaceboLow <- "cadetblue4"  

# Set all Paths
setwd(getwd())
path_in <- "./input/" 
path_out <- "./output/" 

sub_IMT <- "instrumental_motivation/"
sub_PET <- "PET/"
sub_Pulse <- "Pulse/"

if (file.exists(paste(path_out, sub_IMT, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_IMT, sep = ""))}

if (file.exists(paste(path_out, sub_PET, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_PET, sep = ""))}

if (file.exists(paste(path_out, sub_Pulse, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_Pulse, sep = ""))}

###############################################################################
# 2. LOAD DATA
###############################################################################
# Load PET Data 
d_BPs <- read_excel(paste(path_in,"TUE008_NIMG_BPs_ROIs_28_04_25.xlsx", sep = ""))

# Load IMT data
d_IMT <- read_excel(paste(path_in,"IMT_TUE008_all_output.xlsx", sep = ""))

# Load processed IMT data (includes Inv Slope and Maintenance)
d_IMT <- read_csv(paste(path_in,"IMT_TUE008_Exp_Merg+Seg_AggrTrial_20240807.csv", sep= ""), show_col_types = FALSE)
d_IMT <- d_IMT %>% rename(Reward_Money = Rew_Money, Reward_Magnitude = Rew_Mag, Trial = Trial_ID, ID = Subj_ID, Session = Sess_ID)

# IMT during trial data 
d_IMT_trial <- read_csv(paste(path_in,"IMT_TUE008_Exp_Merg+Seg_20240807.csv", sep= ""), show_col_types = FALSE)
d_IMT_trial <- d_IMT_trial %>% rename(Reward_Money = Rew_Money, Reward_Magnitude = Rew_Mag, Trial = Trial_ID, ID = Subj_ID, Session = Sess_ID)

# Load IMT (Neuroimaging) data combined with EAT data (phenotyping) 
d_IMT_comb <- read_excel(paste(path_in,"T_COMB_NIMG_PHEN.xlsx", sep = ""))

# Load Participant data / Conditions 
d_blood <- read.csv(paste(path_in,"blood_preprocessed.csv", sep = ""), header = TRUE)
d_blood[d_blood == "missing"] <- NA
d_ghrelin <- read.csv(paste(path_in,"TUE008_ghrelin_summary_preprocessed.csv", sep = ""), header = TRUE)

# Prep Confounds
d_blood$cAge <- d_blood$Age - mean(d_blood$Age)
d_blood$cBMI <- d_blood$BMI - mean(d_blood$BMI)
d_blood$cSex <- d_blood$Sex_female - mean(d_blood$Sex_female)

d_ghrelin$fTimepoint <- factor(d_ghrelin$fTimepoint, labels = c("T0", "T2","T3"))
d_ghrelin$fID <- factor(d_ghrelin$fID)
d_ghrelin_long <-  d_ghrelin %>% tidyr::pivot_wider(
                names_from = c(fTimepoint), values_from = c(F_AG, F_DG, logF_AG, logF_DG, res_logF_AG, res_logF_DG))

# Merge Condition file with IMT data 
d_IMT <- merge(d_IMT, d_blood, by = c("ID", "Session"))

d_IMT$fID <- factor(d_IMT$ID)
d_IMT$fSession <- factor(d_IMT$Session)
d_IMT$fGhrelin <- factor(d_IMT$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_IMT$cSession <- d_IMT$Session - mean(d_IMT$Session)

d_IMT$Reward_Money  <- factor(d_IMT$Reward_Money, labels = c("Food","Money"))
d_IMT$Reward_Magnitude <- factor(d_IMT$Reward_Magnitude, labels = c("Low","High"))
d_IMT$catTrial   <- cut(d_IMT$Trial, breaks = c(-Inf, 36, Inf), labels = c("Early", "Late"))

# Aggretate data
dAgg <- 
  d_IMT %>% 
  group_by(ID,fGhrelin,Reward_Magnitude,Reward_Money, cSession) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel), M_AbsForce = mean(Force_Abs), M_S_InvSlope = mean(S_InvSlope), M_S_Slope = mean(S_Slope)) %>%
  mutate(M_RelForce,M_AbsForce,M_S_InvSlope,M_S_Slope )

# Check current Sample Size 
length(unique(d_IMT$fID))

# Specify different effort phases 
dINV <- filter(d_IMT_trial, Time_trial <= 1)
dINV <- merge(dINV, d_blood, by = c("ID", "Session"))
head(dINV)
dINV$fID <- factor(dINV$ID)
dINV$fSession <- factor(dINV$Session)
dINV$fGhrelin <- factor(dINV$Ghrelin, labels = c("Placebo", "Ghrelin"))
dINV$cSession <- dINV$Session - mean(dINV$Session)

dINV$Reward_Money  <- factor(dINV$Reward_Money, labels = c("Food","Money"))
dINV$Reward_Magnitude <- factor(dINV$Reward_Magnitude, labels = c("Low","High"))
dINV$catTrial   <- cut(dINV$Trial, breaks = c(-Inf, 36, Inf), labels = c("Early", "Late"))
dINV$fSex <- factor(dINV$Sex_female, labels = c("male", "female"))

d_IMT_trial <- merge(d_IMT_trial, d_blood, by = c("ID", "Session"))
d_IMT_trial$fID <- factor(d_IMT_trial$ID)
d_IMT_trial$fSession <- factor(d_IMT_trial$Session)
d_IMT_trial$fGhrelin <- factor(d_IMT_trial$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_IMT_trial$cSession <- d_IMT_trial$Session - mean(d_IMT_trial$Session)

d_IMT_trial$Reward_Money  <- factor(d_IMT_trial$Reward_Money, labels = c("Food","Money"))
d_IMT_trial$Reward_Magnitude <- factor(d_IMT_trial$Reward_Magnitude, labels = c("Low","High"))
d_IMT_trial$catTrial   <- cut(d_IMT_trial$Trial, breaks = c(-Inf, 36, Inf), labels = c("Early", "Late"))
d_IMT_trial$fSex <- factor(d_IMT_trial$Sex_female, labels = c("male", "female"))
d_IMT_trial <- merge(d_IMT_trial, d_ghrelin_long, by = c("fID", "Session"))


d_AggINV <- dINV %>%
    group_by(ID,fGhrelin,Reward_Magnitude,Reward_Money) %>%
    summarize(MdRelForce = mean(Force_rel), maxdRelForce = max(Force_rel), maxTempDer1 = max(Rel_Dev1), meanTempDer1 = mean(Rel_Dev1))



# Prep PET Data 
d_BPs$fID <- factor(d_BPs$ID)
d_BPs$M_Caudate = (d_BPs$Caudate_r + d_BPs$Caudate_l) / 2
d_BPs$M_Putamen = (d_BPs$Putamen_r + d_BPs$Putamen_l) / 2
d_BPs$M_Accumbens = (d_BPs$Accumbens_r + d_BPs$Accumbens_l) / 2
d_BPs$cSession <- d_BPs$Session - mean(d_BPs$Session)

d_BPs <- merge(d_BPs, d_blood, by = c("ID", "Session"))
d_BPs <- merge(d_BPs, d_ghrelin_long, by = c("ID", "fID", "Session"))

d_BPs_Agg <- merge(d_BPs, dAgg, by = c("ID", "cSession"))


# Merge IMT comb and ghrelin 
Ghrelin_cond <- select(d_blood, ID,Session, fGhrelin)
d_IMT_comb2 <- merge(d_IMT_comb, Ghrelin_cond, by = c("ID", "Session"), all = TRUE)

# BOLD TASK DATA 
d_Task_fMRI <- read_excel(paste(path_in,"betas_extractedROIs.xlsx", sep= ""))
d_Task_fMRI <- merge(d_Task_fMRI, Ghrelin_cond, by = c("ID", "Session"), all = FALSE)

# BOLD Pulses and Matches 
d_pulseDelta <- read_csv(paste(path_out, "Pulse/Deltas_Pulse025.csv", sep= ""), show_col_types = FALSE)
d_pulseMatches <- read_csv(paste(path_out,"Pulse/Deltas_Matches025.csv", sep= ""), show_col_types = FALSE)



#####################################################################
# (3) DESCRIPTIVES TABLE 
#####################################################################

# Copy for descriptives table 
descriptives <- d_blood 
descriptives$fID <- factor(descriptives$ID)
head(descriptives)
descriptives <- 
  descriptives %>% 
  group_by(fID) %>%
  dplyr::summarize(BMI = mean(BMI), Age = mean(Age), T0_Glk = mean(T0_Glk),T0_Ins_p = mean(T0_Ins_p), 
  T0_HOMA_IR = mean(T0_HOMA_IR), T0_TG = mean(T0_TG), 
  # T0_Est = mean(as.numeric(T0_Est), na.rm = TRUE), 
  T0_TyG = mean(T0_TyG),Weight = mean(Weight), Sex = mean(Sex_female),
   Blood_pressure_systol.mmHg._T0 = mean(Blood_pressure_systol.mmHg._T0),Blood_pressure_diastol.mmHg._T0 = mean(Blood_pressure_diastol.mmHg._T0))
descriptives$fSex<- factor(descriptives$Sex, labels = c("male", "female"))

# Create Table with basic sample information for Paper 

  # Now create Table labels 
    table1::label(descriptives$BMI) <- "BMI [kg/m2]"
    table1::label(descriptives$fSex) <- "Sex"
    table1::label(descriptives$Age) <- "Age [years]"
   # table1::label(descriptives$F_AG) <- "Acyl ghrelin [pg/mol])"
   # table1::label(descriptives$F_DG) <- "Des-acyl ghrelin [pg/mol]"
    table1::label(descriptives$T0_Glk) <- "Glucose [mg/dl]"
    #table1::label(descriptives$T0_Est) <- "Estrogen [pmol/l]"
    #table1::label(descriptives$T0_Prg) <- "Progesteron [nmol/l]"
    #table1::label(descriptives$T0_Tst) <- "Testosteron [nmol/l]"
    table1::label(descriptives$T0_Ins_p) <- "Insulin [pmol/l]"
    table1::label(descriptives$T0_HOMA_IR) <- "HOMA IR"
    table1::label(descriptives$T0_TG) <- "Triglycerides [pmol/l]"
    table1::label(descriptives$T0_TyG) <- "Triglyceride Index"
    table1::label(descriptives$Weight) <- "Weight [kg]"
    table1::label(descriptives$Blood_pressure_systol.mmHg._T0) <- "Systolic blood\npressure (baseline)"
    table1::label(descriptives$Blood_pressure_diastol.mmHg._T0) <- "Diastolic blood\npressure (baseline)"
    #table1::label(descriptives$Blood_pressure_systol.mmHg._T1) <- "Systolic blood presure (post)"
    #table1::label(descriptives$Blood_pressure_diastol.mmHg._T1) <- "Diastolic blood presure (post)"

    head(descriptives)
    # Set render such that 2 places after comma 
    my.render.cont <- function(x) {
        with(table1::stats.apply.rounding(table1::stats.default(x), digits=4), 
        {MEAN <- as.numeric(MEAN)
        SD <- as.numeric(SD) 
        c("","Mean (SD)"=sprintf("%0.2f (&plusmn; %0.2f)", MEAN, SD))})
    }


   # Now create the Table, specify Rows, Split by MDD Group status 
    descrp_table <- table1::table1(~   Age + fSex + BMI + Weight + T0_Glk + T0_Ins_p + T0_TG  + T0_HOMA_IR  + T0_TyG  + Blood_pressure_systol.mmHg._T0 + 
    Blood_pressure_diastol.mmHg._T0 | fSex, data = descriptives, 
                render.continuous = my.render.cont)

    # Add stats to descriptives tabl 
    gtsummary::reset_gtsummary_theme()
    gtsummary::theme_gtsummary_mean_sd()
    gtsummary::theme_gtsummary_journal(journal = "jama")
    #> Setting theme `JAMA`

    table_stats <-
    descriptives %>%
        select(fSex, Age, BMI,  T0_Glk, T0_Ins_p , T0_TG , T0_HOMA_IR , T0_TyG  ) %>%
    gtsummary::tbl_summary(by = fSex,  statistic = list(all_continuous() ~ "{mean} (±{sd})"), digits = all_continuous() ~ 1) %>%
     #add_p(test = c(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test.no.correct")) %>% 
    gtsummary::add_overall() %>% 
    gtsummary::bold_labels() %>% 
    #add_stat_label() %>%
    # add a header to the statistic column, which is hidden by default
    # adding the header will also unhide the column
    #modify_header(statistic ~ "**Test Statistic**") %>%
    #modify_fmt_fun(statistic ~ style_sigfig) 
    gtsummary::as_flex_table()%>%
    flextable::save_as_docx( path = "Descriptives.docx")


###############################################################################
# 4. IMT PLOTS 
###############################################################################
p1 <- 
  ggplot(aes(y = Force_rel,x = Reward_Magnitude,fill = fGhrelin),data = d_IMT) +
  geom_bar(fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary',position = position_dodge(width = 0.9), 
           alpha = 0.9) +
  stat_dots(aes(x = Reward_Magnitude,y = M_RelForce,fill = fGhrelin),data = dAgg, position = position_dodge(width = 0.9), 
            side = 'both', alpha = 0.5) +
  scale_fill_manual(guide = guide_legend(title="Reward magnitude"),values = c("coral1","coral4")) +
  #coord_cartesian(ylim = c(30,95)) +
    facet_grid(. ~  Reward_Money) +
  geom_errorbar(aes(group = fGhrelin),data = d_IMT,
                linewidth = 1.0,width = 0.25,fun.data = mean_cl_normal,fun.args = list(conf.int = 0.95),stat = 'summary',position = position_dodge(width = 0.9)) +
  theme(legend.position = 'bottom', text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  xlab(label = 'Reward_Magnitude') +
  ylab(label = 'Relative Force [%] ')

ggsave(paste(path_out, sub_IMT, "IMT_General_Conds.png", sep=""),  
        plot = p1,  height = 10, width = 10, units = "in", dpi = 600, bg = "white")

# Increase relative effort with Ghrelin per Condition 
p2a <- 
  ggplot(aes(x = fGhrelin,y = Force_rel), data = d_IMT) +
  geom_smooth(aes(group=ID, color=ID), size = 2, method = 'rlm', alpha = 0) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) +
            theme(legend.position = 'bottom', text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0), 
        panel.spacing = unit(-10,'lines'),
        strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm"))) +
    facet_grid(. ~  Reward_Money) +
  xlab(label = 'Condition') +
  ylab(label = 'Relative force [%] ') + scale_x_discrete(labels=c("P","G"))
  
 
p2b  <- 
  ggplot(aes(x = fGhrelin,y = Force_rel), data = d_IMT) +
  geom_smooth(aes(group=ID, color=ID), size = 2, method = 'rlm', alpha = 0)+
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) +
   theme(text = element_text(face = 'bold',size = 16.0),axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm"))) +
  xlab(label = 'Condition') +
  ylab(label = 'Relative force [%] ') 


##############################################################
# PLOT EFFORT OVER TRIAL TIME 
# Plot different phases (i.e., invigoration and maintenance)
##############################################################

pDer1 <- 
  ggplot(aes(x = Time_trial,y = Rel_Dev1),data = dINV) +
  geom_smooth(aes(group = fID),data = dINV,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  facet_grid(. ~ Reward_Magnitude) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,1.3), ylim = c(-2,6.5)) +
  theme(legend.position = "bottom",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.3) +
  xlab(label = 'Time [s]') +
  ylab(label = 'Temporal derivative of effort') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceDerivative1_RewMag.png", sep=""),  
        plot = pDer1,  height = 6.5, width = 8, units = "in", dpi = 600, bg = "white")

pDer2 <- 
  ggplot(aes(x = Time_trial,y = Rel_Dev1),data = dINV) +
  geom_smooth(aes(group = fID),data = dINV,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,1.3), ylim = c(-2,6.5)) +
      facet_grid(. ~ Reward_Money) +
theme(legend.position = "bottom",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.3) +
  xlab(label = 'Time [s]') +
  ylab(label = 'Temporal derivative of effort') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceDerivative1_FM.png", sep=""),  
        plot = pDer2,  height = 6.5, width = 8, units = "in", dpi = 600, bg = "white")

d_IMT_trial$Rew_Money
d_IMT_trial$Reward_Money  <- factor(d_IMT_trial$Reward_Money, labels = c("Food","Money"))

pForceTrial_FM <- 
  ggplot(aes(x = Time_trial,y = Force_rel),data = d_IMT_trial) +
  geom_smooth(aes(group = fID),data = d_IMT_trial,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,3.3), ylim = c(-2,100)) +
    facet_grid(. ~ Reward_Money ) +
  theme(legend.position = "bottom",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.3) +
  xlab(label = 'Time [s]') +
  ylab(label = 'Relative Force') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_FoodMoney.png", sep=""),  
        plot = pForceTrial_FM,  height = 5, width = 5, units = "in", dpi = 600, bg = "white")

pForceTrial_RM <- 
  ggplot(aes(x = Time_trial,y = Force_rel),data = d_IMT_trial) +
  geom_smooth(aes(group = fID),data = d_IMT_trial,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,3.3), ylim = c(-2,100)) +
    facet_grid(. ~  Reward_Magnitude) +
  theme(legend.position = "bottom",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.3) +
  xlab(label = 'Time [s]') +
  ylab(label = 'Relative Force') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_Mag.png", sep=""),  
        plot = pForceTrial_RM,  height = 5, width = 5, units = "in", dpi = 600, bg = "white")

ForceTrial_Conditions <- plot_grid(pForceTrial_FM, pForceTrial_RM, labels = "none", label_size = 12, ncol=2)
ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_Panel.png", sep=""),  
        plot = ForceTrial_Conditions,  height = 7, width = 10, units = "in", dpi = 600, bg = "white")

# Check for the Sex effects with Ghrelin 

pForceTrial_Sex <- 
  ggplot(aes(x = Time_trial,y = Force_rel),data = d_IMT_trial) +
  geom_smooth(aes(group = fID),data = d_IMT_trial,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,3.3), ylim = c(-5,110)) +
    facet_grid(. ~  fSex ) +
  theme(legend.position = "bottom",text = element_text(face = 'bold',size = 16.0),
        axis.text = element_text(face = 'plain',size = 16.0),axis.text.x = element_text(size = 16.0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0.3,linetype = "dashed", color = "gray") +
  xlab(label = 'Time [s]') +
  ylab(label = 'Relative Force') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_GhrelinSex.png", sep=""),  
        plot = pForceTrial_Sex,  height = 4, width = 4, units = "in", dpi = 600, bg = "white")

pDer2 <- 
  ggplot(aes(x = Time_trial,y = Rel_Dev1),data = dINV) +
  geom_smooth(aes(group = fID),data = dINV,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,1.3), ylim = c(-2,6.5)) +
      facet_grid(. ~ fSex) +
theme(legend.position = "bottom",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
  geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.3) +
  xlab(label = 'Time [s]') +
  ylab(label = 'Temporal derivative of effort') 

  ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_GhrelinSex_Der.png", sep=""),  
        plot = pDer2,  height = 6.5, width = 8, units = "in", dpi = 600, bg = "white")

###############################################################################
# 5.  IMT STATISTICS 
###############################################################################



# Prep Models 
d_IMT$cForce_rel <- d_IMT$Force_rel - mean(d_IMT$Force_rel) 
d_IMT$cTrial <- d_IMT$Trial - mean(d_IMT$Trial)

contrasts(d_IMT$Reward_Money) <-  contr.treatment(levels(d_IMT$Reward_Money), base = 2)
contrasts(d_IMT$Reward_Magnitude) <-  contr.treatment(levels(d_IMT$Reward_Magnitude), base = 1)

## Model Comparisons: Does Ghrelin add Information? Does ghrelin increase Motivation to work for rewards? 
## Important: for model comparison use REML FALSE 

# Null Model: Without Ghrelin 
fm_0 <- lmer(Force_rel ~ Reward_Magnitude + Reward_Money + cTrial + cSession + cAge + cSex + cBMI + 
            (1  + cTrial + Reward_Magnitude + Reward_Money |fID), d_IMT, REML = FALSE)
summary(fm_0)

# Add Ghrelin Model 
fm_1a <- lmer(Force_rel ~  fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial + fGhrelin  + (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_1a)

# Add Ghrelin Model, and more complex random effect structure 
fm_2a <- lmer(Force_rel ~  fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_2a)

# Model: 3-way Interaction
fm_3a <- lmer(Force_rel ~  fGhrelin * Reward_Magnitude  * Reward_Money + cTrial + cSession + cAge + cSex + cBMI + 
            (1 + cTrial  + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_3a)

# Formel Test of Model Fit

# Test Deviance
# The output shows χ2 statistics representing the difference in deviance between successive models,
# as well as p values based on likelihood ratio test comparisons. 
anova(fm_0, fm_1a, fm_2a, fm_3a)

# Test Model Godness of Fit Criterion 
AIC <- AIC(fm_0, fm_1a , fm_2a, fm_3a)
BIC <- BIC(fm_0, fm_1a , fm_2a,fm_3a)
dev <- c(deviance(fm_0) ,deviance(fm_1a) , deviance(fm_2a), deviance (fm_3a))

# Give Column with Model IDs 
data1 <- cbind(Model = rownames(AIC), AIC)
rownames(data1) <- 1:nrow(data1)
data2 <- cbind(Model = rownames(BIC), BIC)
rownames(data2) <- 1:nrow(data2)

# Calculate Different to Null Model 
AIC_Delta <- data1  %>% mutate(AIC_Delta =  AIC- AIC[Model == "fm_0"])
BIC_Delta <- data2 %>% mutate(BIC_Delta =  BIC- BIC[Model == "fm_0"])
dev_Delta <- dev - dev[1]
dev_Delta2 <- dev - dev[2]

Models <- merge(AIC_Delta,BIC_Delta, by = "Model") %>% select(c(Model, AIC_Delta, BIC_Delta)) 
Models <-  Models %>% tidyr::pivot_longer(cols = "AIC_Delta":"BIC_Delta",  names_to = "Type", values_to = "value")
Models <- subset(Models, Model != "fm_0")

ModelCompP <- ggplot(Models, aes(x = Model,
                  y = Type,
                  fill = value))+geom_tile() + geom_text(aes(label=round(value)),  color = "white", size = 5) + 
                  theme(legend.position = "right",text = element_text(face = 'bold',size = 16.0),
                      axis.text = element_text(face = 'bold',size = 14),
                      axis.text.x = element_text(size = 14, angle = 65, vjust = 1, hjust=1)) +
                  scale_fill_viridis(discrete=FALSE, direction=-1) + 
                  scale_x_discrete(labels = c('ADD Ghrelin' , 'ADD Random\nEffects',"ADD 3-way\nInteraction"))

## WINNING MODEL: 2a! 
# Create Html Output to Word Document with Results 
sjPlot:: tab_model(fm_0, fm_1a, fm_2a, fm_3a,
                   p.val = "satterthwaite",
                   show.re.var=TRUE,
                  dv.labels = c("Base", "ADD Ghrelin", "ADD Random Effects", "3- way Interaction"), 
                  file= paste(path_out, sub_IMT, "IMT_ModelComp", ".doc", sep = ""))

sjPlot:: tab_model( fm_2a,
                   p.val = "satterthwaite",
                   show.re.var=TRUE,
                  dv.labels = c("IMT Winning Model"), 
                  file= paste(path_out, sub_IMT, "IMT_WinningModel", ".doc", sep = ""))

# Plot Winning Model Estimates of Interest

# Extract Random Effects (fixed + random) to plot
re <- coef(fm_2a)

Random_slopes <- re$fID %>%
  tibble::rownames_to_column(var = "fID") %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(c("fGhrelinGhrelin:Reward_MagnitudeHigh",
                           "fGhrelinGhrelin:Reward_MoneyFood")),
    names_to = "Coefficient",
    values_to = "Random_Effect"
  )

# Plot with automatic error bars
slopes_ghrelin_Model <- ggplot(Random_slopes, aes(x = Random_Effect , y = Coefficient)) +
    stat_summary(fun = mean, geom = "bar", fill = "darkgrey", alpha = 0.7, width = 0.6) +  # Mean bars
  geom_jitter(height = .25, size = 4, alpha = 0.6, color = "black") +  # Individual dots
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = "black", width = 0.2, size = 1.5) + 
  geom_vline(xintercept = 0, size = 1, color = "grey") +  # Horizontal line at y = 0
  labs(y = "", x = "Slopes Ghrelin x Reward") +
  theme(legend.position = "none",text = element_text(face = 'bold',size = 18),
                      axis.text = element_text(face = 'bold',size = 16),
                      axis.text.x = element_text(size = 16))+
  scale_y_discrete(labels=c("Magnitude","Type")) 


ggsave(paste(path_out, sub_IMT, "IMT_Slopes_Ghrelin.png", sep=""),  
        plot = slopes_ghrelin_Model,  height = 5, width = 7, units = "in", dpi = 600, bg = "white")



# Get Random Slopes selection for further analysis 

Random_slopes_single <- re$fID %>%
  tibble::rownames_to_column(var = "fID") %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(c("Reward_MoneyFood",
                           "Reward_MagnitudeHigh")),
    names_to = "Coefficient",
    values_to = "Random_Effect"
  ) %>% select(fID, Random_Effect, Coefficient)

Random_slopes_singleGhrelin <- re$fID %>%
  tibble::rownames_to_column(var = "fID") %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(c("fGhrelinGhrelin:Reward_MagnitudeHigh",
                           "fGhrelinGhrelin:Reward_MoneyFood")),
    names_to = "Coefficient",
    values_to = "Random_Effect"
  ) %>% select(fID, Random_Effect, Coefficient)

write.csv(Random_slopes_singleGhrelin, paste(path_out, sub_IMT, "Random_slopes_singleGhrelin.csv", sep =""), row.names=FALSE)

########################################################################
# Now: Sensitivity analyses 

# (a) Sex differences [mostly male mice only used for ghrelin effects]
fm_2b <- lmer(Force_rel ~  fGhrelin* cSex + fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial +  fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_2b)

anova(fm_2a, fm_2b)


# (b) Ghrelin Belief [is it belief ghrelin rather than ghrelin that increases motivation?]
d_IMT$fGhrelin_Belief <- factor(d_IMT$Test_blinding_pp, labels = c("Placebo", "Ghrelin"))

fm_1c <- lmer(Force_rel ~  fGhrelin_Belief * Reward_Magnitude  + fGhrelin_Belief * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial + fGhrelin_Belief * (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_1c)

anova(fm_0, fm_1c)
anova(fm_0, fm_2a)

BIC <- BIC(fm_0, fm_1c , fm_2a)
data2 <- cbind(Model = rownames(BIC), BIC)
rownames(data2) <- 1:nrow(data2)
BIC_Delta <- data2 %>% mutate(BIC_Delta =  BIC- BIC[Model == "fm_0"])

# Check Belief AND Ghrelin 
fm_1d <- lmer(Force_rel ~  (fGhrelin_Belief + fGhrelin)* Reward_Magnitude  + (fGhrelin_Belief+fGhrelin) * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial + (fGhrelin + fGhrelin_Belief) * (Reward_Magnitude + Reward_Money) |fID), d_IMT, REML = FALSE)
summary(fm_1d) ## nearly unidentifable, not good model, inspect visually for S1 and S2 

# Check Plasma Ghrelin levels (pre-scan)
d_IMT_PAG <- merge(d_IMT, d_ghrelin_long, by = c("fID", "Session"))
fm_2a_re <- lmer(Force_rel ~    fGhrelin * (Reward_Magnitude + Reward_Money)  + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial   + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_PAG, REML = FALSE)
summary(fm_2a_re)

# Check Ghrelin blood levels 
fm_2a_plasma <- lmer(Force_rel ~  res_logF_AG_T2 + fGhrelin * (Reward_Magnitude + Reward_Money)  + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial +  res_logF_AG_T2  + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_PAG, REML = FALSE)
summary(fm_2a_plasma)

anova(fm_2a_re, fm_2a_plasma)

# Plot Dependency on Plasma Ghrelin levels 
d_IMT_trial$Med_Split_AG_T2 <- as.numeric(d_IMT_trial$res_logF_AG_T2 > median(na.omit(d_IMT_trial$res_logF_AG_T2)))
d_IMT_trial$Med_Split_AG_T2 <- factor(d_IMT_trial$Med_Split_AG_T2, labels = c("low AG (T2)","high AG (T2)"))


pForceTrial_PlasmaGhrelin <- 
  ggplot(aes(x = Time_trial,y = Force_rel),data = d_IMT_trial) +
  geom_smooth(aes(group = fID),data = d_IMT_trial,color = "darkgray",
              method="gam", formula = y ~ s(x, bs = "cs"), size = 0.2, alpha = 0.1) +
  geom_smooth(aes(color = fGhrelin),
              method="gam", formula = y ~ s(x, bs = "cs"), size = 2) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
  coord_cartesian(xlim = c(0,3.3), ylim = c(-5,110)) +
    facet_grid(. ~  Med_Split_AG_T2 + Reward_Money) +
  theme(legend.position = "bottom",text = element_text(face = 'bold',size = 16.0),
        axis.text = element_text(face = 'plain',size = 16.0),axis.text.x = element_text(size = 16.0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0.3,linetype = "dashed", color = "gray") +
  xlab(label = 'Time [s]') +
  ylab(label = 'Relative Force') 

ggsave(paste(path_out, sub_IMT, "IMT_ForceTime_PlasmaAG.png", sep=""),  
        plot = pForceTrial_PlasmaGhrelin,  height = 4, width = 4, units = "in", dpi = 600, bg = "white")



###############################################################################
# Sensitivity Checks
###############################################################################

# Test min max absolute force values 
d_IMT$MinEffort
d_IMT$MaxEffort
min(d_IMT$Force_Abs)

# Compute subject-level means per condition
effort_summary <- d_IMT %>%
  group_by(fID, fGhrelin) %>%
  mutate(MinEffort_real = min(Force_Abs), 
  MaxEffort_real=max(Force_Abs)) %>%
  summarise(
    MinEffort_mean = max(MinEffort_real, na.rm = TRUE),
    MaxEffort_mean = max(MaxEffort_real, na.rm = TRUE),
    .groups = "drop"
  )

head(effort_summary)

# Wide format for paired tests
effort_wide <- effort_summary %>%
  select(fID, fGhrelin, MinEffort_mean, MaxEffort_mean) %>%
  pivot_wider(
    names_from  = fGhrelin,
    values_from = c(MinEffort_mean, MaxEffort_mean)
  )

# Plot DABESTR Plot 

# MinEffort
Dabestr_MinEff <- effort_summary %>%
  dplyr::select(fID, fGhrelin, MinEffort_mean) %>%
  rename(
    Group      = fGhrelin,
    Measurement = MinEffort_mean
  )

# MaxEffort
Dabestr_MaxEff <- effort_summary %>%
  dplyr::select(fID, fGhrelin, MaxEffort_mean) %>%
  rename(
    Group      = fGhrelin,
    Measurement = MaxEffort_mean
  )

#--------------------------------------------------
# 3) dabest objects (paired, sequential)
#--------------------------------------------------

# MinEffort
dabest_min <- dabestr::load(
  Dabestr_MinEff,
  x      = Group,
  y      = Measurement,
  idx    = c("Placebo", "Ghrelin"),  # or c("Saline", "Ghrelin")
  paired = "sequential",
  id_col = fID
)

dabest_min_diff <- dabestr::mean_diff(dabest_min)
print(dabest_min_diff)

# MaxEffort
dabest_max <- dabestr::load(
  Dabestr_MaxEff,
  x      = Group,
  y      = Measurement,
  idx    = c("Placebo", "Ghrelin"),  # or c("Saline", "Ghrelin")
  paired = "sequential",
  id_col = fID
)

dabest_max_diff <- dabestr::mean_diff(dabest_max)
print(dabest_max_diff)

#--------------------------------------------------
# 4) Plots (one for Min, one for Max)
#--------------------------------------------------

# ---- MinEffort plot ----
dabest_min_plot <- dabestr::dabest_plot(
  dabest_min_diff,
  raw_marker_size   = 1,
  raw_marker_alpha  = 0.4,
  effsize_marker_size = 5
) +
  ggtitle("Absolute Min Effort") +
  labs(
    x = NULL,
    y = "MinEffort (Ghrelin - Placebo)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 18),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 20),
    axis.text.x   = element_text(size = 20),
    legend.position = "none"
  )

# ---- MaxEffort plot ----
dabest_max_plot <- dabestr::dabest_plot(
  dabest_max_diff,
  raw_marker_size   = 1,
  raw_marker_alpha  = 0.4,
  effsize_marker_size = 5
) +
  labs(
    x = NULL,
    y = "MaxEffort (Ghrelin - Placebo)"
  ) +
    ggtitle("Absolute Max Effort") +

  theme_minimal(base_size = 18) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 18),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 20),
    axis.text.x   = element_text(size = 20),
    legend.position = "none"
  )


combined_plot <- cowplot::plot_grid(
  dabest_min_plot,
  dabest_max_plot,
  labels = c("A", "B"),
  ncol = 2,
  align = "h",
  label_size = 18
)

ggsave(paste(path_out, sub_IMT, "IMT_MinMaxForce.png", sep=""),  
        plot = combined_plot,  height = 5, width = 10, units = "in", dpi = 600, bg = "white")

# Same for Sessions 

# Compute subject-level means per condition
effort_summary <- d_IMT %>%
  group_by(fID, Session) %>%
  mutate(MinEffort_real = min(Force_Abs), 
  MaxEffort_real=max(Force_Abs)) %>%
  summarise(
    MinEffort_mean = max(MinEffort_real, na.rm = TRUE),
    MaxEffort_mean = max(MaxEffort_real, na.rm = TRUE),
    .groups = "drop"
  )

head(effort_summary)

# Wide format for paired tests
effort_wide <- effort_summary %>%
  select(fID, Session, MinEffort_mean, MaxEffort_mean) %>%
  pivot_wider(
    names_from  = Session,
    values_from = c(MinEffort_mean, MaxEffort_mean)
  )

# Plot DABESTR Plot 

# MinEffort
Dabestr_MinEff <- effort_summary %>%
  dplyr::select(fID, Session, MinEffort_mean) %>%
  rename(
    Group      = Session,
    Measurement = MinEffort_mean
  )

# MaxEffort
Dabestr_MaxEff <- effort_summary %>%
  dplyr::select(fID, Session, MaxEffort_mean) %>%
  rename(
    Group      = Session,
    Measurement = MaxEffort_mean
  )

#--------------------------------------------------
# 3) dabest objects (paired, sequential)
#--------------------------------------------------

# MinEffort
dabest_min <- dabestr::load(
  Dabestr_MinEff,
  x      = Group,
  y      = Measurement,
  idx    = c(1,2),  # or c("Saline", "Ghrelin")
  paired = "sequential",
  id_col = fID
)

dabest_min_diff <- dabestr::mean_diff(dabest_min)
print(dabest_min_diff)

# MaxEffort
dabest_max <- dabestr::load(
  Dabestr_MaxEff,
  x      = Group,
  y      = Measurement,
  idx    = c(1,2),  # or c("Saline", "Ghrelin")
  paired = "sequential",
  id_col = fID
)

dabest_max_diff <- dabestr::mean_diff(dabest_max)
print(dabest_max_diff)

#--------------------------------------------------
# 4) Plots (one for Min, one for Max)
#--------------------------------------------------

# ---- MinEffort plot ----
dabest_min_plot <- dabestr::dabest_plot(
  dabest_min_diff,
  raw_marker_size   = 1,
  raw_marker_alpha  = 0.4,
  effsize_marker_size = 5
) +
  ggtitle("Absolute Min Effort") +
  labs(
    x = NULL,
    y = "MinEffort (S1 - S2)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 18),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 20),
    axis.text.x   = element_text(size = 20),
    legend.position = "none"
  )

# ---- MaxEffort plot ----
dabest_max_plot <- dabestr::dabest_plot(
  dabest_max_diff,
  raw_marker_size   = 1,
  raw_marker_alpha  = 0.4,
  effsize_marker_size = 5
) +
  labs(
    x = NULL,
    y = "MaxEffort (S1 - S2)"
  ) +
    ggtitle("Absolute Max Effort") +

  theme_minimal(base_size = 18) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 18),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.y   = element_text(size = 20),
    axis.text.x   = element_text(size = 20),
    legend.position = "none"
  )


combined_plot_Session <- cowplot::plot_grid(
  dabest_min_plot,
  dabest_max_plot,
  labels = c("A", "B"),
  ncol = 2,
  align = "h",
  label_size = 18
)

ggsave(paste(path_out, sub_IMT, "IMT_MinMaxForce_Session.png", sep=""),  
        plot = combined_plot_Session,  height = 5, width = 12, units = "in", dpi = 600, bg = "white")


###############################################################################
# 7. Neuroimaging Analysis (SPM Covariate for FC Model)
###############################################################################

# Get 1 estimate per participant (within) for relative effort as a reward senstivtiy index to use as covariate in SPM 
#  High-Low Rewards, Food-Money Average 
#  Food-Money Average, High-Low Rewards  


dAgg_RewIndex <- 
  d_IMT %>% 
  group_by(ID,fGhrelin,cSession, Reward_Magnitude,Reward_Money, cBMI, cAge, cSex) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel)) %>%
  # separate columns for Conditions in order to substract within subject 
  tidyr::pivot_wider(names_from = c(Reward_Money, Reward_Magnitude), values_from = M_RelForce) %>%
  mutate(M_RelForce_Delta4_HL_AvgFoodMoney = (mean(c(Food_High, Money_High)))- (mean(c(Food_Low, Money_Low))))  %>%
  mutate(M_RelForce_Delta5_DiffFoodMoney_AvgHL = mean(c( (Food_High - Money_High), (Food_Low - Money_Low))))

MeanForce <- 
  d_IMT %>% 
  group_by(ID,fGhrelin) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel))

# Write final csv file for further analysis
write.csv(MeanForce, paste(path_out, sub_IMT, "TUE008_IMT_MeanEffort.csv", sep =""), row.names=FALSE)
write.csv(dAgg_RewIndex, paste(path_out, sub_IMT, "TUE008_IMT_RewIndex.csv", sep =""), row.names=FALSE)


###############################################################################
# 9. PET SIMPLE PLOTS 
###############################################################################

# Simple PET Plots q/ Covariates 

d_BPs_Corr <-  d_BPs %>% 
    select(c("fID",  "fGhrelin", "cSex", "cAge", "M_Putamen","M_Caudate","M_Accumbens")) %>% 
     reshape(idvar= "fID",
                    v.names= c("M_Putamen","M_Caudate","M_Accumbens"),
                    timevar= "fGhrelin",
                    direction = "wide")


# PET Correspondence between Sessions 

plot_BP_correlation1 <- 
  ggplot(aes(x = M_Putamen.Placebo ,y = M_Putamen.Ghrelin ),data = d_BPs_Corr) +
  geom_point(size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP Putamen Placebo') +
  ylab(label = 'BP Putamen Ghrelin') +  ggtitle("") +
  smplot2::sm_statCorr() +
  geom_abline(slope=1, intercept = 0, linetype = "dashed")

plot_BP_correlation2 <- 
  ggplot(aes(x = M_Caudate.Placebo ,y =M_Caudate.Ghrelin ),data = d_BPs_Corr) +
  geom_point(size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP Caudate Placebo') +
  ylab(label = 'BP Caudate Ghrelin') +  ggtitle("") +
  smplot2::sm_statCorr()  +
  geom_abline(slope=1, intercept = 0, linetype = "dashed")

plot_BP_correlation3 <- 
  ggplot(aes(x = M_Accumbens.Placebo ,y = M_Accumbens.Ghrelin ),data = d_BPs_Corr) +
  geom_point(size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP Accumbens Placebo') +
  ylab(label = 'BP Accumbens Ghrelin') +  ggtitle("") +
  smplot2::sm_statCorr()  +
  geom_abline(slope=1, intercept = 0, linetype = "dashed")


BP_correlations <- plot_grid(plot_BP_correlation1, plot_BP_correlation2, plot_BP_correlation3, 
labels = "auto", label_size = 12, ncol=3)

ggsave( paste(path_out, sub_PET, "/basic/BP_correlation_GP.png", sep = ""), 
        plot = BP_correlations,  height = 5, width = 10, units = "in", dpi = 600, bg = "white")

# PET and Covariates (Sanity Checks) 

d_BPs_Agg_Wide <-  d_BPs_Agg %>% tidyr::pivot_longer(
                cols = c("M_Putamen", "M_Caudate", "M_Accumbens"), 
                names_to = "ROI",
                values_to = "BP")

d_BPs_Agg_Wide_Merge <- d_BPs_Agg_Wide %>% 
  group_by(fID, ROI, Age) %>%
  summarize(M_RelForce_Mean = mean(M_RelForce), M_BP = mean(BP), M_BMI = mean(BMI))

plot_PET_baseChar <- 
  ggplot(aes(x = Age ,y = M_BP ),data = d_BPs_Agg_Wide_Merge) +
  geom_point( alpha =1)+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "darkgray") +
  stat_cor(label.y = 3.3, label.x = 15) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  ylab(label = 'BP(ND)') +
  xlab(label = 'Age') +  ggtitle("") + 
  facet_wrap(~ROI, labeller = as_labeller(c(M_Accumbens = "Accumbens", M_Caudate = "Caudate", M_Putamen = "Putamen")))

ggsave( paste(path_out, sub_PET, "basic/BP_Age.png", sep = ""), 
        plot = plot_PET_baseChar,  height = 5, width = 6, units = "in", dpi = 600, bg = "white")


plot_PET_baseChar <- 
  ggplot(aes(x = M_BMI ,y = M_BP ),data = d_BPs_Agg_Wide_Merge) +
  geom_point( alpha =1)+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "darkgray") +
  stat_cor(label.y = 3.3, label.x = 15) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  ylab(label = 'BP(ND)') +
  xlab(label = 'BMI') +  ggtitle("") + 
  facet_wrap(~ROI, labeller = as_labeller(c(M_Accumbens = "Accumbens", M_Caudate = "Caudate", M_Putamen = "Putamen")))

ggsave( paste(path_out, sub_PET, "basic/BP_BMI.png", sep = ""), 
        plot = plot_PET_baseChar,  height = 5, width = 6, units = "in", dpi = 600, bg = "white")


###############################################################################
# 10. PET GHRELIN 
###############################################################################


# Make Boxplot for common visualization 
d_BPs_Box <-  d_BPs %>% pivot_longer(
                cols = c("M_Putamen", "M_Caudate", "M_Accumbens"), 
                names_to = "ROI",
                values_to = "BP")

Box_ghrelin <- 
  ggplot(d_BPs_Box, aes(x = ROI, y = BP, fill = fGhrelin)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(aes(color= fGhrelin), alpha = 0.5, 
                 position = position_jitterdodge(jitter.width = 0.3)) +
  scale_color_manual(guide = guide_legend(title="Condition"),values = c( color_Ghrelin, color_Placebo)) +
  scale_fill_manual(guide = "none",values = c(color_Ghrelin, color_Placebo)) +
   theme(text = element_text(face = 'bold',size = 16.0),axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm"))) +
  xlab(label = 'Region of Interest') +
  ylab(label = 'BP[ND] ') + 
  scale_x_discrete(labels=c("Accumbens", "Caudate","Putamen"))

ggsave( paste(path_out, sub_PET, "ghrelin/BP_Boxplot_Ghrelin.png", sep = ""), 
        plot = Box_ghrelin,  height = 5, width = 6, units = "in", dpi = 600, bg = "white")


###############################################################################
# 11. PET STATS 
###############################################################################

d_BPs$fGhrelin <- factor(d_BPs$fGhrelin, levels = c("Placebo", "Ghrelin"))  

# Check PET BP, Ghrelin, and Covariates 

# Single ROIS 
fm_PET_base_all0 <- lmer(M_Caudate ~  cSession  + fGhrelin +(1  |fID) , d_BPs)
summary(fm_PET_base_all0)

fm_PET_base_all1 <- lmer(M_Caudate ~  cSession  + fGhrelin + cSex  +(1 |fID) , d_BPs)
summary(fm_PET_base_all1)

fm_PET_base_all2 <- lmer(M_Caudate ~  cSession  + fGhrelin  + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all2)

fm_PET_base_all3 <- lmer(M_Caudate ~  cSession  + fGhrelin + cBMI  +(1 |fID), d_BPs)
summary(fm_PET_base_all3)

fm_PET_base_all4 <- lmer(M_Caudate ~  cSession  + fGhrelin + cSex + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all4)

fm_PET_base_all5 <- lmer(M_Caudate ~  cSession  + fGhrelin + cSex + cBMI   +(1 |fID), d_BPs)
summary(fm_PET_base_all5)

fm_PET_base_all6 <- lmer(M_Caudate ~  cSession  + fGhrelin + cBMI + cAge   +(1 |fID), d_BPs)
summary(fm_PET_base_all6)

fm_PET_base_all7 <- lmer(M_Caudate ~  cSession  + fGhrelin + cSex + cAge + cBMI +  +(1 |fID), d_BPs)
summary(fm_PET_base_all7)

# res_logHOMA_T0 

fm_PET_base_all7 <- lmer(M_Caudate ~  cSession  + fGhrelin * res_logHOMA_T0  + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all7)

anova(fm_PET_base_all0,fm_PET_base_all1)
anova(fm_PET_base_all0,fm_PET_base_all2)
anova(fm_PET_base_all0,fm_PET_base_all3)
anova(fm_PET_base_all0,fm_PET_base_all4)
anova(fm_PET_base_all0,fm_PET_base_all5)
anova(fm_PET_base_all0,fm_PET_base_all6)
anova(fm_PET_base_all0,fm_PET_base_all7)


fm_PET_base_all0 <- lmer(M_Putamen ~  cSession  + fGhrelin +(1  |fID) , d_BPs)
summary(fm_PET_base_all0)

fm_PET_base_all1 <- lmer(M_Putamen ~  cSession  + fGhrelin + cSex  +(1 |fID) , d_BPs)
summary(fm_PET_base_all1)

fm_PET_base_all2 <- lmer(M_Putamen ~  cSession  + fGhrelin  + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all2)

fm_PET_base_all3 <- lmer(M_Putamen ~  cSession  + fGhrelin + cBMI  +(1 |fID), d_BPs)
summary(fm_PET_base_all3)

fm_PET_base_all4 <- lmer(M_Putamen ~  cSession  + fGhrelin + cSex + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all4)

fm_PET_base_all5 <- lmer(M_Putamen ~  cSession  + fGhrelin + cSex + cBMI   +(1 |fID), d_BPs)
summary(fm_PET_base_all5)

fm_PET_base_all6 <- lmer(M_Putamen ~  cSession  + fGhrelin + cBMI + cAge   +(1 |fID), d_BPs)
summary(fm_PET_base_all6)

fm_PET_base_all7 <- lmer(M_Putamen ~  cSession  + fGhrelin + cSex + cAge + cBMI +  +(1 |fID), d_BPs)
summary(fm_PET_base_all7)

# res_logHOMA_T0 

anova(fm_PET_base_all0,fm_PET_base_all1)
anova(fm_PET_base_all0,fm_PET_base_all2)
anova(fm_PET_base_all0,fm_PET_base_all3)
anova(fm_PET_base_all0,fm_PET_base_all4)
anova(fm_PET_base_all0,fm_PET_base_all5)
anova(fm_PET_base_all0,fm_PET_base_all6)
anova(fm_PET_base_all0,fm_PET_base_all7)

fm_PET_base_all_W1 <- lmer(M_Putamen ~  cSession  + fGhrelin * cBMI + cAge + cBMI + cAge +  +(1 |fID), d_BPs)
summary(fm_PET_base_all_W1)

fm_PET_base_all_W1 <- lmer(M_Putamen ~  cSession  + fGhrelin * res_logHOMA_T0 + cBMI + cAge + cBMI + cAge +  +(1 |fID), d_BPs)
summary(fm_PET_base_all_W1)

## 
fm_PET_base_all0 <- lmer(M_Accumbens ~  cSession  + fGhrelin +(1  |fID) , d_BPs)
summary(fm_PET_base_all0)

fm_PET_base_all1 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cSex  +(1 |fID) , d_BPs)
summary(fm_PET_base_all1)

fm_PET_base_all2 <- lmer(M_Accumbens ~  cSession  + fGhrelin  + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all2)

fm_PET_base_all3 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cBMI  +(1 |fID), d_BPs)
summary(fm_PET_base_all3)

fm_PET_base_all4 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cSex + cAge  +(1 |fID), d_BPs)
summary(fm_PET_base_all4)

fm_PET_base_all5 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cSex + cBMI   +(1 |fID), d_BPs)
summary(fm_PET_base_all5)

fm_PET_base_all6 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cBMI + cAge   +(1 |fID), d_BPs)
summary(fm_PET_base_all6)

fm_PET_base_all7 <- lmer(M_Accumbens ~  cSession  + fGhrelin + cSex + cAge + cBMI +  +(1 |fID), d_BPs)
summary(fm_PET_base_all7)

# res_logHOMA_T0 

anova(fm_PET_base_all0,fm_PET_base_all1)
anova(fm_PET_base_all0,fm_PET_base_all2)
anova(fm_PET_base_all0,fm_PET_base_all3)
anova(fm_PET_base_all0,fm_PET_base_all4)
anova(fm_PET_base_all0,fm_PET_base_all5)
anova(fm_PET_base_all0,fm_PET_base_all6)
anova(fm_PET_base_all0,fm_PET_base_all7) 

##############################################################################
# 12. Spike + IMT  
##############################################################################

d_pulseDelta_Filtered <- d_pulseDelta %>% filter(fBaseCond == "TaskFree")

d_IMT_spikes <- merge(d_IMT, d_pulseDelta_Filtered, by = c("fID", "fGhrelin"))
d_IMT_spikes$cPulse_Delta <- d_IMT_spikes$Pulse_Delta - mean(d_IMT_spikes$Pulse_Delta)

head(d_IMT_spikes)

# (0.) Redo Force Model with Ghrelin Effect aka WINNING MODEL 
fm_2b_pulse_base <- lmer(Force_rel ~   fGhrelin * (Reward_Magnitude + Reward_Money) + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial +  fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikes, REML = FALSE)
summary(fm_2b_pulse_base)

# (1). Add PULSE DELTA + Interaction with Ghrelin
fm_2b_pulse <- lmer(Force_rel ~  fGhrelin * cPulse_Delta + fGhrelin * (Reward_Magnitude + Reward_Money) + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial  + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikes, REML = FALSE)
summary(fm_2b_pulse)

# (2). Add PULSE DELTA  + Interaction with Reward Types/Mag
fm_2b_pulse2 <- lmer(Force_rel ~   cPulse_Delta * (Reward_Magnitude + Reward_Money) + fGhrelin * (Reward_Magnitude + Reward_Money) + cTrial + cSession + cBMI + cSex + cAge  + 
            (1 +  cTrial + fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikes, REML = FALSE)
summary(fm_2b_pulse2)

# Model comparison
anova(fm_2b_pulse_base, fm_2b_pulse)
anova(fm_2b_pulse_base, fm_2b_pulse2)

# Create Html Output to Word Document with Results 
sjPlot:: tab_model(fm_2b_pulse2,
                   p.val = "satterthwaite",
                   show.re.var=TRUE,
                   show.se= TRUE, 
                  dv.labels = c("Effort - Pulses"), 
                  file= paste(path_out, sub_Pulse, "Model_IMT_Pulse.doc", sep = ""))


color_other <- "gray40"
color_food <- "#2E8B57"   
color_low     <- "#B22222"  

# Food vs. Money interaction (Reward_Money) at representative Pulse levels
# Build grid of EMMs for Pulse × RewardType × Ghrelin
dat <- d_IMT_spikes
pulse_seq <- seq(min(dat$cPulse_Delta, na.rm = TRUE),
                 max(dat$cPulse_Delta, na.rm = TRUE),
                 length.out = 40)

emm_effort_food <- emmeans(
  fm_2b_pulse2,
  specs = ~ cPulse_Delta * Reward_Money ,
  at = list(cPulse_Delta = pulse_seq),
  cov.reduce = mean,
  type = "response"
) %>% as.data.frame()

# Plot with ribbons
p_effort_food <- ggplot(emm_effort_food,
                        aes(x = cPulse_Delta, y = emmean,
                            color = Reward_Money, fill = Reward_Money)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Food" = color_food, "Money" = color_other)) +
  scale_fill_manual(values  = c("Food" = color_food, "Money" = color_other)) +
  labs(x = "Δ Pulse count (TaskFree, NAcc)",
       y = "Relative effort",
       color = "Reward type", fill = "Reward type") +
        theme(legend.position = "none",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 

emm_effort_RewMag <- emmeans(
  fm_2b_pulse2,
  specs = ~ cPulse_Delta * Reward_Magnitude ,
  at = list(cPulse_Delta = pulse_seq),
  cov.reduce = mean,
  type = "response"
) %>% as.data.frame()

# Plot with ribbons
p_effort_magntiude <- ggplot(emm_effort_RewMag,
                        aes(x = cPulse_Delta, y = emmean,
                            color = Reward_Magnitude, fill = Reward_Magnitude)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("High" = color_other, "Low" = color_low)) +
  scale_fill_manual(values  = c("High" = color_other, "Low" = color_low)) +
  labs(x = "Δ Pulse count (TaskFree, NAcc)",
       y = "Relative effort",
       color = "Reward magnitude", fill = "Reward magnitude") + 
               theme(legend.position = "none",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


Sens_Pulses <-  cowplot::plot_grid( p_effort_food, p_effort_magntiude,
          labels=c("", ""),
          ncol = 1,
          nrow = 2,
          rel_widths= c(1),
          label_size = 18)

ggsave( paste(path_out, sub_Pulse, "IMT_Pulses_Sensitivities.png", sep = ""), 
         height = 4, width = 4, units = "in", dpi = 600, bg = "white")


##############################################################################
# 12. Spike Matches + IMT  
##############################################################################

# Important: Only use TaskFree ( or Task, not both!)
d_pulseMatches_Filtered <- d_pulseMatches %>% filter(fBaseCond == "TaskFree")

d_IMT_spikesMatches <- merge(d_IMT, d_pulseMatches_Filtered, by = c("fID", "fGhrelin"))
head(d_IMT_spikesMatches)
d_IMT_spikesMatches$cMatches_Delta <- d_IMT_spikesMatches$Matches_Delta - mean(d_IMT_spikesMatches$Matches_Delta)

fm_2b_matches_base <- lmer(Force_rel ~   fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge + 
            (1 +  cTrial +  fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikesMatches, REML = FALSE)
summary(fm_2b_matches_base)

fm_2b_matches <- lmer(Force_rel ~  fGhrelin* cMatches_Delta + fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge +  
            (1 +  cTrial +  fGhrelin * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikesMatches, REML = FALSE)
summary(fm_2b_matches)

fm_2b_matches2 <- lmer(Force_rel ~  cMatches_Delta * Reward_Magnitude + cMatches_Delta * Reward_Money + fGhrelin * Reward_Magnitude + fGhrelin * Reward_Money + cTrial + cSession + cBMI + cSex + cAge +  
            (1 +  cTrial + (fGhrelin) * (Reward_Magnitude + Reward_Money) |fID), d_IMT_spikesMatches, REML = FALSE)
summary(fm_2b_matches2)


# Model comparison
anova(fm_2b_matches_base, fm_2b_matches )
anova(fm_2b_matches_base, fm_2b_matches2 )

# Create Html Output to Word Document with Results 
sjPlot:: tab_model(fm_2b_matches2,
                   p.val = "satterthwaite",
                   show.re.var=TRUE,
                    show.se= TRUE, 
                  dv.labels = c("Effort - Matches"), 
                  file= paste(path_out, sub_Pulse, "Model_IMT_Matches.doc", sep = ""))


# Food vs. Money interaction (Reward_Money) at representative Pulse levels
# Build grid of EMMs for Pulse × RewardType × Ghrelin
dat <- d_IMT_spikesMatches
pulse_seq <- seq(min(dat$cMatches_Delta, na.rm = TRUE),
                 max(dat$cMatches_Delta, na.rm = TRUE),
                 length.out = 40)

emm_effort_food <- emmeans(
  fm_2b_matches2,
  specs = ~ cMatches_Delta * Reward_Money ,
  at = list(cMatches_Delta = pulse_seq),
  cov.reduce = mean,
  type = "response"
) %>% as.data.frame()

# Plot with ribbons
p_effort_food_M <- ggplot(emm_effort_food,
                        aes(x = cMatches_Delta, y = emmean,
                            color = Reward_Money, fill = Reward_Money)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Food" = color_food, "Money" = color_other)) +
  scale_fill_manual(values  = c("Food" = color_food, "Money" = color_other)) +
        guides(colour=guide_legend(ncol=1,nrow=2)) + 
  labs(x = "Δ Matches (TaskFree, Hypo-NAcc)",
       y = "Relative effort",
       color = "Reward\ntype", fill = "Reward\ntype") +
        theme(legend.position = "right",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 

emm_effort_RewMag <- emmeans(
  fm_2b_matches2,
  specs = ~ cMatches_Delta * Reward_Magnitude ,
  at = list(cMatches_Delta = pulse_seq),
  cov.reduce = mean,
  type = "response"
) %>% as.data.frame()

# Plot with ribbons
p_effort_magntiude_M <- ggplot(emm_effort_RewMag,
                        aes(x = cMatches_Delta, y = emmean,
                            color = Reward_Magnitude, fill = Reward_Magnitude)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("High" = color_other, "Low" = color_low)) +
  scale_fill_manual(values  = c("High" = color_other, "Low" = color_low)) +
    guides(colour=guide_legend(ncol=1,nrow=2)) + 
  labs(x = "Δ Matches (TaskFree, Hypo-NAcc)",
       y = "Relative effort",
       color = "Reward\nmagnitude", fill = "Reward\nmagnitude") + 
               theme(legend.position = "right",text = element_text(face = 'bold',size = 14.0),
        axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


Sens_Matches <- cowplot::plot_grid( p_effort_food_M, p_effort_magntiude_M,
          labels=c("", ""),
          ncol = 1,
          nrow = 2,
          rel_widths= c(1),
          label_size = 18)

ggsave( paste(path_out, sub_Pulse, "IMT_Matches_Sensitivities.png", sep = ""), 
         height = 4, width = 5, units = "in", dpi = 600, bg = "white")

cowplot::plot_grid( Sens_Pulses, Sens_Matches,
          labels=c("", ""),
          ncol = 2,
          nrow = 1,
          rel_widths= c(1,1.3),
          label_size = 18)

ggsave( paste(path_out, sub_Pulse, "IMT_ALL_Sensitivities.png", sep = ""), 
         height = 4, width = 10, units = "in", dpi = 600, bg = "white")



###################################################################
# STATS for PET and IMT General Effect 
###################################################################

d_BPs_Agg_Wide_Caudate <- d_BPs_Agg_Wide %>%  filter(ROI == "M_Caudate")
d_BPs_Agg_Wide_NAC <- d_BPs_Agg_Wide %>%  filter(ROI == "M_Accumbens")
d_BPs_Agg_Wide_Putamen <- d_BPs_Agg_Wide %>%  filter(ROI == "M_Putamen")

d_BPs_Agg_Wide_Caudate$cBP <- d_BPs_Agg_Wide_Caudate$BP - mean(d_BPs_Agg_Wide_Caudate$BP)
d_BPs_Agg_Wide_NAC$cBP <- d_BPs_Agg_Wide_NAC$BP - mean(d_BPs_Agg_Wide_NAC$BP)
d_BPs_Agg_Wide_Putamen$cBP <- d_BPs_Agg_Wide_Putamen$BP - mean(d_BPs_Agg_Wide_Putamen$BP)

# Test DA and REWARD SENSITIVITY 
# This is simplified to direct diff in conditions (no fit issues)

head(d_BPs)
d_BPs_Sensitivity <- merge(d_BPs, dAgg_RewIndex, by = c("ID", "fGhrelin"))
d_BPs_Sensitivity$M_RelForce_Delta5_DiffFoodMoney_AvgHL # Food Sensitivity 
d_BPs_Sensitivity$M_RelForce_Delta4_HL_AvgFoodMoney # Reward Sensitivity 

d_BPs_Agg_merge <- d_BPs_Sensitivity %>% 
  group_by(fID, fGhrelin, cSession.x, cSex.x, cAge.x, cBMI.x) %>%
  summarize(FoodSen_Mean = mean(M_RelForce_Delta5_DiffFoodMoney_AvgHL),
  RewSens_Mean = mean(M_RelForce_Delta4_HL_AvgFoodMoney) ,M_Putamen_Mean = mean(M_Putamen), 
  M_Caudate_Mean = mean(M_Caudate),M_Accumbens_Mean = mean(M_Accumbens)) %>% 
                tidyr::pivot_longer(
                cols = c("M_Putamen_Mean", "M_Caudate_Mean", "M_Accumbens_Mean"), 
                names_to = "ROI",
                values_to = "BP")


plot_PET_Force1 <- 
  ggplot(aes(x = BP ,y = FoodSen_Mean ),data = d_BPs_Agg_merge) +
  geom_point(aes(color= fGhrelin), alpha =1, size = 2.5 )+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "black") +
  stat_cor(label.y = 7, label.x = 1.5)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP(ND) ') +
  ylab(label = 'Food sensitivity ') +  ggtitle("") + 
  facet_wrap(~ ROI , labeller = as_labeller(c(M_Accumbens_Mean = "Accumbens", M_Caudate_Mean = "Caudate", M_Putamen_Mean = "Putamen")))

ggsave( paste(path_out, sub_PET, "BP_IMT_FoodSens.png", sep = ""), 
        plot = plot_PET_Force1,  height = 3.5, width = 7, units = "in", dpi = 600, bg = "white")



plot_PET_Force1 <- 
  ggplot(aes(x = BP ,y = RewSens_Mean ),data = d_BPs_Agg_merge) +
  geom_point(aes(color= fGhrelin), alpha =1, size = 2.5 )+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "black") +
  stat_cor( label.y = 60, label.x = 1.5)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP(ND) ') +
  ylab(label = 'Reward Sensitivity ') +  ggtitle("") + 
  facet_wrap(~ ROI , labeller = as_labeller(c(M_Accumbens_Mean = "Accumbens", M_Caudate_Mean = "Caudate", M_Putamen_Mean = "Putamen")))

ggsave( paste(path_out, sub_PET, "BP_IMT_RewardSens.png", sep = ""), 
        plot = plot_PET_Force1,  height = 3.5, width = 7, units = "in", dpi = 600, bg = "white")

d_BPs_Agg_merge_AVG <- d_BPs_Sensitivity %>% 
  group_by(fID, cSex.x, cAge.x) %>%
  summarize(FoodSen_Mean = mean(M_RelForce_Delta5_DiffFoodMoney_AvgHL),
  RewSens_Mean = mean(M_RelForce_Delta4_HL_AvgFoodMoney) ,M_Putamen_Mean = mean(M_Putamen), 
  M_Caudate_Mean = mean(M_Caudate),M_Accumbens_Mean = mean(M_Accumbens)) %>% tidyr::pivot_longer(
                cols = c("M_Putamen_Mean", "M_Caudate_Mean", "M_Accumbens_Mean"), 
                names_to = "ROI",
                values_to = "BP")

plot_PET_Force1 <- 
  ggplot(aes(x = BP ,y = RewSens_Mean ),data = d_BPs_Agg_merge_AVG) +
  geom_point( alpha =1, size = 2.5 )+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "black") +
  stat_cor()+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP(ND) ') +
  ylab(label = 'Reward Sensitivity ') +  ggtitle("") + 
  facet_wrap(~ ROI , scale = "free_x", labeller = as_labeller(c(M_Accumbens_Mean = "Accumbens", M_Caudate_Mean = "Caudate", M_Putamen_Mean = "Putamen")))

ggsave( paste(path_out, sub_PET, "BP_IMT_RewardSens_AVG.png", sep = ""), 
        plot = plot_PET_Force1,  height = 3.5, width = 7, units = "in", dpi = 600, bg = "white")

plot_PET_Force_FS <- 
  ggplot(aes(x = BP ,y = FoodSen_Mean ),data = d_BPs_Agg_merge_AVG) +
  geom_point( alpha =1, size = 2.5 )+
   geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1, color = "black") +
  stat_cor()+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP(ND) ') +
  ylab(label = 'Food Sensitivity ') +  ggtitle("") + 
  facet_wrap(~ ROI , scale = "free_x", labeller = as_labeller(c(M_Accumbens_Mean = "Accumbens", M_Caudate_Mean = "Caudate", M_Putamen_Mean = "Putamen")))

BP_IMT_Index <- plot_grid(plot_PET_Force1, plot_PET_Force_FS, 
             label_size = 12, ncol=1)

ggsave( paste(path_out, sub_PET, "BP_IMT_FS_RS.png", sep = ""), 
        plot = BP_IMT_Index,  height = 7, width = 7, units = "in", dpi = 600, bg = "white")


# Test multivariate  
d_BPs_Agg_merge_MV <- d_BPs_Agg_merge %>% tidyr::pivot_wider(
                names_from = "ROI", values_from = "BP" )

d_BPs_Agg_merge_MV$cRewSens_Mean <- d_BPs_Agg_merge_MV$RewSens_Mean - mean(d_BPs_Agg_merge_MV$RewSens_Mean) 
d_BPs_Agg_merge_MV$cFoodSen_Mean <- d_BPs_Agg_merge_MV$FoodSen_Mean - mean(d_BPs_Agg_merge_MV$FoodSen_Mean) 

fm_PET_IMT1 <- manova(cbind(M_Putamen_Mean, M_Caudate_Mean, M_Accumbens_Mean) ~    fGhrelin + (cRewSens_Mean ) + cSession.x  + cSex.x + cAge.x + cBMI.x , d_BPs_Agg_merge_MV)
summary(fm_PET_IMT1)

fm_PET_IMT1 <- manova(cbind(M_Putamen_Mean, M_Caudate_Mean, M_Accumbens_Mean) ~    fGhrelin + (cFoodSen_Mean ) + cSession.x  + cSex.x + cAge.x + cBMI.x , d_BPs_Agg_merge_MV)
summary(fm_PET_IMT1)

dAgg_MV <- 
  d_IMT %>% 
  group_by(fID,fGhrelin,Session, cSession, cBMI, cAge, cSex) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel)) 

dAgg_MV <- merge(d_BPs_Agg_merge_MV,dAgg_MV, by = c("fID", "fGhrelin") )
dAgg_MV$cRelEff <- dAgg_MV$M_RelForce - mean(dAgg_MV$M_RelForce) 

fm_PET_IMT1 <- manova(cbind(M_Putamen_Mean, M_Caudate_Mean, M_Accumbens_Mean) ~    fGhrelin + (cRelEff ) + cSession.x  + cSex.x + cAge.x + cBMI.x , dAgg_MV)
summary(fm_PET_IMT1)

####################################################################
# IMT AND BOLD DATA 
####################################################################

# Average Force over Session 
dAgg <- 
  d_IMT %>% 
  #filter(Reward_Money == "Food") %>%
  group_by(ID,fGhrelin,Session, cAge, cBMI, cSex) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel), M_AbsForce = mean(Force_Abs))

dAgg_RewIndex <- 
  d_IMT %>% 
  group_by(ID,fGhrelin,Session, cSession, Reward_Magnitude,Reward_Money, cBMI, cAge, cSex) %>%
  dplyr::summarize(M_RelForce = mean(Force_rel)) %>%
  # separate columns for Conditions in order to substract within subject 
  tidyr::pivot_wider(names_from = c(Reward_Money, Reward_Magnitude), values_from = M_RelForce) %>%
  mutate(RewSens = (mean(c(Food_High, Money_High)))- (mean(c(Food_Low, Money_Low))))  %>%
  mutate(FoodSens = mean(c( (Food_High - Money_High), (Food_Low - Money_Low))))


# Merge Force with extracted Betas 
d_Task_fMRI_Force <- merge(dAgg, d_Task_fMRI, by = c("ID", "Session", "fGhrelin"), all = FALSE)
head(d_Task_fMRI_Force)

d_Task_fMRI_Force2 <- merge(dAgg_RewIndex, d_Task_fMRI, by = c("ID", "Session", "fGhrelin"), all = FALSE)
head(d_Task_fMRI_Force2)

# Filter for Contrast 
# 27: Workblock food + Money 
# 21: CUE food + Money > Feedback foood + money 
# 9: workblock FOOD 
# 8: workblock MONEY 
# 7: anticipation food + money 

d_Task_fMRI_Force_filter = d_Task_fMRI_Force[d_Task_fMRI_Force$Contrast == 27, ]

d_Task_fMRI_Force_filter$fID <- factor(d_Task_fMRI_Force_filter$ID)
d_Task_fMRI_Force_filter$cSession <- d_Task_fMRI_Force_filter$Session - mean(d_Task_fMRI_Force_filter$Session)

# DA Midbrain 
d_Task_fMRI_DA_Midbrain <- d_Task_fMRI_Force_filter %>% filter(ROI == 144 | ROI == 143 | ROI == 135 | ROI == 136 | ROI == 139 | ROI ==140) %>% 
    group_by(fID, ROI, cAge, cBMI, cSex, fGhrelin,cSession, M_RelForce) %>% 
    summarise(betas = mean(betas)) %>% 
    tidyr::pivot_wider(names_from = c(ROI), values_from = betas)  

names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "144"] <- "VTA1"
names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "143"] <- "VTA2"
names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "135"] <- "SNC1"
names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "136"] <- "SNC2"
names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "139"] <- "SNr1"
names(d_Task_fMRI_DA_Midbrain)[names(d_Task_fMRI_DA_Midbrain) == "140"] <- "SNr2"

d_Task_fMRI_DA_Midbrain$M_BETAS = (d_Task_fMRI_DA_Midbrain$VTA1 + d_Task_fMRI_DA_Midbrain$VTA2) + (d_Task_fMRI_DA_Midbrain$SNC1 + d_Task_fMRI_DA_Midbrain$SNC2) + (d_Task_fMRI_DA_Midbrain$SNr1 + d_Task_fMRI_DA_Midbrain$SNr2) / 6

d_Task_fMRI_DA_Midbrain$cM_RelForce <- d_Task_fMRI_DA_Midbrain$M_RelForce - mean(d_Task_fMRI_DA_Midbrain$M_RelForce)


fm_DA <- lmer(M_BETAS ~  cM_RelForce + fGhrelin + (1  |fID)  , d_Task_fMRI_DA_Midbrain)
summary(fm_DA)

fm_DA <- lmer(cM_RelForce ~  M_BETAS * fGhrelin + (1  |fID)  , d_Task_fMRI_DA_Midbrain)
summary(fm_DA)

plot_REWARD_DA  <- 
  ggplot(aes(x = M_RelForce, y = M_BETAS),data = d_Task_fMRI_DA_Midbrain) +
    geom_point(aes(color = fGhrelin), size = 3.5, alpha = .6) +
  geom_smooth(  method = 'rlm', alpha = 0.5, linewidth = 2 , color = "black") +
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
    theme(legend.position = 'bottom', text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0)) + 
   ylab(label = 'Betas midbrain \n[workblock money]') +
  xlab(label = 'Relative Effort (money)') +  ggtitle("") 

ggsave(paste(path_out, sub_IMT, "Betas_Workblock_Money.png", sep=""),  
        plot = plot_REWARD_DA,  height = 4, width = 5, units = "in", dpi = 600, bg = "white")

# FOR CINGULATE CUE PHASE 

d_Task_fMRI_CC <- d_Task_fMRI_Force_filter %>% filter(ROI == 33) %>% 
    group_by(fID, ROI, cAge, cBMI, cSex, fGhrelin,cSession, M_RelForce) %>% 
    summarise(betas = mean(betas)) %>% 
    tidyr::pivot_wider(names_from = c(ROI), values_from = betas)  

names(d_Task_fMRI_CC)[names(d_Task_fMRI_CC) == "33"] <- "Postcentral_right"

d_Task_fMRI_CC$M_BETAS = d_Task_fMRI_CC$Postcentral_right

d_Task_fMRI_CC$cM_RelForce <- d_Task_fMRI_CC$M_RelForce - mean(d_Task_fMRI_CC$M_RelForce)

fm_CC <- lm(M_BETAS ~    fGhrelin + cM_RelForce , d_Task_fMRI_CC)
summary(fm_CC)


plot_REWARD_DA  <- 
  ggplot(aes(x = M_RelForce, y = M_BETAS),data = d_Task_fMRI_CC) +
    geom_point(aes(color = fGhrelin), size = 3.5, alpha = .6) +
geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 2, color= "black" ) +
#sm_statCorr() + 
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
    theme(legend.position = 'bottom', text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0)) + 
   ylab(label = 'Betas right postcentral gyrus \n[workblock food+money]') +
  xlab(label = 'Relative Effort (food+money)') +  ggtitle("") 

ggsave(paste(path_out, sub_IMT, "Betas_Effort_Post_Central.png", sep=""),  
        plot = plot_REWARD_DA,  height = 5, width = 7, units = "in", dpi = 600, bg = "white")




# Striatum 
d_Task_fMRI_STRIATUM <- d_Task_fMRI_Force_filter %>% filter(ROI == 105 | ROI == 104 | ROI == 94 | ROI == 95 | ROI == 96 | ROI ==97) %>% 
    group_by(fID, ROI, cAge, cBMI, cSex, fGhrelin,cSession, M_RelForce) %>% 
    summarise(betas = mean(betas)) %>% 
    tidyr::pivot_wider(names_from = c(ROI), values_from = betas)  

names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "94"] <- "Caudate_l"
names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "95"] <- "Caudate_r"
names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "105"] <- "Accumbens_l"
names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "104"] <- "Accumbens_r"
names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "96"] <- "Putamen_l"
names(d_Task_fMRI_STRIATUM)[names(d_Task_fMRI_STRIATUM) == "97"] <- "Putamen_r"

d_Task_fMRI_STRIATUM$M_BETAS = (d_Task_fMRI_STRIATUM$Caudate_r + d_Task_fMRI_STRIATUM$Caudate_l) + (d_Task_fMRI_STRIATUM$Putamen_r + d_Task_fMRI_STRIATUM$Putamen_l) + (d_Task_fMRI_STRIATUM$Accumbens_r + d_Task_fMRI_STRIATUM$Accumbens_l) / 6
d_Task_fMRI_STRIATUM$M_BETAS =   (d_Task_fMRI_STRIATUM$Accumbens_r + d_Task_fMRI_STRIATUM$Accumbens_l) / 2

d_Task_fMRI_STRIATUM$cM_RelForce <- d_Task_fMRI_STRIATUM$M_RelForce - mean(d_Task_fMRI_STRIATUM$M_RelForce)


fm_STR <- lmer(M_BETAS ~  cM_RelForce  + fGhrelin  + cSession + cBMI + cSex + cAge + 
            (1  |fID), d_Task_fMRI_STRIATUM)
summary(fm_STR)


plot_REWARD_STRIATUM  <- 
  ggplot(aes(x = M_BETAS , y = M_RelForce),data = d_Task_fMRI_STRIATUM) +
    geom_point(aes(color = fGhrelin), size = 1.5, alpha = .6) +
  geom_smooth( aes(color = fGhrelin), method = 'rlm', alpha = 0.5, linewidth = 1) +
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin )) +
   theme(legend.position = 'bottom', text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0)) + 
  geom_hline(yintercept = 0) +
  xlab(label = 'betas [workblock food +money]') +
  ylab(label = 'Relative effort (%)') +  ggtitle("") 
























