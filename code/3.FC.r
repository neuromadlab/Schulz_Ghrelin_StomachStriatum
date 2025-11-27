###############################################################################
# TUE008 Neuroimaging Study: Plotting FC during PET-MR 
# Corinna Schulz, June 2024 
###############################################################################

# Mainly: Double Check SPM Results by doing ROI-analyses (should yield sim. effects)
# Investigate temporal evolution of functional connectivity/visualize 
# check correspondences FC with PET, IMT

###############################################################################
# (1) SET UP
###############################################################################

# Install and load all required libraries
if (!require("librarian")) install.packages("librarian")
librarian::shelf(readxl, lme4, lmerTest, foreign, MASS, 
                 sjPlot, stargazer, table1, RColorBrewer, wordcloud, wordcloud2, tm, MASS, readr,
                 ggplot2, ggpubr, cowplot, viridis, tidybayes, dplyr, httpgd, 
                 languageserver, dabestr, smplot2, tidyr, dplyr, ggridges 

)

# install.packages("reshape")
# R.matlab, reshape

#install.packages("devtools")
#devtools::install_github("GRousselet/rogme")
library(rogme)


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

sub_FC <- "FC/"
sub_PET <- "PET/"

# Load Condition file, with Matching of Conn and Study ID 
d_conds <- read_excel(paste(path_in,"T_COV_Tabelle.xlsx", sep = ""))

# Important: ID is here the CONN Subject ID
length(unique(d_conds$ID))

# Read in FC Values 
d_FC <- read_excel(paste(path_in,"FC_CONN_extractedROIs.xlsx", sep = ""))

length(unique(d_FC$ID_Conn))

# Prepare FC DF in longformat 

# d_FC_long <- data.table::melt(d_FC, measure.vars = 2:9, variable.name = "PhaseSession")
d_FC_long  <- d_FC %>%
  pivot_longer(cols = c(2:9),
               names_to = c("PhaseSession"))

# Add Columns indicating PHASE and SESSION
d_FC_long <- d_FC_long %>% 
  mutate(Session = case_when( PhaseSession == "Phase1_S1" | PhaseSession == "Phase2_S1" | PhaseSession == "Phase3_S1" | PhaseSession == "Phase4_S1"   ~ 1,
                                  PhaseSession == "Phase1_S2" | PhaseSession == "Phase2_S2" | PhaseSession == "Phase3_S2" | PhaseSession == "Phase4_S2"   ~ 2)) %>% 
  mutate(Phase = case_when( PhaseSession == "Phase1_S1" | PhaseSession == "Phase1_S2"  ~ 1, 
                                  PhaseSession == "Phase2_S1" | PhaseSession == "Phase2_S2"  ~ 2,
                                  PhaseSession == "Phase3_S1" | PhaseSession == "Phase3_S2"  ~ 3,
                                  PhaseSession == "Phase4_S1" | PhaseSession == "Phase4_S2"  ~ 4)) %>% dplyr::rename(ID = ID_Conn)

# Load Atlas labels file 
d_labels <- read_excel(paste(path_in,"Atlas.xlsx", sep = ""), col_names = FALSE)

# Merge Condition file with FC data and Covariates 
# ID here == CONN ID
d_FC_long_cond <- merge(d_FC_long, d_conds, by = c("ID", "Session","Phase"))
# Write FC sumamry file before overwriting STUDY ID WITH CONN ID 
# write.csv(d_FC_long_cond, "/mnt/TUE_general/SynologyDrive/10_analysis/fMRI/FC/Hypothalamus_FC/Extracted_ROIs/ExtractedFC_preprocessed.csv", row.names=FALSE)
length(unique(d_FC_long_cond$ID))

############################
# Do the same for NAcc Seed 

# Read in FC Values 
d_FC_NAcc <- read_excel(paste(path_in,"FC_CONN_extractedROIs_NAcc.xlsx", sep = ""))

length(unique(d_FC_NAcc$ID_Conn))

# Prepare FC DF in longformat 

# d_FC_long <- data.table::melt(d_FC, measure.vars = 2:9, variable.name = "PhaseSession")
d_FC_long_NAcc  <- d_FC_NAcc %>%
  pivot_longer(cols = c(2:9),
               names_to = c("PhaseSession"))

# Add Columns indicating PHASE and SESSION
d_FC_long_NAcc <- d_FC_long_NAcc %>% 
  mutate(Session = case_when( PhaseSession == "Phase1_S1" | PhaseSession == "Phase2_S1" | PhaseSession == "Phase3_S1" | PhaseSession == "Phase4_S1"   ~ 1,
                                  PhaseSession == "Phase1_S2" | PhaseSession == "Phase2_S2" | PhaseSession == "Phase3_S2" | PhaseSession == "Phase4_S2"   ~ 2)) %>% 
  mutate(Phase = case_when( PhaseSession == "Phase1_S1" | PhaseSession == "Phase1_S2"  ~ 1, 
                                  PhaseSession == "Phase2_S1" | PhaseSession == "Phase2_S2"  ~ 2,
                                  PhaseSession == "Phase3_S1" | PhaseSession == "Phase3_S2"  ~ 3,
                                  PhaseSession == "Phase4_S1" | PhaseSession == "Phase4_S2"  ~ 4)) %>% dplyr::rename(ID = ID_Conn)


# Merge Condition file with FC data and Covariates 
# ID here == CONN ID
d_FC_long_cond_NAcc <- merge(d_FC_long_NAcc, d_conds, by = c("ID", "Session","Phase"))
length(unique(d_FC_long_cond_NAcc$ID))

############################################################################################
# Load VAS data
d_vas <- read_excel(paste(path_in,"VASstate_TUE008_NIMG_all_output.xlsx", sep = ""))

# Label VAS data
d_vas$fItem <- factor(d_vas$Item, labels = c(
  "hungry", "thirsty", "tired", "full", "active", "distressed",
  "interested", "excited", "upset", "strong", "guilty", "scared",
  "hostile", "inspired", "proud", "irritable", "enthusiastic",
  "ashamed", "alert", "nervous", "determined",
  "attentive", "jittery", "afraid"))

d_vas$fTimepoint <- factor(d_vas$Timepoint, labels = c("T0", "T1", "T2","T3"))
d_vas$fID <- factor(d_vas$ID) # Actual Study ID 

# Metabolic State Items 
d_hunger <- filter(d_vas, Item == 1)
d_full <- filter(d_vas, Item == 4)
d_hunger$MetState <- (d_hunger$Rating - d_full$Rating) 
d_hunger <- d_hunger %>% select( fID,Session, fTimepoint, MetState) %>% filter(fTimepoint== "T3")

length(unique(d_hunger$ID))
length(unique(d_FC_long_cond$ID)) # CONN ID!! 
length(unique(d_FC_long_cond$fID)) # STUDY ID 

d_FC_long_cond$fSession <- factor(d_FC_long_cond$Session)
d_FC_long_cond$fGhrelin <- factor(d_FC_long_cond$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_FC_long_cond$fPhase <- factor(d_FC_long_cond$Phase)
d_FC_long_cond$cSession <- d_FC_long_cond$Session - mean(d_FC_long_cond$Session)
d_FC_long_cond$fROI <- factor(d_FC_long_cond$ROI)
d_FC_long_cond$cAge <- d_FC_long_cond$Age - mean(d_FC_long_cond$Age)
d_FC_long_cond$cBMI <- d_FC_long_cond$BMI - mean(d_FC_long_cond$BMI)
d_FC_long_cond$cSex <- d_FC_long_cond$Sex - mean(d_FC_long_cond$Sex)

d_FC_long_cond_NAcc$fSession <- factor(d_FC_long_cond_NAcc$Session)
d_FC_long_cond_NAcc$fGhrelin <- factor(d_FC_long_cond_NAcc$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_FC_long_cond_NAcc$fPhase <- factor(d_FC_long_cond_NAcc$Phase)
d_FC_long_cond_NAcc$cSession <- d_FC_long_cond_NAcc$Session - mean(d_FC_long_cond_NAcc$Session)
d_FC_long_cond_NAcc$fROI <- factor(d_FC_long_cond_NAcc$ROI)
d_FC_long_cond_NAcc$cAge <- d_FC_long_cond_NAcc$Age - mean(d_FC_long_cond_NAcc$Age)
d_FC_long_cond_NAcc$cBMI <- d_FC_long_cond_NAcc$BMI - mean(d_FC_long_cond_NAcc$BMI)
d_FC_long_cond_NAcc$cSex <- d_FC_long_cond_NAcc$Sex - mean(d_FC_long_cond_NAcc$Sex)

d_FC_long_cond <- merge(d_FC_long_cond, d_hunger, by = c("fID", "Session"))
length(unique(d_FC_long_cond$fID)) # STUDY ID 

d_FC_long_cond_NAcc <- merge(d_FC_long_cond_NAcc, d_hunger, by = c("fID", "Session"))
length(unique(d_FC_long_cond_NAcc$fID)) # STUDY ID 

length(unique(d_FC_long_cond$fROI))
length(unique(d_FC_long_cond$fID))
length(unique(d_FC_long_cond$ID))
length(unique(d_FC_long_cond$fTimepoint))
head(d_FC_long_cond)

# Collapse over Session 

d_FC_long_cond_COLLAPSED <- d_FC_long_cond %>% 
    group_by(fID, cSession , fGhrelin, Session, fROI )%>% 
  summarize(Mean_FC = mean(value)) 


## LOAD IMT DATA 

d_IMT_INDEX <- read_csv(paste(path_in,"TUE008_IMT_IndexRewardEffort.csv", sep = ""), show_col_types = FALSE)

d_IMT_INDEX$fID <- factor(d_IMT_INDEX$ID) # Actual Study ID 
d_FC_long_cond_IMT <- merge(d_FC_long_cond, d_IMT_INDEX, by = c("fID", "fGhrelin"))
d_FC_long_cond_NAcc_IMT <- merge(d_FC_long_cond_NAcc, d_IMT_INDEX, by = c("fID", "fGhrelin"))

length(unique(d_FC_long_cond_IMT$fID)) # STUDY ID 

# M_RelForce_Delta4_HL_AvgFoodMoney 
# M_RelForce_Delta5_DiffFoodMoney_AvgHL 
head(d_FC_long_cond_IMT)
head(d_FC_long_cond)
head(d_IMT_INDEX)

## FOR SPM COVARIATE ANALYSIS 
# Save T_COV with additional column: IMT Reward Sensitivity 
d_IMT_INDEX_SPM <- d_IMT_INDEX %>% select(c("ID", "Session", "M_RelForce_Delta4_HL_AvgFoodMoney",  "M_RelForce_Delta5_DiffFoodMoney_AvgHL"))
d_cond_SPM <- merge(d_conds, d_IMT_INDEX_SPM, by = c("ID", "Session"))            
#write.csv(d_cond_SPM, "/mnt/TUE_general/SynologyDrive/10_analysis/fMRI/FC/Hypothalamus_FC/Extracted_ROIs/T_COV_IMT.csv", row.names=FALSE)

d_VAS_Hunger <- read_csv(paste(path_in,"Delta_Hunger_Ghrelin.csv", sep = ""), col_names = TRUE)
d_VAS_Hunger <- d_VAS_Hunger %>% select(fID, Session, Rating_Delta)

colnames(d_FC_long_cond_IMT)[colnames(d_FC_long_cond_IMT) == 'Session.x'] <- 'Session'
d_FC_long_cond_IMT <- merge(d_FC_long_cond_IMT, d_VAS_Hunger, by = c("fID", "Session"))

colnames(d_FC_long_cond_NAcc_IMT)[colnames(d_FC_long_cond_NAcc_IMT) == 'Session.x'] <- 'Session'
d_FC_long_cond_NAcc_IMT <- merge(d_FC_long_cond_NAcc_IMT, d_VAS_Hunger, by = c("fID", "Session"))

head(d_FC_long_cond_IMT)
# Load HOMA_IR Data (preprocessed)
d_HOMA <- read_csv(paste(path_in,"T_COV_Insulin.csv", sep = ""),  show_col_types = TRUE)
d_FC_long_cond_IMT <- merge(d_FC_long_cond_IMT, d_HOMA, by = c("fID", "Session"))
d_FC_long_cond_NAcc_IMT <- merge(d_FC_long_cond_NAcc_IMT, d_HOMA, by = c("fID", "Session"))

# Load Ghrelin Plasma Data 
d_ghrelin <- read.csv(paste(path_in,"TUE008_ghrelin_summary_preprocessed.csv", sep = ""), header = TRUE)

d_ghrelin$fTimepoint <- factor(d_ghrelin$fTimepoint, labels = c("T0", "T2","T3"))
d_ghrelin$fID <- factor(d_ghrelin$fID)
d_ghrelin_long <-  d_ghrelin %>% tidyr::pivot_wider(
                names_from = c(fTimepoint), values_from = c(F_AG, F_DG, logF_AG, logF_DG, res_logF_AG, res_logF_DG))

# Load PET DATA 
d_BPs <- read_excel(paste(path_in,"TUE008_NIMG_BPs_ROIs_28_04_25.xlsx", sep = ""))
d_BPs$fID <- factor(d_BPs$ID)
d_FC_long_cond_IMT <- merge(d_FC_long_cond_IMT, d_BPs, by = c("fID", "Session"), all = TRUE)
length(unique(d_FC_long_cond_IMT$fID)) # STUDY ID 
d_FC_long_cond_NAcc_IMT <- merge(d_FC_long_cond_NAcc_IMT, d_BPs, by = c("fID", "Session"), all = TRUE)

r = 153
# Merge Plasma data with FC data 
d_plasma_change_Hyp <- merge(d_FC_long_cond_IMT, d_ghrelin_long, by = c("fID","Session"))
d_plasma_change_Hyp <-d_plasma_change_Hyp  %>%   select(c("fID",  "fGhrelin", "res_logF_AG_T2", "res_logF_AG_T3", "res_logF_DG_T2", "res_logF_DG_T3")) %>% 
          mutate(Infusion_AG = as.numeric(res_logF_AG_T3) - as.numeric(res_logF_AG_T2), 
          Infusion_DG = as.numeric(res_logF_DG_T3) - as.numeric(res_logF_DG_T2) )
d_plasma_change_Hyp <- merge(d_FC_long_cond_IMT, d_plasma_change_Hyp, by = c("fID","fGhrelin"))

d_plasma_change_NAcc <- merge(d_FC_long_cond_NAcc_IMT, d_ghrelin_long, by = c("fID","Session"))
d_plasma_change_NAcc <-d_plasma_change_NAcc  %>%   select(c("fID",  "fGhrelin", "res_logF_AG_T2", "res_logF_AG_T3", "res_logF_DG_T2", "res_logF_DG_T3")) %>% 
          mutate(Infusion_AG = as.numeric(res_logF_AG_T3) - as.numeric(res_logF_AG_T2) , 
                    Infusion_DG = as.numeric(res_logF_DG_T3) - as.numeric(res_logF_DG_T2) )
d_plasma_change_NAcc <- merge(d_FC_long_cond_NAcc_IMT, d_plasma_change_NAcc, by = c("fID","fGhrelin"))

# -------------------------------------------------------------------------

#     Striatum: (atlas ~= 94 & atlas ~= 95 & atlas ~= 96 & atlas ~= 97& atlas ~= 104 & atlas ~= 105 )  ; 

# For Seed: Hypothalamus, ROI: Striatum 
# Only take Caudate + NACc (as those were sign. on SPM analsysis to check here)
FC_Hypo_Striatum <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == "94" | d_FC_long_cond_IMT$fROI == "95" | d_FC_long_cond_IMT$fROI == "104" | d_FC_long_cond_IMT$fROI == "105" ,]

FC_Hypo_Striatum_COLLAPSED <- FC_Hypo_Striatum %>% 
        group_by(fID, cAge, cSex.x, cBMI, fPhase, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value)) %>% # filter(fPhase == "1")
        group_by(fID, cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin)%>%
        summarise(Mean_FC = mean(Mean_FC))

fm_1 <- lmer(Mean_FC  ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x + (1 |fID)  , FC_Hypo_Striatum_COLLAPSED)
summary(fm_1)

## Results corresponds with SPM analyses, ghrelin sign. increases FC 

FC_Hypo_Striatum_COLLAPSED$cHunger_delta <-FC_Hypo_Striatum_COLLAPSED$Rating_Delta - mean(FC_Hypo_Striatum_COLLAPSED$Rating_Delta)
fm_VAS <- lmer(Mean_FC  ~  fGhrelin  + cHunger_delta + cSession.x + cAge + cBMI  + cSex.x + (1| fID)  , FC_Hypo_Striatum_COLLAPSED)
summary(fm_VAS)

p_FC_FS <- 
      ggplot(aes(x = Rating_Delta ,y = Mean_FC  ),data = FC_Hypo_Striatum_COLLAPSED) +
      geom_point(aes(size = 3.5,  color = fGhrelin)) +
      geom_smooth( aes(group = fGhrelin, color = fGhrelin), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "Functional Connectivity")) )+
      xlab(label = 'Delta Hunger') +
      stat_cor(aes(group = fGhrelin))  +
      ylim(c(-0.2, 0.2))

# For Seed: NAcc, ROI: Striatum 

FC_NAcc_Striatum <- d_FC_long_cond_NAcc_IMT[d_FC_long_cond_NAcc_IMT$fROI == "95" | d_FC_long_cond_NAcc_IMT$fROI == "96" | d_FC_long_cond_NAcc_IMT$fROI == "97" | d_FC_long_cond_NAcc_IMT$fROI == "98"  ,]

FC_NAcc_Striatum_COLLAPSED <- FC_NAcc_Striatum %>% 
        group_by(fID, cAge, cSex.x, fPhase, cBMI, MetState, Rating_Delta, M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value))   %>% #filter(fPhase == "1")
        group_by(fID, cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin)%>%
        summarise(Mean_FC = mean(Mean_FC))

fm_1 <- lmer(Mean_FC  ~ fGhrelin +  cSex.x + cAge + cBMI + cSession.x + (1  |fID)  , FC_NAcc_Striatum_COLLAPSED)
summary(fm_1)

# Again corresponds to SPM Analyses 

FC_NAcc_Striatum_COLLAPSED$cHunger_delta <-FC_NAcc_Striatum_COLLAPSED$Rating_Delta - mean(FC_NAcc_Striatum_COLLAPSED$Rating_Delta)
fm_VAS <- lmer(Mean_FC  ~  fGhrelin + cHunger_delta + cSession.x + cAge + cBMI  + cSex.x + (1| fID)  , FC_NAcc_Striatum_COLLAPSED)
summary(fm_VAS)

p_FC_FS2 <- 
      ggplot(aes(x = Rating_Delta ,y = Mean_FC  ),data = FC_NAcc_Striatum_COLLAPSED) +
      geom_point(aes(size = 3.5,  color = fGhrelin)) +
      geom_smooth( aes(group = fGhrelin, color = fGhrelin), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "Functional Connectivity")) )+
      xlab(label = 'Delta Hunger') +
      stat_cor(aes(group = fGhrelin))  +
      ylim(c(-0.2, 0.2))


G <- cowplot::plot_grid(p_FC_FS, p_FC_FS2,
                  ncol = 2, nrow = 1, labels = "AUTO", label_size = 18, 
                  scale = .95, hjust = 0, align = "h")

ggsave(paste(path_out, "FC/Hypo_Striatum_DeltaHunger_Bolus.png", sep=""), 
  plot = G, height = 5, width = 8, units = "in", dpi = 600, bg = "white")

# Correlation with Plasma Ghrelin levels 

FC_Hypo_Striatum <- d_plasma_change_Hyp[d_plasma_change_Hyp$fROI == "95" | d_plasma_change_Hyp$fROI == "96" | d_plasma_change_Hyp$fROI == "94" | d_plasma_change_Hyp$fROI == "95" |   d_plasma_change_Hyp$fROI == "104" | d_plasma_change_Hyp$fROI == "105" ,]

FC_Hypo_Striatum_COLLAPSED <- FC_Hypo_Striatum %>% 
        group_by(fID,Infusion_AG, Infusion_DG,  cAge, cSex.x, cBMI, fPhase, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value)) %>% #filter(fPhase == "1")
        group_by(fID,Infusion_AG, Infusion_DG,  cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(Mean_FC))

FC_Hypo_Striatum_COLLAPSED$cInfusion_AG <- FC_Hypo_Striatum_COLLAPSED$Infusion_AG - mean(FC_Hypo_Striatum_COLLAPSED$Infusion_AG)

fm_1 <- lmer(Mean_FC  ~  cInfusion_AG  + cSex.x + cAge + cBMI + cSession.x + (1 |fID)  , FC_Hypo_Striatum_COLLAPSED)
summary(fm_1)

# Striatum
#FC_NAcc_Striatum <- d_plasma_change_NAcc[d_plasma_change_NAcc$fROI == "94" | d_plasma_change_NAcc$fROI == "95" | d_plasma_change_NAcc$fROI == "96" | d_plasma_change_NAcc$fROI == "97"  ,]

# Just Putamen (whole brain sign.) 
FC_NAcc_Striatum <- d_plasma_change_NAcc[d_plasma_change_NAcc$fROI == "96" | d_plasma_change_NAcc$fROI == "97"  ,]
FC_NAcc_NTS <- d_plasma_change_NAcc[d_plasma_change_NAcc$fROI == "153"  ,]
FC_NAcc_PaCiG <- d_plasma_change_NAcc[d_plasma_change_NAcc$fROI == "54"  | d_plasma_change_NAcc$fROI == "53"  ,]

FC_NAcc_Striatum_COLLAPSED <- FC_NAcc_Striatum %>% 
        group_by(fID, Infusion_AG,Infusion_DG,  cAge, cSex.x, fPhase, cBMI, MetState, Rating_Delta, M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value))   %>% #filter(fPhase == "1")
        group_by(fID,Infusion_AG, Infusion_DG,  cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(Mean_FC))

FC_NAcc_NTS_COLLAPSED <- FC_NAcc_NTS %>% 
        group_by(fID, Infusion_AG,Infusion_DG,  cAge, cSex.x, fPhase, cBMI, MetState, Rating_Delta, M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value))  %>% #filter(fPhase == "1")
        group_by(fID,Infusion_AG, Infusion_DG,  cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(Mean_FC))

FC_NAcc_PaCiG_COLLAPSED <- FC_NAcc_PaCiG %>% 
        group_by(fID, Infusion_AG,Infusion_DG,  cAge, cSex.x, fPhase, cBMI, MetState, Rating_Delta, M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value))  %>% #filter(fPhase == "1")
        group_by(fID,Infusion_AG, Infusion_DG,  cAge, cSex.x, cBMI, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(Mean_FC))

FC_NAcc_Striatum_COLLAPSED$cInfusion_AG <- FC_NAcc_Striatum_COLLAPSED$Infusion_AG - mean(FC_NAcc_Striatum_COLLAPSED$Infusion_AG)
fm_1 <- lmer(Mean_FC  ~ cInfusion_AG +  cSex.x + cAge + cBMI + cSession.x + (1  |fID)  , FC_NAcc_Striatum_COLLAPSED)
summary(fm_1)

fm_1 <- lmer(Mean_FC  ~ Infusion_AG +  cSex.x + cAge + cBMI + cSession.x + (1  |fID)  , FC_NAcc_NTS_COLLAPSED)
summary(fm_1)

fm_1 <- lmer(Mean_FC  ~ Infusion_AG +  cSex.x + cAge + cBMI + cSession.x + (1  |fID)  , FC_NAcc_PaCiG_COLLAPSED)
summary(fm_1)

# Plot Functional Connectivity Changes and Ghrelin Plasma Changes 
p_FC_PlasmaGhrelin <- 
      ggplot(aes(x = Infusion_AG ,y = Mean_FC  ),data = FC_NAcc_Striatum_COLLAPSED) +
      geom_point(aes(size = 3.5,  color = fGhrelin)) +
      geom_smooth( aes(), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "FC (NAcc - Putamen)")) )+
      xlab(label = 'Infusion-induced changes in Plasma Acyl Ghrelin') +
      stat_cor(aes(), label.y = 0.1) 

ggsave(paste(path_out, "FC/NAcc_Putamen_DeltaPlasmaGhrelin.png", sep=""), 
  height = 4, width = 6, units = "in", dpi = 600, bg = "white")

p_FC_PlasmaGhrelin <- 
      ggplot(aes(x = Infusion_AG ,y = Mean_FC  ),data = FC_NAcc_PaCiG_COLLAPSED) +
      geom_point(aes(size = 3.5,  color = fGhrelin)) +
      geom_smooth( aes(), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "FC (NAcc - PaCiG)")) )+
      xlab(label = 'Infusion-induced changes in Plasma Acyl Ghrelin') +
      stat_cor(aes(), label.y = 0.1) 

ggsave(paste(path_out, "FC/NAcc_PaCig_DeltaPlasmaGhrelin.png", sep=""), 
  height = 4, width = 6, units = "in", dpi = 600, bg = "white")

p_FC_PlasmaGhrelin <- 
      ggplot(aes(x = Infusion_AG ,y = Mean_FC  ),data = FC_Hypo_Striatum_COLLAPSED) +
      geom_point(aes(size = 3.5,  color = fGhrelin)) +
      geom_smooth( aes(), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "FC (Hypo - Striatum)")) )+
      xlab(label = 'Infusion-induced changes in Plasma Acyl Ghrelin') +
      stat_cor(aes(), label.y = 0.1) 

ggsave(paste(path_out, "FC/Hypo_Striatum_DeltaPlasmaGhrelin.png", sep=""), 
  height = 5, width = 6, units = "in", dpi = 600, bg = "white")

# Test Delta Delta (i.e., Delta pre post AND ghrelin vs. saline)
FC_NAcc_PaCiG_COLLAPSED_DELTA <- FC_NAcc_PaCiG_COLLAPSED  %>% ungroup() %>%  select(c("fID", "Mean_FC",  "fGhrelin", "Infusion_AG")) %>% 
              pivot_wider(names_from = fGhrelin, values_from = c(Infusion_AG, Mean_FC)) %>%
               mutate(
                    Ghrelin_Infusion_AG = Infusion_AG_Ghrelin - Infusion_AG_Placebo,
                    Ghrelin_FC_Mean     = Mean_FC_Ghrelin - Mean_FC_Placebo)

p_FC_PlasmaGhrelin <- 
      ggplot(aes(x = Ghrelin_Infusion_AG ,y = Ghrelin_FC_Mean  ),data = FC_NAcc_PaCiG_COLLAPSED_DELTA) +
      geom_point(aes(size = 3.5)) +
      geom_smooth( aes(), size = 1.5, method = 'rlm', alpha = 0.5) +
      scale_color_manual(guide = guide_legend(title="Group"),values = c( color_Placebo, color_Ghrelin)) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 16.0),
            axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) +
      ylab(label = expression(bold(Delta ~ "FC (Hypo - Striatum)")) )+
      xlab(label = 'Ghrelin-induced changes in Plasma Acyl Ghrelin (vs. saline)') +
      stat_cor(aes()) 

# Correlation Hypo and NAcc FC 

FC_Hypo_Striatum_COLLAPSED <- FC_Hypo_Striatum %>% 
        group_by(fID, cAge, cSex.x, cBMI, MetState, Rating_Delta,M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value)) 

FC_NAcc_Striatum_COLLAPSED <- FC_NAcc_Striatum %>% 
        group_by(fID, cAge, cSex.x, cBMI, MetState, Rating_Delta, M_RelForce_Delta4_HL_AvgFoodMoney, M_RelForce_Delta5_DiffFoodMoney_AvgHL, cSession.x, fGhrelin) %>% 
        summarise(Mean_FC = mean(value))


cor.test(FC_Hypo_Striatum_COLLAPSED$Mean_FC, FC_NAcc_Striatum_COLLAPSED$Mean_FC)

Both_FC <- merge(FC_Hypo_Striatum_COLLAPSED, FC_NAcc_Striatum_COLLAPSED, by = c("fID", "fGhrelin"))
Both_FC$cHunger_delta <- Both_FC$Rating_Delta.x - mean(Both_FC$Rating_Delta.x)
Both_FC$cMean_FC_Hypo <- Both_FC$Mean_FC.x - mean(Both_FC$Mean_FC.x)
Both_FC$cMean_FC_NAcc <- Both_FC$Mean_FC.y  - mean(Both_FC$Mean_FC.y )


Both_FC_MV_V2 <- lmer(cHunger_delta  ~ fGhrelin +  cMean_FC_Hypo + cMean_FC_NAcc   +cSession.x.x + cAge.x + cBMI.x  + cSex.x.x + (1 | fID)  , Both_FC)
summary(Both_FC_MV_V2)



###############################################################################################
## FC AND PET 

# Filter for PET ROIS 
head(d_FC_long_cond_IMT)
FC_Caudate_l <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 95,]
FC_Caudate_r <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 94,]
FC_Accumbens_l <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 105,]
FC_Accumbens_r <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 104,]
FC_Putamen_l <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 97,]
FC_Putamen_r <- d_FC_long_cond_IMT[d_FC_long_cond_IMT$fROI == 96,]

library(tidyr)
d_FC_wide <- d_FC_long_cond_IMT %>% filter(fROI == "95" | fROI == "94" | fROI == "105" | fROI == "104" | fROI == "97" | fROI == "96") %>% 
    group_by(fID, fROI, cAge, cBMI, cSex.x, fGhrelin,cSession.x, Session) %>% 
    summarise(Mean_FC = mean(value)) %>% 
    tidyr::pivot_wider(names_from = c(fROI), values_from = Mean_FC)  

names(d_FC_wide)[names(d_FC_wide) == "95"] <- "Caudate_l"
names(d_FC_wide)[names(d_FC_wide) == "94"] <- "Caudate_r"
names(d_FC_wide)[names(d_FC_wide) == "105"] <- "Accumbens_l"
names(d_FC_wide)[names(d_FC_wide) == "104"] <- "Accumbens_r"
names(d_FC_wide)[names(d_FC_wide) == "97"] <- "Putamen_l"
names(d_FC_wide)[names(d_FC_wide) == "96"] <- "Putamen_r"

d_FC_wide$M_Caudate_FC = (d_FC_wide$Caudate_r + d_FC_wide$Caudate_l) / 2
d_FC_wide$M_Putamen_FC = (d_FC_wide$Putamen_r + d_FC_wide$Putamen_l) / 2
d_FC_wide$M_Accumbens_FC = (d_FC_wide$Accumbens_r + d_FC_wide$Accumbens_l) / 2

# Combine PET ROIS longformat
d_BPs$fID <- factor(d_BPs$ID)
d_BPs$M_Caudate = (d_BPs$Caudate_r + d_BPs$Caudate_l) / 2
d_BPs$M_Putamen = (d_BPs$Putamen_r + d_BPs$Putamen_l) / 2
d_BPs$M_Accumbens = (d_BPs$Accumbens_r + d_BPs$Accumbens_l) / 2
d_BPs$cSession <- d_BPs$Session - mean(d_BPs$Session)


d_BPs_long <-  d_BPs %>% pivot_longer(
                cols = c("M_Caudate", "M_Putamen", "M_Accumbens"), 
                names_to = "ROI",
                values_to = "BP")

FC_PET <- merge(d_FC_wide, d_BPs_long, by = c("fID", "Session"))
length(unique(FC_PET$fID))

FC_Caudate <- FC_PET[FC_PET$ROI == "M_Caudate",] 
FC_Caudate$cBP <- FC_Caudate$BP - mean(FC_Caudate$BP)

fm_PET_0 <- lm(M_Caudate_FC ~ fGhrelin+ cSession.x + cSex.x + cAge + cBMI , FC_Caudate)
summary(fm_PET_0)
fm_PET_1 <- lm(M_Caudate_FC ~ fGhrelin * cBP + cSex.x  + cAge + cBMI + cSession.x   , FC_Caudate)
summary(fm_PET_1)

# Test whether adding BP is signficiant 
anova(fm_PET_0, fm_PET_1)

FC_Putamen <- FC_PET[FC_PET$ROI == "M_Putamen",] 
FC_Putamen$cBP <- FC_Putamen$BP - mean(FC_Putamen$BP)

fm_PET_0 <- lm(M_Putamen_FC ~ fGhrelin+ cSex.x + cAge + cBMI + cSession.x  , FC_Putamen)
summary(fm_PET_0)

fm_PET <- lm(M_Putamen_FC ~ fGhrelin * cBP + cSex.x + cAge + cBMI + cSession.x   ,FC_Putamen)
summary(fm_PET)
anova(fm_PET_0, fm_PET)


FC_NAcc <- FC_PET[FC_PET$ROI == "M_Accumbens",] 
FC_NAcc$cBP <- FC_NAcc$BP - mean(FC_NAcc$BP)

fm_PET_0 <- lm(M_Accumbens_FC ~ fGhrelin+ cSex.x + cAge + cBMI + cSession.x  , FC_NAcc)
summary(fm_PET_0)

fm_PET <- lm(M_Accumbens_FC ~ fGhrelin * cBP + cSex.x + cAge + cBMI + cSession.x   , FC_NAcc)
summary(fm_PET)
anova(fm_PET_0, fm_PET)







# ----- same for NACC- Striatum --- 

head(d_FC_long_cond_NAcc_IMT)
FC_Caudate_l <- d_FC_long_cond_NAcc_IMT[d_FC_long_cond_NAcc_IMT$fROI == 95,]
FC_Caudate_r <- d_FC_long_cond_NAcc_IMT[d_FC_long_cond_NAcc_IMT$fROI == 94,]
FC_Putamen_l <- d_FC_long_cond_NAcc_IMT[d_FC_long_cond_NAcc_IMT$fROI == 97,]
FC_Putamen_r <- d_FC_long_cond_NAcc_IMT[d_FC_long_cond_NAcc_IMT$fROI == 96,]

library(tidyr)
d_FC_wide <- d_FC_long_cond_NAcc_IMT %>% filter(fROI == "95" | fROI == "94" | fROI == "97" | fROI == "96") %>% 
    group_by(fID, fROI, cAge, cBMI, cSex.x, fGhrelin,cSession.x, Session) %>% 
    summarise(Mean_FC = mean(value)) %>% 
    tidyr::pivot_wider(names_from = c(fROI), values_from = Mean_FC)  

names(d_FC_wide)[names(d_FC_wide) == "95"] <- "Caudate_l"
names(d_FC_wide)[names(d_FC_wide) == "94"] <- "Caudate_r"
names(d_FC_wide)[names(d_FC_wide) == "97"] <- "Putamen_l"
names(d_FC_wide)[names(d_FC_wide) == "96"] <- "Putamen_r"

d_FC_wide$M_Caudate_FC = (d_FC_wide$Caudate_r + d_FC_wide$Caudate_l) / 2
d_FC_wide$M_Putamen_FC = (d_FC_wide$Putamen_r + d_FC_wide$Putamen_l) / 2

# Combine PET ROIS longformat
d_BPs$fID <- factor(d_BPs$ID)
d_BPs$M_Caudate = (d_BPs$Caudate_r + d_BPs$Caudate_l) / 2
d_BPs$M_Putamen = (d_BPs$Putamen_r + d_BPs$Putamen_l) / 2
d_BPs$cSession <- d_BPs$Session - mean(d_BPs$Session)


d_BPs_long <-  d_BPs %>% pivot_longer(
                cols = c("M_Caudate", "M_Putamen","M_Accumbens"), 
                names_to = "ROI",
                values_to = "BP")

FC_PET <- merge(d_FC_wide, d_BPs_long, by = c("fID", "Session"))
head(FC_PET)
FC_Caudate <- FC_PET[FC_PET$ROI == "M_Accumbens",] 
FC_Caudate$cBP <- FC_Caudate$BP - mean(FC_Caudate$BP)

fm_PET_0 <- lm(M_Caudate_FC ~ fGhrelin+ cSex.x + cAge + cBMI + cSession.x  , FC_Caudate)
summary(fm_PET_0)
fm_PET_1 <- lm(M_Caudate_FC ~ fGhrelin * cBP + cSex.x  + cAge + cBMI + cSession.x   , FC_Caudate)
summary(fm_PET_1)

# Test whether adding BP is signficiant 
anova(fm_PET_0, fm_PET_1)

FC_Putamen <- FC_PET[FC_PET$ROI == "M_Accumbens",] 
FC_Putamen$cBP <- FC_Putamen$BP - mean(FC_Putamen$BP)

fm_PET_0 <- lm(M_Putamen_FC ~ fGhrelin+ cSex.x + cAge + cBMI + cSession.x  , FC_Putamen)
fm_PET <- lm(M_Putamen_FC ~ fGhrelin * cBP + cSex.x + cAge + cBMI + cSession.x   ,FC_Putamen)
summary(fm_PET)
anova(fm_PET_0, fm_PET)

plot_FC_PET <- 
  ggplot(aes(x = M_Putamen_FC ,y = BP), data = FC_Putamen) +
  geom_point(aes( color = fGhrelin), size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP (Putamen)') +
  ylab(label = 'Delta FC (NAcc-Putamen)') +  ggtitle("") +
  smplot2::sm_statCorr()  

ggsave( paste(path_out, sub_PET, "FC_PET_Putamen.png", sep = ""), 
        plot = plot_FC_PET,  height = 5, width = 5, units = "in", dpi = 600, bg = "white")

plot_FC_PET <- 
  ggplot(aes(x = M_Caudate_FC ,y = BP), data = FC_Caudate) +
  geom_point(aes( color = fGhrelin), size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP (Caudate)') +
  ylab(label = 'Delta FC (NAcc-Caudate)') +  ggtitle("") +
  smplot2::sm_statCorr()  

ggsave( paste(path_out, sub_PET, "FC_PET_Caudate.png", sep = ""), 
        plot = plot_FC_PET,  height = 5, width = 5, units = "in", dpi = 600, bg = "white")


FC_Caudate_BP_NAcc <- FC_PET[FC_PET$ROI == "M_Accumbens",] 
FC_Caudate_BP_NAcc$cBP <- FC_Caudate_BP_NAcc$BP - mean(FC_Caudate_BP_NAcc$BP)

plot_FC_PET <- 
  ggplot(aes(x = M_Putamen_FC ,y = BP), data = FC_Caudate_BP_NAcc) +
  geom_point(aes( color = fGhrelin), size = 3.5, alpha =1)+
  scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo, color_Ghrelin)) + 
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'BP (NAcc)') +
  ylab(label = 'Delta FC (NAcc-Putamen)') +  ggtitle("") +
  smplot2::sm_statCorr()  

ggsave( paste(path_out, sub_PET, "FC_PET_PutamenNacc.png", sep = ""), 
        plot = plot_FC_PET,  height = 5, width = 5, units = "in", dpi = 600, bg = "white")

##################################################################
# FC & IMT DATA 
################################################################## 

d_FC_wide <- d_FC_long_cond_IMT %>% filter(fROI == "143" | fROI == "144" |  fROI == "135" | fROI == "136" |   fROI == "95" | fROI == "94" | fROI == "105" | fROI == "104" | fROI == "97" | fROI == "96") %>% 
    group_by(fID, fROI, cAge, cBMI, cSex.x, fGhrelin,cSession.x) %>% 
    summarise(Mean_FC = mean(value)) %>% 
    tidyr::pivot_wider(names_from = c(fROI), values_from = Mean_FC)  

names(d_FC_wide)[names(d_FC_wide) == "95"] <- "Caudate_l"
names(d_FC_wide)[names(d_FC_wide) == "94"] <- "Caudate_r"
names(d_FC_wide)[names(d_FC_wide) == "105"] <- "Accumbens_l"
names(d_FC_wide)[names(d_FC_wide) == "104"] <- "Accumbens_r"
names(d_FC_wide)[names(d_FC_wide) == "97"] <- "Putamen_l"
names(d_FC_wide)[names(d_FC_wide) == "96"] <- "Putamen_r"

names(d_FC_wide)[names(d_FC_wide) == "143"] <- "VTA_r"
names(d_FC_wide)[names(d_FC_wide) == "144"] <- "VTA_l"
names(d_FC_wide)[names(d_FC_wide) == "135"] <- "SN_c_r"
names(d_FC_wide)[names(d_FC_wide) == "136"] <- "SN_c_l"

d_FC_wide$M_Caudate_FC = (d_FC_wide$Caudate_r + d_FC_wide$Caudate_l) / 2
d_FC_wide$M_Putamen_FC = (d_FC_wide$Putamen_r + d_FC_wide$Putamen_l) / 2
d_FC_wide$M_Accumbens_FC = (d_FC_wide$Accumbens_r + d_FC_wide$Accumbens_l) / 2

d_FC_wide$M_VTA_FC = (d_FC_wide$VTA_r + d_FC_wide$VTA_l) / 2
d_FC_wide$M_SNc_FC = (d_FC_wide$SN_c_r + d_FC_wide$SN_c_l) / 2

# Plot FC BOXPLOT analogous to PET Boxplot for same ROIS for comparability 

d_FC_wide_Box <-  d_FC_wide %>% pivot_longer(
                cols = c("M_Putamen_FC", "M_Caudate_FC", "M_Accumbens_FC","M_VTA_FC", "M_SNc_FC"), 
                names_to = "ROI",
                values_to = "FC")


fm_FC <- lm(cbind(M_Putamen_FC, M_Caudate_FC,M_Accumbens_FC)  ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x  , d_FC_wide)
summary(fm_FC)
anova(fm_FC)
car::Anova(fm_FC)
vcov(fm_FC)

fm_FC <- lm(cbind(M_VTA_FC, M_SNc_FC,M_Accumbens_FC)  ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x  , d_FC_wide)
summary(fm_FC)
car::Anova(fm_FC)

fm_FC <- lm(FC ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x  , d_FC_wide_Box[d_FC_wide_Box$ROI == "M_Putamen_FC",])
summary(fm_FC)

fm_FC <- lm(FC ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x   , d_FC_wide_Box[d_FC_wide_Box$ROI == "M_Caudate_FC",])
summary(fm_FC)

fm_FC <- lm(FC ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x  , d_FC_wide_Box[d_FC_wide_Box$ROI == "M_Accumbens_FC",])
summary(fm_FC)

fm_FC <- lm(FC ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x   , d_FC_wide_Box[d_FC_wide_Box$ROI == "M_VTA_FC",])
summary(fm_FC)

fm_FC <- lm(FC ~ fGhrelin  + cSex.x + cAge + cBMI + cSession.x   , d_FC_wide_Box[d_FC_wide_Box$ROI == "M_SNc_FC",])
summary(fm_FC)


Box_ghrelin <- 
  ggplot(d_FC_wide_Box, aes(x = ROI, y = FC, fill = fGhrelin)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_point(aes(color= fGhrelin), alpha = 0.5, 
                 position = position_jitterdodge(jitter.width = 0.3)) +
  scale_color_manual(guide = guide_legend(title="Condition"),values = c( color_Placebo, color_Ghrelin)) +
  scale_fill_manual(guide = "none",values = c(color_Placebo, color_Ghrelin)) +
   theme(text = element_text(face = 'bold',size = 16.0),axis.text = element_text(face = 'plain',size = 16.0),
        axis.text.x = element_text(size = 16.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm"))) +
  xlab(label = 'Region of Interest') +
  ylab(label = 'Delta FC ') + 
  scale_x_discrete(labels=c("Accumbens", "Caudate","Putamen","VTA", "SNc"))

ggsave( paste(path_out, sub_FC, "BP_Boxplot_Ghrelin.png", sep = ""), 
        plot = Box_ghrelin,  height = 5, width = 8, units = "in", dpi = 600, bg = "white")



########################################################################
# Functional Connectivity Temporal Evoluition 
########################################################################

# Test whether in STRIATUM functional connectivity shifts 
d_Hyp_Striatum <- read_excel(paste(path_in,"Hypothalamus_tval_GP_Striatum.xlsx", sep = ""))
d_NAcc_Striatum <- read_excel(paste(path_in,"Accumbens_tval_GP_Striatum.xlsx", sep = ""))

# Convert wide format to long format
d_Hyp_Striatum <- d_Hyp_Striatum %>%
  tidyr::pivot_longer(cols = c( Bolus, Infusion1, Infusion2, Task), 
               names_to = "Condition", 
               values_to = "T_values") %>%
               mutate(Region = "Hypothalamus-Striatum")

d_NAcc_Striatum <- d_NAcc_Striatum %>%
  tidyr::pivot_longer(cols = c(Bolus, Infusion1, Infusion2, Task), 
               names_to = "Condition", 
               values_to = "T_values") %>%
               mutate(Region = "Accumbens-Striatum")

# Combine both datasets
df_combined <- bind_rows(d_Hyp_Striatum, d_NAcc_Striatum)
df_combined$Region <- factor(df_combined$Region, levels = c("Hypothalamus-Striatum", "Accumbens-Striatum"))  

# Ensure "Condition" is a factor with correct order
df_combined$Condition <- factor(df_combined$Condition, levels = c( "Bolus", "Infusion1", "Infusion2", "Task"))

library(ggridges)
theme_set(theme_cowplot(font_size = 12))

ggplot(df_combined, aes(x = T_values, y = Condition, fill = Condition) ) +
  geom_density_ridges(scale = 6, rel_min_height = 0.001,  
  alpha = .8, color = "white") +
  scale_fill_manual(values = c(
                                "Bolus" = "#56B4E9",   
                             "Infusion1" = "#0072B2",  
                             "Infusion2" = "#CC79A7",  
                             "Task" = "#D55E00")) +  # Smooth color transition from Bolus to Task
  labs(title = "",
       x = "T-values [ghrelin > saline]",
       y = "") +
  facet_wrap(~Region, ncol = 1) + 
  geom_vline(xintercept = 0, 
            color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 2, 
            color = "grey", linewidth = 1, linetype = "dashed") +
  theme_ridges(font_size = 13, grid = FALSE) +
  theme(legend.position = "none", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),
        axis.text.x = element_text(family="sans", size = 12.0), axis.title.x = element_text(hjust = 0.5),  # Center x-axis label
    axis.title.y = element_text(hjust = 0.5)) +
  scale_y_discrete(labels=c(expression(Delta ~ "Bolus"), expression(Delta ~ "early Infusion"), 
                 expression(Delta ~ "late Infusion") , expression(Delta ~ "Task")))

ggsave(paste(path_out, "FC/GGRidge_Striatum_G>P_Tvals.png", sep=""), 
height = 5, width = 6, units = "in", dpi = 600, bg = "white")

# STATISTICS FOR TEMP EVOLUTION 

###################################################################################
# SHIFT FUNCTION USING ROGME PACKAGE 

# For shift function domnt use Delta FC but absolut FC (incl. baseline) 
d_Hyp_Striatum <- read_excel(paste(path_in,"Hypothalamus_tval_GP_Striatum_absolut.xlsx", sep = ""))
# Convert wide format to long format
d_Hyp_Striatum <- d_Hyp_Striatum %>%
  tidyr::pivot_longer(cols = c( Baseline, Bolus, Infusion1, Infusion2, Task), 
               names_to = "Condition", 
               values_to = "T_values") %>%
               mutate(Region = "Hypothalamus-Striatum")


#> ----------------------------------------------------
#> compute shift function
set.seed(7)
g2 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Baseline"]
g1 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Bolus"]
df <- rogme::mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf0 <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf0 <- add_sf_lab(psf0, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf0[[1]] <- psf0[[1]] +  labs(x = "Bolus quantiles of \nt-values (ghrelin > saline)",
                             y = "Bolus - Baseline \nquantile differences (a.u.)")

g2 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Bolus"]
g1 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Infusion1"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

fBasics::ks2Test(g2, g1)


#> plot shift function
psf <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf <- add_sf_lab(psf, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf[[1]] <- psf[[1]] +  labs(x = "Early Infusion quantiles of \nt-values (ghrelin > saline)",
                             y = "Early Infusion - Bolus \nquantile differences (a.u.)")


## It shows a non-uniform shift between the marginal distributions, 
# with overall trend of growing differences as we progress towards the right tails 
# of the distributions. In other words, among larger observations, observations in 
# Infusion1 tend to be higher than in Bolus.

# ADD other phases 
g2 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Infusion1"]
g1 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Infusion2"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf2 <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf2 <- add_sf_lab(psf2, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf2[[1]] <- psf2[[1]] +  labs(x = "Late Infusion quantiles of \nt-values (ghrelin > saline)",
                             y = "Late - early Infusion \nquantile differences (a.u.)")


g2 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Infusion2"]
g1 <- d_Hyp_Striatum$T_values[d_Hyp_Striatum$Condition=="Task"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf3 <- plot_sf(sf, plot_theme = 2)
psf3 <- add_sf_lab(psf3, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf3[[1]] <- psf3[[1]] +  labs(x = "Task quantiles of \nt-values (ghrelin > saline)",
                             y = "Task - Late Infusion \nquantile differences (a.u.)")


#> combine plots
cowplot::plot_grid(psf0[[1]], psf[[1]], psf2[[1]], psf3[[1]],
          labels=c("A", "B", "C","D"),
          ncol = 4,
          nrow = 1,
          rel_heights = c(1, 1, 1,1),
          label_size = 18,
          hjust = -1,
          scale = .95,
          align ="v")

# save figure
ggsave(paste(path_out, sub_FC, "Shift_functions_Hypothalamus_Striatum_absolutFC.png", sep=""), 
  height = 5, width = 17, units = "in", dpi = 600, bg = "white")



#############################
# Same for NAcc-Striatum
#> --------------------------

# For shift function domnt use Delta FC but absolut FC (incl. baseline) 
d_NAcc_Striatum <- read_excel(paste(path_in,"Accumbens_tval_GP_Striatum_absolut.xlsx", sep = ""))

# Convert wide format to long format
d_NAcc_Striatum <- d_NAcc_Striatum %>%
  tidyr::pivot_longer(cols = c( Baseline, Bolus, Infusion1, Infusion2, Task), 
               names_to = "Condition", 
               values_to = "T_values") %>%
               mutate(Region = "Hypothalamus-Striatum")

#> --------------------------
#> compute shift function
set.seed(7)
g2 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Baseline"]
g1 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Bolus"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf40 <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf40 <- add_sf_lab(psf40, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf40[[1]] <- psf40[[1]] +  labs(x = "Bolus quantiles of \nt-values (ghrelin > saline)",
                             y = "Bolus - Baseline \nquantile differences (a.u.)")


g2 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Bolus"]
g1 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Infusion1"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

fBasics::ks2Test(g2, g1)

#> plot shift function
psf4 <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf4 <- add_sf_lab(psf4, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf4[[1]] <- psf4[[1]] +  labs(x = "Early Infusion quantiles of \nt-values (ghrelin > saline)",
                             y = "Early Infusion - Bolus \nquantile differences (a.u.)")


## It shows a non-uniform shift between the marginal distributions, 
# with overall trend of growing differences as we progress towards the right tails 
# of the distributions. In other words, among larger observations, observations in 
# Infusion1 tend to be higher than in Bolus.

# ADD other phases 
g2 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Infusion1"]
g1 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Infusion2"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf5 <- plot_sf(sf, plot_theme = 2)
#> add labels for deciles 1 & 9
psf5 <- add_sf_lab(psf5, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf5[[1]] <- psf5[[1]] +  labs(x = "Late Infusion quantiles of \nt-values (ghrelin > saline)",
                             y = "Late - early Infusion \nquantile differences (a.u.)")


g2 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Infusion2"]
g1 <- d_NAcc_Striatum$T_values[d_NAcc_Striatum$Condition=="Task"]
df <- mkt2(g1,g2)
sf <- shiftdhd(data = df, formula = obs ~ gr, nboot = 1000)

#> plot shift function
psf6 <- plot_sf(sf, plot_theme = 2)
psf6 <- add_sf_lab(psf6, sf, 
                  y_lab_nudge = .1, 
                  text_size = 4)
#> change axis labels
psf6[[1]] <- psf6[[1]] +  labs(x = "Task quantiles of \nt-values (ghrelin > saline)",
                             y = "Task - Late Infusion \nquantile differences (a.u.)")


#> combine plots
cowplot::plot_grid(psf0[[1]], psf[[1]], psf2[[1]], psf3[[1]],psf40[[1]], psf4[[1]],  psf5[[1]], psf6[[1]],
          labels=c("A","","", "", "B","","",""),
          ncol = 4,
          nrow = 2,
          rel_heights = c(1, 1, 1,1),
          label_size = 18,
          hjust = -4,
          scale = .85,
          align ="v")

# save figure
ggsave(paste(path_out, sub_FC, "Shift_functions_Hypo_NAcc_Abs.png", sep=""), 
  height = 10, width = 20, units = "in", dpi = 600, bg = "white")

