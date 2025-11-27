
######################################################
# TUE008 PET-MR Study 
# Corinna Schulz, corinna.schulz96@gmail.com 
# 2024 
#######################################################

# (1) SET UP

# Install and load all required libraries
if (!require("librarian")) install.packages("librarian")
librarian::shelf(readxl, lme4, lmerTest, foreign, MASS, 
                 sjPlot, stargazer, table1,
                 ggplot2, ggpubr, cowplot, viridis, tidybayes, dplyr, httpgd, 
                 languageserver, dabestr, svglite,emmeans)


#Set all themes
theme_set(theme_cowplot(font_size = 12))

# Set all Paths
setwd(getwd())
path_in <- "./input/" 
path_out <- "./output/" 

sub_VAS <-  "metabolic_state/"
sub_IMT <- "instrumental_motivation/"
sub_PET <- "PET/"

if (file.exists(paste(path_out, sub_VAS, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_VAS, sep = ""))}
if (file.exists(paste(path_out, sub_IMT, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_IMT, sep = ""))}
if (file.exists(paste(path_out, sub_PET, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_PET, sep = ""))}

# Set Colors
color_Placebo <- "darkblue"  
color_Ghrelin <- "darkgoldenrod"  

# Load VAS data
d_vas <- read_excel(paste(path_in, "VASstate_TUE008_NIMG_all_output.xlsx", sep = ""))

# Load Participant data (Blood data and Conditions)
d_blood <- read.csv(paste(path_in,"blood_preprocessed.csv", sep = ""), header = TRUE)
d_ghrelin <- read.csv(paste(path_in,"TUE008_ghrelin_summary_preprocessed.csv", sep = ""), header = TRUE)

# Load PET Data 
d_BPs <- read_excel(paste(path_in,"TUE008_NIMG_BPs_ROIs_28_04_25.xlsx", sep = ""))

# Label VAS data
d_vas$fItem <- factor(d_vas$Item, labels = c(
  "hungry", "thirsty", "tired", "full", "active", "distressed",
  "interested", "excited", "upset", "strong", "guilty", "scared",
  "hostile", "inspired", "proud", "irritable", "enthusiastic",
  "ashamed", "alert", "nervous", "determined",
  "attentive", "jittery", "afraid"))

d_vas$fTimepoint <- factor(d_vas$Timepoint, labels = c("T0", "T1", "T2","T3"))

length(unique(d_vas$ID))

# Prep Confounds
d_blood$cAge <- d_blood$Age - mean(d_blood$Age)
d_blood$cBMI <- d_blood$BMI - mean(d_blood$BMI)
d_blood$cSex <- d_blood$Sex_female - mean(d_blood$Sex_female)
d_blood$fSex <- factor(d_blood$Sex_female, labels = c("male", "female"))

# Merge Condition file with VAS data 
d_NIMG <- merge(d_vas, d_blood, by = c("ID", "Session"))
d_NIMG$fID <- factor(d_NIMG$ID)
d_NIMG$fGhrelin <- factor(d_NIMG$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_NIMG$cGhrelin <- d_NIMG$Ghrelin - mean(d_NIMG$Ghrelin)
d_NIMG$cSession <- d_NIMG$Session - mean(d_NIMG$Session)

d_NIMG$Item_cat <- factor(d_NIMG$Item_cat, labels = c("HState","Mood"))
d_NIMG$Item_rec <- factor(d_NIMG$Item_rec, labels = c("-Mood","HState","+Mood"))

#d_NIMG <- merge(d_NIMG, d_BPs, by = c("ID","Session"), all = TRUE)
length(unique(d_NIMG$ID))

d_ghrelin$fTimepoint <- factor(d_ghrelin$fTimepoint, labels = c("T0", "T2","T3"))
d_ghrelin$fID <- factor(d_ghrelin$fID)

d_NIMG <- d_NIMG %>%
  left_join(d_ghrelin, by = c("fID" = "fID", "Session" = "Session", "fTimepoint" = "fTimepoint"))

# Metabolic State Items 
d_hunger <- filter(d_NIMG, Item == 1)
d_full <- filter(d_NIMG, Item == 4)
d_hunger$MetState <- (d_hunger$Rating - d_full$Rating)

# Positive and Negative Affect Scales 
d_PA <- d_NIMG %>%
  filter((Item %in% c(5,7,8,10,14,15,17,19,21,22))) %>% 
  group_by(fID, fGhrelin,cGhrelin, Session, fTimepoint , cAge, cBMI, cSession, cSex, res_logF_AG )  %>% 
  summarize(VAS_PA_s  = mean(Rating)) 

d_NA <- d_NIMG %>%
  filter((Item %in% c(6,9,11,12,13,16,18,20,23,24))) %>% 
  group_by(fID, fGhrelin,Session,  cGhrelin, fTimepoint, cAge, cBMI, cSession , cSex, res_logF_AG  )  %>% 
  summarize(VAS_NA_s = mean(Rating)) 

d_PA_NA <- d_NIMG %>%
  group_by(fID, fGhrelin,Session,  cGhrelin, fTimepoint, cAge, cBMI, cSex, cSession, res_logF_AG)  %>% 
  summarize(VAS_PANA =  mean(Rating[Item %in% c(5,7,8,10,14,15,17,19,21,22)] - mean(Rating[Item %in% c(6,9,11,12,13,16,18,20,23,24)])) )

d_PANA <- merge(d_PA, d_NA, by = c("fID","fTimepoint", "fGhrelin"))

# Delta Values 

d_PA_Delta <-  d_PA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3") %>%
  group_by(fID, fGhrelin) %>%
  mutate(RatingValue_Delta = VAS_PA_s - VAS_PA_s[fTimepoint == "T2"]) %>%
  filter(fTimepoint == "T3")

d_NA_Delta <-  d_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3") %>%
  group_by(fID, fGhrelin) %>%
  mutate(RatingValue_Delta = VAS_NA_s - VAS_NA_s[fTimepoint == "T2"]) %>%
  filter(fTimepoint == "T3")

d_PA_NA_Delta <- d_PA_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3") %>%
  group_by(fID, fGhrelin) %>%
  mutate(RatingValue_Delta = VAS_PANA - VAS_PANA[fTimepoint == "T2"]) %>%
  filter(fTimepoint == "T3")

d_PANA_Intervention <- d_PA_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")

d_NA_Intervention <- d_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")

d_PA_Intervention <- d_PA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")

d_Hunger_Intervention <- d_hunger %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")

d_Full_Intervention <- d_full %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")

# FOR NIMG SPM ANALYSIS 
d_Hunger_Intervention_Delta <- d_Hunger_Intervention %>%
   group_by(fID, fGhrelin) %>%
   mutate(Rating_Delta = Rating - Rating[fTimepoint == "T2"]) %>%
   filter(fTimepoint == "T3") %>%
   select(fID, Session, Rating_Delta)

#write.csv(d_Hunger_Intervention_Delta, paste(path_in, "Delta_Hunger_Ghrelin.xlsx", sep = ""), row.names=FALSE)

#######################################################################################
# METABOLIC STATE 
#######################################################################################

# Hunger across time
pVAS1 <- ggplot(aes(x = fTimepoint, y = Rating, color = fGhrelin), data = d_hunger) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Hunger Across Time",
       x = "Timepoint",
       y = "Hunger (VAS)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  xlab(label = "Timepoint") +
  ylab(label = "Hunger (VAS)") +
  theme(legend.position = "bottom") 

  ggsave(paste(path_out, sub_VAS, "MetStates_Ghrelin.png", sep = ""),
        height = 4, width = 5, units = "in", dpi = 600, bg = "transparent")

pVAS2 <- ggplot(aes(x = fTimepoint, y = Rating, color = fGhrelin), data = d_full) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Fullness Across Time",
       x = "Timepoint",
       y = "Fullness (VAS)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  xlab(label = "Timepoint") +
  ylab(label = "Fullness (VAS)") +
  theme(legend.position = "bottom") 

# Metabolic state
pVAS3 <- ggplot(aes(x = fTimepoint, y = MetState, color = fGhrelin), data = d_hunger) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "MetState Across Time",
       x = "Timepoint",
       y = "MetState (VAS)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  xlab(label = "Timepoint") +
  ylab(label = "MetState (VAS)") +
  theme(legend.position = "bottom") 

pMetStates <- cowplot::plot_grid(pVAS1, pVAS2, label_size = 12, ncol=2, rel_widths = c(1, 1))
print(pMetStates)
ggsave(paste(path_out, sub_VAS, "MetStates_Ghrelin.png", sep = ""),
       plot = pMetStates,  height = 4, width = 10, units = "in", dpi = 600, bg = "white")



# --------------------------------------
# Simplified Figures for the Hunger Effect 

# Hunger 
head(d_hunger)
d_MetState_plot <- d_hunger %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")%>%
  mutate(PRE_POST = ifelse(fTimepoint %in% c("T2"), "PRE", "POST"))

d_MetState_plot2 <- d_MetState_plot %>%
  group_by(fID, fGhrelin, PRE_POST) %>%
  summarize(VAS_MetState_M =  mean(Rating)) %>% 
  group_by(fID, fGhrelin ) %>% 
  mutate(VAS_MetState_Delta = VAS_MetState_M - VAS_MetState_M[PRE_POST=="PRE"]) %>% 
  filter(PRE_POST == "POST") 

ggplot(d_MetState_plot2, aes(x = fGhrelin, y = VAS_MetState_Delta, fill = fGhrelin, color = fGhrelin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  # Mean bars with CI
  stat_summary(fun.data = mean_cl_boot, geom = "bar", 
               width = 0.5, alpha = 0.6, color = NA) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
             width = 0.2) +

  # Individual participant dots
  # geom_jitter(width = 0.1, size = 2, alpha = 0.7) +

  # Colors
  scale_fill_manual(values = c(color_Placebo, color_Ghrelin),
                    labels = c("Placebo", "Ghrelin"),
                    name = "Condition") +
  scale_color_manual(values = c(color_Placebo, color_Ghrelin),
                     guide = "none") +

  labs(title = "Hunger Change (Post – Pre)",
       x = "Condition",
       y = "Δ Hunger (Post - Pre)") +
  theme(legend.position = "none") + 
  theme( legend.position = "none", legend.justification = c("center"), text = element_text(face = 'bold',size = 12.0),
        axis.text = element_text(face = 'plain',size = 18),
        axis.title=element_text(size=18), axis.text.x = element_text(size = 18), 
        title = element_text(size = 20)) 

ggsave(paste(path_out, sub_VAS, "Hunger_Ghrelin_Delta.png", sep = ""),
        height = 5, width = 5, units = "in", dpi = 600, bg = "white")


# Compute effect size for this simple difference:.groups

# Step 1: Reshape to wide format
df_wide_hunger <- d_MetState_plot2 %>%
  select(fID, fGhrelin, VAS_MetState_Delta) %>%
  tidyr::pivot_wider(names_from = fGhrelin, values_from = VAS_MetState_Delta)

# Step 2: Compute difference (Ghrelin - Placebo)
df_wide_hunger <- df_wide_hunger %>%
  mutate(diff = Ghrelin - Placebo)

# Step 3: Compute Cohen's dz
dz <- mean(df_wide_hunger$diff, na.rm = TRUE) / sd(df_wide_hunger$diff, na.rm = TRUE)
dz

# ---- Hunger and Ghrelin -------
cor.test(d_hunger$res_logF_AG, d_hunger$Rating)

#cor_ghr_hunger <- d_hunger %>% group_by(fTimepoint) %>% summarize(r=cor(Rating, res_logF_AG))

p4<-
  ggplot(aes(x = Rating,y = res_logF_AG, fill = fGhrelin),data = d_hunger) +
  #smplot2::sm_statCorr()  +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  #geom_smooth(method = 'rlm', alpha = 0.1) +
  geom_smooth(aes(group = fGhrelin, color = fGhrelin), method = 'rlm', size = 1.5, alpha = 0.5) +
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Correlation of hunger and ghrelin levels") +
    facet_wrap(~fTimepoint, scale = "free") +
  ylab(label = 'Fasting levels of acyl ghrelin (res log)') +
  xlab(label = 'VAS Hunger ratings')

ggsave(paste(path_out, sub_VAS , "Hunger_AG.png", sep=""), 
              plot = p4,  height =6, width = 10, units = "in", dpi = 600, bg = "white")



#######################################################################################
# MOOD STATE 
#######################################################################################

pPA <- ggplot(aes(x = fTimepoint, y = VAS_PA_s, color = fGhrelin), data = d_PA) +
  #geom_line(aes(group = ID), color = "darkgray", size = 1, alpha = 0.25) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Positive Affect",
       x = "Timepoint",
       y = "Positive Affect (mean)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  theme(legend.position = "bottom")

pNA <- ggplot(aes(x = fTimepoint, y = VAS_NA_s, color = fGhrelin), data = d_NA) +
  #geom_line(aes(group = ID), color = "darkgray", size = 1, alpha = 0.25) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Negative Affect",
       x = "Timepoint",
       y = "Negative Affect (mean)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  theme(legend.position = "bottom")

pNAPA <- ggplot(aes(x = fTimepoint, y = VAS_PANA, color = fGhrelin), data = d_PA_NA) +
  #geom_line(aes(group = ID), color = "darkgray", size = 1, alpha = 0.25) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Mood state",
       x = "Timepoint",
       y = "Mood state (pos - neg affect)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  theme(legend.position = "bottom")

pPANAS <- cowplot::plot_grid(pPA, pNA, pNAPA, labels = "auto", label_size = 12, ncol=3, rel_widths = c(1, 1, 1))
print(pPANAS)
ggsave(paste(path_out, sub_VAS, "PANAS_Ghrelin.png", sep = ""),
       plot = pPANAS,  height = 5, width = 12, units = "in", dpi = 600, bg = "white")

# --------------------------------------
# Simplified Figures for the Mood Effect 

# Mood 
d_PA_NA_plot <- d_PA_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")%>%
  mutate(PRE_POST = ifelse(fTimepoint %in% c("T2"), "PRE", "POST"))

d_PA_NA_plot2 <- d_PA_NA_plot %>%
  group_by(fID, fGhrelin, PRE_POST) %>%
  summarize(VAS_PANA_M =  mean(VAS_PANA)) %>% 
  group_by(fID, fGhrelin ) %>% 
  mutate(VAS_PANA_Delta = VAS_PANA_M - VAS_PANA_M[PRE_POST=="PRE"]) %>% 
  filter(PRE_POST == "POST") 

pre_avg <- d_PA_NA_plot %>%
  filter(PRE_POST == "PRE") %>%
  group_by(fID, fGhrelin) %>%
  summarize(pre_mean = mean(VAS_PANA))

d_PA_NA_delta <- d_PA_NA_plot %>%
  filter(PRE_POST == "POST") %>%
  left_join(pre_avg, by = c("fID", "fGhrelin")) %>%
  mutate(VAS_PANA_Delta = VAS_PANA - pre_mean)


ggplot(d_PA_NA_delta, aes(x = fGhrelin, y = VAS_PANA_Delta, fill = fGhrelin, color = fGhrelin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  # Mean bars with CI
  stat_summary(fun.data = mean_cl_boot, geom = "bar", 
               width = 0.5, alpha = 0.6, color = NA) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               width = 0.2) +

  # Individual participant dots
  #geom_jitter(width = 0.1, size = 2, alpha = 0.7) +

  # Colors
  scale_fill_manual(values = c(color_Placebo, color_Ghrelin),
                    labels = c("Placebo", "Ghrelin"),
                    name = "Condition") +
  scale_color_manual(values = c(color_Placebo, color_Ghrelin),
                     guide = "none") +

  labs(title = "Mood Change (Post – Pre)",
       x = "Condition",
       y = "Δ Mood (Post - Pre)") +
  theme(legend.position = "none") + 
  theme( legend.position = "none", legend.justification = c("center"), text = element_text(face = 'bold',size = 12.0),
        axis.text = element_text(face = 'plain',size = 18),
        axis.title=element_text(size=18), axis.text.x = element_text(size = 18), 
        title = element_text(size = 20)) 

ggsave(paste(path_out, sub_VAS, "PANAS_Ghrelin_Delta.png", sep = ""),
        height = 5, width = 5, units = "in", dpi = 600, bg = "white")


# Compute effect size for this simple difference:.groups

# Step 1: Reshape to wide format
df_wide_mood <- d_PA_NA_delta %>% ungroup() %>%
  dplyr::select(fID, fGhrelin, VAS_PANA_Delta) %>%
  tidyr::pivot_wider(names_from = fGhrelin, values_from = VAS_PANA_Delta)

# Step 2: Compute difference (Ghrelin - Placebo)
df_wide_mood <- df_wide_mood %>%
  mutate(diff = Ghrelin - Placebo)

# Step 3: Compute Cohen's dz
dz <- mean(df_wide_mood$diff, na.rm = TRUE) / sd(df_wide_mood$diff, na.rm = TRUE)
dz


# -- Hunger and Mood changes correlated? 
# Correlated Delta (Ghrelin Saline) Hunger with Delta (Ghrelin Saline) Mood 

cor.test(df_wide_hunger$diff, df_wide_mood$diff)
deltas_hunger_mood <- merge(df_wide_hunger,df_wide_mood, by = "fID")

plot_hunger_mood_relation  <- 
  ggplot(aes(x = diff.x, y = diff.y),data = deltas_hunger_mood) +
   geom_point(size = 3.5) +
 geom_smooth( method = 'rlm', alpha = 0.5, linewidth = 1) +
  theme(legend.position = "bottom", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) +
  xlab(label = 'Ghrelin-induced changes in hunger') +
  ylab(label = 'Ghrelin-induced changes in mood') +  ggtitle("") 

# Change in Ghrelin & Change in Hunger/ MOOD 

d_ghrelin_change <- d_PA_NA %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")%>%
  mutate(PRE_POST = ifelse(fTimepoint %in% c("T2"), "PRE", "POST"))

d_ghrelin_change2 <- d_hunger %>%
  filter(fTimepoint == "T2" | fTimepoint == "T3")%>%
  mutate(PRE_POST = ifelse(fTimepoint %in% c("T2"), "PRE", "POST"))

d_ghrelin_change_HM <- merge(d_ghrelin_change, d_ghrelin_change2, by = c("fID", "fGhrelin","PRE_POST","res_logF_AG"))

d_ghrelin_change_HM <- d_ghrelin_change_HM %>%
  group_by(fID, fGhrelin) %>%
  mutate(
    res_logF_AG_Delta = res_logF_AG - res_logF_AG[PRE_POST=="PRE"], 
        res_logF_DG_Delta = res_logF_DG - res_logF_DG[PRE_POST=="PRE"], 
    VAS_PANA_Delta = VAS_PANA - VAS_PANA[PRE_POST=="PRE"],
    Hunger_Delta = Rating - Rating[PRE_POST=="PRE"]) %>% 
  filter(PRE_POST == "POST") 

M_NA <- lmer(VAS_PANA_Delta ~  res_logF_DG_Delta * fGhrelin + cSession.x +  cSex.x + cAge.x + cBMI.x   + (1 |fID), d_ghrelin_change_HM)
summary(M_NA)

M_NA <- lmer(Hunger_Delta ~  res_logF_AG_Delta   + cSession.x +  cSex.x + cAge.x + cBMI.x  + (1|fID), d_ghrelin_change_HM)
summary(M_NA)

plot_AG_Hunger <-
  ggplot(aes(x = res_logF_AG_Delta,y = Hunger_Delta),data = d_ghrelin_change_HM) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  geom_smooth(method = 'rlm', alpha = 0.1) +
  #smplot2::sm_statCorr() +
  geom_smooth(aes(), method = 'rlm', size = 1.5, alpha = 0.5) +
 scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  xlab(label = 'Change in fasting levels of acyl ghrelin (res log) (T2-T3)') +
  ylab(label = 'Change in hunger ratings (T2-T3)') +
  stat_cor(aes(), label.y = 75) 


ggsave(paste(path_out, sub_VAS , "AG_Hunger_Change.png", sep=""), 
              plot = plot_AG_Hunger,  height =4, width = 5, units = "in", dpi = 600, bg = "white")

plot_AG_Mood <-
  ggplot(aes(x = res_logF_AG_Delta,y = VAS_PANA_Delta, color = fGhrelin),data = d_ghrelin_change_HM) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  #geom_smooth(method = 'rlm', alpha = 0.1) +
  #geom_smooth(aes(group = fGhrelin, color = fGhrelin), method = 'rlm', size = 1.5, alpha = 0.5) +
scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  xlab(label = 'Change in fasting levels of acyl ghrelin (res log) (T2-T3)') +
  ylab(label = 'Change in Mood ratings (T2-T3)')+
   smplot2::sm_statCorr()  

ggsave(paste(path_out, sub_VAS , "AG_Mood_Change.png", sep=""), 
              plot = plot_AG_Mood,  height =6, width = 8, units = "in", dpi = 600, bg = "white")


# Single VAS Items exploration 
# Idea: Ghrelin might be more related to arousal than valence? 
# Cleaveland dot plot 

d_items_intervention <- d_NIMG %>%
   group_by(fID, fGhrelin, fItem, Item_rec) %>%
   mutate(Rating_Delta = Rating - Rating[fTimepoint == "T2"]) %>%
      filter(fTimepoint == "T3", Item_rec != "HState") %>%
   #filter(fTimepoint == "T3") %>%
   select(fID, fGhrelin, Session, fItem, Rating_Delta, Item_rec)%>%
   group_by(fGhrelin, fItem, Item_rec) %>%
  summarize(M_Rating_Delta = mean(Rating_Delta)) %>%
     group_by( fItem,Item_rec) %>%
  mutate(M_Rating_Delta_Delta = M_Rating_Delta - M_Rating_Delta[fGhrelin == "Placebo"])

head(d_items_intervention)

Names <- c( `-Mood` = "Negative affect",
            `+Mood` = "Positive affect", 
            `HState` = "Physiological state")

Names <- c( `-Mood` = "Negative affect",
            `+Mood` = "Positive affect")

Items_Infusion <- ggplot(d_items_intervention, aes(x = M_Rating_Delta , y = reorder(fItem,abs(M_Rating_Delta_Delta) )), color = fGhrelin) +
    geom_line(aes(group = fItem), alpha = .5, size = 1) +
facet_grid( ~ Item_rec, labeller = as_labeller(Names)) +
  geom_point(aes(color = fGhrelin), size =3, alpha = .9) +
  geom_vline(xintercept = 0, color = "darkgray", size = 1, alpha = 0.25) +
  scale_color_manual(guide = guide_legend(title="Condition"),values = c( color_Placebo, color_Ghrelin)) +
    labs(x = "Rating change (pre- to post- infusion)", y = "Item") +
     theme(legend.position = "none", legend.justification = c("center"),text = element_text(family="sans", face = 'bold',size = 16.0),
        axis.text = element_text(family="sans", face = 'plain',size = 12.0),axis.text.x = element_text(family="sans", size = 12.0)) 
  
ggsave(paste(path_out, sub_VAS , "Items_Infusion_Change.png", sep=""), 
              plot = Items_Infusion,  height =6, width = 8, units = "in", dpi = 600, bg = "white")


#################################################################
# Plasma Ghrelin ------------------------------------------------
#################################################################

# Plasma Ghrelin across time

# Raw values 
pVAS_raw <- ggplot(aes(x = fTimepoint, y = F_AG, color = fGhrelin), data = d_hunger) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Plasma Ghrelin Across Time",
       x = "Timepoint",
       y = "Plasma Ghrelin (raw)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  xlab(label = "Timepoint") +
  ylab(label = "Plasma Ghrelin (raw)") +
  theme(legend.position = "bottom") 






# Res log 
pVAS1 <- ggplot(aes(x = fTimepoint, y = res_logF_AG, color = fGhrelin), data = d_hunger) +
  geom_line(aes(group = fGhrelin, color = fGhrelin), stat = "summary", size = 2) +
  stat_summary(aes(group = fGhrelin, color = fGhrelin),
               fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.1)) +
  scale_color_manual(guide = guide_legend(title = "Group"),
                     values = c(color_Placebo, color_Ghrelin),
                     labels = c("Placebo", "Ghrelin")) +
  labs(title = "Plasma Ghrelin Across Time",
       x = "Timepoint",
       y = "Plasma Ghrelin (res log)",
       color = "Group") +
  theme(text = element_text(face = "bold", size = 12.0),
        axis.text = element_text(face = "plain", size = 12.0),
        axis.text.x = element_text(size = 12.0),
        strip.text.x = element_text(margin = margin(0.15, 0, 0.15, 0, "cm"))) +
  xlab(label = "Timepoint") +
  ylab(label = "Plasma Ghrelin (res log)") +
  theme(legend.position = "bottom") 

# Change in Acyl and Change in Desacyl 
plot_AG_DG <-
  ggplot(aes(x = res_logF_AG_Delta,y = res_logF_DG_Delta, color = fGhrelin),data = d_ghrelin_change_HM) +
  geom_point(aes(size = BMI, color = fGhrelin)) +
  #geom_smooth(method = 'rlm', alpha = 0.1) +
  #geom_smooth(aes(group = fGhrelin, color = fGhrelin), method = 'rlm', size = 1.5, alpha = 0.5) +
scale_color_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  scale_fill_manual(guide = guide_legend(title="Group"),values = c(color_Placebo,color_Ghrelin)) +
  theme(text = element_text(face = 'bold',size = 12.0),axis.text = element_text(face = 'plain',size = 12.0),
        axis.text.x = element_text(size = 12.0), strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
        plot.title = element_text(hjust = 0.5)) +
  xlab(label = 'Change in levels of acyl ghrelin (res log) (T2-T3)') +
  ylab(label = 'Change in levels of desacyl ghrelin (res log) (T2-T3)')+
   smplot2::sm_statCorr()  

ggsave(paste(path_out, sub_VAS , "AG_DG_Change.png", sep=""), 
              plot = plot_AG_DG,  height =6, width = 8, units = "in", dpi = 600, bg = "white")





# ------------------------------------------------------------
# Raincloud plots (3 timepoints, lines per subject, 2 conditions)
# ------------------------------------------------------------

#remotes::install_github('jorvlan/raincloudplots')
library(raincloudplots)

pal <- c(color_Placebo, color_Ghrelin, color_Placebo,
         color_Ghrelin, color_Placebo, color_Ghrelin)

df_2x3 <- data_2x2(
  array_1 = d_hunger$F_AG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Placebo"],
  array_2 = d_hunger$F_AG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Ghrelin"],
  array_3 = d_hunger$F_AG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Placebo"], 
  array_4 = d_hunger$F_AG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Ghrelin"],
  array_5 = d_hunger$F_AG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Placebo"], 
  array_6 = d_hunger$F_AG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Ghrelin"],
  labels = (c('Placebo','Ghrelin')),
  jit_distance = .05,
  jit_seed = 321) 

raincloud_2x3_vertical <- raincloud_2x3_repmes(
  data = df_2x3,
  colors = pal,
  fills = pal,
  size = 2,
  alpha = .6,
  ort = 'v') +

scale_x_continuous(breaks=c(1,2,3), labels=c("T-1", "T-2", "T-3"), limits=c(0, 4)) +
  xlab("Timepoint") +
  ylab("Acyl ghrelin plasma (pg/ml)") +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "plain", size = 12),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )


raincloud_2x3_vertical
ggsave(paste(path_out, sub_VAS , "AG_Timepoints_Raw.png", sep=""), 
              plot = raincloud_2x3_vertical,  height =20, width = 8, units = "in", dpi = 600, bg = "white")

# Plot Log with labels pg/mol 
means_raw <- d_hunger %>%
  filter(fTimepoint %in% c("T0","T2","T3"),
         fGhrelin %in% c("Placebo","Ghrelin")) %>%
  group_by(fGhrelin, fTimepoint) %>%
  summarise(
    mean_log  = mean(res_logF_AG, na.rm = TRUE),   # y-position on log-scale
    mean_raw  = mean(F_AG,    na.rm = TRUE)    # label text in pg/ml
  , .groups = "drop") %>%
  mutate(
    # map T0/T2/T3 to 1/2/3 and nudge by group to sit beside the raincloud halves
    x = dplyr::recode(fTimepoint, T0 = 1, T2 = 2, T3 = 3),
    #y = ifelse(fGhrelin == "Placebo", x + 0.1, x - 0.1),
    label = ifelse(fTimepoint == "T3", sprintf("%.1f pg/ml", mean_raw), sprintf("%.1f", mean_raw)), 
    y_pos = ifelse(fGhrelin == "Placebo", mean_log - 0.23, mean_log + 0.23)  # +0.05 log-units upwards
)

df_2x3 <- data_2x2(
  array_1 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Placebo"],
  array_2 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Ghrelin"],
  array_3 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Placebo"], 
  array_4 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Ghrelin"],
  array_5 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Placebo"], 
  array_6 = d_hunger$res_logF_AG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Ghrelin"],
  labels = (c('Placebo','Ghrelin')),
  jit_distance = .05,
  jit_seed = 321) 

raincloud_2x3_vertical <- raincloud_2x3_repmes(
  data = df_2x3,
  colors = pal,
  fills = pal,
  size = 2,
  alpha = .5,
  ort = 'v') +

scale_x_continuous(breaks=c(1,2,3), labels=c("T0\nFasted", "T2\nPost-meal", "T3\nPost-infusion"), limits=c(0, 4)) +
  xlab("Timepoint") +
  ylab("") +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "plain", size = 12),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(2, 1, 2, 1)
  ) +
    coord_cartesian(ylim = c(-3, 4), xlim = c(0.7, 4), clip = "off") +  
  geom_label(
    data = means_raw,
    aes(x = x + 0.45,             
        y = y_pos,
        label = label,
        fill = fGhrelin),         
    inherit.aes = FALSE,
    hjust = 0,                   
    fontface = "bold",
    size = 4,
    color = "white",              
    label.size = 0,               
    alpha = 1,                 
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal) # use your group palette


#ggsave(paste(path_out, sub_VAS , "AG_Timepoints_Log.png", sep=""), 
 #             plot = raincloud_2x3_vertical,  height =5.5, width = 9.5, units = "in", dpi = 600, bg = "white")

# Combine with Data from Phenotyping 
# --- 1) Read phenotyping data ---
TUE008_Phen_data <- read.csv(file.path(path_in, "TUE008_PHEN_ghrelin_summary_preprocessed.csv"))

# If fTimepoint is missing, set all to "T0" (fasted baseline)
if (!"fTimepoint" %in% names(TUE008_Phen_data)) {
  TUE008_Phen_data$fTimepoint <- "T0"
}

# Keep only Fasted (T0) for phenotyping baseline; drop NAs
phen_fasted <- TUE008_Phen_data %>%
  dplyr::select(ID, fTimepoint, fMDD, F_AG, res_logF_AG) %>%
  dplyr::filter(is.finite(res_logF_AG), is.finite(F_AG))

# Means for label boxes (log for y-position, raw for text)
means_phen <- phen_fasted %>%
  group_by(fMDD) %>%
  summarise(
    mean_log = mean(res_logF_AG, na.rm = TRUE),
    mean_raw = mean(F_AG,        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x = ifelse(fMDD == "MDD", 1, 1),
    label = sprintf("%.1f pg/ml", mean_raw),
    y_pos = ifelse(fMDD == "HCP", mean_log + 0.18, mean_log - 0.18)
  ) 

library(ggrain)


p_phen <- ggplot(phen_fasted, aes(1, res_logF_AG , fill = fMDD, color = fMDD)) +
  geom_rain(alpha = .5, size = 2,    violin.args  = list(linewidth = 0.25, alpha = 0.5),    # thinner density outline    
            boxplot.args = list(color = "black", outlier.shape = NA)) +
 theme_cowplot(font_size = 12) +
 scale_x_continuous(breaks=c(1), labels=c("T0\nFasted")) +
  xlab("Timepoint") +
  ylab("Plasma acyl ghrelin, res log") +
 
  theme(
    legend.position = "none",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "plain", size = 12),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(2, 1, 2, 1)
  ) +
  scale_fill_manual(values=c("darkorange", "dodgerblue")) +
  scale_color_manual(values=c("darkorange", "dodgerblue")) +
    coord_cartesian(ylim = c(-3, 4), clip = "off") +   # tighter; still room for labels

  # guides(fill = 'none', color = 'none') +
    geom_label(
    data = means_phen,
    aes(x = x + 0.2,             # push to the right of the cloud
        y = y_pos,
        label = label,
        fill = fMDD),         # box color by condition
    inherit.aes = FALSE,
    hjust = 0,                    # left-align text
    fontface = "bold",
    size = 4,
    color = "white",              # text always black
    label.size = 0,               # no border
    alpha = 1,                  # transparency of the fill box
    show.legend = FALSE
  ) 

# Assemble Plots together: 

# Get N 
n_phen <- phen_fasted %>%
  filter(!is.na(res_logF_AG)) %>%
  summarise(N = n_distinct(ID)) %>%
  pull(N)

n_nimg  <- n_distinct(d_hunger$fID)  # or d_hunger$fID

combo_core <- cowplot::plot_grid(
  p_phen,
  raincloud_2x3_vertical ,
  ncol = 2,
  rel_widths = c(0.5, 1.10), 
  labels = c(sprintf("Phenotyping (N = %d)", n_phen),
             sprintf("Neuroimaging (N = %d)", n_nimg)),
  label_size = 14,
  label_fontface = "bold",
  # left-align labels within each panel
  hjust = 0, vjust = 1,
  label_x = c(0.25, 0.25),
  label_y = c(0.99, 0.99)) 


ggsave(paste(path_out, sub_VAS, "AG_Timepoints_Phen.png", sep = ""),
      plot = combo_core, width = 8, height = 4, units = "in", dpi = 600, bg = "white")

# Percentage Change? 

library(emmeans)

# Make sure T0 is the reference
d_hunger$fTimepoint <- relevel(d_hunger$fTimepoint, ref = "T0")
mod <- lmer(logF_AG ~ fGhrelin * fTimepoint + cBMI + cAge + cSex + (1 |fID), data = d_hunger)
summary(mod)
# Estimated marginal means on the log scale, then contrasts vs T0 within each condition
emm <- emmeans(mod, ~ fTimepoint | fGhrelin)

# % change vs T0, within each condition (GMR - 1)*100; type="response" back-transforms via exp
pct_vs_T0 <- emmeans::contrast(emm, "trt.vs.ctrl", ref = "T0") %>%
  summary(type = "response")  %>% # does back to geometric mean!
  transform(pct_change = 100 * (estimate-1),  SE_pct = 100 * SE )  # transform to perc. change
pct_vs_T0



# Descriptive ghrelin value +/- std. deb. 

sum_tp_cond <- d_hunger %>%
  filter(fTimepoint %in% c("T0","T2","T3")) %>%
  group_by(fGhrelin, fTimepoint) %>%
  summarise(
    N         = sum(is.finite(F_AG)),
    mean_pgml = mean(F_AG, na.rm = TRUE),
    sd_pgml   = sd(F_AG,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(fGhrelin, fTimepoint)

sum_tp_cond

# --------------------------------# 
# Same for Desacyl 
# --------------------------------# 

df_2x3 <- data_2x2(
  array_1 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Placebo"],
  array_2 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T0" & d_hunger$fGhrelin == "Ghrelin"],
  array_3 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Placebo"], 
  array_4 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T2" & d_hunger$fGhrelin == "Ghrelin"],
  array_5 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Placebo"], 
  array_6 = d_hunger$res_logF_DG[d_hunger$fTimepoint == "T3" & d_hunger$fGhrelin == "Ghrelin"],
  labels = (c('Placebo','Ghrelin')),
  jit_distance = .05,
  jit_seed = 321) 

# Plot Log with labels pg/mol 
means_raw <- d_hunger %>%
  filter(fTimepoint %in% c("T0","T2","T3"),
         fGhrelin %in% c("Placebo","Ghrelin")) %>%
  group_by(fGhrelin, fTimepoint) %>%
  summarise(
    mean_log  = mean(res_logF_DG, na.rm = TRUE),   # y-position on log-scale
    mean_raw  = mean(F_DG,    na.rm = TRUE)    # label text in pg/ml
  , .groups = "drop") %>%
  mutate(
    # map T0/T2/T3 to 1/2/3 and nudge by group to sit beside the raincloud halves
    x = dplyr::recode(fTimepoint, T0 = 1, T2 = 2, T3 = 3),
    #y = ifelse(fGhrelin == "Placebo", x + 0.1, x - 0.1),
    label = sprintf("%.1f pg/ml", mean_raw), 
    y_pos = ifelse(fGhrelin == "Placebo", mean_log - 0.2, mean_log + 0.2)  # +0.05 log-units upwards
)

raincloud_2x3_vertical <- raincloud_2x3_repmes(
  data = df_2x3,
  colors = pal,
  fills = pal,
  size = 2,
  alpha = .5,
  ort = 'v') +
scale_x_continuous(breaks=c(1,2,3), labels=c("Fasted", "Post-Meal", "Post-Infusion"), limits=c(0, 4)) +
  xlab("Timepoint") +
  ylab("") +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "plain", size = 12),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(2, 1, 2, 1)
  ) +
    coord_cartesian(ylim = c(-3, 4), xlim = c(0.7, 4), clip = "off") +   # tighter; still room for labels
  geom_label(
    data = means_raw,
    aes(x = x + 0.45,             # push to the right of the cloud
        y = y_pos,
        label = label,
        fill = fGhrelin),         # box color by condition
    inherit.aes = FALSE,
    hjust = 0,                    # left-align text
    fontface = "bold",
    size = 4,
    color = "white",              # text always black
    label.size = 0,               # no border
    alpha = 1,                  # transparency of the fill box
    show.legend = TRUE
  ) +
  scale_fill_manual(values = pal) # use your group palette

phen_fasted <- TUE008_Phen_data %>%
  dplyr::select(ID, fTimepoint, fMDD, F_DG, res_logF_DG) %>%
  dplyr::filter(is.finite(res_logF_DG), is.finite(F_DG))

means_phen <- phen_fasted %>%
  group_by(fMDD) %>%
  summarise(
    mean_log = mean(res_logF_DG, na.rm = TRUE),
    mean_raw = mean(F_DG,        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x = ifelse(fMDD == "MDD", 1, 1),
    label = sprintf("%.1f pg/ml", mean_raw),
    y_pos = ifelse(fMDD == "HCP", mean_log + 0.2, mean_log - 0.2)
  ) 


p_phen <- ggplot(phen_fasted, aes(1, res_logF_DG , fill = fMDD, color = fMDD)) +
  geom_rain(alpha = .5, size = 2,    violin.args  = list(linewidth = 0.25, alpha = 0.5),    # thinner density outline    
            boxplot.args = list(color = "black", outlier.shape = NA)) +
 theme_cowplot(font_size = 12) +
 scale_x_continuous(breaks=c(1), labels=c("Fasted")) +
  xlab("Timepoint") +
  ylab("Plasma des-acyl ghrelin, res log") +
 
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "plain", size = 12),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(margin = margin(0.15,0,0.15,0, "cm")),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(2, 1, 2, 1)
  ) +
  scale_fill_manual(values=c("darkorange", "dodgerblue")) +
  scale_color_manual(values=c("darkorange", "dodgerblue")) +
    coord_cartesian(ylim = c(-3, 4), clip = "off") +   # tighter; still room for labels

  # guides(fill = 'none', color = 'none') +
    geom_label(
    data = means_phen,
    aes(x = x + 0.2,             # push to the right of the cloud
        y = y_pos,
        label = label,
        fill = fMDD),         # box color by condition
    inherit.aes = FALSE,
    hjust = 0,                    # left-align text
    fontface = "bold",
    size = 4,
    color = "white",              # text always black
    label.size = 0,               # no border
    alpha = 1,                  # transparency of the fill box
    show.legend = TRUE
  ) 

# Assemble Plots together: 

# Get N 
n_phen <- phen_fasted %>%
  filter(!is.na(res_logF_DG)) %>%
  summarise(N = n_distinct(ID)) %>%
  pull(N)

n_nimg  <- n_distinct(d_hunger$fID)  # or d_hunger$fID

combo_core_des <- cowplot::plot_grid(
  p_phen,
  raincloud_2x3_vertical ,
  ncol = 2,
  rel_widths = c(0.5, 1.10), 
  #labels = c(sprintf("Phenotyping (N = %d)", n_phen),
             #sprintf("Neuroimaging (N = %d)", n_nimg)),
  #label_size = 14,
  #label_fontface = "bold",
  # left-align labels within each panel
  hjust = 0, vjust = 1,
  label_x = c(0.25, 0.25),
  label_y = c(0.99, 0.99)) 

ggsave(paste(path_out, sub_VAS, "DG_Timepoints_Phen.png", sep = ""),
      plot = combo_core_des, width = 8, height = 4, units = "in", dpi = 600, bg = "white")

#

# ----- Ghrelin and Insulin sensitivity
M_NA <- lmer(res_logF_AG_Delta  ~  res_logHOMA_T0 * Hunger_Delta   + cSession.x +  cSex.x + cAge.x + cBMI.x  + (1|fID), d_ghrelin_change_HM)
summary(M_NA)

M_NA <- lmer(res_logF_AG  ~  res_logTyG_T0 * fTimepoint    + cSession +  cSex + cAge + cBMI  + (1|fID), d_hunger)
summary(M_NA)

M_NA <- lmer(res_logF_AG_Delta  ~  res_logTyG_T0 + Hunger_Delta   + cSession.x +  cSex.x + cAge.x + cBMI.x  + (1|fID), d_ghrelin_change_HM)
summary(M_NA)

#######################################################################################
# TEST FOR SUBJECTIVE RATINGS 
####################################################################################### 

contrasts(d_NA$fTimepoint) <-  contr.treatment(levels(d_NA$fTimepoint), base = 1)
contrasts(d_PA$fTimepoint) <-  contr.treatment(levels(d_PA$fTimepoint), base = 1)
contrasts(d_PA_NA$fTimepoint) <-  contr.treatment(levels(d_PA_NA$fTimepoint), base = 1)
contrasts(d_PA_NA$fTimepoint) <-  contr.treatment(levels(d_PA_NA$fTimepoint), base = 4)
contrasts(d_hunger$fTimepoint) <-  contr.treatment(levels(d_hunger$fTimepoint), base = 1)

# Mood 
M_NA <- lmer(VAS_NA_s ~  fGhrelin * fTimepoint + cSession +  cSex + cAge + cBMI  + (1+ fTimepoint +cGhrelin|fID), d_NA)
summary(M_NA)

M_PA <- lmer(VAS_PA_s ~ fGhrelin  * fTimepoint  + cSession  + cAge + cBMI + cSex + (1+ fTimepoint+ cGhrelin|fID), d_PA)
summary(M_PA)

M_VAS_NAPA <- lmer(VAS_PANA ~  fGhrelin * fTimepoint   + cSession + cAge + cBMI + cSex + (1 + fGhrelin |fID), d_PA_NA)
summary(M_VAS_NAPA)

# Post hoc contrasts, only interested in T2-T3 Ghrelin: Hypothesis: increase mood & hunger 
emms <-  emmeans(M_VAS_NAPA, ~ fGhrelin*fTimepoint)
IC_st <- contrast(emms, interaction = "pairwise", adjust = "none", side= ">")
emms <-  emmeans(M_PA, ~ fGhrelin*fTimepoint)
IC_st <- contrast(emms, interaction = "pairwise", adjust = "none", side= ">")
emms <-  emmeans(M_NA, ~ fGhrelin*fTimepoint)
IC_st <- contrast(emms, interaction = "pairwise", adjust = "none", side= "<")

# Now try for next project to contast T0,1,2 vs. T3 
P1 = c(-1/3, -1/3,  -1/3,  1, 0 , 0 ,0 ,0 )

G1 = c(0 , 0 ,0 ,0, -1/3, -1/3, -1/3, 1)

Interaction = G1 - P1

emms1 <-  emmeans(M_VAS_NAPA, ~  fTimepoint * fGhrelin )
contrast(emms1, method = list("G - P (T0,1,2 - T3)" =   Interaction  ), adjust = "none", side= ">")
eff_size(emms1, method = list("G - P (T0,1,2 - T3)" =   Interaction  ), side = ">", sigma = sigma(M_VAS_NAPA), edf = df.residual(M_VAS_NAPA))

emms1 <-  emmeans(M_PA, ~  fTimepoint * fGhrelin )
contrast(emms1, method = list("G - P (T0,1,2 - T3)" =   Interaction  ), adjust = "none", side= ">")


# Metabolic 

M_VAS_Hunger <- lmer(Rating ~ fGhrelin  * fTimepoint  + cSession + cAge + cBMI + cSex + (1 |fID), d_hunger)
summary(M_VAS_Hunger)

M_VAS_Full <- lmer(Rating ~ fGhrelin  * fTimepoint + cSession + cAge + cBMI + cSex + (1 |fID), d_full)
summary(M_VAS_Full)


emms <-  emmeans(M_VAS_Hunger, ~ fGhrelin*fTimepoint)
IC_st <- contrast(emms, interaction = "pairwise", adjust = "none", side= ">")

emms2 <-  emmeans(M_VAS_Full, ~ fGhrelin*fTimepoint)
IC_st2 <- contrast(emms, interaction = "pairwise", adjust = "none", side= "<")

# Write down NAcc BP Model + NAcc BP + Matches 
sjPlot:: tab_model(M_VAS_Hunger, M_VAS_Full,
                   p.val = "satterthwaite",
                   show.re.var=TRUE,
                  dv.labels = c("Hunger", "Fullness"), 
                  file= paste(path_out, sub_VAS, "VAS_Ghrelin", ".doc", sep = ""))

# Tidy contrast tables with CIs and p-values
tab_h <- summary(IC_st, infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>%
  mutate(
    across(where(is.numeric) & !matches("p.value"), ~ round(.x, 2)),
    p.value = sprintf("%.3f", p.value)
  )
# Keep useful columns (adjust names if your factors/columns differ)
pick_cols <- intersect(
  c("fGhrelin_pairwise","fTimepoint_pairwise","estimate","SE","df","t.ratio","p.value","lower.CL","upper.CL"),
  names(tab_h)
)

ft_h <- flextable::regulartable(tab_h %>% select(all_of(pick_cols))) %>% flextable::autofit()

# Write to Word document alongside headings
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "EMMeans contrasts – VAS Hunger", style = "heading 1")
doc <- flextable::body_add_flextable(doc, ft_h)
print(doc, target = file.path(path_out, paste0(sub_VAS, "VAS_Ghrelin_contrasts.docx")))




