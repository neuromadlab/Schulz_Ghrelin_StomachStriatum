###############################################################################
# TUE008 Neuroimaging Study: BOLD Spikes during PET-MR 
# Corinna Schulz, October 205
###############################################################################

# Stats + Plots for Matches/Pulses + Ghrelin (density plots) 
# Sensititivty: checks whether original (.025% results hold true for different settings)

###############################################################################
# (1) SET UP
###############################################################################

# Install and load all required libraries
if (!require("librarian")) install.packages("librarian")
librarian::shelf(readxl, lme4, lmerTest, 
                 ggplot2, scales, ggpubr, emmeans, cowplot, viridis, tidybayes, dplyr, 
                 tidyr, dplyr, readr) 

#Set all themes
theme_set(theme_cowplot(font_size = 12))

# Set Colors
color_Placebo <- "darkblue"  
color_Ghrelin <- "darkgoldenrod"  

# Set all Paths
setwd(getwd())
path_in <- "./input/" 
path_out <- "./output/" 

sub_Pulse <- "Pulse/"

# Load Condition file, with Matching of Conn and Study ID 
d_conds <- read_excel(paste(path_in,"T_COV_Tabelle.xlsx", sep = ""))

# Load Pulse Detection Results 
# ATTENTION: Subject ID is CONN ID, Phase too, don't just merge! 
PulseCounts <- read_excel(paste(path_in,"PulseCount_allROIs.xlsx", sep = ""))
PulseCooccur <-  read_excel(paste(path_in,"PulseCooccurance_Clean_Hypothalamus-NAcc.xlsx", sep = ""))

PulseCooccur$fBaseCond <- factor(PulseCooccur$BaseCond)

# Load PET DATA 
d_BPs <- read_excel(paste(path_in,"TUE008_NIMG_BPs_ROIs_28_04_25.xlsx", sep = ""))
d_BPs$fID <- factor(d_BPs$ID)

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


PulseCounts_Thr22 <- read_excel(paste(path_in,"PulseCount_allROIs_Thresh22.xlsx", sep = ""))
PulseCounts_Thr28 <- read_excel(paste(path_in,"PulseCount_allROIs_Thresh28.xlsx", sep = ""))
PulseCounts_Thr30 <- read_excel(paste(path_in,"PulseCount_allROIs_Thresh3.xlsx", sep = ""))
PulseCounts_Thr32 <- read_excel(paste(path_in,"PulseCount_allROIs_Thresh32.xlsx", sep = ""))
PulseCounts_Thr34 <- read_excel(paste(path_in,"PulseCount_allROIs_Thresh34.xlsx", sep = ""))

# Load Mean Effort for Pulses - BP - Effort 
TUE008_IMT_MeanEffort <- read_csv(paste(path_out, "instrumental_motivation/TUE008_IMT_MeanEffort.csv", sep = ""))

####################################
# Do Analysis: Ghrelin effect + IMT effect 

# Take 025 as real one, sensitivity for some other thresholds to see robustness of results
idx_thres = c("025","028","030")

# Save all relevant STATS
models_list_PulseCount <- list()  
models_list_MatchCount <- list()  
models_list_PulseCount_ACROSS <- list()  
models_list_MatchCount_ARCOSS <- list()  



count = 0 
for (idx in idx_thres) {

    count = count+1 

    if(count ==1){ 
      PulseCounts <- read_excel(paste(path_in,"PulseCount_allROIs.xlsx", sep = ""))
      PulseCooccur <-  read_excel(paste(path_in,"PulseCooccurance_Clean_Hypothalamus-NAcc.xlsx", sep = ""))

    }else if (count == 2) {
       PulseCounts <- PulseCounts_Thr28
        PulseCooccur <-  read_excel(paste(path_in,"PulseCooccurance_Clean_28_Hypothalamus-NAcc.xlsx", sep = ""))

    }else if (count == 3) {
      PulseCounts <- PulseCounts_Thr30
      PulseCooccur <-  read_excel(paste(path_in,"PulseCooccurance_Clean_30_Hypothalamus-NAcc.xlsx", sep = ""))
    }


# Prepare Pulse Count 

    PulseCounts <- PulseCounts %>% rename(ID = Subject, Ghrelin = Condition) 
    PulseCounts <- PulseCounts %>% mutate(ID = as.numeric(ID))
    PulseCounts$Ghrelin <- factor(PulseCounts$Ghrelin) 

    d_conds_Pulses <- d_conds %>% filter(Phase == 1) %>% select("ID","fID", "Ghrelin","Session", "Age","Sex","BMI","Hunger_delta","Food_Sensitivity","Reward_Sensitivity")
    d_conds_Pulses$Ghrelin <- factor(d_conds_Pulses$Ghrelin, labels = c("Placebo", "Ghrelin")) 

    PulseCounts_Cov <- merge(PulseCounts, d_conds_Pulses, by = c("ID", "Ghrelin")) 
    PulseCounts_Cov$cAge <- PulseCounts_Cov$Age - mean(PulseCounts_Cov$Age)
    PulseCounts_Cov$cSex <- PulseCounts_Cov$Sex - mean(PulseCounts_Cov$Sex)
    PulseCounts_Cov$cBMI <- PulseCounts_Cov$BMI - mean(PulseCounts_Cov$BMI)
    PulseCounts_Cov$cSession <- PulseCounts_Cov$Session - mean(PulseCounts_Cov$Session)
    PulseCounts_Cov$fGhrelin <- relevel(factor(PulseCounts_Cov$Ghrelin), ref = "Placebo")
    PulseCounts_Cov$fBaseCond <- relevel(factor(PulseCounts_Cov$BaseCond), ref = "Baseline")
    PulseCounts_Cov$fROI <- (factor(PulseCounts_Cov$ROI))

    # Filter For NAcc (most plausible given animal data)
    PulseCounts_Cov_ROI <- PulseCounts_Cov %>% filter(fROI == "NAcc")

    PulseCount_NAcc_Phases_Cov <- merge(PulseCounts_Cov_ROI,d_BPs_long, by= c("fID", "cSession" )) 
    PulseCount_NAcc_Phases_Cov <- PulseCount_NAcc_Phases_Cov %>% filter(ROI.y == "M_Accumbens")
    PulseCount_NAcc_Phases_Cov$cBP <- PulseCount_NAcc_Phases_Cov$BP - mean(PulseCount_NAcc_Phases_Cov$BP)
    PulseCount_NAcc_Phases_Cov$fBaseCond <- relevel(factor(PulseCount_NAcc_Phases_Cov$fBaseCond), ref = "Baseline")
    PulseCount_NAcc_Phases_Cov$fBaseCond <- factor(PulseCount_NAcc_Phases_Cov$fBaseCond, levels = c("Baseline","TaskFree","Task"))

# Prepare Coocurance 

  # Match with Covariates (Attention CONN ID!) 
  PulseCooccur$fBaseCond <- factor(PulseCooccur$BaseCond)
  PulseCooccur <- PulseCooccur %>% rename(ID = Subject, Ghrelin = Condition) 
  PulseCooccur <-    PulseCooccur %>% mutate(ID = as.numeric(ID))
  PulseCooccur$Ghrelin <- factor(PulseCooccur$Ghrelin) 

  d_conds_Pulses <- d_conds %>% filter(Phase == 1) %>% select("ID","fID", "Ghrelin","Session", "Age","Sex","BMI","Hunger_delta")
  d_conds_Pulses$Ghrelin <- factor(d_conds_Pulses$Ghrelin, labels = c("Placebo", "Ghrelin")) 
  d_conds_Pulses$Ghrelin <- factor(d_conds_Pulses$Ghrelin, labels = c("Placebo", "Ghrelin")) 

  PulseCooccur_Cov_corrected <- merge(PulseCooccur, d_conds_Pulses, by = c("ID", "Ghrelin")) 
  PulseCooccur_Cov_corrected$cAge <- PulseCooccur_Cov_corrected$Age - mean(PulseCooccur_Cov_corrected$Age)
  PulseCooccur_Cov_corrected$cSex <- PulseCooccur_Cov_corrected$Sex - mean(PulseCooccur_Cov_corrected$Sex)
  PulseCooccur_Cov_corrected$cBMI <- PulseCooccur_Cov_corrected$BMI - mean(PulseCooccur_Cov_corrected$BMI)
  PulseCooccur_Cov_corrected$cSession <- PulseCooccur_Cov_corrected$Session - mean(PulseCooccur_Cov_corrected$Session)
  PulseCooccur_Cov_corrected$fGhrelin <- relevel(factor(PulseCooccur_Cov_corrected$Ghrelin), ref = "Placebo")

  d_BPs_NAcc <-  d_BPs_long  %>% filter(ROI == "M_Accumbens") # fID Is real ID

  PulseCooccur_Cov_corrected_PET <- merge(PulseCooccur_Cov_corrected, d_BPs_NAcc, by = c("fID", "Session") )
  PulseCooccur_Cov_corrected_PET$cBP <- PulseCooccur_Cov_corrected_PET$BP - mean(PulseCooccur_Cov_corrected_PET$BP)


  #####################################
  # Descriptives Plots/ Summary of Counts 
  # Test PET 
  # Test Effort 
  # Ghrelin Effect 
  #####################################

  color_baseline <- "gray40"
  color_taskfree <- "#2E8B57"   # darkseagreen4 / dark green
  color_task     <- "#B22222"   # firebrick / dark red

  # Table 
  pulse_summary <- PulseCount_NAcc_Phases_Cov %>%
    group_by(fGhrelin, fBaseCond) %>%      # or use Pair, fGhrelin if you prefer
    summarise(
      mean_count = mean(PulseCount, na.rm = TRUE),
      var_count  = var(PulseCount, na.rm = TRUE),
      range_min  = min(PulseCount, na.rm = TRUE),
      range_max  = max(PulseCount, na.rm = TRUE),
      n          = n()
    ) 

  # CLEAN LABEL 
  PulseCount_NAcc_Phases_Cov <- PulseCount_NAcc_Phases_Cov %>%
    mutate(fBaseCond = factor(fBaseCond,
                              levels = c("Baseline", "TaskFree", "Task")), 
          fGhrelin = factor(fGhrelin, 
                              levels = c("Placebo","Ghrelin")))


  base_tbl <- PulseCount_NAcc_Phases_Cov %>%
    filter(fBaseCond == "Baseline") %>%
    group_by(fID, fGhrelin, BP) %>%
    summarise(baseline_val = sum(PulseCount, na.rm = TRUE), .groups = "drop")

  PulseCounts_Deltas <- PulseCount_NAcc_Phases_Cov %>%
    left_join(base_tbl, by = c("fID","fGhrelin","BP")) %>%
    mutate(Pulse_Delta = PulseCount - baseline_val) %>% filter(fBaseCond != "Baseline")%>%
    group_by(fID, fGhrelin, BP, fBaseCond)  %>% 
    summarise(Pulse_Delta = sum(Pulse_Delta, na.rm = TRUE))

  write.csv(PulseCounts_Deltas, paste(path_out, sub_Pulse, "Deltas_Pulse", idx, ".csv", sep = ""), row.names=FALSE)

  # Test Pulse Count and PET, and Ghrelin 

  m_rate_PET <- glmer(
    PulseCount ~ fGhrelin * fBaseCond + cBP * fBaseCond  + cSession + cAge + cSex + cBMI + (1 + fGhrelin| fID),
    data = PulseCount_NAcc_Phases_Cov, 
    family = poisson(link = "log")) 
  summary(m_rate_PET)  
  print(paste("Tested Pulse Count:", idx))

  # save model 
  models_list_PulseCount[count] <- list(m_rate_PET) 

  # Test 3-way interactions (=> does not improve, stick with BP general effect)
  m_rate_PET_Int <- glmer(
    PulseCount ~ fGhrelin * fBaseCond * cBP   + cSession + cAge + cSex + cBMI + (1 + fGhrelin| fID),
    data = PulseCount_NAcc_Phases_Cov, 
    family = poisson(link = "log")) 
  summary(m_rate_PET_Int)  

  anova(m_rate_PET,m_rate_PET_Int)

  # Test Overdispersion; OKAY
  #library(DHARMa)
  #sim <- simulateResiduals(m_rate_PET, n = 1000)  
  #plot(sim)
  #testDispersion(sim)
  #simulationOutput <- simulateResiduals(fittedModel = m_rate_PET)
  #plot(simulationOutput)


  ## Test one across the phases (most sensible, BP is only across entire session) 
  PulseCount_NAcc_Phases_Cov_ACROSS <- PulseCount_NAcc_Phases_Cov %>% 
        mutate(fPhase_dicho = case_when( Phase == 1 ~ "Intervention",
                                        Phase == 2 ~ "Baseline",
                                        Phase == 3 ~ "Intervention",
                                        Phase == 4 ~ "Intervention", 
                                        Phase == 5 ~ "Intervention"))

  PulseCount_NAcc_Phases_Cov_ACROSS$fPhase_dicho <- relevel(factor(PulseCount_NAcc_Phases_Cov_ACROSS$fPhase_dicho), ref = "Baseline") # 2 is Baseline in Conn lingo 

    m_rate_PET_ACROSS <- glmer(
    PulseCount ~ fGhrelin * fPhase_dicho + cBP * fPhase_dicho  + cSession + cAge + cSex + cBMI + (1 + fGhrelin| fID),
    data = PulseCount_NAcc_Phases_Cov_ACROSS, 
    family = poisson(link = "log")) 
  summary(m_rate_PET_ACROSS)  
 
  # get the p-value for the interaction fPhase_dichoIntervention:cBP
  s <- summary(m_rate_PET_ACROSS)
  ct <- coef(s)  
  pval_pulse_BP <- round(ct["fPhase_dichoIntervention:cBP", grep("^Pr\\(", colnames(ct), value = TRUE)],3)

  print(paste("Tested Pulse Count ACROSS:", idx))

  # save model 
  models_list_PulseCount_ACROSS[count] <- list(m_rate_PET_ACROSS) 




  # Density Plots 
  Density_pulses <- ggplot(PulseCount_NAcc_Phases_Cov,
        aes(x = PulseCount, fill = fBaseCond, color = fBaseCond)) +
    facet_wrap(~fGhrelin) + 

    geom_density(alpha = 0.3, linewidth = 1.0, adjust = 1.2) +
  scale_fill_manual(
      name = "Phase",
      values = c(color_baseline, color_taskfree, color_task)
    ) +  
    scale_color_manual(
      name = "Phase",
      values = c(color_baseline, color_taskfree, color_task)
    ) +
    labs(
      x = "Pulse count (NAcc)",
      y = "Density"
    ) +
  theme(legend.position = "none",text = element_text(face = 'bold',size = 14),
              axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 

  Density_pulses_SplitGhrelin <- ggplot(PulseCount_NAcc_Phases_Cov,
        aes(x = PulseCount, fill = fGhrelin, color = fGhrelin)) +
    facet_wrap(~fBaseCond) + 

    geom_density(alpha = 0.3, linewidth = 1.0, adjust = 1.2) +
  scale_fill_manual(
      name = "Condition",
      values = c(color_Placebo, color_Ghrelin)
    ) +  
    scale_color_manual(
      name = "Condition",
      values = c(color_Placebo, color_Ghrelin)
    ) +
    labs(
      x = "Pulse count (NAcc)",
      y = "Density"
    ) +
  theme(legend.position = "none",text = element_text(face = 'bold',size = 14),
              axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


  # Density Delta Plots 
  PulseCounts_Deltas <- PulseCounts_Deltas %>%
    mutate(fBaseCond = factor(fBaseCond,
                              levels = c( "TaskFree", "Task"))) 

  Density_pulsesDelta <- ggplot(PulseCounts_Deltas,
        aes(x = Pulse_Delta, fill = fGhrelin, color = fGhrelin)) +
    facet_wrap(~fBaseCond) + 
    geom_density(alpha = 0.3, linewidth = 1.0, adjust = 1) +
      scale_fill_manual(
      name = "Condition",
      values = c(color_Placebo, color_Ghrelin)
    ) +  
    scale_color_manual(
      name = "Condition",
      values = c( color_Placebo, color_Ghrelin)
    ) +
    labs(
      x = "Δ Pulse count (NAcc)",
      y = "Density"
    ) +
      theme(legend.position = "none",text = element_text(face = 'bold',size = 14),
              axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0))
      

  ## FINAL PLOT EMMEANS 
  # Plot Ratio estimates 

      dat <- PulseCount_NAcc_Phases_Cov   # <- data used to fit m_rate_PET
      dat$fBaseCond <- factor(dat$fBaseCond, levels = c("Baseline","TaskFree","Task"))

      # 1) Build a modest grid over cBP (fewer points -> cleaner x-axis)
      cBP_seq <- seq(min(dat$cBP, na.rm = TRUE),
                  max(dat$cBP, na.rm = TRUE),
                  length.out = 20)

      # 2) EMMs of phase by cBP (averaged over other covariates)
      summary(m_rate_PET)
      emm <- emmeans(
      m_rate_PET,
      specs = ~ fBaseCond | cBP,
      at    = list(cBP = cBP_seq),
      type  = "response",
      cov.reduce = mean
      )

      # 3) Phase-vs-Baseline rate ratios at each cBP
      rr <- contrast(
      emm,
      method = list(
          "Task Free / Baseline" = c(-1, 1, 0), 
          "Task / Baseline"      = c(-1, 0, 1)),
      by   = "cBP",
      type = "response"
      )

      # 4) Get CIs, harmonize column names
      rr_df <- summary(rr, infer = TRUE, level = 0.95) %>% as.data.frame()


      f_to_num <- function(x) {
      if (is.factor(x)) as.numeric(levels(x))[x] else as.numeric(x)
      }

      rr_df <- rr_df %>%
      mutate(
          cBP      = f_to_num(cBP),                 # or BP_raw if you created it
          ratio    = f_to_num(get(if ("ratio" %in% names(rr_df)) "ratio" else "estimate")),
          lower.CL = f_to_num(if ("lower.CL" %in% names(rr_df)) lower.CL else asymp.LCL),
          upper.CL = f_to_num(if ("upper.CL" %in% names(rr_df)) upper.CL else asymp.UCL),
          contrast = factor(contrast, levels = c( "Task Free / Baseline", "Task / Baseline"))
      )

      p_rr_bp <- ggplot(rr_df, aes(x = cBP, y = ratio)) +
      geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
      geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
                  alpha = 0.2, fill = "gray70", color = NA) +
      geom_line(linewidth = 1.2, color = "black") +
      facet_wrap(~ contrast, nrow = 1) +
      scale_x_continuous(breaks = pretty_breaks(n = 5),
                          labels = label_number(accuracy = 0.05)) +
      scale_y_continuous(trans = "log10",
                          breaks = c(0, 0.5, 0.75, 1, 1.5, 2, 2.5, 3,4),
                          labels = label_number(accuracy = 0.1)) +
      labs(x = "NAcc BP (centered)",
          y = "Count ratio (log scale)",
          title = sprintf("Model-adjusted effects, p = %.3f",pval_pulse_BP)) + 
      theme(legend.position = "none",text = element_text(face = 'bold',size = 14.0),
                  axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


      # Print Stats as well     
      # 1) Simple slopes of BP (cBP) within each phase on the log (link) scale
      sl <- emtrends(
      m_rate_PET,
      ~ fBaseCond,
      var = "cBP",
      type = "link",        # slope on log scale (β per 1-unit cBP)
      cov.reduce = mean
      )

      # 2) BP × Phase interaction = differences in BP slopes between phases
      #    (Task vs Baseline, TaskFree vs Baseline, Task vs TaskFree)
      int_tests <- contrast(
      sl,
      method = list(
          "TaskFree vs Baseline"   = c(-1, 1, 0), 
          "Task vs Baseline"       = c(-1, 0, 1)), 
      infer = TRUE
      ) %>%
      as.data.frame() %>%
      mutate(
          int_label = sprintf("β = %.3f (p = %.3f)",
                              estimate, p.value))


      # Choose x/y positions from your rr_df used to build p_rr_bp
      # If you plotted against BP_raw, use that; else cBP.
      x_var <- if ("BP_raw" %in% names(rr_df)) "BP_raw" else "cBP"
      y_hi  <- if ("upper.CL" %in% names(rr_df)) "upper.CL" else if ("asymp.UCL" %in% names(rr_df)) "asymp.UCL" else NA

      pos_df <- rr_df %>%
          group_by(contrast) %>%
          summarise(
          x = min(.data[[x_var]], na.rm = TRUE) + 0.02 * diff(range(.data[[x_var]], na.rm = TRUE)),
          y = max(.data[[y_hi]],   na.rm = TRUE) * 0.95,
          .groups = "drop"
          ) 

      # Build labels for Task / Baseline and Task Free / Baseline facets
      lab_map <- int_tests %>%
          mutate(
          contrast = recode(contrast,
                              "TaskFree vs Baseline" = "Task Free / Baseline",
                              "Task vs Baseline" = "Task / Baseline"),                
          label = int_label
          )%>%
          mutate(contrast = factor(contrast, levels = c( "Task Free / Baseline", "Task / Baseline")))


      ann_df <- left_join(pos_df, lab_map, by = "contrast")

      
      p_rr_bp_stats <- p_rr_bp +
          geom_text(data = ann_df,
                  aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.8)

  ###############################################
  ### PULSE MATCHES 
  ###############################################

  # Summarize pulse counts per ROI
  PulseCooccur %>%
    group_by(Ghrelin, fBaseCond) %>%  # or ROI, fGhrelin, etc.
    summarise(
      mean_count = mean(Matches, na.rm = TRUE),
      var_count  = var(Matches, na.rm = TRUE),
      range_min  = min(Matches, na.rm = TRUE),
      range_max  = max(Matches, na.rm = TRUE),
      n = n()
    ) 

  PulseCooccur <- PulseCooccur %>%
    mutate(fBaseCond = factor(fBaseCond,
                              levels = c("Baseline", "TaskFree", "Task")), 
          fGhrelin = factor(Ghrelin, 
                              levels = c("Placebo","Ghrelin")))


  base_tbl <- PulseCooccur_Cov_corrected_PET %>%
    filter(fBaseCond == "Baseline") %>%
    group_by(fID, fGhrelin, BP) %>%
    summarise(baseline_val = sum(Matches, na.rm = TRUE), .groups = "drop")

  PulseMatches_Deltas <- PulseCooccur_Cov_corrected_PET %>%
    left_join(base_tbl, by = c("fID","fGhrelin","BP")) %>%
    mutate(Matches_Delta = Matches - baseline_val) %>% filter(fBaseCond != "Baseline")%>%
    group_by(fID, fGhrelin, fBaseCond, BP)  %>% 
    summarise(Matches_Delta = sum(Matches_Delta, na.rm = TRUE), .groups = "drop")

  write.csv(PulseMatches_Deltas, paste(path_out, sub_Pulse, "Deltas_Matches" ,idx, ".csv", sep=""), row.names=FALSE)

  # Test 
  PulseCooccur_Cov_corrected_PET$fBaseCond <- factor(PulseCooccur_Cov_corrected_PET$fBaseCond, levels = c("Baseline","TaskFree","Task"))
  m_rate_matches <- glmer(
    Matches ~   fGhrelin * fBaseCond + cBP * fBaseCond  + cSession.x + cAge + cSex + cBMI + (1 + fGhrelin| fID),
    data = PulseCooccur_Cov_corrected_PET, 
    family = poisson(link = "log")
  )
  summary(m_rate_matches)
  print(paste("Tested Match Count:", idx))

  # Save stats 
  models_list_MatchCount[count] <- list(m_rate_matches) 



  ## Test one across the phases (most sensible, BP is only across entire session) 
  PulseCooccur_Cov_corrected_PET_ACROSS <- PulseCooccur_Cov_corrected_PET %>% 
        mutate(fPhase_dicho = case_when( Phase == 1 ~ "Intervention",
                                        Phase == 2 ~ "Baseline",
                                        Phase == 3 ~ "Intervention",
                                        Phase == 4 ~ "Intervention", 
                                        Phase == 5 ~ "Intervention"))

  PulseCooccur_Cov_corrected_PET_ACROSS$fPhase_dicho <- relevel(factor(PulseCooccur_Cov_corrected_PET_ACROSS$fPhase_dicho), ref = "Baseline") # 2 is Baseline in Conn lingo 

    m_rate_matches_ACROSS <- glmer(
    Matches ~   fGhrelin * fPhase_dicho + cBP * fPhase_dicho  + cSession.x + cAge + cSex + cBMI + (1 + fGhrelin| fID),
    data = PulseCooccur_Cov_corrected_PET_ACROSS, 
    family = poisson(link = "log")) 
  summary(m_rate_matches_ACROSS)  
 
  # get the p-value for the interaction fPhase_dichoIntervention:cBP
  s <- summary(m_rate_matches_ACROSS)
  ct <- coef(s)  
  pval_match_BP <- round(ct["fPhase_dichoIntervention:cBP", grep("^Pr\\(", colnames(ct), value = TRUE)],3)

  print(paste("Tested Matches  ACROSS:", idx))

  # save model 
  models_list_MatchCount_ARCOSS[count] <- list(m_rate_matches_ACROSS) 


  Density_matches <- ggplot(PulseCooccur,
        aes(x = Matches, fill = fBaseCond, color = fBaseCond)) +
        geom_density(alpha = 0.3, linewidth = 1.0, adjust = 1.2) +
          facet_wrap(~fGhrelin) + 
        scale_fill_manual(
      name = "Phase",
      values = c(color_baseline, color_taskfree, color_task)) + 
      scale_color_manual(
      name = "Phase",
      values = c(color_baseline, color_taskfree, color_task)
    ) +
        labs(
      x = "Pulse matches (Hypo-NAcc)",
      y = "Density") +
  theme(legend.position = "right",text = element_text(face = 'bold',size = 14),
              axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 

  PulseMatches_Deltas <- PulseMatches_Deltas %>%
    mutate(fBaseCond = factor(fBaseCond,
                              levels = c( "TaskFree", "Task")))

  Density_matches_Delta <- ggplot(PulseMatches_Deltas,
        aes(x = Matches_Delta, fill = fGhrelin, color = fGhrelin)) +
        geom_density(alpha = 0.3, linewidth = 1.0, adjust = 1) +
      facet_wrap(~ fBaseCond) +
      scale_fill_manual(
        name = "Condition",
        values = c( color_Placebo, color_Ghrelin), 
        breaks = c("Placebo","Ghrelin"),
        labels = c("Saline","Ghrelin")) + 
      scale_color_manual(
        name = "Condition",
        values = c( color_Placebo, color_Ghrelin), 
        breaks = c("Placebo","Ghrelin"),
        labels = c("Saline","Ghrelin")) + 
        labs(
          x = "Δ Pulse matches (Hypo-NAcc)",
          y = "Density") +
          theme(legend.position = "right",text = element_text(face = 'bold',size = 14),
              axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


  #> combine plots
  cowplot::plot_grid(Density_pulses, Density_matches,
            labels=c("", ""),
            ncol = 2,
            nrow = 1,
            rel_widths= c(1, 1.3),
            label_size = 18)

  # save figure
  ggsave(paste(path_out,sub_Pulse, "PulseMatches_Descriptives", idx, ".png", sep=""), 
    height = 3, width = 10, units = "in", dpi = 600, bg = "white")

  #> combine plots
  cowplot::plot_grid(Density_pulsesDelta, Density_matches_Delta,
            labels=c("", ""),
            ncol = 2,
            nrow = 1,
            rel_widths= c(1, 1.3),
            label_size = 18)

  # save figure
  ggsave(paste(path_out, sub_Pulse,"PulseMatches_Delta_Descriptives", idx, ".png", sep=""), 
    height = 2.5, width = 10, units = "in", dpi = 600, bg = "white")


  ## FINAL PLOT EMMEANS 
  # Plot Ratio estimates 

      dat <- PulseCooccur_Cov_corrected_PET   # <- data used to fit m_rate_matches
      dat$fBaseCond <- factor(dat$fBaseCond, levels = c("Baseline","TaskFree","Task"))

      # 1) Build a modest grid over cBP (fewer points -> cleaner x-axis)
      cBP_seq <- seq(min(dat$cBP, na.rm = TRUE),
                  max(dat$cBP, na.rm = TRUE),
                  length.out = 20)

      # 2) EMMs of phase by cBP (averaged over other covariates)
      summary(m_rate_matches)
      emm <- emmeans(
      m_rate_matches,
      specs = ~ fBaseCond | cBP,
      at    = list(cBP = cBP_seq),
      type  = "response",
      cov.reduce = mean
      )

      # 3) Phase-vs-Baseline rate ratios at each cBP
      rr <- contrast(
      emm,
      method = list(
          "Task Free / Baseline" = c(-1, 1, 0), 
          "Task / Baseline"      = c(-1, 0, 1)),
      by   = "cBP",
      type = "response"
      )

      # 4) Get CIs, harmonize column names
      rr_df <- summary(rr, infer = TRUE, level = 0.95) %>% as.data.frame()


      f_to_num <- function(x) {
      if (is.factor(x)) as.numeric(levels(x))[x] else as.numeric(x)
      }

      rr_df <- rr_df %>%
      mutate(
          cBP      = f_to_num(cBP),                 
          ratio    = f_to_num(get(if ("ratio" %in% names(rr_df)) "ratio" else "estimate")),
          lower.CL = f_to_num(if ("lower.CL" %in% names(rr_df)) lower.CL else asymp.LCL),
          upper.CL = f_to_num(if ("upper.CL" %in% names(rr_df)) upper.CL else asymp.UCL),
          contrast = factor(contrast, levels = c( "Task Free / Baseline", "Task / Baseline"))
      )

      p_rr_bp_matches <- ggplot(rr_df, aes(x = cBP, y = ratio)) +
      geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
      geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
                  alpha = 0.2, fill = "gray70", color = NA) +
      geom_line(linewidth = 1.2, color = "black") +
      facet_wrap(~ contrast, nrow = 1) +
      scale_x_continuous(breaks = pretty_breaks(n = 5),
                          labels = label_number(accuracy = 0.05)) +
      scale_y_continuous(trans = "log10",
                          breaks = c(0, 0.5, 0.75, 1, 1.5, 2, 2.5, 3,4),
                          labels = label_number(accuracy = 0.1)) +
      labs(x = "NAcc BP (centered)",
          y = "Matches ratio (log scale)",
          title = sprintf("Model-adjusted effects, p = %.3f",pval_match_BP)) + 
      theme(legend.position = "none",text = element_text(face = 'bold',size = 14.0),
                  axis.text = element_text(face = 'plain',size = 12.0),axis.text.x = element_text(size = 12.0)) 


      # Print Stats as well     
      # 1) Simple slopes of BP (cBP) within each phase on the log (link) scale
      sl <- emtrends(
      m_rate_matches,
      ~ fBaseCond,
      var = "cBP",
      type = "link",        # slope on log scale (β per 1-unit cBP)
      cov.reduce = mean
      )

      # 2) BP × Phase interaction = differences in BP slopes between phases
      #    (Task vs Baseline, TaskFree vs Baseline, Task vs TaskFree)
      int_tests <- contrast(
      sl,
      method = list(
          "TaskFree vs Baseline"   = c(-1, 1, 0), 
          "Task vs Baseline"       = c(-1, 0, 1)), 
      infer = TRUE
      ) %>%
      as.data.frame() %>%
      mutate(
          int_label = sprintf("β = %.3f (p = %.3f)",
                              estimate, p.value))


      # Choose x/y positions from your rr_df used to build p_rr_bp
      # If you plotted against BP_raw, use that; else cBP.
      x_var <- if ("BP_raw" %in% names(rr_df)) "BP_raw" else "cBP"
      y_hi  <- if ("upper.CL" %in% names(rr_df)) "upper.CL" else if ("asymp.UCL" %in% names(rr_df)) "asymp.UCL" else NA

      pos_df <- rr_df %>%
          group_by(contrast) %>%
          summarise(
          x = min(.data[[x_var]], na.rm = TRUE) + 0.02 * diff(range(.data[[x_var]], na.rm = TRUE)),
          y = max(.data[[y_hi]],   na.rm = TRUE) * 0.95,
          .groups = "drop"
          ) 

      # Build labels for Task / Baseline and Task Free / Baseline facets
      lab_map <- int_tests %>%
          mutate(
          contrast = recode(contrast,
                              "TaskFree vs Baseline" = "Task Free / Baseline",
                              "Task vs Baseline" = "Task / Baseline"),                
          label = int_label
          )%>%
          mutate(contrast = factor(contrast, levels = c( "Task Free / Baseline", "Task / Baseline")))


      ann_df <- left_join(pos_df, lab_map, by = "contrast")

      
      p_rr_bp_matches_stats <- p_rr_bp_matches +
          geom_text(data = ann_df,
                  aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.8)


  cowplot::plot_grid( p_rr_bp_stats, p_rr_bp_matches_stats,
            
            ncol = 2,
            nrow = 1,
            label_size = 18,
            rel_widths = c(1,1.3),
            hjust = -1,
            scale = .95,
            align ="v")

  ggsave( paste(path_out, sub_Pulse, "Pulses_Cooc_BP_EMMs_", idx, ".png", sep = ""), 
          height = 4, width = 10, units = "in", dpi = 600, bg = "white")
  
}

# Write all relevant models (loop doesnt work with tab_model)

# Pulses 
sjPlot::tab_model(models_list_PulseCount[1], 
                   show.re.var=TRUE,
                  dv.labels = c("Pulses"), 
                   use.viewer = TRUE, 
                    show.se = TRUE,      
                  file= paste(path_out, sub_Pulse, "Model_Pulse_", idx_thres[1], ".html", sep = ""))

# Pulses 
sjPlot::tab_model(models_list_PulseCount_ACROSS[1], 
                   show.re.var=TRUE,
                  dv.labels = c("Pulses"), 
                  use.viewer = TRUE, 
                    show.se = TRUE,      
                  file= paste(path_out, sub_Pulse, "Model_Pulse_", idx_thres[1], "_ACROSS.html", sep = ""))


sjPlot::tab_model(models_list_PulseCount[2], 
                   show.re.var=TRUE,
                  dv.labels = c("Pulses"), 
                  use.viewer = TRUE, 
                    show.se = TRUE,      
                  file= paste(path_out, sub_Pulse, "Model_Pulse_", idx_thres[2], ".html", sep = ""))

sjPlot::tab_model(models_list_PulseCount[3], 
                   show.re.var=TRUE,
                  dv.labels = c("Pulses"), 
                   use.viewer = TRUE, 
                    show.se = TRUE,      
                  file= paste(path_out, sub_Pulse, "Model_Pulse_", idx_thres[3], ".html", sep = ""))





# Matches 
sjPlot::tab_model(models_list_MatchCount[1], 
                   show.re.var=TRUE,
                  dv.labels = c("Matches"), 
                  use.viewer = FALSE, 
                  show.se = TRUE,      
                  file= paste(path_out, sub_Pulse, "Model_Match_", idx_thres[1], ".html", sep = ""))

sjPlot::tab_model(models_list_MatchCount_ARCOSS[1], 
                   show.re.var=TRUE,
                  dv.labels = c("Matches"), 
                  use.viewer = FALSE, 
                  show.se = TRUE,      
                   file= paste(path_out, sub_Pulse, "Model_Match_", idx_thres[1], "_ACROSS.html", sep = ""))



sjPlot::tab_model(models_list_MatchCount[2], 
                   show.re.var=TRUE,
                  dv.labels = c("Matches"), 
                  use.viewer = FALSE, 
                  show.se = TRUE, 
                  file= paste(path_out, sub_Pulse, "Model_Match_", idx_thres[2], ".html", sep = ""))

sjPlot::tab_model(models_list_MatchCount[3], 
                   show.re.var=TRUE,
                  dv.labels = c("Matches"), 
                  use.viewer = FALSE, 
                  show.se = TRUE, 
                  file= paste(path_out, sub_Pulse, "Model_Match_", idx_thres[3], ".html", sep = ""))

