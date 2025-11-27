
######################################################
# TUE008 PET-MR Study 
# Corinna Schulz, corinna.schulz96@gmail.com 
# 2024 
# Subjective Sideeffect reports from ghrelin infusion
#######################################################
## Following tutorial here: https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a 


if (!require("librarian")) install.packages("librarian")
librarian::shelf(RColorBrewer, wordcloud, wordcloud2, tm, MASS, readr,
ggplot2, ggpubr, cowplot, gridExtra, viridis, tidybayes, dplyr, httpgd, 
languageserver, dabestr, readxl, vcd, quanteda, fastai,quanteda.textplots, quanteda.textstats,
simstudy, tidyr, ggalluvial, likert,patchwork,gt, gtsummary,flextable)

#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
#install.packages("ggwordcloud")  # if not installed
library(ggwordcloud)
library(dplyr)
library(tm)
library(tidyr)
library(ggplot2)

# Set Colors
color_Placebo <- "darkblue"  
color_Ghrelin <- "darkgoldenrod"  

# Set all Paths
setwd(getwd())
path_in <- "./input/" 
path_out <- "./output/" 

sub_effects <-  "sideeffects/"

if (file.exists(paste(path_out, sub_effects, sep = "")) == FALSE){
  dir.create(paste(path_out, sub_effects, sep = ""))}


# Load Participant data 
d_words <- read.csv(paste(path_in,"blood_preprocessed.csv", sep = ""), header = TRUE)

d_words_ghrelin <- filter(d_words, Ghrelin == 1)
d_words_saline <- filter(d_words, Ghrelin == 0)

#Create a vector containing only the text
text <- d_words_saline$Words_coded
# Create a corpus  
docs <- Corpus(VectorSource(text))

docs <- docs %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
docs <- tm_map(docs, content_transformer(tolower))
docs <- tm_map(docs, removeWords, stopwords("english"))

dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)

set.seed(2345) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1,
max.words=15, random.order=FALSE, rot.per=0.15, colors=brewer.pal(2, "Dark2"))

wordcloud2(data=df, size=0.8, color='random-dark')

# -- compare two wordclouds --- 

text_ghrelin <- d_words %>% filter(Ghrelin == 1) %>% pull(Words_coded)
text_saline  <- d_words %>% filter(Ghrelin == 0) %>% pull(Words_coded)

# --- Create corpus-cleaning function ---
clean_text <- function(text_vec) {
  docs <- Corpus(VectorSource(text_vec))
  docs <- docs %>%
    tm_map(removeNumbers) %>%
    tm_map(removePunctuation) %>%
    tm_map(stripWhitespace) %>%
    tm_map(content_transformer(tolower)) %>%
    tm_map(removeWords, stopwords("english"))
  
  dtm <- TermDocumentMatrix(docs)
  matrix <- as.matrix(dtm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  return(word_freq)
}

# --- Get frequency vectors ---
freq_ghrelin <- clean_text(text_ghrelin)
freq_saline  <- clean_text(text_saline)

# --- Combine frequencies into one dataframe ---
all_words <- union(names(freq_ghrelin), names(freq_saline))

df_compare <- data.frame(
  word = all_words,
  ghrelin = ifelse(all_words %in% names(freq_ghrelin), freq_ghrelin[all_words], 0),
  saline  = ifelse(all_words %in% names(freq_saline),  freq_saline[all_words],  0)
)

# 

# --- Function to clean text and return word frequencies ---
get_word_freq <- function(text_vec) {
  docs <- Corpus(VectorSource(text_vec)) %>%
    tm_map(removeNumbers) %>%
    tm_map(removePunctuation) %>%
    tm_map(stripWhitespace) %>%
    tm_map(content_transformer(tolower)) %>%
    tm_map(removeWords, stopwords("english"))
  
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  freq <- sort(rowSums(m), decreasing = TRUE)
  tibble(word = names(freq), freq = freq)
}

# --- Apply to both conditions ---
freq_ghrelin <- d_words %>% filter(Ghrelin == 1) %>% pull(Words_coded) %>% get_word_freq()
freq_saline  <- d_words %>% filter(Ghrelin == 0) %>% pull(Words_coded) %>% get_word_freq()

# --- Add condition labels and combine ---
freq_ghrelin$Condition <- "Ghrelin"
freq_saline$Condition  <- "Saline"

df_words <- bind_rows(freq_ghrelin, freq_saline)
# --- Set consistent word colors across both facets ---
unique_words <- unique(df_words$word)
n_colors <- length(unique_words)

# Set color palette (extend if needed)
word_colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(n_colors), unique_words)

# --- Plot ---
wordcloud_plot <- ggplot(df_words, aes(label = word, size = freq, color = word)) +
  geom_text_wordcloud_area(rm_outside = TRUE, shape = "square", grid_size = 4) +
  scale_size_area(max_size = 50) +
  scale_color_manual(values = word_colors, guide = "none") +
  facet_wrap(~Condition) +
  theme_void() +  # removes all background, ticks, and labels
  theme(
    strip.text = element_text(face = "bold", size = 20, hjust = 0.5),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(paste(path_out, sub_effects, "wordcloud_conditions.png", sep = ""),
       plot = wordcloud_plot,
       width = 6, height = 6, units = "in", dpi = 600)

#############
# Quantedo library Wordcloud 

# Prepare CORPUS from DF 
Words_prep <- corpus(d_words, text_field = "Words_coded")
docid <- paste(d_words$fGhrelin)
docnames(Words_prep) <- docid
print(Words_prep)

# Create a dfm grouped by president
dfmat_pres <- tokens(Words_prep, remove_punct = TRUE) |>
  tokens_group(groups = fGhrelin) |>
  dfm()
tstat_keyness <- textstat_keyness(dfmat_pres, target = "Ghrelin")
# Plot estimated word keyness
textplot_keyness(tstat_keyness, , min_count = 2,   color = c(color_Ghrelin, color_Placebo), 
  labelcolor = "black",
  labelsize = 4,
  font = NULL, show_legend = FALSE, show_reference = TRUE)

ggsave(paste(path_out, sub_effects, "ChiSquare_Sideeffects.png", sep = ""), 
        height = 4, width = 7, units = "in", dpi = 600, bg = "white")


# Plot frequency of Sideeffects in Cleaveland Dot Plot 
tstat_freq_inaug <- textstat_frequency(dfmat_pres, groups = fGhrelin)

ggplot(tstat_freq_inaug, aes(x = frequency, y = reorder(feature, frequency)), color = fGhrelin) +
    geom_line(aes(group = feature)) +
    geom_point(aes(color = group)) + 
  scale_color_manual(guide = guide_legend(title="Condition"),values = c(color_Ghrelin, color_Placebo)) +
    labs(x = "Frequency", y = "Feature")

ggsave(paste(path_out, sub_effects, "Sideffects_Cleaveland.png", sep = ""), 
        height = 6, width = 4, units = "in", dpi = 600, bg = "white")


# Prepare CORPUS from DF 
Words_prep <- corpus(d_words, text_field = "Words_agg")
docid <- paste(d_words$fGhrelin)
docnames(Words_prep) <- docid
# Create a dfm grouped
dfmat_pres <- tokens(Words_prep, remove_punct = TRUE) |>
  tokens_group(groups = fGhrelin) |>
  dfm()

tstat_keyness <- textstat_keyness(dfmat_pres, target = "Ghrelin",measure = "chi2", correction = "default" )
textplot_keyness(tstat_keyness, , min_count = 2,   color = c(color_Ghrelin, color_Placebo), 
  labelcolor = "black",
  labelsize = 4,
  font = NULL, show_legend = FALSE, show_reference = TRUE)
tstat_freq_inaug <- textstat_frequency(dfmat_pres, groups = fGhrelin)

ggsave(paste(path_out, sub_effects, "ChiSquare_Sideeffects2.png", sep = ""), 
        height = 3, width = 13, units = "in", dpi = 600, bg = "white")

ggplot(tstat_freq_inaug, aes(x = frequency, y = reorder(feature, frequency)), color = fGhrelin) +
    geom_line(aes(group = feature), size = 3, colour = '#D0D0D0') +
    geom_point(aes(color = group), size = 4) + 
  scale_color_manual(guide = guide_legend(title="Condition"),values = c(color_Ghrelin, color_Placebo)) +
    labs(x = "Frequency", y = "Feature") +
    scale_x_continuous(breaks=c(2,4,6,8,10,12))

ggsave(paste(path_out, sub_effects, "Sideffects_Cleaveland_Agg.png", sep = ""), 
        height = 4, width = 4, units = "in", dpi = 600, bg = "white")

##############
##################

# Plot Frequency Plot using VCD package 

d_words <- d_words %>% 
            mutate(Hungry_Ghrelin = ifelse(Yes_hungry == 1 & Ghrelin== 1,1, ifelse(Yes_hungry == 1 & Ghrelin== 0,2, 0)))

d_words$hungry <- factor(d_words$Yes_hungry, labels = c("No", "Yes"))
d_words$condition <- factor(d_words$Ghrelin, labels = c("Placebo", "Ghrelin"))
d_words$belief<- factor(d_words$Test_blinding_pp, labels = c("Placebo", "Ghrelin"))
d_words$importance_IMT <- factor(d_words$IMT_RewType_importance, levels = c("money", "equally", "food"))

d_words$hungry <- ordered(d_words$hungry, 
                              levels=c("No", "Yes"))

# Plot Frequency Plot: SUBEJCTIVE HUNGER REPORTED RECALL) END OF SESSION
artH <- xtabs(~condition +hungry, data = d_words)
mosaic(artH, gp = shading_max, direction = c("v", "h"), main="Reported side effect: hunger")
summary(artH)

chisq.test(artH, correct = FALSE)

ggsave(paste(path_out, sub_effects, "Hunger_Frequency.png", sep = ""), 
        height = 4, width = 10, units = "in", dpi = 600, bg = "white")

# Plot Frequency Plot: BLINDING PARTICIPANTS
art <- xtabs(~belief + condition, data = d_words)
mosaic(art, gp = shading_max, split_vertical = FALSE, main="Successful blinding of participants")
summary(art)
chisq.test(art, correct = FALSE)

art1 <- xtabs(~belief + condition, data = filter(d_words, Session == 1))
art2 <- xtabs(~belief + condition, data = filter(d_words, Session == 2))
mosaic(art1, gp = shading_max, split_vertical = FALSE, main="Session 1")
mosaic(art2, gp = shading_max, split_vertical = FALSE, main="Session 2")

chisq.test(art1, correct = FALSE)
chisq.test(art2, correct = FALSE)

# Plot Frequency Plot: IMT IMPORTANCE 
art_IMT <- xtabs(~importance_IMT + condition, data = d_words)
mosaic(art_IMT, boot = FALSE, gp = shading_max, split_vertical = FALSE , main="Subjective importance of \nRewards (Debriefing IMT)")
# which corresponds to the mosaic plot p-value.
max_test <- summary(art_IMT, type = "maximum")
fisher.test(art_IMT, simulate.p.value = TRUE, B = 1e6)


mosaic(art, gp = shading_max, split_vertical = FALSE, main="Successful blinding of participants")
summary(art)


# IMT 

max_prop <- d_words %>% 
  group_by(fGhrelin) %>%
  count(importance_IMT) %>%
  mutate(freq = n / sum(n)) %>% 
  .$freq %>%
  max

max_prop <- plyr::round_any(max_prop, 0.05, f = ceiling)


p2 <- d_words %>% 
    select(ID, fGhrelin, IMT_RewType_importance) %>% 
  tidyr::pivot_wider(names_from = fGhrelin, values_from = IMT_RewType_importance) %>% 
  make_long(Placebo, Ghrelin) %>% 
  mutate(node = factor(node, levels = c("money", "equally", "food")),
         next_node = factor(next_node, levels = c("money", "equally", "food"))) %>% 
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node))) +
  geom_sankey(alpha = 0.7,
              node.color = 'black') +
  geom_sankey_label(aes(label = node), alpha = 0.75,
                    size = 5, color = "black", fill = "gray80") +
  scale_x_discrete(expand = c(0.05,0.05)) +
  scale_fill_viridis_d(option = "E", alpha = 0.2) +
  theme_sankey(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(reverse = T, nrow = 1)) +
  labs(title = "",
       fill = "'One reward \nmore important..'",
       x = "")

ggsave(paste(path_out, sub_effects, "IMT_Reports_Ghrelin.png", sep = ""), 
        height = 8, width = 6, units = "in", dpi = 600, bg = "white")



