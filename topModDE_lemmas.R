# TOPIC MODELING - FEMINADE - term = lemma
# loading libraries
library(tm)
library(NLP)
library(tidyverse)
library(tidytext)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(wordcloud)
library(RColorBrewer)
library(stopwords)
library(lda)
library(topicmodels)

# Step 1: Importing data
# 1) Create a list of the files to read in with list.files()
files <- list.files(path = "femina_corpus/feminaDE_tagged",
                    pattern = "*.txt",
                    full.names = TRUE)

# 2) Create an empty data frame and 3) bind_rows() with each file in the for loop
text_files <- data.frame()
for (i in 1:length(files)){
  this_file <- read_tsv(files[i], col_names = FALSE)
  this_file$filename <- gsub("\\_tagged\\.txt|feminaDE_tagged\\/", "", files[i])
  text_files <- bind_rows(text_files, 
                          this_file)
}

# 3) Change column names with rename()
text_files <- text_files %>%
  rename(token = X1,
         tag = X2,
         lemma = X3)

# Step 2: converting data frame into tibble
text_files_df <- tibble(text_files)

# Step 3: Preprocessing: filtering tags NN, NE, ADJA, ADJD, VV*
tags <- c("ADJA", "ADJD", "NN", "NE", "VMFIN", "VMINF", "VMPP", "VVFIN", 
          "VVIMP", "VVINF", "VVIZU", "VVPP")
feminaDE_anv <- text_files_df %>%
  filter(tag %in% tags)

# creating custom stop word
custom_stop_word <- tribble(
  ~lemma,      ~lexicon,
  "<unknown>", "CUSTOM"
  )

feminaDE_anv_clean <- feminaDE_anv %>%
  anti_join(custom_stop_word, by = "lemma")

feminaDE_anv_clean %>%
  filter(lemma == "<unknow>")


# Step 4: counting lemmas and transforming dataframe in DTM
# The topicmodels package in R takes a dtm as input and produces a model that can be tidied by tidytext
# cast_dtm function from tidytext package - 
feminaDE_dtm <- feminaDE_anv_clean %>%
  count(filename, lemma, sort = TRUE) %>%
  cast_dtm(document=filename, term=lemma, value=n)

as.matrix(feminaDE_dtm)
feminaDE_dtm

# Step 5: topic modeling - Parameters from Cel McCracken
## 10 Topics 
feminaDE_mod = LDA(feminaDE_dtm, k=10, method="Gibbs",
                     control=list(alpha=0.001, seed=10005, burnin= 500,
                                  delta= 0.1, iter=4000,  thin= 100))

          
top_terms_mod <-tidy(feminaDE_mod, matrix = "beta") %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  slice(seq_len(15)) %>%
  arrange(topic, desc(beta)) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  mutate(topic = paste("Topic", str_pad(topic, width = 2, pad = "0"), sep = " "))

write_csv(top_terms_mod, "outputs_feminaDE/topTermsMod15.csv")

# Step 6: Plot the topics
top_terms_mod %>%
  mutate(topic = factor(topic),
         term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(term, beta, fill = log(beta))) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "grey20", size = 0.2) +
  scale_x_reordered() +
  facet_wrap(~topic, scales = "free_y", ncol = 3) +
  coord_flip() +
  theme_minimal() +
  scale_fill_distiller(palette = "Dark2") +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
  labs(title= "FeminaDE: Strongest Words by Topic", y = NULL, x = NULL)

ggsave("outputs_feminaDE/topterms_feminaDE15.png")

# Step 7: Top terms by documments
topTerms_docDE <-tidy(feminaDE_mod, matrix = "gamma") %>%
  group_by(topic) %>%
  arrange(desc(gamma)) %>%
  slice(seq_len(15)) %>%
  arrange(topic, desc(gamma)) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  mutate(topic = paste("Topic", str_pad(topic, width = 2, pad = "0"), sep = " "))

write_csv(top_terms_mod, "outputs_feminaDE/topTermsDocDE.csv")

topTerms_docDE %>%
  mutate(topic = factor(topic),
         document = reorder_within(document, gamma, topic)) %>%
  ggplot(aes(document, gamma, fill = log(gamma))) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "grey20", size = 0.2) +
  scale_x_reordered() +
  facet_wrap(~topic, scales = "free_y", ncol = 3) +
  coord_flip() +
  theme_minimal() +
  scale_fill_distiller(palette = "Dark2") +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
  labs(title= "FeminaDE: Topic by documents", y = NULL, x = NULL)

ggsave("outputs_feminaDE/toptermsDoc_feminaDE.png")

