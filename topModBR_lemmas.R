# TOPIC MODELING - FEMINABR - term = lemma
# loading libraires
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
files <- list.files(path = "femina_corpus/feminaBR_tagged",
                    pattern = "*.txt",
                    full.names = TRUE)

# 2) Create an empty data frame and 3) bind_rows() with each file in the for loop
text_files <- data.frame()
for (i in 1:length(files)){
  this_file <- read_tsv(files[i], col_names = FALSE)
  this_file$filename <- gsub("\\_tagged\\.txt|feminaBR_tagged\\/", "", files[i])
  text_files <- bind_rows(text_files, 
                          this_file)
}

# 3) Change column names with rename()
text_files <- text_files %>%
  rename(token = X1,
         tag = X2,
         lemma = X3)

write_csv(text_files, "outputs_feminaBR/text_filesBR.csv")

# Step 2: converting data frame into tibble
text_files_df <- tibble(text_files)

# Step 3: Preprocessing: filtering tags NN, NE, ADJA, ADJD, VV*
tags <- c("AO0", "AQ0", "AQC", "AQD", "AQS", "NCCN", "NCCP", "NCCS", "NCFN", 
          "NCFP", "NCFS", "NCMN", "NCMP", "NCMS", "NP0", "VMG", "VMI", "VMM", 
          "VMN", "VMP", "VMS")
feminaBR_anv <- text_files_df %>%
  filter(tag %in% tags)

# creating custom stop word
custom_stop_word <- tribble(
  ~lemma,      ~lexicon,
  "<unknown>", "CUSTOM",
  "ir",        "CUSTOM",
  "no",        "CUSTOM"
)

feminaBR_anv_clean <- feminaBR_anv %>%
  anti_join(custom_stop_word, by = "lemma")

write_csv(feminaBR_anv_clean, "outputs_feminaBR/feminaBR_anv_clean.csv")

# Step 4: counting lemmas and transforming dataframe in DTM
# The topicmodels package in R takes a dtm as input and produces a model that can be tidied by tidytext
# cast_dtm function from tidytext package - 
feminaBR_dtm <- feminaBR_anv_clean %>%
  count(filename, lemma, sort = TRUE) %>%
  cast_dtm(document=filename, term=lemma, value=n)

# Step 5: Topic modeling - Parameters from Cel McCracken
## 10 Topics
feminaBR_mod = LDA(feminaBR_dtm, k=10, method="Gibbs",
                   control=list(alpha=0.001, seed=10005, burnin= 500,
                                delta= 0.1, iter=4000,  thin= 100))

top_terms_modBR <-tidy(feminaBR_mod, matrix = "beta") %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  slice(seq_len(15)) %>%
  arrange(topic, desc(beta)) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  mutate(topic = paste("Topic", str_pad(topic, width = 2, pad = "0"), sep = " "))

write_csv(top_terms_modBR, "outputs_feminaBR/topTermsModBR15.csv")


# Step 6: Plot the topics
top_terms_modBR %>%
  mutate(topic = factor(topic),
         term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(term, beta, fill = log(beta))) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "grey20", size = 0.2) +
  scale_x_reordered() +
  facet_wrap(~topic, scales = "free", ncol = 3) +
  coord_flip() +
  theme_minimal() +
  scale_fill_distiller(palette = "Dark2") +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.text.x = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")) +
  labs(title= "FeminaBR_mod: Strongest Words by Topic", y = NULL, x = NULL)

ggsave("outputs_feminaBR/topterms_feminaBR15.png")

# Step 7: Top terms - by Documents
topTerms_docBR <-tidy(feminaBR_mod, matrix = "gamma") %>%
  group_by(topic) %>%
  arrange(desc(gamma)) %>%
  slice(seq_len(15)) %>%
  arrange(topic, desc(gamma)) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  mutate(topic = paste("Topic", str_pad(topic, width = 2, pad = "0"), sep = " "))

write_csv(top_terms_mod, "outputs_feminaBR/topTermsDocBR.csv")

topTerms_docBR %>%
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
  labs(title= "FeminaBR: Topic by documents", y = NULL, x = NULL)

ggsave("outputs_feminaBR/toptermsDoc_feminaBR.png")
