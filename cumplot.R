# load libraries
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
# or just library(tidyverse)

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

positions_path <- args[1]
output_pdf <- gsub('.csv', '_cumplot.pdf', positions_path)
output_png <- gsub('.csv', '_cumplot.png', positions_path)

# read in positions file
positions <- read_csv(positions_path, show_col_types = FALSE)

# cumulative plot
cum_plot <- positions |>
  # arrange by position
  arrange(potential_mut_site) |>
  # group by sequence and context
  group_by(seq_name, context) |>
  # find cumulative potential and match sites for each group
  mutate(cum_potential = cumsum(prop_context),
         cum_match = cumsum(mut_match),
         context = factor(str_to_sentence(context), levels = c('Primary', 'Control'))) |>
  # plot cumulative potential and match sites for each group
  ggplot(aes(x = cum_potential,
             y = cum_match,
             col = seq_name,
             linetype = context)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Cumulative number of potential sites',
       y = 'Cumulative number of matches',
       linetype = 'Context',
       col = 'Sequence')

ggsave(output_pdf, cum_plot, width = 10, height = 7)
ggsave(output_png, cum_plot, width = 10, height = 7)