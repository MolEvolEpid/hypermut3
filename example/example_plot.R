# load libraries
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
# or just library(tidyverse)

# read in positions file
positions_strict <- read_csv('example-strict-positions.csv', show_col_types = FALSE)
positions_partial <- read_csv('example-partial-positions.csv', show_col_types = FALSE)

# cumulative plot (for primary)
cum_matches_plot <- bind_rows(positions_strict %>% mutate(mode = 'Strict'),
                              positions_partial %>% mutate(mode = 'Partial')) %>%
  # arrange by position
  arrange(potential_mut_site) |>
  # group by sequence, context, and match mode
  group_by(seq_name, context, mode) |>
  # find cumulative potential and match sites for each group
  mutate(cum_potential = cumsum(prop_context),
         cum_match = cumsum(mut_match),
         context = factor(str_to_sentence(context), levels = c('Primary', 'Control')),
         mode = factor(mode, levels = c('Strict', 'Partial'))) |>
  # plot cumulative potential and match sites for each group
  ggplot(aes(x = cum_potential,
             y = cum_match,
             col = seq_name,
             linetype = context)) +
  facet_grid(~mode) +
  geom_line() +
  guides(color=guide_legend(ncol=3)) +
  theme_classic() +
  labs(x = 'Cumulative number of potential sites',
       y = 'Cumulative number of matches',
       linetype = 'Context',
       col = 'Sequence') +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        legend.position = 'bottom', legend.direction = 'vertical',
  )

ggsave('example.png', cum_matches_plot, width = 7, height = 5)
