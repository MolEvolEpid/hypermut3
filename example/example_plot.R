# load library
library(tidyverse)

# read in positions file
positions_strict <- read_csv('example-strict-positions.csv', comment = '#')
positions_partial <- read_csv('example-partial-positions.csv', comment = '#')

# cumulative plot (for primary)
cum_matches_plot <- bind_rows(positions_strict %>% mutate(type = 'Strict'),
                              positions_partial %>% mutate(type = 'Partial')) %>%
  filter(control == 0) %>%
  arrange(potential_mut_site) %>%
  group_by(seq_name, type) %>%
  mutate(cum_potential = cumsum(prop_control),
         cum_match = cumsum(mut_match),
         type = factor(type, levels = c('Strict', 'Partial'))) %>%
  ggplot(aes(x = cum_potential, y = cum_match, col = seq_name)) +
  facet_grid(~type) +
  geom_line() +
  guides(color=guide_legend(ncol=2)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        legend.position = 'bottom', legend.direction = 'vertical',
  ) +
  labs(x = 'Cumulative number of potential sites', y = 'Cumulative number of matches', col = '')

ggsave('example.png', cum_matches_plot, width = 7, height = 5)
