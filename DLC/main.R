# Load Packages -----------------------------------------------------------

# data loading
library(readxl); library(data.table)

# data manipulation
library(tidyverse); library(dplyr); library(tidyr); library(rlang); library(stringr); 
library(purrr); library(forcats); library(broom); library(glue); library(lubridate)

# analysis & visualization
library(ggplot2); library(nortest); library(hrbrthemes); library(gtools);
library(directlabels)

# Load Necessary Datasets -------------------------------------------------
Pupil_x = read_excel("Pupil_m131_WG04182023DLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                     sheet = "x")
Pupil_y = read_excel("Pupil_m131_WG04182023DLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                     sheet = "y")
Pupil_likelihood = read_excel("Pupil_m131_WG04182023DLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                              sheet = "likelihood")

rhythm_x = read_excel("rythem at beginning m666_20230905_50ms_80dbDLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                     sheet = "x")
rhythm_y = read_excel("rythem at beginning m666_20230905_50ms_80dbDLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                     sheet = "y")
rhythm_likelihood = read_excel("rythem at beginning m666_20230905_50ms_80dbDLC_Resnet50_GRANT mouse facesAug31shuffle1_snapshot_200.xlsx",
                              sheet = "likelihood")
  

Pupil_y = 
  Pupil_y %>%
  mutate(eye_spread = eye_peak - eye_lower_lid,
         pupil_spread = pupil_peak - pupil_lower,
         whisker_deflection_left = whisker_top_left - whisker_base_left,
         whisker_deflection_right = whisker_top_right - whisker_base_right,
         cheek_spread = cheek_top_left - cheek_bottom_right)

rhythm_y = 
  rhythm_y %>%
  mutate(eye_spread = eye_peak - eye_lower_lid,
         whisker_deflection_left = whisker_top_left - whisker_base_left,
         whisker_deflection_right = whisker_top_right - whisker_base_right,
         cheek_spread = cheek_top_left - cheek_bottom_right)


graph =
Pupil_y %>%
  rename(Frame = bodyparts) %>%
  select(Frame, eye_spread, pupil_spread, cheek_spread, cheek_middle, cheek_top_right, jaw, whisker_base_left) %>%
  gather(key = "bodypart", value = "y", -Frame) %>% 
  ggplot(aes(x = Frame, y = y, color = bodypart)) +
    geom_line(linewidth = 1) +
    labs(x = "Frame",
         y = "y position",
         color = "bodypart") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(250, 1500)) +
  # scale_y_continuous(expand = c(0, 0.01)) +
  facet_wrap(~ bodypart, scales = "free_y", ncol = 1, strip.position = "left") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    legend.position = "none"
    # legend.position = c(.9,.85)
  )

print(graph)

ggsave(file="rythm1.svg", plot=graph, width=10, height=8)

rhythm_y %>%
  rename(Frame = bodyparts) %>%
  # select() %>%
  gather(key = "bodypart", value = "y", -Frame) %>% 
  ggplot(aes(x = Frame, y = y, color = bodypart)) +
  geom_line() +
  labs(x = "Frame",
       y = "y position",
       color = "bodypart") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(50, 600)) +
  # scale_y_continuous(expand = c(0, 0.01)) +
  # facet_wrap(~ bodypart, scales = "free_y") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = rgb(235, 235, 235, 255, maxColorValue = 255)),
    # legend.position = "none"
    # legend.position = c(.9,.85)
  )
