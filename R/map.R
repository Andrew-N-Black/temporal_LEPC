library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(patchwork)
library(ggplot2)
library(ggspatial)
library(svglite)

# --- 1. Load data ---------------------------------------------------------
df <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity_unrel.xlsx")

# --- 2. Robust GPS parsing -------------------------------------------------
# Grabs the first decimal number (latitude) and the trailing negative
# decimal number (longitude), regardless of whether GPS uses a comma,
# space, or other separator between them.
df <- df %>%
  mutate(
    GPS = as.character(GPS),
    lat = as.numeric(str_extract(GPS, "^-?\\d+\\.?\\d*")),
    lon = as.numeric(str_extract(GPS, "-\\d+\\.?\\d*\\s*$"))
  )

bad <- df %>% filter(is.na(lat) | is.na(lon))
if (nrow(bad) > 0) {
  message("Could not parse GPS for these IDs — check raw values:")
  print(bad %>% select(ID, GPS))
}

df <- df %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(between(lat, -90, 90), between(lon, -180, 180))

# --- 3. Basemap: lower-48 US states only -----------------------------------
us_states <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!name %in% c("Alaska", "Hawaii"))

# --- 4. Padded box around the NM sample points -----------------------------
pad <- 7   # degrees of padding — increase/decrease to zoom the box in/out

box <- tibble(
  xmin = min(df$lon, na.rm = TRUE) - pad,
  xmax = max(df$lon, na.rm = TRUE) + pad,
  ymin = min(df$lat, na.rm = TRUE) - pad,
  ymax = max(df$lat, na.rm = TRUE) + pad
)

grp_colors <- c("Past" = "cadetblue", "Present" = "black")

# --- 5. Overview map: contiguous US with box highlighting NM sample area ---
overview_map <- ggplot() +
  geom_sf(data = us_states, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_rect(data = box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", linewidth = 1) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Sample Location - New Mexico", x = NULL, y = NULL)

# --- 6. Zoomed-in panel: sample points colored by GRP ----------------------
zoom_map <- ggplot() +
  geom_sf(data = us_states, fill = "grey95", color = "grey50", linewidth = 0.3) +
  geom_point(data = df, aes(x = lon, y = lat, color = GRP), size = 1, alpha = 0.85) +
  scale_color_manual(values = grp_colors) +
  coord_sf(xlim = c(box$xmin, box$xmax), ylim = c(box$ymin, box$ymax), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "red", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 6),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- 7. Combine into a single figure (zoomed panel inset on the overview) --
combined <- overview_map +
  inset_element(zoom_map, left = 0.45, bottom = 0.32, right = 0.99, top = 0.98)

ggsave("/Users/andrewblack/claude/aging/nm_sample_map.png", combined,
       width = 10, height = 7, dpi = 300, bg = "white")

ggsave("/Users/andrewblack/claude/aging/nm_sample_map.svg", combined,
       width = 10, height = 7, device = svglite::svglite, bg = "white")


# --- Standalone nautical-style north arrow, matching the B&W map style ---
north_arrow_plot <- ggplot() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  annotation_north_arrow(
    which_north = "true",
    location = "center",
    style = north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "black",
      text_col = "black"
    ),
    height = unit(1, "npc"),
    width  = unit(1, "npc")
  )

ggsave("/Users/andrewblack/claude/aging/north_arrow.svg", north_arrow_plot,
       width = 3, height = 3, device = svglite::svglite, bg = "transparent")
