library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(leaflet)
library(RColorBrewer)
library(ggfortify)

# Load the data
water_data <- read_excel("DEP PFAS in Fish Database 09_12_23_final.xlsx", sheet = "W_Data")

# Filter to real samples only (not blanks, dups, etc.) not QCs either
water_clean <- water_data %>%
  filter(`Sample Type (R, B, D, FB, EB)` == "R", `QC sample (Y/N)` == "N")

# there were two samples at each site, except for Falls Pond and Hopedale Pond which had two
# I would ideally have liked to have seen a min of 3 samples per site.
print(water_clean %>%
  group_by(Waterbody) %>%
  summarise(n_samples = n_distinct(FieldSampNum)) %>%
  arrange(desc(n_samples)), n=52)

# Impute <MDL only if compound was detected somewhere
detected_compounds <- water_clean %>%
  group_by(`Analyte/Abbrv`) %>%
  summarise(any_detected = any(!is.na(Result_Num)), .groups = "drop")

# Join with main data
water_estimated <- water_clean %>%
  left_join(detected_compounds, by = "Analyte/Abbrv") %>%
  mutate(
    Result_Est = case_when(
      !is.na(Result_Num) ~ Result_Num,
      grepl("<", Result) & any_detected ~ as.numeric(MDL) / 2,
      TRUE ~ NA_real_
    )
  )

# Average  PFAS values per Waterbody + compound for the two sites that have two samples
pfas_water_avg <- water_estimated %>%
  group_by(Waterbody, `Analyte/Abbrv`) %>%
  summarise(mean_value = mean(Result_Est, na.rm = TRUE), .groups = "drop")

# Pivot wider so each Waterbody has 1 row, compounds as columns
pfas_water_wide <- pfas_water_avg %>%
  pivot_wider(names_from = `Analyte/Abbrv`, values_from = mean_value)

# Save spatial info and remove it from the clustering matrix
pfas_geo <- water_estimated %>%
  select(Waterbody, `Sampling Location Lat`, `Sampling Location Long`) %>%
  distinct()

pfas_matrix <- pfas_water_wide %>%
  select(-Waterbody) %>%
  select(where(~ !all(is.na(.))))

# there are some stray NaNs in there, let's clean that up
# Get median MDL per compound (ignoring non-numeric or missing)

mdl_lookup <- water_estimated %>%
  filter(!is.na(MDL)) %>%
  group_by(`Analyte/Abbrv`) %>%
  summarise(median_mdl = median(as.numeric(MDL), na.rm = TRUE), .groups = "drop") %>%
  deframe()  # named vector: names = PFAS acronyms

# Replace NaN in each compound column with ½ the MDL
pfas_matrix <- pfas_matrix %>%
  mutate(across(
    everything(),
    ~ ifelse(is.na(.), mdl_lookup[cur_column()] / 2, .)
  ))

# Scale the data and perform clustering
pfas_scaled <- scale(pfas_matrix)
dist_matrix <- dist(pfas_scaled)
pfas_hclust <- hclust(dist_matrix, method = "ward.D2")

# Determine optimal number of clusters visually with scree plot
fviz_nbclust(pfas_scaled, kmeans, method = "wss")

# Cut tree into clusters (choose k = 5 based on scree plot)
k <- 5
pfas_geo$cluster <- as.factor(cutree(pfas_hclust, k = k))

# Map the clusters
pal <- colorFactor(brewer.pal(k, "Dark2"), domain = pfas_geo$cluster)

leaflet(pfas_geo) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    lng = ~`Sampling Location Long`,
    lat = ~`Sampling Location Lat`,
    color = ~pal(cluster),
    radius = 6,
    popup = ~paste("Sample:", Waterbody, "<br>Cluster:", cluster),
    stroke = TRUE, fillOpacity = 0.8
  ) %>%
  addLegend("bottomright", pal = pal, values = ~cluster, title = "Cluster")
