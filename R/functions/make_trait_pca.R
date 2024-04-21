## TRAITS (PCA)
mean2_trait_pca <- make_trait_pca(trait_mean2)
make_trait_pca <- function(trait_mean){

  set.seed(32)

  # make wide trait table
  cwm_fat <- trait_mean2 %>%
    # remove nutrient ratio traits
    filter(!trait %in% c("CN_ratio", "NP_ratio")) |>
    select(trait, turfID, origSiteID, warming, grazing, Namount_kg_ha_y, Nitrogen_log, mean) %>%
    #select(Gradient:mean, Mean_elevation, GS, SoilMoisture, SoilTemperature) %>%
    pivot_wider(names_from = "trait", values_from = "mean") %>%
    ungroup()

  # environemntal variables
   #env_alpine <- cwm_fat %>%
  #   filter(origSiteID == "Alpine") %>%
   #  select(warming, grazing, Nitrogen_log)

   #env_subalpine <- cwm_fat %>%
  #   filter(origSiteID == "Sub-alpine") %>%
   #  select(warming, grazing, Nitrogen_log)

  #pca_output_alpine <- cwm_fat %>%
   # filter(origSiteID == "Alpine") %>%
    #select(-(turfID:Nitrogen_log)) %>%
    #rda(scale = TRUE, center = TRUE)

  #pca_output_subalpine <- cwm_fat %>%
   # filter(origSiteID == "Sub-alpine") %>%
  #  select(-(turfID:Nitrogen_log)) %>%
  #  rda(scale = TRUE, center = TRUE)

  pca_output <- cwm_fat %>%
    select(-(turfID:Nitrogen_log)) %>%
    rda(scale = TRUE, center = TRUE)

   #env_fit <- envfit(pca_output_alpine, env_subalpine)
   #env_fit
   #env_out <- scores(env_fit, "factors") |>
    # as_tibble() |>
    # mutate(label = c("Ambient", "Warming",
    #                  "Control", "Natural"),
    #        class = "Environment",
    #        figure_names = c("Ambient", "Warming",
    #                         "Control", "Natural"))

  pca_sites <- bind_cols(
    cwm_fat %>%
      select(turfID:Nitrogen_log),
    fortify(pca_output, display = "sites")
  )

  # arrows
  pca_traits <- fortify(pca_output, display = "species") %>%
    mutate(trait_trans = label) %>%
    fancy_trait_name_dictionary() #|>
    # add environmental variables
    #bind_rows(env_out) |>
    # mutate(class = as.character(class),
    #        class = factor(class, levels = c("Size", "Leaf economics", "Isotopes", "Environment")))

  # permutation test
  # traits
   raw <- cwm_fat %>% select(-(turfID:Nitrogen_log))
  # # meta data
   meta <- cwm_fat %>% select(turfID:Nitrogen_log) %>%
     mutate(turfID = factor(turfID))
  #
  # # adonis test
   if(meta %>% distinct(grazing) %>% count() == 2){
     adonis_result <- adonis2(raw ~ warming * grazing + warming * origSiteID + grazing * origSiteID, data = meta, permutations = 999, method = "euclidean")
   } else {
     adonis_result <- adonis2(raw ~ warming * Nitrogen_log + warming * origSiteID + Nitrogen_log * origSiteID, data = meta, permutations = 999, method = "euclidean")
   }

  outputList <- list(pca_sites, pca_traits, pca_output, adonis_result)

  return(outputList)
}



make_pca_plot <- function(g_trait_pca, n_trait_pca, col2){

  trait_pca <- g_trait_pca
  trait_pca <- trait_mean2_pca

  # both gradients
  e_B1 <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  pcg1 <- trait_pca[[1]] %>%
    as.data.frame() %>%
    #mutate(site_grazing = paste0(substr(origSiteID, 1, 1), "_", substr(grazing, 1, 1)),
    #       site_grazing = factor(site_grazing, levels = c("A_C", "A_N", "S_C", "S_N"))) |>
    ggplot(aes(x = PC1, y = PC2, colour = warming, shape = origSiteID)) +
    geom_point(size = 3) +
    coord_equal() +
    #stat_ellipse(aes(colour = Mean_elevation)) +
    scale_colour_manual(name = "", values = c("grey40", col_palette[2])) +
    scale_shape_manual(name = "Site and origin", values = c(16, 1)) +
    #scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1, name = "Elevation m a.s.l.", limits = c(range[1], range[2])) +
    labs(x = glue("PCA1 ({round(e_B1[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B1[2] * 100, 1)}%)"),
         tag = "a)") +
    theme_bw()
pcg1

trait_pca[[2]] <- trait_pca[[2]] %>%
  mutate(figure_names = case_match(trait_trans,
                                   "plant_height_cm" ~ "Height~(cm)",
                                   "dry_mass_g" ~ "Dry~mass~(g)",
                                   "leaf_area_cm2" ~ "Area~(cm^2)",
                                   "leaf_thickness_mm" ~ "Thickness~(mm)",
                                   "ldmc" ~ "LDMC",
                                   "sla_cm2_g" ~ "SLA~(cm^2*g^{-1})"))
  arrowsg1 <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]] |>
                mutate(PC2 = case_when(label == "leaf_area_cm2" ~ 0.1,
                                       label == "dry_mass_g" ~ -0.6,
                                       label == "plant_height_cm" ~ -1.0,
                                       label == "leaf_thickness_mm" ~ -0.5,
                                       TRUE ~ PC2),
                       PC1 = case_when(label == "dry_mass_g" ~ -1.7,
                                       label == "sla_cm2_g" ~ -0.5,
                                       TRUE ~ PC1)),
              aes(x = PC1+0.3, y = PC2, label = figure_names),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE) +
    coord_equal() +
    labs(x = "PCA1", y = "PCA2",
         tag = "b)") +
    scale_x_continuous(expand = c(.2, 0))  +
    #scale_linetype_manual(name = "", values = c("solid", "dashed", "solid", "solid")) +
    #scale_colour_manual(name = "", values = c("black", "grey40", "grey70", "cornflowerblue")) +
    theme_bw()

arrowsg1

  trait_pca <- n_trait_pca

  # both gradients
  e_B2 <- eigenvals(trait_pca[[3]])/sum(eigenvals(trait_pca[[3]]))

  pcan2 <- trait_pca[[1]] %>%
    mutate(site_warming = paste0(substr(origSiteID, 1, 1), "_", substr(warming, 1, 1)),
           site_warming = factor(site_warming, levels = c("A_A", "A_W", "S_A", "S_W"))) |>
    ggplot(aes(x = PC1, y = PC2, colour = warming, shape = origSiteID, alpha = Nitrogen_log)) +
    geom_point(size = 2) +
    coord_equal() +
    scale_shape_manual(name = "", values = c(17, 16)) +
    scale_colour_manual(name = "", values = c("grey", col_palette[2])) +
    scale_alpha(name = "log(nitrogen)", range = c(0.3, 1)) +
    #scale_colour_viridis_c(end = 0.8, option = "inferno", direction = -1) +
    labs(x = glue("PCA1 ({round(e_B2[1] * 100, 1)}%)"),
         y = glue("PCA2 ({round(e_B2[2] * 100, 1)}%)"),
         tag = "b)") +
    theme_bw()

  arrowsn2 <- trait_pca[[1]] %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(data = trait_pca[[2]],
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 inherit.aes = FALSE) +
    geom_text(data = trait_pca[[2]], aes(x = PC1+0.3, y = PC2, label = trait_fancy),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE)
    geom_text(data = trait_pca[[2]] |>
                mutate(PC2 = case_when(label == "leaf_area_cm2_log" ~ 0.05,
                                       label == "dry_mass_g_log" ~ 0.2,
                                       label == "plant_height_cm_log" ~ -0.15,
                                       TRUE ~ PC2)),
              aes(x = PC1+0.3, y = PC2, label = figure_names),
              size = 2.5,
              inherit.aes = FALSE,
              show.legend = FALSE, parse = TRUE) +
    coord_equal() +
    labs(x = "PCA1", y = "PCA2",
         tag = "d)") +
    scale_x_continuous(expand = c(.2, 0)) +
    theme_bw()

  (pcg1 + pcan2) / (arrowsg1 + arrowsn2) + plot_layout(guides = "collect")
  pcg1 + arrowsg1 + plot_layout(guides = "collect")

}
