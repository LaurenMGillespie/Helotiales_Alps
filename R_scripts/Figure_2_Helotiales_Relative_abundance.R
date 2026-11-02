## Figure 2. Helotiales relative abundances across the altitudinal gradient 

#(a) Relative abundance of all Helotiales within the entire soil fungal community and relative abundance of target plant-beneficial 
#lineages within just the soil Helotiales community. The plot number indicates the plot altitude (m). ERM lineages are indicated in 
#shades of blue and non-ERM beneficial Helotiales in red-purple shades. Plots are sorted by increasing altitude, the habitat type of 
#the plot is indicated below. (b) Correlation between Helotiales relative abundance in soil and altitude. (c) Correlation between 
#Glomeromycota relative abundance in soil and altitude. R2 values and p-values from LMM analyses with altitude as an explanatory 
#variable and the plot as a random variable are indicated.

# Required packages -------------------------------------
library("dplyr")
library("ggplot2")
library("forcats")
library("ggpmisc")
library("lme4")
library("lmerTest")
library("patchwork")

# Load data ---------------------------------------------

orchamp_metabar_blast_agg_years <- readRDS("data/orchamp_metabar_blast_agg_years_20260327.rds")


# Figure 2a ---------------------------------------------

palette_Fig2 = c(
  "O._maius_OM19"                       = "#2D5DABFF",
  "O._maius_OMA_B"                      = "#A7E6D7FF",
  "O._maius_Zn_Oidma1"                  = "#0C1F33FF",
  "P._ericae_AC2_IMI_392567"            = "#3873B4FF",
  "P._ericae_EriJKE06"                  = "#7ED7D9FF",
  "P._ericae_LivJKC02"                  = "#0277BCFF",
  "C._variabilis_VPRI_42388_T"          = "#418CBFFF",
  "G._debralockiae_T6G9_T"              = "#5EBFD3FF",
  "H._hepaticicola_UAMH7357"            = "#2546A1FF",
  "M._bicolor_E_Melbi2"                 = "#4CA5CAFF",
  "M._variabilis_T_Melva1"              = "#11304CFF",
  "A._macrosclerotiorum_Acema1"         = "#AA0088FF",
  "Ascomycota_sp._Da_Bac_I4"            = "#C83771FF",
  "Ascomycota_sp._Da_PoT1_C7"           = "#FF56D6FF",
  "C._finlandica_PRF15"                 = "#D42AFFFF",
  "C._rhizophila_CBS109839"            = "#5500D4FF",
  "Cadophora_sp._Fe06"                  = "#FFAACCFF",
  "Calycina_sp._Bisps1"                 = "#7B009FFF",
  "Grass_root_mycorrhizal_PPO-1"        = "#FF0000FF",
  "Grass_root_mycorrhizal_PPO-2"        = "#D35F5FFF",
  "H._ericae_UAMH5828"                  = "#D45500FF",
  "H._varicosporoides_Hymvar1"          = "#C87137FF",
  "Helotiales_sp_1-39"                  = "#FDAF91FF",
  "Hyaloscyphaceae_sp._1M05"            = "#CCABFFFF",
  "L._orchidicola_UAMH5422"             = "#FF6600FF",
  "L._orchidicola_UAMH8152"             = "#54278FFF",
  "Leptodontidium_sp._EF-M"             = "#803300FF",
  "Leptodontidium_sp._MPI-SDFR-AT-0119"  = "#C83737FF",
  "M._variabilis_UAMH8861"              = "#FF728CFF",
  "M._vraolstadiae_G1"                  = "#E9C6AFFF",
  "P._finlandica_CBS444.86"             = "#FFD5D5FF",
  "P._fortinii_CBS179.46"               = "#DE8787FF",
  "P._fortinii_PFO_F"                   = "#916F7CFF",
  "P._scopiformis_Phisc1"               = "#E9AFC6FF",
  "P._subalpina_UAMH-11012"             = "#860042FF",
  "Rhexocercosporidium_sp._Rhesp1"      = "#A65628FF",
  "T._maxilliforme_1-40"                = "#D81E00FF",
  "Tapesia_sp._Da_MaP_I9"               = "#AC93A7FF",
  "Helotiales_other"                    = "#FFFFFFFF"
)


### (a) Relative abundance of all Helotiales within the entire soil fungal community ####
fungi_stackedbarplot = orchamp_metabar_blast_agg_years_2 %>%
  filter(code_analysis == "fung_relative_abundance") %>% 
  filter(reads_mean >0) %>% 
  select(codeplot,  altitude, isolate, code_analysis_4, value_mean) %>%
  group_by(codeplot, altitude, code_analysis_4) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  filter(!grepl("Glomero", code_analysis_4)) %>%
  filter(!code_analysis_4 == "Helotiales_all") %>%
  mutate(codeplot =  fct_reorder(codeplot, altitude),
         code_analysis_4 = factor(code_analysis_4, levels=c("Helotiales_other",   "A._macrosclerotiorum_Acema1" ,  "Ascomycota_sp._Da_Bac_I4" ,  "Ascomycota_sp._Da_PoT1_C7" ,  "C._finlandica_PRF15" ,  "C._rhizophila_CBS109839" ,  "Cadophora_sp._Fe06" ,  "Calycina_sp._Bisps1" ,  "Grass_root_mycorrhizal_PPO-1" ,  "Grass_root_mycorrhizal_PPO-2" ,  "H._ericae_UAMH5828" ,  "H._varicosporoides_Hymvar1" ,  "Helotiales_sp_1-39" ,  "Hyaloscyphaceae_sp._1M05" ,  "L._orchidicola_UAMH5422" ,  "L._orchidicola_UAMH8152" ,  "Leptodontidium_sp._EF-M" ,  "Leptodontidium_sp._MPI-SDFR-AT-0119" ,  "M._variabilis_UAMH8861" ,  "M._vraolstadiae_G1" ,  "P._finlandica_CBS444.86" ,  "P._fortinii_CBS179.46" ,  "P._fortinii_PFO_F" ,  "P._scopiformis_Phisc1" ,  "P._subalpina_UAMH-11012" ,  "Rhexocercosporidium_sp._Rhesp1" ,  "T._maxilliforme_1-40" ,  "Tapesia_sp._Da_MaP_I9", "C._variabilis_VPRI_42388_T" ,  "G._debralockiae_T6G9_T" ,  "H._hepaticicola_UAMH7357" ,  "M._bicolor_E_Melbi2" ,  "M._variabilis_T_Melva1" ,  "O._maius_OM19" ,  "O._maius_OMA_B" ,  "O._maius_Zn_Oidma1" ,  "P._ericae_AC2_IMI_392567" ,  "P._ericae_EriJKE06" ,  "P._ericae_LivJKC02" ))) %>% 
  ggplot(., aes(x = codeplot, y = (value_mean*100), fill= code_analysis_4)) +
  geom_bar(position="stack", stat="identity", color = "black") + 
  scale_y_continuous(expand = c(0, 0), lim = c(0, 38),
                     breaks = seq(0.00, 38, by = 10)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values=palette_Fig2,
                    limits = c("Helotiales_other",    "C._variabilis_VPRI_42388_T" ,  "G._debralockiae_T6G9_T" ,  "H._hepaticicola_UAMH7357" ,  "M._bicolor_E_Melbi2" ,  "M._variabilis_T_Melva1" ,  "O._maius_OM19" ,  "O._maius_OMA_B" ,  "O._maius_Zn_Oidma1" ,  "P._ericae_AC2_IMI_392567" ,  "P._ericae_EriJKE06" ,  "P._ericae_LivJKC02" ,  "A._macrosclerotiorum_Acema1" ,  "Ascomycota_sp._Da_Bac_I4" ,  "Ascomycota_sp._Da_PoT1_C7" ,  "C._finlandica_PRF15" ,  "C._rhizophila_CBS109839" ,  "Cadophora_sp._Fe06" ,  "Calycina_sp._Bisps1" ,  "Grass_root_mycorrhizal_PPO-1" ,  "Grass_root_mycorrhizal_PPO-2" ,  "H._ericae_UAMH5828" ,  "H._varicosporoides_Hymvar1" ,  "Helotiales_sp_1-39" ,  "Hyaloscyphaceae_sp._1M05" ,  "L._orchidicola_UAMH5422" ,  "L._orchidicola_UAMH8152" ,  "Leptodontidium_sp._EF-M" ,  "Leptodontidium_sp._MPI-SDFR-AT-0119" ,  "M._variabilis_UAMH8861" ,  "M._vraolstadiae_G1" ,  "P._finlandica_CBS444.86" ,  "P._fortinii_CBS179.46" ,  "P._fortinii_PFO_F" ,  "P._scopiformis_Phisc1" ,  "P._subalpina_UAMH-11012" ,  "Rhexocercosporidium_sp._Rhesp1" ,  "T._maxilliforme_1-40" ,  "Tapesia_sp._Da_MaP_I9"),
                    labels=c("Helotiales other",    "C. variabilis VPRI 42388 T" ,  "G. debralockiae T6G9 T" ,  "H. hepaticicola UAMH7357" ,  "M. bicolor E Melbi2" ,  "M. variabilis T Melva1" ,  "O. maius OM19" ,  "O. maius OMA B" ,  "O. maius Zn Oidma1" ,  "P. ericae AC2 IMI 392567" ,  "P. ericae EriJKE06" ,  "P. ericae LivJKC02" ,  "A. macrosclerotiorum Acema1" ,  "Ascomycota sp. Da Bac I4" ,  "Ascomycota sp. Da PoT1 C7" ,  "C. finlandica PRF15" ,  "C. rhizophila CBS109839" ,  "Cadophora sp. Fe06" ,  "Calycina sp. Bisps1" ,  "Grass root mycorrhizal PPO-1" ,  "Grass root mycorrhizal PPO-2" ,  "H. ericae UAMH5828" ,  "H. varicosporoides Hymvar1" ,  "Helotiales sp 1-39" ,  "Hyaloscyphaceae sp. 1M05" ,  "L. orchidicola UAMH5422" ,  "L. orchidicola UAMH8152" ,  "Leptodontidium sp. EF-M" ,  "Leptodontidium sp. MPI-SDFR-AT-0119" ,  "M. variabilis UAMH8861" ,  "M. vraolstadiae G1" ,  "P. finlandica CBS444.86" ,  "P. fortinii CBS179.46" ,  "P. fortinii PFO F" ,  "P. scopiformis Phisc1" ,  "P. subalpina UAMH-11012" ,  "Rhexocercosporidium sp. Rhesp1" ,  "T. maxilliforme 1-40" ,  "Tapesia sp. Da MaP I9")) + 
  guides(fill=guide_legend(title="Lineages", position = "top", ncol = 7)) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), clip = "off") +
  ylab("Relative abundance within the\nentire soil fungal community (%)") +  
  theme_bw() +
  theme(axis.title.y = element_text(size = 13,  face = "bold"),  
        axis.text.y = element_text(size = 13),
        plot.title = element_text(hjust = 0),
        legend.title = element_text(angle = 90,  hjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  ) 

### (a) Relative abundance of target plant-beneficial lineages within just the soil Helotiales community ####
Helo_stackedbarplot = orchamp_metabar_blast_agg_years_2 %>%
  filter(code_analysis == "fung_helotiales_relative_abundance") %>%
  group_by(codeplot, altitude, code_analysis_4) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!grepl("Glomero", code_analysis_4)) %>%
  filter(!code_analysis_4 == "Helotiales_all") %>%
  mutate(codeplot =  fct_reorder(codeplot, altitude),
         code_analysis_4 = factor(code_analysis_4, levels=c("Helotiales_other",   "A._macrosclerotiorum_Acema1" ,  "Ascomycota_sp._Da_Bac_I4" ,  "Ascomycota_sp._Da_PoT1_C7" ,  "C._finlandica_PRF15" ,  "C._rhizophila_CBS109839" ,  "Cadophora_sp._Fe06" ,  "Calycina_sp._Bisps1" ,  "Grass_root_mycorrhizal_PPO-1" ,  "Grass_root_mycorrhizal_PPO-2" ,  "H._ericae_UAMH5828" ,  "H._varicosporoides_Hymvar1" ,  "Helotiales_sp_1-39" ,  "Hyaloscyphaceae_sp._1M05" ,  "L._orchidicola_UAMH5422" ,  "L._orchidicola_UAMH8152" ,  "Leptodontidium_sp._EF-M" ,  "Leptodontidium_sp._MPI-SDFR-AT-0119" ,  "M._variabilis_UAMH8861" ,  "M._vraolstadiae_G1" ,  "P._finlandica_CBS444.86" ,  "P._fortinii_CBS179.46" ,  "P._fortinii_PFO_F" ,  "P._scopiformis_Phisc1" ,  "P._subalpina_UAMH-11012" ,  "Rhexocercosporidium_sp._Rhesp1" ,  "T._maxilliforme_1-40" ,  "Tapesia_sp._Da_MaP_I9", "C._variabilis_VPRI_42388_T" ,  "G._debralockiae_T6G9_T" ,  "H._hepaticicola_UAMH7357" ,  "M._bicolor_E_Melbi2" ,  "M._variabilis_T_Melva1" ,  "O._maius_OM19" ,  "O._maius_OMA_B" ,  "O._maius_Zn_Oidma1" ,  "P._ericae_AC2_IMI_392567" ,  "P._ericae_EriJKE06" ,  "P._ericae_LivJKC02" ))) %>% 
  ggplot(., aes(x = codeplot, y = (value_mean*100), fill= code_analysis_4)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100.000000001),
                     breaks = seq(0.0, 100, by = 20)) +
  scale_x_discrete(expand = c(0, 0),
                   position = "top") +
  scale_fill_manual(values=palette_JA,
                    limits = c("Helotiales_other",    "C._variabilis_VPRI_42388_T" ,  "G._debralockiae_T6G9_T" ,  "H._hepaticicola_UAMH7357" ,  "M._bicolor_E_Melbi2" ,  "M._variabilis_T_Melva1" ,  "O._maius_OM19" ,  "O._maius_OMA_B" ,  "O._maius_Zn_Oidma1" ,  "P._ericae_AC2_IMI_392567" ,  "P._ericae_EriJKE06" ,  "P._ericae_LivJKC02" ,  "A._macrosclerotiorum_Acema1" ,  "Ascomycota_sp._Da_Bac_I4" ,  "Ascomycota_sp._Da_PoT1_C7" ,  "C._finlandica_PRF15" ,  "C._rhizophila_CBS109839" ,  "Cadophora_sp._Fe06" ,  "Calycina_sp._Bisps1" ,  "Grass_root_mycorrhizal_PPO-1" ,  "Grass_root_mycorrhizal_PPO-2" ,  "H._ericae_UAMH5828" ,  "H._varicosporoides_Hymvar1" ,  "Helotiales_sp_1-39" ,  "Hyaloscyphaceae_sp._1M05" ,  "L._orchidicola_UAMH5422" ,  "L._orchidicola_UAMH8152" ,  "Leptodontidium_sp._EF-M" ,  "Leptodontidium_sp._MPI-SDFR-AT-0119" ,  "M._variabilis_UAMH8861" ,  "M._vraolstadiae_G1" ,  "P._finlandica_CBS444.86" ,  "P._fortinii_CBS179.46" ,  "P._fortinii_PFO_F" ,  "P._scopiformis_Phisc1" ,  "P._subalpina_UAMH-11012" ,  "Rhexocercosporidium_sp._Rhesp1" ,  "T._maxilliforme_1-40" ,  "Tapesia_sp._Da_MaP_I9"),
                    labels=c("Helotiales other",    "C. variabilis VPRI 42388 T" ,  "G. debralockiae T6G9 T" ,  "H. hepaticicola UAMH7357" ,  "M. bicolor E Melbi2" ,  "M. variabilis T Melva1" ,  "O. maius OM19" ,  "O. maius OMA B" ,  "O. maius Zn Oidma1" ,  "P. ericae AC2 IMI 392567" ,  "P. ericae EriJKE06" ,  "P. ericae LivJKC02" ,  "A. macrosclerotiorum Acema1" ,  "Ascomycota sp. Da Bac I4" ,  "Ascomycota sp. Da PoT1 C7" ,  "C. finlandica PRF15" ,  "C. rhizophila CBS109839" ,  "Cadophora sp. Fe06" ,  "Calycina sp. Bisps1" ,  "Grass root mycorrhizal PPO-1" ,  "Grass root mycorrhizal PPO-2" ,  "H. ericae UAMH5828" ,  "H. varicosporoides Hymvar1" ,  "Helotiales sp 1-39" ,  "Hyaloscyphaceae sp. 1M05" ,  "L. orchidicola UAMH5422" ,  "L. orchidicola UAMH8152" ,  "Leptodontidium sp. EF-M" ,  "Leptodontidium sp. MPI-SDFR-AT-0119" ,  "M. variabilis UAMH8861" ,  "M. vraolstadiae G1" ,  "P. finlandica CBS444.86" ,  "P. fortinii CBS179.46" ,  "P. fortinii PFO F" ,  "P. scopiformis Phisc1" ,  "P. subalpina UAMH-11012" ,  "Rhexocercosporidium sp. Rhesp1" ,  "T. maxilliforme 1-40" ,  "Tapesia sp. Da MaP I9")) + 
  guides(fill=guide_legend(title="Lineages", position = "top", ncol = 7)) +
  ylab("Relative abundance within the\nHelotiales soil community (%)") +  
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), clip = "off") +
  annotate("text", x = -3, y = 98, label = "Plot code", angle = 90, hjust=-0.25, size = 5) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13,  face = "bold"),  
        axis.text.y = element_text(size = 13),
        plot.title = element_text(hjust = 0),
        legend.title = element_text(angle = 90,  hjust=0.5),
        plot.subtitle = element_text(angle = 90, hjust = 1.35),
        axis.title.x = element_blank(),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust=0))




# Habitat heatmap
habitat_heatmap <- orchamp_metabar_blast_agg_years_2 %>%
  filter(code_analysis == "fung_relative_abundance") %>%
  filter(!grepl("Glomero", code_analysis_4)) %>%
  filter(!code_analysis_4 == "Helotiales_all") %>%
  select(codeplot, code_analysis, altitude, CorineLandCover_simplifie, CorineLandCover_simplifie_2, habitat_type) %>%
  distinct() %>% 
  mutate(CorineLandCover_simplifie_2 = case_when(
           is.na(CorineLandCover_simplifie_2) & habitat_type == "Milieu forestier" ~ "Forest",
           is.na(CorineLandCover_simplifie_2) & habitat_type == "Milieu ouvert" ~ "Grassland",
           TRUE ~ CorineLandCover_simplifie_2),
         CorineLandCover_simplifie =  case_when(
           CorineLandCover_simplifie == "Forets de coniferes" ~ "Evergreen forest",
           CorineLandCover_simplifie == "Forets de feuillus" ~ "Deciduous forest",
           CorineLandCover_simplifie == "Forets mixtes" ~ "Mixed forest",
           CorineLandCover_simplifie == "Pelouses et patures" ~ "Grassland and pastures",
           CorineLandCover_simplifie == "Pelouses ecorches" ~ "Stripped grassland",
           CorineLandCover_simplifie == "Landes" ~ "Heathland",
           TRUE ~ "NA")) %>% 
  mutate(codeplot =  fct_reorder(codeplot, altitude)) %>%
  ggplot(., aes(codeplot, code_analysis)) +                           
  geom_tile(aes(fill = CorineLandCover_simplifie), color = "black") +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#56B4E9", "#008000ff","#55d400ff", "#808000FF",  "gray"),
                      breaks=c("Stripped grassland", "Grassland and pastures", "Heathland", "Evergreen forest", "Mixed forest", "Deciduous forest", "NA")) +
  coord_cartesian(xlim = c(0, NA), clip = "off") +
  scale_x_discrete(expand = c(0, 0),
                   position = "top") +
  scale_y_discrete(labels= c("Habitat\ntype"),
    expand = c(0, 0))  +
  labs(fill="Habitat type") +
  guides(fill=guide_legend(position = "bottom", ncol = 6)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 13,  face = "bold"),  
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0, "lines"))


### (b) Correlation between Helotiales relative abundance in soil and altitude.  ####
scatterplot_Helotiales_all = orchamp_metabar_blast_agg_years %>%
  filter(code_analysis == "fung_relative_abundance") %>% 
  select(codeplot,  idtiers, altitude, isolate, code_analysis_4, value_mean) %>% 
  group_by(codeplot, idtiers, altitude, code_analysis_4) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)) %>% 
  ungroup() %>%
  filter(!grepl("Glomero", code_analysis_4)) %>%
  filter(!grepl("Helotiales_all", code_analysis_4)) %>%
  group_by(codeplot, idtiers, altitude) %>%  
  summarize( value_mean_sum = sum(value_mean, na.rm = TRUE)) %>% 
  ungroup() %>%
  ggplot(., aes(x=altitude, y=(value_mean_sum*100))) + 
  geom_point(color = "black") +
  stat_poly_line(se = FALSE) +
  annotate("text", x = 500, y = 90,
           hjust = 0, 
           label = "~R^2~'c = 0.63; '~R^2~'m = 0.14'",
           parse = TRUE) +
  annotate("text", x = 500, y = 85,
           hjust = 0, 
           label = "~italic(p)~'< 0.001'",
           parse = TRUE) +
  annotate("text", x = 500, y = 80,
           hjust = 0, 
           label = "~n == 517",
           parse = TRUE) +
  coord_cartesian(xlim = c(200, 3025), expand=FALSE, clip = "off") +
  labs(x = "Altitude (m)", 
       y = "Helotiales relative abundance within\nthe entire soil fungal community (%)") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100),
                     breaks = seq(0.0, 100, by = 20)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(500, 3000, by = 1000)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13,  face = "bold"),  
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13,  face = "bold"),  
        axis.text.x = element_text(size = 13),
        plot.margin = margin(0.15,0.15,0.15,0.8, "cm"))

scatterplot_Helotiales_all

### (c) Correlation between Glomeromycota relative abundance in soil and altitude. ####
scatterplot_Glomero_all = orchamp_metabar_blast_agg_years %>%
  filter(code_analysis == "fung_relative_abundance") %>% 
  select(codeplot,  idtiers, altitude, isolate, code_analysis_4, value_mean) %>% #
  group_by(codeplot, idtiers, altitude, code_analysis_4) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  filter(grepl("Glomero", code_analysis_4)) %>%
  filter(!grepl("Helotiales_all", code_analysis_4)) %>%
  group_by(codeplot, idtiers, altitude) %>%  
  summarize( value_mean_sum = sum(value_mean, na.rm = TRUE)) %>% 
  ungroup() %>%
  ggplot(., aes(x=altitude, y=(value_mean_sum*100))) +
  geom_point(color = "black") +
  stat_poly_line(se = FALSE) +
  annotate("text", x = 500, y = 90,
           hjust = 0, 
           label = "~R^2~'c = 0.43; '~R^2~'m = 0.01'",
           parse = TRUE) +
  annotate("text", x = 500, y = 85,
           hjust = 0, 
           label = "~italic(p)~'= 0.06'",
           parse = TRUE) +
  annotate("text", x = 500, y = 80,
           hjust = 0, 
           label = "~n == 517",
           parse = TRUE) +
  coord_cartesian(xlim = c(200, 3025), expand=FALSE, clip = "off") +
  labs(x = "Altitude (m)", 
       y = "Glomeromycota relative abundance within\nthe entire soil fungal community (%)") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 100),
                     breaks = seq(0.0, 100, by = 20)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = seq(500, 3000, by = 1000)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 13,  face = "bold"),  
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13,  face = "bold"),  
        axis.text.x = element_text(size = 13),
        plot.margin = margin(0.15,0.15,0.15,0.8, "cm"))

scatterplot_Glomero_all


### Increasing altitude triangle ####
increasing_altitude = data.frame(x =c(1:174), y = c(1,2,3)) %>%
  ggplot(., aes(x = x, y = y)) +
  geom_point(color = "white") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), clip = "off") +
  annotate("segment", x = 3, y = 2, xend = 11, yend = 2, 
           linewidth = 1.25, linejoin = "mitre",
           arrow = arrow(type = "closed", length = unit(0.2, "npc"))) +
  annotate("text", x=0, y=2, label = "Increasing altitude", color = "black", 
           angle = 0, size = 5, fontface = "bold", hjust = 1) +
  annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "darkgray")+
  theme_void()+
  theme(plot.margin = margin(0.5,0.15,0.15,0.5, "cm"))

### LMM Helotiales relative abundance vs altitude ###
only_target_helotiales = orchamp_metabar_blast_agg_years %>%
  filter(code_analysis == "fung_relative_abundance") %>% 
  select(codeplot,  idtiers, altitude, isolate, code_analysis_4, value_mean) %>% 
  group_by(codeplot, idtiers, altitude, code_analysis_4) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  filter(!grepl("Glomero", code_analysis_4)) %>%
  filter(!grepl("Helotiales_all|Helotiales_other", code_analysis_4)) %>% 
  group_by(codeplot, idtiers, altitude) %>%  
  summarize( value_mean_sum = sum(value_mean, na.rm = TRUE)) %>% 
  ungroup() 
m0 = lmer(value_mean_sum ~ altitude + (1|codeplot), data = only_target_helotiales, REML = F, control=lmerControl(optimizer="bobyqa"))
summary(m0)

r2_nakagawa(update(m0, REML = T), by_group = FALSE, tolerance = 1e-05)


### Combine figures ###
fungi_stackedbarplot + scatterplot_Helotiales_all + plot_spacer() + plot_spacer() + Helo_stackedbarplot + scatterplot_Glomero_all + plot_spacer() + plot_spacer() + increasing_altitude + plot_spacer() + habitat_heatmap + plot_spacer() +
  plot_layout(nrow = 6, ncol = 2,
    guides = "collect", 
    heights = c(4, -0.7, 4, -0.80, 0.60, 0.30), 
    widths = unit(c(4, 0.70), c('null', 'null', 'cm'))) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 18),
                                legend.position = 'top'))



ggsave("Figures/Fig2.pdf", width = 70, height = 30, units = "cm")
