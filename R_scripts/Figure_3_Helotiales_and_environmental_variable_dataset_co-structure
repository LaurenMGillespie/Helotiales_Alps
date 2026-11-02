# Figure 3

# Co-inertia analyses (CIA) between Helotiales lineage relative abundance and environmental variables including 
# (a) climate variables, (b) surface soil variables, (c) vegetation indices, and (d) plant family relative abundance. 
# Fungal and plant relative abundances were normalized before analysis (log transformation). Lineage names are listed 
# in Table S2 and environmental variable details are in Table S3. Abbreviations: LW = long wave; SW = short wave; 
# CWS = climatic water stress; AG = of α-Glucosidase; BG = β-Glucosidase; CB = β-D-cellubiosidase; LAP = Leucine aminopeptidase; 
# NAG = N-acetyl-β-Glucosaminidase; PHOS = Phosphatase; XYL = β-Xylosidase; tolal_EEA = extracellular enzyme activity; 
# EEN_EEP_ratio = ratio of exoenzymes targeting nitrogen and phosphorus ((LAP+NAG)/PHOS); EEN_relative = relative activity 
# of exoenzymes targeting nitrogen ((LAP+NAG)/total EEA); EEP_relative = relative activity of exoenzymes targeting phosphorus 
# (PHOS/total_EEA); SOM = soil organic matter; H’ = Shannon diversity index.


# Required packages ------------------------------------------------------------------------
libarary("dplyr")
library("data.table")
library("ade4")
library("adegraphics")
library("vegan")


# import data ------------------------------------------------------------------------------
orchamp_metabar_blast_agg_years <- readRDS("data/orchamp_metabar_blast_agg_years_20260327.rds")  


# color palette ----------------------------------------------------------------------------
colors_pca = orchamp_metabar_blast_agg_years %>%
  filter(grepl("Helotiales_all|Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3)) %>%
  dplyr::select(code_analysis_3, code_analysis_4, root_endophyte_status, mycorrhizal_status, lifestyle, lifestyle_JA) %>%
distinct() %>% 
  mutate(
  color_lifestyle_JA = case_when(
    lifestyle_JA == "Endophyte" ~ "#E41A1C",
    lifestyle_JA == "Ericoid_Mycorrhizal" ~ "#377EB8",
    lifestyle_JA == "Other_or_Unknown" ~ "#4DAF4A",
    lifestyle_JA == "Mixed" ~ "black"
  ))


# Figure 3a. Helotiales and Meteorological variable dataset co-structure ------------------

# create meteo data subset
meteo_variables_agg = orchamp_metabar_blast_agg_years_2 %>%
  filter(file_name %in% c("SafranAlpes.RData", "CrocusAlpes.RData")) %>%
  select(code_analysis_2) %>%
  distinct() %>%  
  dplyr::pull(., code_analysis_2)

orchamp_wide_meteo_agg = orchamp_metabar_blast_agg_years %>%
  group_by(codeplot, code_analysis, idtiers, code_analysis_4 , file_name, altitude, slope, X_L93, Y_L93) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  select(code_analysis, code_analysis_4 , idtiers, codeplot, value_mean, altitude, slope, X_L93, Y_L93) %>%  
  filter(grepl("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1", code_analysis_4) &
           code_analysis == "fung_relative_abundance") %>%
  mutate(value_mean = log10((value_mean + (0.0000113/2)))) %>% 
  tidyr::pivot_wider(names_from = code_analysis_4 , values_from = value_mean) %>%
  mutate_at(c('altitude', 'slope', 'X_L93', 'Y_L93'), as.numeric) %>% 
  mutate_at(vars(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")), ~replace(., is.na(.), 0)) %>% 
  tidyr::drop_na() %>%
  select(-code_analysis)

orchamp_wide_meteo_agg_sub = orchamp_metabar_blast_agg_years %>%
  select(codeplot, code_analysis_3, value_mean) %>%
  filter(code_analysis_3 %in% meteo_variables_agg) %>%
  mutate(code_analysis_3 =  case_when(
    #code_analysis_3 == "In situ soil subsurface temperature in degree Celsius" ~ "HOBO_soil_temp",
    code_analysis_3 == "TG1" ~ "Soil_temp",
    code_analysis_3 == "TG4" ~ "Soil_temp4",
    code_analysis_3 == "Rainf" ~ "Rainfall",
    code_analysis_3 == "Snowf" ~ "Snowfall",
    code_analysis_3 == "DSN_T_ISBA" ~ "Snow_depth",
    code_analysis_3 == "ETP" ~ "Evapotranspiration",
    code_analysis_3 == "DIR_SWdown" ~ "Direct_SW_radiation",
    code_analysis_3 == "SCA_SWdown" ~ "Diffuse_SW_radiation",
    code_analysis_3 == "LWdown" ~ "LW_radiation",
    code_analysis_3 == "VPD" ~ "Water_vapor_pressure",
    code_analysis_3 == "Wind" ~ "Wind_speed",
    code_analysis_3 == "Humrel" ~ "Relative_humidity",
    code_analysis_3 == "TairMax" ~ "Max_air_temp",
    code_analysis_3 == "TairMin" ~ "Min_air_temp",
    code_analysis_3 == "TairMean" ~ "Mean_air_temp",
    code_analysis_3 == "net_solar_radiation" ~ "Net_solar_radiation",
    code_analysis_3 == "net_precipitation" ~ "Net_precipitation",
    code_analysis_3 == "CWS" ~ "CWS",
    code_analysis_3 == "SW_solar_radiation" ~ "SW_solar_radiation",
    TRUE ~ code_analysis_3)) %>% 
  tidyr:: pivot_wider(names_from = code_analysis_3, values_from = value_mean) %>%
  tidyr::drop_na()

setDT(orchamp_wide_meteo_agg)[orchamp_wide_meteo_agg_sub, c("GDD", "FDD", "Direct_SW_radiation", "Snow_depth", "Evapotranspiration", "Relative_humidity", #"HOBO_soil_temp", 
                                                            "LW_radiation", "Nebulosity", "Rainfall", "Diffuse_SW_radiation",  "Snowfall", "Soil_temp", 
                                                            "Soil_temp4", "Max_air_temp", "Mean_air_temp", "Min_air_temp", "Water_vapor_pressure",  "Wind_speed",
                                                            "Net_solar_radiation", "Net_precipitation", "CWS", "SW_solar_radiation" ) := 
                                .(GDD, FDD, Direct_SW_radiation, Snow_depth, Evapotranspiration, Relative_humidity, 
                                  LW_radiation, Nebulosity, Rainfall, Diffuse_SW_radiation,  Snowfall, Soil_temp, 
                                  Soil_temp4, Max_air_temp, Mean_air_temp, Min_air_temp, Water_vapor_pressure,  Wind_speed,
                                  Net_solar_radiation, Net_precipitation, CWS, SW_solar_radiation), on=c("codeplot")]


orchamp_wide_meteo_agg = orchamp_wide_meteo_agg %>%
  tidyr::drop_na() %>%
  rename(Slope = slope, Altitude = altitude, Longitude = X_L93, Latitude = Y_L93
         ) %>%
  filter(!codeplot == "TAN_2300")

#Indirect ordination 
pca.fungi <- orchamp_wide_meteo_agg %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

pca.meteo <- orchamp_wide_meteo_agg %>%
  select(matches("Altitude|Slope|Longitude|Latitude|Net_solar_radiation|CWS|SW_solar_radiation|GDD|FDD|Direct_SW_radiation|Snow_depth|Evapotranspiration|Relative_humidity|LW_radiation|Nebulosity|Rainfall|Diffuse_SW_radiation|Snowfall|Soil_temp|Mean_air_temp|Water_vapor_pressure|Wind_speed")) %>% 
  select(-Soil_temp4) %>% 
  dudi.pca(., scale = TRUE, scannf = FALSE, nf = 2)

# Co-inertia analyses 
coia.fungimeteo <- coinertia(pca.fungi, pca.meteo, scannf=FALSE)

colors_fungmeteo = coia.fungimeteo$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)

meteo_coia = s.arrow(coia.fungimeteo$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), plot = FALSE, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
fung_coia = s.arrow(coia.fungimeteo$c1,  plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungmeteo, plines = list(col = colors_fungmeteo), xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
pca_fungmeteo = superpose(meteo_coia, fung_coia, plot = TRUE)


# Mantel test 
orchamp_wide_meteo_agg_mantel = orchamp_metabar_blast_agg_years_2 %>%
  group_by(codeplot, code_analysis, code_analysis_2, code_analysis_3, file_name, altitude, slope, X_L93, Y_L93) %>% #isolate, code_analysis, code_analysis_2, 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  filter(grepl("Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
           code_analysis == "fung_relative_abundance") %>% 
  select(code_analysis, code_analysis_2, codeplot, value_mean, altitude, slope, X_L93, Y_L93) %>%  
  mutate(value_mean = log10((value_mean + (0.0000113/2)))) %>% 
  tidyr:: pivot_wider(names_from = code_analysis_2, values_from = value_mean) %>%
  mutate_at(c('altitude', 'slope', 'X_L93', 'Y_L93'), as.numeric) %>% 
  mutate_at(vars(-matches("codeplot|altitude|slope|X_L93|Y_L93")), ~replace(., is.na(.), 0)) %>% 
  tidyr::drop_na() %>%
  select(-code_analysis)

setDT(orchamp_wide_meteo_agg_mantel)[orchamp_wide_meteo_agg_sub, c("GDD", "FDD", "Direct_SW_radiation", "Snow_depth", "Evapotranspiration", "Relative_humidity",  
                                                            "LW_radiation", "Nebulosity", "Rainfall", "Diffuse_SW_radiation",  "Snowfall", "Soil_temp", 
                                                            "Soil_temp4", "Max_air_temp", "Mean_air_temp", "Min_air_temp", "Water_vapor_pressure",  "Wind_speed",
                                                            "Net_solar_radiation", "Net_precipitation", "CWS", "SW_solar_radiation" ) := 
                                .(GDD, FDD, Direct_SW_radiation, Snow_depth, Evapotranspiration, Relative_humidity, 
                                  LW_radiation, Nebulosity, Rainfall, Diffuse_SW_radiation,  Snowfall, Soil_temp, 
                                  Soil_temp4, Max_air_temp, Mean_air_temp, Min_air_temp, Water_vapor_pressure,  Wind_speed,
                                  Net_solar_radiation, Net_precipitation, CWS, SW_solar_radiation), on=c("codeplot")]


orchamp_wide_meteo_agg_mantel = orchamp_wide_meteo_agg_mantel %>%
  tidyr::drop_na() %>%
  filter(!codeplot == "TAN_2300")

fung_rel_abun_mantel = orchamp_wide_meteo_agg_mantel %>% 
  select(-names(orchamp_wide_meteo_agg_sub), -codeplot, -altitude, -slope, -X_L93, -Y_L93) %>%
  mutate_if(is.numeric, ~(.x+5.247952)) %>%
  vegan::vegdist(., method = "bray")

meteo_rel_abun_mantel = orchamp_wide_meteo_agg_mantel %>%
  select(names(orchamp_wide_meteo_agg_sub), altitude, slope, X_L93,-Y_L93) %>%
  select(-codeplot) %>%
  vegan::vegdist(., method = "bray")

fung_meteo_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  meteo_rel_abun_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_meteo_mantel



# Figure 3b. Helotiales and Surface Soil variable dataset co-structure ------------------

# create surface soil variables subset
orchamp_wide_enz_agg = orchamp_metabar_blast_agg_years %>%
  select(codeplot, value_mean, code_analysis, code_analysis_4, idtiers) %>% 
  mutate(code_analysis_4 = case_when(
    code_analysis_4 == "CN_2mm" ~ "CN",
    code_analysis_4 == "Carbon_prct_2mm" ~ "C",
    code_analysis_4 == "EEN_EEP_ratio" ~ "EEN_EEP_ratio",
    code_analysis_4 == "EEN_relative" ~ "EEN_relative",
    code_analysis_4 == "EEP_relative" ~ "EEP_relative",
    code_analysis_4 == "Enzyme_AG" ~ "Enzyme_AG",
    code_analysis_4 == "Enzyme_BG" ~ "Enzyme_BG",
    code_analysis_4 == "Enzyme_CB" ~ "Enzyme_CB",
    code_analysis_4 == "Enzyme_LAP" ~ "Enzyme_LAP",
    code_analysis_4 == "Enzyme_NAG" ~ "Enzyme_NAG",
    code_analysis_4 == "Enzyme_PHOS" ~ "Enzyme_PHOS",
    code_analysis_4 == "Enzyme_XYL" ~ "Enzyme_XYL",
    code_analysis_4 == "Nitrogen_prct_2mm" ~ "N",
    code_analysis_4 == "Phosphorus" ~ "P",
    code_analysis_4 == "SOM_2mm" ~ "SOM_2mm",
    code_analysis_4 == "SOM_5.6mm" ~ "SOM_5.6mm",
    code_analysis_4 == "pH - H2O 1:4_2mm" ~ "pH - H2O 1:4_2mm",
    code_analysis_4 == "pH - H2O 1:5_2mm" ~ "pH - H2O 1:5_2mm",
    code_analysis_4 == "total_EEA" ~ "total_EEA",
    TRUE ~ code_analysis_4
  )) %>%
  filter(grepl("^CN$|^C$|total_EEA|EEN_EEP_ratio|EEN_relative|EEP_relative|Enzyme_LAP|Enzyme_NAG|Enzyme_PHOS|Enzyme_AG|Enzyme_CB|Enzyme_BG|Enzyme_XYL|SOM_2mm|SOM_5.6mm|^N$|^P$|pH - H2O 1:4_2mm|pH - H2O 1:5_2mm", code_analysis_4)|
        (grepl("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1", code_analysis_4) &
         code_analysis == "fung_relative_abundance")) %>% 
  group_by(codeplot, idtiers,  code_analysis_4) %>%  
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & code_analysis_4 %in% c("Helotiales_sp_1-39", "T._maxilliforme_1-40", "L._orchidicola_UAMH8152", "P._finlandica_CBS444.86", "M._vraolstadiae_G1", "H._ericae_UAMH5828", "C._rhizophila_CBS109839", "Grass_root_mycorrhizal_PPO-1", "Grass_root_mycorrhizal_PPO-2", "M._variabilis_UAMH8861", "A._macrosclerotiorum_Acema1", "Helotiales_other", "Calycina_sp._Bisps1", "Leptodontidium_sp._EF-M", "C._finlandica_PRF15", "P._fortinii_PFO_F", "P._fortinii_CBS179.46", "P._ericae_AC2_IMI_392567", "Ascomycota_sp._Da_Bac_I4", "Ascomycota_sp._Da_PoT1_C7", "Tapesia_sp._Da_MaP_I9", "H._varicosporoides_Hymvar1", "Leptodontidium_sp._MPI-SDFR-AT-0119", "O._maius_OMA_B", "C._variabilis_VPRI_42388_T", "P._ericae_EriJKE06", "P._ericae_LivJKC02", "O._maius_OM19", "Cadophora_sp._Fe06", "Hyaloscyphaceae_sp._1M05", "M._bicolor_E_Melbi2", "M._variabilis_T_Melva1", "L._orchidicola_UAMH5422", "H._hepaticicola_UAMH7357", "O._maius_Zn_Oidma1", "P._subalpina_UAMH-11012", "G._debralockiae_T6G9_T", "P._scopiformis_Phisc1", "Rhexocercosporidium_sp._Rhesp1") ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
    )) %>%
  tidyr::pivot_wider(names_from = code_analysis_4, values_from = value_mean) %>%
  mutate_at(vars(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")), ~replace(., is.na(.), 0)) %>%
  mutate(pH = case_when(
    (!is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ (`pH - H2O 1:4_2mm` + `pH - H2O 1:5_2mm`)/2,
    (is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:5_2mm`,
    (!is.na(`pH - H2O 1:4_2mm`) & is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:4_2mm`,
    TRUE ~ NA
  ), 
  SOM = ((SOM_2mm +`SOM_5.6mm`)/2) 
  ) %>%
  select(-`pH - H2O 1:4_2mm`, -`pH - H2O 1:5_2mm`, -SOM_2mm, -`SOM_5.6mm`
         ) %>% 
  tidyr::drop_na() 

#Indirect ordination 
pca.fungi <- orchamp_wide_enz_agg %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

pca.enz <- orchamp_wide_enz_agg %>%
  select(matches("^CN$|^C$|total_EEA|EEN_EEP_ratio|EEN_relative|EEP_relative|Enzyme_LAP|Enzyme_NAG|Enzyme_PHOS|Enzyme_AG|Enzyme_CB|Enzyme_BG|Enzyme_XYL|^SOM$|^N$|^P$|^pH$")) %>%
  dudi.pca(., scale = TRUE, scannf = FALSE, nf = 2)

#Coinertia Analysis
coia.fungienz <- coinertia(pca.fungi, pca.enz, scannf=FALSE)

colors_fungenz = coia.fungienz$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)

enz_coia = s.arrow(coia.fungienz$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.55))
fung_coia = s.arrow(coia.fungienz$c1, plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungenz, plines = list(col = colors_fungenz), xlim = c(-0.55, 0.55), ylim = c(-0.55, 0.55)) 
superpose(enz_coia, fung_coia, plot = TRUE)


# Mantel test 
orchamp_wide_enz_agg_mantel = orchamp_metabar_blast_agg_years %>%
  select(codeplot, value_mean, code_analysis, code_analysis_2, code_analysis_3, idtiers) %>% 
  mutate(code_analysis_2 = case_when(
    code_analysis_2 == "CN_2mm" ~ "CN",
    code_analysis_2 == "Carbon_prct_2mm" ~ "C",
    code_analysis_2 == "EEN_EEP_ratio" ~ "EEN_EEP_ratio",
    code_analysis_2 == "EEN_relative" ~ "EEN_relative",
    code_analysis_2 == "EEP_relative" ~ "EEP_relative",
    code_analysis_2 == "Enzyme_AG" ~ "Enzyme_AG",
    code_analysis_2 == "Enzyme_BG" ~ "Enzyme_BG",
    code_analysis_2 == "Enzyme_CB" ~ "Enzyme_CB",
    code_analysis_2 == "Enzyme_LAP" ~ "Enzyme_LAP",
    code_analysis_2 == "Enzyme_NAG" ~ "Enzyme_NAG",
    code_analysis_2 == "Enzyme_PHOS" ~ "Enzyme_PHOS",
    code_analysis_2 == "Enzyme_XYL" ~ "Enzyme_XYL",
    code_analysis_2 == "Nitrogen_prct_2mm" ~ "N",
    code_analysis_2 == "Phosphorus" ~ "P",
    code_analysis_2 == "SOM_2mm" ~ "SOM_2mm",
    code_analysis_2 == "SOM_5.6mm" ~ "SOM_5.6mm",
    code_analysis_2 == "pH - H2O 1:4_2mm" ~ "pH - H2O 1:4_2mm",
    code_analysis_2 == "pH - H2O 1:5_2mm" ~ "pH - H2O 1:5_2mm",
    code_analysis_2 == "total_EEA" ~ "total_EEA",
    TRUE ~ code_analysis_2
  )) %>%
  filter(grepl("^CN$|^C$|total_EEA|EEN_EEP_ratio|EEN_relative|EEP_relative|Enzyme_LAP|Enzyme_NAG|Enzyme_PHOS|Enzyme_AG|Enzyme_CB|Enzyme_BG|Enzyme_XYL|SOM_2mm|SOM_5.6mm|^N$|^P$|pH - H2O 1:4_2mm|pH - H2O 1:5_2mm", code_analysis_2)|
           (grepl("Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
              code_analysis == "fung_relative_abundance")
  ) %>% 
  group_by(codeplot, idtiers,  code_analysis_2, code_analysis_3) %>% # c'est pour faire la somme de Helotiales_other et Glomero 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & code_analysis_3 %in% c("Helotiales_sp_1-39", "T._maxilliforme_1-40", "L._orchidicola_UAMH8152", "P._finlandica_CBS444.86", "M._vraolstadiae_G1", "H._ericae_UAMH5828", "C._rhizophila_CBS109839", "Grass_root_mycorrhizal_PPO-1", "Grass_root_mycorrhizal_PPO-2", "M._variabilis_UAMH8861", "A._macrosclerotiorum_Acema1", "Helotiales_other", "Calycina_sp._Bisps1", "Leptodontidium_sp._EF-M", "C._finlandica_PRF15", "P._fortinii_PFO_F", "P._fortinii_CBS179.46", "P._ericae_AC2_IMI_392567", "Ascomycota_sp._Da_Bac_I4", "Ascomycota_sp._Da_PoT1_C7", "Tapesia_sp._Da_MaP_I9", "H._varicosporoides_Hymvar1", "Leptodontidium_sp._MPI-SDFR-AT-0119", "O._maius_OMA_B", "C._variabilis_VPRI_42388_T", "P._ericae_EriJKE06", "P._ericae_LivJKC02", "O._maius_OM19", "Cadophora_sp._Fe06", "Hyaloscyphaceae_sp._1M05", "M._bicolor_E_Melbi2", "M._variabilis_T_Melva1", "L._orchidicola_UAMH5422", "H._hepaticicola_UAMH7357", "O._maius_Zn_Oidma1", "P._subalpina_UAMH-11012", "G._debralockiae_T6G9_T", "P._scopiformis_Phisc1", "Rhexocercosporidium_sp._Rhesp1") ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
  )) %>%
  select(-code_analysis_3) %>%
  tidyr::pivot_wider(names_from = code_analysis_2, values_from = value_mean) %>% 
  mutate_at(vars(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ]))), ~replace(., is.na(.), 0)) %>%
  mutate(pH = case_when(
    (!is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ (`pH - H2O 1:4_2mm` + `pH - H2O 1:5_2mm`)/2,
    (is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:5_2mm`,
    (!is.na(`pH - H2O 1:4_2mm`) & is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:4_2mm`,
    TRUE ~ NA
  ), 
  SOM = ((SOM_2mm +`SOM_5.6mm`)/2) 
  ) %>%
  select(-`pH - H2O 1:4_2mm`, -`pH - H2O 1:5_2mm`, -SOM_2mm, -`SOM_5.6mm`
  ) %>% 
  tidyr::drop_na() 

library(vegan)
fung_rel_abun_mantel = orchamp_wide_enz_agg_mantel  %>% 
  select(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  mutate_if(is.numeric, ~(.x+5.247952)) %>% 
  vegan::vegdist(., method = "bray")

enz_rel_abun_mantel = orchamp_wide_enz_agg_mantel  %>%
  select(-matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  select(-codeplot, -idtiers) %>% 
  vegan::vegdist(., method = "bray")


fung_enz_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  enz_rel_abun_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_enz_mantel



# Figure 3c. Helotiales and Vegetation indices dataset co-structure -----------------------

# create vegetatiojn indices subset

veg_variables_agg = orchamp_metabar_blast_agg_years %>%
  filter(code_analysis_3 %in% c("vegetation_Ericaceae_percent",
                                "vegetation_nonAMplant_percent", 
                                "total_vegetation_cover", 
                            "vegetation_absolute_abundance", 
                            "vegetation_shannon_diversity", 
                            "vegetation_species_richness") ) %>%
  select(code_analysis_3) %>%
  distinct() %>%  
  dplyr::pull(., code_analysis_3)


library(data.table)
orchamp_wide_veg_variables_agg = orchamp_metabar_blast_agg_years %>%
  select(code_analysis_4, codeplot, idtiers, value_mean) %>%
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  filter(grepl("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1", code_analysis_4)) %>%
  mutate(value_mean =  log10((value_mean + (0.0000113/2)))) %>%  
  tidyr::pivot_wider(names_from = code_analysis_4, values_from = value_mean) %>%
  mutate_at(vars(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")), ~replace(., is.na(.), 0)) %>% 
  tidyr::drop_na()

orchamp_wide_veg_variables_agg_sub = orchamp_metabar_blast_agg_years %>%
  select(codeplot, code_analysis_2, idtiers, value_mean) %>%
  filter(code_analysis_2 %in% veg_variables_agg) %>%
  tidyr:: pivot_wider(names_from = code_analysis_2, values_from = value_mean) %>%
  tidyr::drop_na()

setDT(orchamp_wide_veg_variables_agg)[orchamp_wide_veg_variables_agg_sub, c("vegetation_Ericaceae_percent",
                                                        "vegetation_nonAMplant_percent",
                                                        "total_vegetation_cover", 
                                                        "vegetation_absolute_abundance", 
                                                        "vegetation_shannon_diversity", 
                                                        "vegetation_species_richness") := 
                                .(vegetation_Ericaceae_percent,
                                  vegetation_nonAMplant_percent,
                                    total_vegetation_cover, 
                                  vegetation_absolute_abundance, 
                                  vegetation_shannon_diversity, 
                                  vegetation_species_richness), on=c("codeplot")]

orchamp_wide_veg_variables_agg = orchamp_wide_veg_variables_agg %>%
  tidyr::drop_na()

#Indirect ordination 
pca.fungi <- orchamp_wide_veg_variables_agg %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

pca.veg_variables <- orchamp_wide_veg_variables_agg %>% 
  select(matches("vegetation_Ericaceae_percent|vegetation_nonAMplant_percent|total_vegetation_cover|vegetation_shannon_diversity|vegetation_species_richness" )) %>% 
  rename(`H'` = vegetation_shannon_diversity) %>%
  dudi.pca(., scale = TRUE, scannf = FALSE, nf = 2)

#Coinertia Analysis
coia.fungivegvariables <- coinertia(pca.fungi, pca.veg_variables, scannf=FALSE)


colors_fungvegvariables = coia.fungivegvariables$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)

vegvariable_coia = s.arrow(coia.fungivegvariables$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), xlim = c(-0.9, 0.9), ylim = c(-0.9, 0.9))
fung_coia = s.arrow(coia.fungivegvariables$c1, plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungvegvariables, plines = list(col = colors_fungvegvariables), xlim = c(-0.9, 0.9), ylim = c(-0.9, 0.9))
superpose(vegvariable_coia, fung_coia, plot = TRUE)

#### Mantel test ####
orchamp_wide_vegindice_fung_mantel = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_other|plant_relative_abundance|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")|
           code_analysis %in% c("vegetation_Ericaceae_percent", "vegetation_nonAMplant_percent", "vegetation_species_richness")
         ) %>% 
  group_by(codeplot, idtiers,  code_analysis_3) %>%
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & !code_analysis_3 %in% c("vegetation_Ericaceae_percent", "vegetation_nonAMplant_percent", "vegetation_species_richness") ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
  )) %>%
  filter(!codeplot %in% c("CER_1865", "CER_1920", "CER_2075", "CER_2215", "CER_2290", "CER_2315", "DEN_1500", "DEN_1720", "DEN_1870", "DEN_2100", "DEN_2340", "DES_1500", "DES_1730", 
                          "DES_1970", "DES_2130", "DES_2330", "RAC_0520", "SEYS_0460", "SEYS_0615", "SEYS_0840", "SEYS_1035", "SEYS_1350", "SEYS_1500", "SEYS_1700", "SEYS_1870", "VAUL_0830",
                          "VAUL_1100", "VAUL_1290")) %>% # no vegetation sampling at these plots
  tidyr::pivot_wider(names_from = code_analysis_3, values_from = value_mean) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

library(permute)
library(vegan)
fung_rel_abun_mantel = orchamp_wide_vegindice_fung_mantel %>% 
  select(vegetation_Ericaceae_percent, vegetation_nonAMplant_percent, vegetation_species_richness) %>%
  mutate_if(is.numeric, ~(.x+5.247952)) %>%
  vegan::vegdist(., method = "bray")

vegindice_mantel = orchamp_wide_vegindice_fung_mantel %>%
  select(-codeplot, -idtiers, -vegetation_Ericaceae_percent, -vegetation_nonAMplant_percent, -vegetation_species_richness) %>% 
  mutate_if(is.numeric, ~(.x+5.247952)) %>% 
  vegan::vegdist(., method = "bray")


fung_vegindice_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  vegindice_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_vegindice_mantel


# Figure 3d. Helotiales and Vegetation family dataset co-structure -------------------------

# create Vegetation family variables subset
orchamp_wide_vegfam_fung = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_all|Helotiales_other|plant_relative_abundance|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")|
           code_analysis == "plant_relative_abundance") %>% 
  mutate(lb_nom = case_when(
    is.na(lb_nom) & !is.na(isolate) ~ code_analysis_4,
    TRUE ~ lb_nom
  )) %>% 
  mutate(plant_family = case_when(
    is.na(plant_family) & !is.na(isolate) ~ code_analysis_4,
    TRUE ~ plant_family
  )) %>%
  filter(!lb_nom == "aaa - En attente de détermination") %>%
  group_by(codeplot, idtiers,  code_analysis_4, plant_family, lb_nom) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  group_by(codeplot, idtiers,  code_analysis_4, plant_family
  ) %>%
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & !code_analysis_4 == "plant_relative_abundance" ~ log10((value_mean + (0.0000113/2))),
    !is.na(value_mean) & code_analysis_4 == "plant_relative_abundance" ~ log10((value_mean + (0.001006036/2))),
    TRUE ~ value_mean
  )) %>%
  select(-code_analysis_4) %>%
  filter(!codeplot %in% c("CER_1865", "CER_1920", "CER_2075", "CER_2215", "CER_2290", "CER_2315", "DEN_1500", "DEN_1720", "DEN_1870", "DEN_2100", "DEN_2340", "DES_1500", "DES_1730", 
                          "DES_1970", "DES_2130", "DES_2330", "RAC_0520", "SEYS_0460", "SEYS_0615", "SEYS_0840", "SEYS_1035", "SEYS_1350", "SEYS_1500", "SEYS_1700", "SEYS_1870", "VAUL_0830",
                          "VAUL_1100", "VAUL_1290")) %>% # pas d'echantionnage pour la vegetation
  tidyr::pivot_wider(names_from = plant_family, values_from = value_mean) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))


#Indirect ordination 
pca.fungi <- orchamp_wide_vegfam_fung %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

pca.vegfam <- orchamp_wide_vegfam_fung %>%
  select(-matches("codeplot|idtiers|Glomeromycota|Helotiales_all|Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)


#Coinertia Analysis
(coia.fungivegfam <- coinertia(pca.fungi, pca.vegfam, scannf=FALSE))
plot(coia.fungivegfam)
summary(coia.fungivegfam)

colors_fungveg = coia.fungivegfam$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)

veg_coia = s.arrow(coia.fungivegfam$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))
fung_coia = s.arrow(coia.fungivegfam$c1, plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungveg, plines = list(col = colors_fungveg), xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.65))
superpose(veg_coia, fung_coia, plot = TRUE)


# Mantel test

orchamp_wide_vegfam_fung_mantel = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_other|plant_relative_abundance|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")|
           code_analysis == "plant_relative_abundance") %>% 
  mutate(lb_nom = case_when(
    is.na(lb_nom) & !is.na(isolate) ~ code_analysis_2,
    TRUE ~ lb_nom
  )) %>% 
  mutate(plant_family = case_when(
    is.na(plant_family) & !is.na(isolate) ~ code_analysis_2,
    TRUE ~ plant_family
  )) %>%
  filter(!lb_nom == "aaa - En attente de détermination") %>%
  group_by(codeplot, idtiers,  code_analysis_2, plant_family) %>%
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & !code_analysis_2 == "plant_relative_abundance" ~ log10((value_mean + (0.0000113/2))),
    !is.na(value_mean) & code_analysis_2 == "plant_relative_abundance" ~ log10((value_mean + (0.001006036/2))),
    TRUE ~ value_mean
  )) %>%
  select(-code_analysis_2) %>%
  filter(!codeplot %in% c("CER_1865", "CER_1920", "CER_2075", "CER_2215", "CER_2290", "CER_2315", "DEN_1500", "DEN_1720", "DEN_1870", "DEN_2100", "DEN_2340", "DES_1500", "DES_1730", 
                          "DES_1970", "DES_2130", "DES_2330", "RAC_0520", "SEYS_0460", "SEYS_0615", "SEYS_0840", "SEYS_1035", "SEYS_1350", "SEYS_1500", "SEYS_1700", "SEYS_1870", "VAUL_0830",
                          "VAUL_1100", "VAUL_1290")) %>% # pas d'echantionnage pour la vegetation
  tidyr::pivot_wider(names_from = plant_family, values_from = value_mean) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

fung_rel_abun_mantel = orchamp_wide_vegfam_fung_mantel %>% 
  select(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  mutate_if(is.numeric, ~(.x+5.247952)) %>%
  vegan::vegdist(., method = "bray")

vegfam_rel_abun_mantel = orchamp_wide_vegfam_fung_mantel %>%
  select(-matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  select(-codeplot, -idtiers) %>% 
  mutate_if(is.numeric, ~(.x+5.247952)) %>% 
  vegan::vegdist(., method = "bray")


fung_vegfam_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  vegfam_rel_abun_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_vegfam_mantel
