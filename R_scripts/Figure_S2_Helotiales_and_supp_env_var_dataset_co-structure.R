# Figure S2

# Co-inertia analyses between log transformed Helotiales strain relative abundance and (a) the soil profile (i.e. soil pit 
# variables) and (b) log transformed plant species relative abundance. Strain names are listed in Table S2 and environmental 
# variable details are in Table S3.


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

# Figure 2a. Helotiales and soil profile (i.e. soil pit variables) variable dataset co-structure ------------------

# create soil profile data subset
pedo_variables_agg = orchamp_metabar_blast_agg_years_2 %>%
  filter(file_name %in% c("PhysicoChimie_20240711.json" )) %>%
  filter(!code_analysis %in% c("Calcaire (CaCO3) total", "Densité apparente", "Humidité prélèvement", "Pierrosité", "Traitement Carbone", "Aluminium (Al)", "Fer (Fe)", "Silicium (Si)")) %>% # not enough plots for these variables
  select(code_analysis) %>%
  distinct() %>%  
  dplyr::pull(., code_analysis)

orchamp_wide_pedo_agg = orchamp_metabar_blast_agg_years_2 %>%
  group_by(codeplot, code_analysis_4, code_analysis, file_name) %>% 
  summarize(value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%   
  mutate(code_analysis_4 = case_when(
    file_name == "PhysicoChimie_20240711.json" ~ code_analysis,
    TRUE ~ code_analysis_4
  )) %>% 
  select(codeplot, value_mean, code_analysis, code_analysis_4, file_name) %>% 
  mutate(code_analysis_4 = case_when(
    grepl("Aluminium \\(Al\\) CEC", code_analysis_4) ~ gsub("Aluminium \\(Al\\) CEC", "Exchangeable_Al", code_analysis_4),
    grepl("Argile \\(< 2 µm\\)", code_analysis_4) ~ gsub("Argile \\(< 2 µm\\)", "Clay", code_analysis_4),
    grepl("Azote \\(N\\) total", code_analysis_4) ~ gsub("Azote \\(N\\) total", "N_total", code_analysis_4),
    grepl("Broyage \\(< 250 µm\\)", code_analysis_4) ~ gsub("Broyage \\(< 250 µm\\)", "Broyage_0_250µm", code_analysis_4),
    grepl("C\\/N", code_analysis_4) ~ sub("C\\/N", "CN", code_analysis_4),
    grepl("CEC cobaltihexamine", code_analysis_4) ~ gsub("CEC cobaltihexamine", "CEC_cobaltihexamine", code_analysis_4),
    grepl("Calcaire \\(CaCO3\\) total", code_analysis_4) ~ gsub("Calcaire \\(CaCO3\\) total", "CaCO3_total", code_analysis_4),
    grepl("Calcium \\(Ca\\)", code_analysis_4) ~ gsub("Calcium \\(Ca\\)", "Exchangeable_Ca", code_analysis_4),
    grepl("Carbone \\(C\\) organique", code_analysis_4) ~ gsub("Carbone \\(C\\) organique", "Organic_C", code_analysis_4),
    grepl("Densité apparente", code_analysis_4) ~ gsub("Densité apparente", "Bulk_density", code_analysis_4),
    grepl("Fer \\(Fe\\) CEC", code_analysis_4) ~ gsub("Fer \\(Fe\\) CEC", "Exchangeable_Fe", code_analysis_4),
    grepl("Humidité prélèvement", code_analysis_4) ~ gsub("Humidité prélèvement", "Humidite_prelevement", code_analysis_4),
    grepl("Limons fins \\(2\\/20 µm\\)", code_analysis_4) ~ gsub("Limons fins \\(2\\/20 µm\\)", "Fine_silt", code_analysis_4),
    grepl("Limons grossiers \\(20\\/50 µm\\)", code_analysis_4) ~ gsub("Limons grossiers \\(20\\/50 µm\\)", "Coarse_silt", code_analysis_4),
    grepl("Magnésium \\(Mg\\)", code_analysis_4) ~ gsub("Magnésium \\(Mg\\)", "Exchangeable_Mg", code_analysis_4),
    grepl("Manganèse \\(Mn\\)", code_analysis_4) ~ gsub("Manganèse \\(Mn\\)", "Exchangeable_Mn", code_analysis_4),
    grepl("Matiere organique", code_analysis_4) ~ gsub("Matiere organique", "SOM", code_analysis_4),
    grepl("Phosphore \\(P2O5\\)", code_analysis_4) ~ gsub("Phosphore \\(P2O5\\)", "P", code_analysis_4),
    grepl("Pierrosité", code_analysis_4) ~ gsub("Pierrosité", "Stoniness", code_analysis_4),
    grepl("Potassium \\(K\\)", code_analysis_4) ~ gsub("Potassium \\(K\\)", "Exchangeable_K", code_analysis_4),
    grepl("Préparation", code_analysis_4) ~ gsub("Préparation", "Preparation", code_analysis_4),
    grepl("Sables fins \\(50\\/200 µm\\)", code_analysis_4) ~ gsub("Sables fins \\(50\\/200 µm\\)", "Fine_sand", code_analysis_4),
    grepl("Sables grossiers \\(200\\/2000 µm\\)", code_analysis_4) ~ gsub("Sables grossiers \\(200\\/2000 µm\\)", "Coarse_sand", code_analysis_4),
    grepl("Sodium \\(Na\\)", code_analysis_4) ~ gsub("Sodium \\(Na\\)", "Exchangeable_Na", code_analysis_4),
    grepl("Traitement Carbone", code_analysis_4) ~ gsub("Traitement Carbone", "C_treatment", code_analysis_4),
    grepl("DSN_T_ISBA", code_analysis_4) ~ gsub("DSN_T_ISBA", "DSN_T_ISBA", code_analysis_4),
    grepl("In situ soil subsurface temperature in degree Celsius", code_analysis_4) ~ gsub("In situ soil subsurface temperature in degree Celsius", "HOBO_temp_C", code_analysis_4),
    grepl("TG1", code_analysis_4) ~ gsub("TG1", "Crocus_soil_temp1", code_analysis_4),
    grepl("TG4", code_analysis_4) ~ gsub("TG4", "Crocus_soil_temp4", code_analysis_4),
    grepl("Aluminium \\(Al\\)", code_analysis_4) ~ gsub("Aluminium \\(Al\\)", "Al", code_analysis_4),
    grepl("Fer \\(Fe\\)", code_analysis_4) ~ gsub("Fer \\(Fe\\)", "Fe", code_analysis_4),
    grepl("Silicium \\(Si\\)", code_analysis_4) ~ gsub("Silicium \\(Si\\)", "Si", code_analysis_4),
    TRUE ~ code_analysis_4
  )) %>%
  filter(grepl("Exchangeable_Al|Clay|N_total|Broyage_0_250µm|^CN$|CEC_cobaltihexamine|CaCO3_total|^Exchangeable_Ca$|Organic_C|Bulk_density|Exchangeable_Fe|Humidite_prelevement|Fine_silt|Coarse_silt|Exchangeable_Mg$|Exchangeable_Mn$|^SOM$|^P$|Stoniness|^Exchangeable_K$|Preparation|Fine_sand|Coarse_sand|Exchangeable_Na|C_treatment|pH$|Al$|Fe$|Si$", code_analysis_4) |
           (grepl("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1", code_analysis_4) &
           code_analysis == "fung_relative_abundance")) %>%

  filter(!grepl("CaCO3_total|Bulk_density|Humidite_prelevement|Stoniness|Preparation|Broyage_0_250µm|C_treatment", code_analysis_4)) %>%
  filter(!grepl("^Al|^Fe|^Si", code_analysis_4)) %>% # not enough data
  group_by(codeplot, code_analysis_4, ) %>%
  summarize( value_mean = mean(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%  
  mutate(value_mean =  case_when(
    !is.na(value_mean) & grepl("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1", code_analysis_4) ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
  )) %>%  
  tidyr::pivot_wider(names_from = code_analysis_4, values_from = value_mean) %>%
  mutate_at(vars(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1")), ~replace(., is.na(.), 0)) %>% 
  tidyr::drop_na()

df1 = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_all|Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")) %>%
  select(codeplot) %>% distinct()
df2 = orchamp_wide_pedo_agg
setdiff(df1$codeplot, df2$codeplot)


# Indirect ordination 
pca.fungi <- orchamp_wide_pedo_agg %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)
pca.pedo <- orchamp_wide_pedo_agg %>%
select(matches("Exchangeable_Al|Clay|N_total|Broyage_0_250µm|CN$|CEC_cobaltihexamine|CaCO3_total|Exchangeable_Ca$|Organic_C|Bulk_density|Exchangeable_Fe|Humidite_prelevement|Fine_silt|Coarse_silt|Exchangeable_Mg$|Exchangeable_Mn$|SOM$|P$|Stoniness|Exchangeable_K$|Fine_sand|Coarse_sand|Exchangeable_Na|C_treatment|pH$|Al$|Fe$|Si$")) %>% 
dudi.pca(., scannf = FALSE, nf = 2)

# Coinertia Analysis
(coia.fungipedo <- coinertia(pca.fungi, pca.pedo, scannf=FALSE))
(coia.fungisubpedo <- coinertia(pca.fungi_sub, pca.pedo, scannf=FALSE))
summary(coia.fungipedo)

colors_fungipedo = coia.fungipedo$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)


pedo_coia = s.arrow(coia.fungipedo$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), 
                    xlim = c(-0.65, 0.65), ylim = c(-0.65, 0.65) )
fung_coia = s.arrow(coia.fungipedo$c1, plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungipedo, plines = list(col = colors_fungipedo), 
                    xlim = c(-0.65, 0.65), ylim = c(-0.65, 0.65) )
superpose(pedo_coia, fung_coia, plot = TRUE)


# Mantel test 
orchamp_wide_pedo_agg_mantel = orchamp_metabar_blast_agg_years_2 %>%
  group_by(codeplot, code_analysis_2, code_analysis_3, code_analysis, file_name) %>%
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%   
  mutate(code_analysis_2 = case_when(
    file_name == "PhysicoChimie_20240711.json" ~ code_analysis,
    TRUE ~ code_analysis_2
  )) %>% 
  select(codeplot, value_mean, code_analysis, code_analysis_2, file_name) %>% 
  mutate(code_analysis_2 = case_when(
    grepl("Aluminium \\(Al\\) CEC", code_analysis_2) ~ gsub("Aluminium \\(Al\\) CEC", "Exchangeable_Al", code_analysis_2),
    grepl("Argile \\(< 2 µm\\)", code_analysis_2) ~ gsub("Argile \\(< 2 µm\\)", "Clay", code_analysis_2), 
    grepl("Azote \\(N\\) total", code_analysis_2) ~ gsub("Azote \\(N\\) total", "N_total", code_analysis_2),
    grepl("Broyage \\(< 250 µm\\)", code_analysis_2) ~ gsub("Broyage \\(< 250 µm\\)", "Broyage_0_250µm", code_analysis_2),
    grepl("C\\/N", code_analysis_2) ~ sub("C\\/N", "CN", code_analysis_2),
    grepl("CEC cobaltihexamine", code_analysis_2) ~ gsub("CEC cobaltihexamine", "CEC_cobaltihexamine", code_analysis_2),
    grepl("Calcaire \\(CaCO3\\) total", code_analysis_2) ~ gsub("Calcaire \\(CaCO3\\) total", "CaCO3_total", code_analysis_2),
    grepl("Calcium \\(Ca\\)", code_analysis_2) ~ gsub("Calcium \\(Ca\\)", "Exchangeable_Ca", code_analysis_2),
    grepl("Carbone \\(C\\) organique", code_analysis_2) ~ gsub("Carbone \\(C\\) organique", "Organic_C", code_analysis_2),
    grepl("Densité apparente", code_analysis_2) ~ gsub("Densité apparente", "Bulk_density", code_analysis_2),
    grepl("Fer \\(Fe\\) CEC", code_analysis_2) ~ gsub("Fer \\(Fe\\) CEC", "Exchangeable_Fe", code_analysis_2),
    grepl("Humidité prélèvement", code_analysis_2) ~ gsub("Humidité prélèvement", "Humidite_prelevement", code_analysis_2),
    grepl("Limons fins \\(2\\/20 µm\\)", code_analysis_2) ~ gsub("Limons fins \\(2\\/20 µm\\)", "Fine_silt", code_analysis_2),
    grepl("Limons grossiers \\(20\\/50 µm\\)", code_analysis_2) ~ gsub("Limons grossiers \\(20\\/50 µm\\)", "Coarse_silt", code_analysis_2),
    grepl("Magnésium \\(Mg\\)", code_analysis_2) ~ gsub("Magnésium \\(Mg\\)", "Exchangeable_Mg", code_analysis_2),
    grepl("Manganèse \\(Mn\\)", code_analysis_2) ~ gsub("Manganèse \\(Mn\\)", "Exchangeable_Mn", code_analysis_2),
    grepl("Matiere organique", code_analysis_2) ~ gsub("Matiere organique", "SOM", code_analysis_2),
    grepl("Phosphore \\(P2O5\\)", code_analysis_2) ~ gsub("Phosphore \\(P2O5\\)", "P", code_analysis_2),
    grepl("Pierrosité", code_analysis_2) ~ gsub("Pierrosité", "Stoniness", code_analysis_2),
    grepl("Potassium \\(K\\)", code_analysis_2) ~ gsub("Potassium \\(K\\)", "Exchangeable_K", code_analysis_2),
    grepl("Préparation", code_analysis_2) ~ gsub("Préparation", "Preparation", code_analysis_2),
    grepl("Sables fins \\(50\\/200 µm\\)", code_analysis_2) ~ gsub("Sables fins \\(50\\/200 µm\\)", "Fine_sand", code_analysis_2),
    grepl("Sables grossiers \\(200\\/2000 µm\\)", code_analysis_2) ~ gsub("Sables grossiers \\(200\\/2000 µm\\)", "Coarse_sand", code_analysis_2),
    grepl("Sodium \\(Na\\)", code_analysis_2) ~ gsub("Sodium \\(Na\\)", "Exchangeable_Na", code_analysis_2),
    grepl("Traitement Carbone", code_analysis_2) ~ gsub("Traitement Carbone", "C_treatment", code_analysis_2),
    grepl("DSN_T_ISBA", code_analysis_2) ~ gsub("DSN_T_ISBA", "DSN_T_ISBA", code_analysis_2),
    grepl("In situ soil subsurface temperature in degree Celsius", code_analysis_2) ~ gsub("In situ soil subsurface temperature in degree Celsius", "HOBO_temp_C", code_analysis_2),
    grepl("TG1", code_analysis_2) ~ gsub("TG1", "Crocus_soil_temp1", code_analysis_2),
    grepl("TG4", code_analysis_2) ~ gsub("TG4", "Crocus_soil_temp4", code_analysis_2),
    grepl("Aluminium \\(Al\\)", code_analysis_2) ~ gsub("Aluminium \\(Al\\)", "Al", code_analysis_2),
    grepl("Fer \\(Fe\\)", code_analysis_2) ~ gsub("Fer \\(Fe\\)", "Fe", code_analysis_2),
    grepl("Silicium \\(Si\\)", code_analysis_2) ~ gsub("Silicium \\(Si\\)", "Si", code_analysis_2),
    TRUE ~ code_analysis_2
  )) %>%
  filter(
    grepl("Al_CEC|Clay|N_total|Broyage_0_250µm|^CN$|CEC_cobaltihexamine|CaCO3_total|^Exchangeable_Ca$|Organic_C|Bulk_density|Fe_CEC|Humidite_prelevement|Fine_silt|Coarse_silt|Exchangeable_Mg$|Exchangeable_Mn$|^SOM$|^P$|Stoniness|^Exchangeable_K$|Preparation|Fine_sand|Coarse_sand|Exchangeable_Na|C_treatment|pH$|Al$|Fe$|Si$", code_analysis_2) |
      (code_analysis_2 %in% unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ]) &
         code_analysis == "fung_relative_abundance")) %>%
  filter(!grepl("CaCO3_total|Bulk_density|Humidite_prelevement|Stoniness|Preparation|Broyage_0_250µm|C_treatment", code_analysis_2)) %>%
  filter(!grepl("Al$|Fe$|Si$", code_analysis_2)) %>% 
  group_by(codeplot, code_analysis_2) %>%
  summarize( value_mean = mean(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%  
  mutate(value_mean =  case_when(
    !is.na(value_mean) & (code_analysis_2 %in% unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ])) ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
  )) %>%
  tidyr::pivot_wider(names_from = code_analysis_2, values_from = value_mean) %>%
  mutate_at(vars(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ]))), ~replace(., is.na(.), 0)) %>% 
  tidyr::drop_na()


fung_rel_abun_mantel = orchamp_wide_pedo_agg_mantel  %>% 
  select(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ]))) %>%
  mutate_if(is.numeric, ~(.x+5.247952)) %>%
  vegan::vegdist(., method = "bray")

pedo_rel_abun_mantel = orchamp_wide_pedo_agg_mantel  %>%
  select(-matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", code_analysis_2 ]))) %>%
  select(-codeplot) %>% 
  vegan::vegdist(., method = "bray")

fung_pedo_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  pedo_rel_abun_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_pedo_mantel 


# Figure 2b. Helotiales and Vegetation species dataset co-structure ------------------------

# create Vegetation species data subset
orchamp_wide_veg_fung = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_all|Helotiales_other|plant_relative_abundance|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")|
           code_analysis == "plant_relative_abundance") %>% 
  mutate(lb_nom = case_when(
    is.na(lb_nom) & !is.na(isolate) ~ code_analysis_4,
    TRUE ~ lb_nom
  )) %>% 
  filter(!lb_nom == "aaa - En attente de détermination") %>%
  group_by(codeplot, idtiers,  code_analysis_4, lb_nom) %>% 
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
  tidyr::pivot_wider(names_from = lb_nom, values_from = value_mean) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

# Indirect ordination 
pca.fungi <- orchamp_wide_veg_fung %>%
  select(matches("Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

pca.veg <- orchamp_wide_veg_fung %>%
  select(-matches("codeplot|idtiers|Helotiales_all|Glomeromycota|Helotiales_sp_1-39|T._maxilliforme_1-40|L._orchidicola_UAMH8152|P._finlandica_CBS444.86|M._vraolstadiae_G1|H._ericae_UAMH5828|C._rhizophila_CBS109839|Grass_root_mycorrhizal_PPO-1|Grass_root_mycorrhizal_PPO-2|M._variabilis_UAMH8861|A._macrosclerotiorum_Acema1|Helotiales_other|Calycina_sp._Bisps1|Leptodontidium_sp._EF-M|C._finlandica_PRF15|P._fortinii_PFO_F|P._fortinii_CBS179.46|P._ericae_AC2_IMI_392567|Ascomycota_sp._Da_Bac_I4|Ascomycota_sp._Da_PoT1_C7|Tapesia_sp._Da_MaP_I9|H._varicosporoides_Hymvar1|Leptodontidium_sp._MPI-SDFR-AT-0119|O._maius_OMA_B|C._variabilis_VPRI_42388_T|P._ericae_EriJKE06|P._ericae_LivJKC02|O._maius_OM19|Cadophora_sp._Fe06|Hyaloscyphaceae_sp._1M05|M._bicolor_E_Melbi2|M._variabilis_T_Melva1|L._orchidicola_UAMH5422|H._hepaticicola_UAMH7357|O._maius_Zn_Oidma1|P._subalpina_UAMH-11012|G._debralockiae_T6G9_T|P._scopiformis_Phisc1|Rhexocercosporidium_sp._Rhesp1|0")) %>%
  dudi.pca(., scale = FALSE, scannf = FALSE, nf = 2)

# Coinertia Analysis
coia.fungiveg <- coinertia(pca.fungi, pca.veg, scannf=FALSE)

colors_fungveg = coia.fungiveg$c1 %>%
  tibble::rownames_to_column(., "code_analysis_4") %>%
  merge(., colors_pca[, c("code_analysis_4", "color_endo", "color_myco", "color_lifestyle", "color_lifestyle_JA")], by="code_analysis_4", all.x=TRUE, sort=FALSE) %>%
  pull(color_lifestyle_JA)

veg_coia = s.arrow(coia.fungiveg$l1, plabels = list(box = list(draw = FALSE), optim = TRUE), xlim = c(-0.45, 0.45), ylim = c(-0.45, 0.45))
fung_coia = s.arrow(coia.fungiveg$c1, plabels = list(box = list(draw = FALSE), optim = TRUE), plabels.col = colors_fungveg, plines = list(col = colors_fungveg), xlim = c(-0.45, 0.45), ylim = c(-0.45, 0.45))
superpose(veg_coia, fung_coia, plot = TRUE)

# Mantel test
orchamp_wide_veg_fung_mantel = orchamp_metabar_blast_agg_years %>%
  filter((grepl("Helotiales_other|plant_relative_abundance|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
            code_analysis == "fung_relative_abundance")|
           code_analysis == "plant_relative_abundance") %>% 
  mutate(lb_nom = case_when(
    is.na(lb_nom) & !is.na(isolate) ~ code_analysis_3,
    TRUE ~ lb_nom
  )) %>% 
  filter(!lb_nom == "aaa - En attente de détermination") %>%
  group_by(codeplot, idtiers,  code_analysis_3, lb_nom) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>% 
  ungroup() %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & !code_analysis_3 == "plant_relative_abundance" ~ log10((value_mean + (0.0000113/2))),
    !is.na(value_mean) & code_analysis_3 == "plant_relative_abundance" ~ log10((value_mean + (0.001006036/2))),
    TRUE ~ value_mean
  )) %>%
  select(-code_analysis_3) %>%
  filter(!codeplot %in% c("CER_1865", "CER_1920", "CER_2075", "CER_2215", "CER_2290", "CER_2315", "DEN_1500", "DEN_1720", "DEN_1870", "DEN_2100", "DEN_2340", "DES_1500", "DES_1730", 
                          "DES_1970", "DES_2130", "DES_2330", "RAC_0520", "SEYS_0460", "SEYS_0615", "SEYS_0840", "SEYS_1035", "SEYS_1350", "SEYS_1500", "SEYS_1700", "SEYS_1870", "VAUL_0830",
                          "VAUL_1100", "VAUL_1290")) %>% # pas d'echantionnage pour la vegetation
  tidyr::pivot_wider(names_from = lb_nom, values_from = value_mean) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))


library(vegan)
fung_rel_abun_mantel = orchamp_wide_veg_fung_mantel %>% 
  select(matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  mutate_if(is.numeric, ~(.x+5.247953)) %>% 
  vegan::vegdist(., method = "bray")

veg_rel_abun_mantel = orchamp_wide_veg_fung_mantel %>%
  select(-matches(unique(orchamp_metabar_blast_agg_years_2[orchamp_metabar_blast_agg_years_2$group == "Helotiales_targeted" |orchamp_metabar_blast_agg_years_2$group == "Helotiales_other", isolate_code ]))) %>%
  select(-codeplot, -idtiers) %>%
  mutate_if(is.numeric, ~(.x+3.298417)) %>% 
  vegan::vegdist(., method = "bray")


filter(!grepl("Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", lb_nom)) %>% View()


fung_veg_mantel <- vegan::mantel(
  fung_rel_abun_mantel,
  veg_rel_abun_mantel,
  method = "pearson",
  permutations = 999,
  na.rm = TRUE)

fung_veg_mantel
