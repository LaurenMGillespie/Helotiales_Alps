# Required packages -------------------------------------
library("dplyr")
require("data.table")
library("bestNormalize")
library("MuMIn")
library("lme4")
library("lmerTest")
library("forcats")
library("RColorBrewer")
library("ggtext")
library("grid")
library("ggh4x")
library("ggimage")
library("legendry")
library("scales")
library("patchwork")

# Load data ---------------------------------------------

orchamp_metabar_blast_agg_years <- readRDS("data/orchamp_metabar_blast_agg_years_20260327.rds")
orchamp_metabar_blast_agg_years_2 <- readRDS("orchamp_metabar_blast_agg_years_idtiers_20260327.rds")

# curves

# Run dredges--------------------------------------------

## format dataframe
orchamp_agg_all_wide = orchamp_metabar_blast_agg_years %>%
  dplyr::select(id_gradient, codeplot, idtiers, value_mean, code_analysis, code_analysis_3) %>%
  filter(!code_analysis_3 =="Helotialies_all") %>%
  mutate(code_analysis_3 =  case_when(
    code_analysis_3 == "CN_2mm" ~ "CN",
    code_analysis_3 == "Phosphorus" ~ "P",
    code_analysis_3 == "pH - H2O 1:4_2mm" ~ "pH - H2O 1:4_2mm",
    code_analysis_3 == "pH - H2O 1:5_2mm" ~ "pH - H2O 1:5_2mm",
    TRUE ~ code_analysis_3)) %>%
  filter(grepl("vegetation_Ericaceae_percent|vegetation_nonAMplant_percent|vegetation_species_richness|^CN$|^P$|pH - H2O 1:4_2mm|pH - H2O 1:5_2mm|total_vegetation_cover|vegetation_shannon_diversity|vegetation_nonAMplant_percent|vegetation_Ericaceae_percent", code_analysis_3)| 
           (grepl("Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1", code_analysis_3) &
              code_analysis == "fung_relative_abundance")
  ) %>%
  group_by(id_gradient, codeplot, idtiers, code_analysis_3) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & code_analysis_3 %in% c("Helotiales_other", "isolate_1-39", "isolate_1-40", "AF214576.1", "AF486119.1", "AJ292199.1", "AJ319077.1", "AY078144.1", "AY394921.1", "AY599235.1", "AY599236.1", "AY762619.1", "AY762620.1", "Acema1", "Bisps1", "Cadsp1", "DQ148411.1", "DQ485204.1", "EF446149.1", "EF446150.1", "EU938674.1", "EU938675.1", "FM180476.1", "FM180477.1", "FM207641.1", "FN178469.1", "FN178471.1", "Glomeromycota", "Hymvar1", "JAHLEZ010000069.1", "JF304942.1", "JN882306.1", "KF836358.1", "KT827742.1", "KT835080.1", "KU382495.1", "KU886589.1", "LT821494.1", "MF943108.1", "MH860824.1", "MT276008.1", "Melbi2", "Melva1", "OL742649.1", "OM238143.1", "Oidma1", "PRJEB12348", "PRJNA407395", "Rhesp1", "AY176753.1", "Phisc1") ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean
  )) %>%
  tidyr::pivot_wider(names_from = code_analysis_3, values_from = value_mean) %>%
  mutate_at(vars(matches("Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1")), ~replace(., is.na(.), 0)) %>%
  mutate(pH = case_when(
    (!is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ (`pH - H2O 1:4_2mm` + `pH - H2O 1:5_2mm`)/2,
    (is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:5_2mm`,
    (!is.na(`pH - H2O 1:4_2mm`) & is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:4_2mm`,
    TRUE ~ NA
  )) %>%
  dplyr::select(-`pH - H2O 1:4_2mm`, -`pH - H2O 1:5_2mm`) %>%
  tidyr::drop_na()

## add meteo variables 
orchamp_agg_all_wide_meteo = orchamp_metabar_blast_agg_years %>%
  dplyr::select(id_gradient, codeplot, value_mean, code_analysis, code_analysis_3, 
  ) %>%
  filter(!code_analysis_3 =="Helotialies_all") %>%
  mutate(code_analysis_3 =  case_when(
    code_analysis_3 == "TairMean" ~ "Mean_air_temp",
    code_analysis_3 == "net_solar_radiation" ~ "Net_solar_radiation",
    code_analysis_3 == "CWS" ~ "CWS",
    TRUE ~ code_analysis_3)) %>%
  filter(grepl("Mean_air_temp|Net_solar_radiation|CWS", code_analysis_3)) %>%
  group_by(id_gradient, codeplot, code_analysis_3) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  tidyr::pivot_wider(names_from = code_analysis_3, values_from = value_mean) %>%
  tidyr::drop_na()

setDT(orchamp_agg_all_wide)[orchamp_agg_all_wide_meteo, c("Mean_air_temp", "Net_solar_radiation", "CWS") :=
                              .(Mean_air_temp, Net_solar_radiation, CWS), on=c("id_gradient", "codeplot")]



## conduct the tansformations
library(bestNormalize)
orchamp_agg_all_wide_trans = orchamp_agg_all_wide %>%
 mutate(Mean_air_temp_log = log(Mean_air_temp),
        Mean_air_temp_poly = poly(Mean_air_temp, 2, raw=T)[,2],
        Net_solar_radiation_orderNorm = predict(orderNorm(Net_solar_radiation, n_logit_fit = min(length(Net_solar_radiation), 10000), warn = TRUE), newdata = NULL, inverse = FALSE, warn = TRUE),
        Net_solar_radiation_log = log(Net_solar_radiation),
        Net_solar_radiation_poly = poly(Net_solar_radiation, 2, raw=T)[,2],
        CWS_sqrt = sqrt(CWS),
        CWS_log = log(CWS),
        CWS_poly = poly(CWS, 2, raw=T)[,2],
        CN_yeojohnson = predict(yeojohnson(CN, eps = 0.001, standardize = TRUE), newdata = NULL, inverse = FALSE),
        CN_log = log(CN),
        CN_poly = poly(CN, 2, raw=T)[,2],
        P_Log_b = predict(log_x(P, a = NULL, b = 10, standardize = TRUE, eps = 0.001, warn = TRUE), newdata = NULL, inverse = FALSE),
        P_log = log(P),
        P_poly = poly(P, 2, raw=T)[,2],
        pH_asinh = asinh(pH),
        pH_log = log(pH),
        pH_poly = poly(pH, 2, raw=T)[,2],
        total_vegetation_cover_log = log(total_vegetation_cover),
        total_vegetation_cover_poly = poly(total_vegetation_cover, 2, raw=T)[,2],
        veg_shannon_div_boxcox = predict(boxcox(vegetation_shannon_diversity+0.0000000001, standardize = TRUE), newdata = NULL, inverse = FALSE),
        veg_shannon_div_log = log(vegetation_shannon_diversity+0.0000000001),
        veg_shannon_div_poly = poly(vegetation_shannon_diversity, 2, raw=T)[,2],
        veg_species_richness_boxcox = predict(boxcox(vegetation_species_richness+0.0000000001, standardize = TRUE), newdata = NULL, inverse = FALSE),
        veg_species_richness_log = log(vegetation_species_richness),
        veg_species_richness_poly = poly(vegetation_species_richness, 2, raw=T)[,2],
        veg_Ericaceae_per_orderNorm = predict(orderNorm(vegetation_Ericaceae_percent, n_logit_fit = min(length(vegetation_Ericaceae_percent), 10000), warn = TRUE), newdata = NULL, inverse = FALSE, warn = TRUE),
        veg_Ericaceae_per_log = log(vegetation_Ericaceae_percent+0.0000000001),
        veg_Ericaceae_per_poly = poly(vegetation_Ericaceae_percent, 2, raw=T)[,2],
        veg_nonAMplant_per_sqrt = sqrt(vegetation_nonAMplant_percent),
        veg_nonAMplant_per_log = log(vegetation_nonAMplant_percent+0.0000000001),
        veg_nonAMplant_per_poly = poly(vegetation_nonAMplant_percent, 2, raw=T)[,2]
 )

orchamp_agg_all_wide_trans_scaled = orchamp_agg_all_wide_trans %>%
  mutate_at(vars(-matches("id_gradient|codeplot|idtiers|Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1")), ~scale(., center = TRUE, scale = TRUE)[,1])
  

souches <- orchamp_metabar_blast_agg_years %>%
  filter(code_analysis == "fung_relative_abundance") %>%
  filter(!code_analysis_3 == "Helotiales_all") %>%
  select(code_analysis_3) %>%
  distinct() %>%
  pull(code_analysis_3)


## run dredges for each lineage of Helotiales fungi
dredges <- list()

for(i in souches) {
  options(na.action = "na.omit")
  try(
    { 
  globalmodel_indv = lmer(data=orchamp_agg_all_wide_trans_scaled,
                     get(i) ~ Mean_air_temp + Mean_air_temp_log + Mean_air_temp_poly +
                       Net_solar_radiation_orderNorm + Net_solar_radiation_log + Net_solar_radiation_poly +
                       CWS_sqrt + CWS_log + CWS_poly +
                       CN_yeojohnson + CN_log + CN_poly +
                       P_Log_b + P_log + P_poly +
                       pH_asinh + pH_log + pH_poly +
                       veg_species_richness_boxcox + veg_species_richness_log + veg_species_richness_poly +
                       veg_Ericaceae_per_orderNorm + veg_Ericaceae_per_log + veg_Ericaceae_per_poly +
                       veg_nonAMplant_per_sqrt + veg_nonAMplant_per_log + veg_nonAMplant_per_poly +
                       (1|codeplot) + (1|id_gradient),
                     REML = F, na.action=na.fail
  )


  nullmodel <- MuMIn:::.nullFitRE(globalmodel_indv)
  options(na.action <-  "na.fail")
  dredges[[i]] <- dredge(globalmodel_indv,
                         m.lim = c(0,5),
                         rank = "AIC", 
                         extra =list(R2 = function(x) {
                           r.squaredGLMM(x, null = nullmodel)}),
                         subset =
                           !({Mean_air_temp} && {Mean_air_temp_log} && {Mean_air_temp_poly}) &&
                           !({Mean_air_temp} && {Mean_air_temp_log}) &&
                           !({Mean_air_temp}&& {Mean_air_temp_poly}) &&
                           !({Mean_air_temp_log} && {Mean_air_temp_poly}) &&

                           !({Net_solar_radiation_orderNorm} && {Net_solar_radiation_log} && {Net_solar_radiation_poly}) &&
                           !({Net_solar_radiation_orderNorm} && {Net_solar_radiation_log}) &&
                           !({Net_solar_radiation_orderNorm}&& {Net_solar_radiation_poly}) &&
                           !({Net_solar_radiation_log} && {Net_solar_radiation_poly}) &&

                           !({CWS_sqrt} && {CWS_log} && {CWS_poly}) &&
                           !({CWS_sqrt} && {CWS_log}) &&
                           !({CWS_sqrt}&& {CWS_poly}) &&
                           !({CWS_log} && {CWS_poly}) &&

                           !({CN_yeojohnson} && {CN_log} && {CN_poly}) &&
                           !({CN_yeojohnson} && {CN_log}) &&
                           !({CN_yeojohnson}&& {CN_poly}) &&
                           !({CN_log} && {CN_poly}) &&

                           !({P_Log_b} && {P_log} && {P_poly}) &&
                           !({P_Log_b} && {P_log}) &&
                           !({P_Log_b}&& {P_poly}) &&
                           !({P_log} && {P_poly}) &&

                           !({pH_asinh} && {pH_log} && {pH_poly}) &&
                           !({pH_asinh} && {pH_log}) &&
                           !({pH_asinh}&& {pH_poly}) &&
                           !({pH_log} && {pH_poly}) &&
                           
                           !({veg_species_richness_boxcox} && {veg_species_richness_log} && {veg_species_richness_poly}) &&
                           !({veg_species_richness_boxcox} && {veg_species_richness_log}) &&
                           !({veg_species_richness_boxcox}&& {veg_species_richness_poly}) &&
                           !({veg_species_richness_log} && {veg_species_richness_poly}) &&

                           !({veg_Ericaceae_per_orderNorm} && {veg_Ericaceae_per_log} && {veg_Ericaceae_per_poly}) &&
                           !({veg_Ericaceae_per_orderNorm} && {veg_Ericaceae_per_log}) &&
                           !({veg_Ericaceae_per_orderNorm}&& {veg_Ericaceae_per_poly}) &&
                           !({veg_Ericaceae_per_log} && {veg_Ericaceae_per_poly}) &&

                           !({veg_nonAMplant_per_sqrt} && {veg_nonAMplant_per_log} && {veg_nonAMplant_per_poly}) &&
                           !({veg_nonAMplant_per_sqrt} && {veg_nonAMplant_per_log}) &&
                           !({veg_nonAMplant_per_sqrt}&& {veg_nonAMplant_per_poly}) &&
                           !({veg_nonAMplant_per_log} && {veg_nonAMplant_per_poly})
                         )
              }
  )           
  }

saveRDS(dredges, file = "data/dredges.rds")



# For the dredge on Helotiales_all
## format dataframe

orchamp_agg_all_wide_global = orchamp_metabar_blast_agg_years %>%
  select(id_gradient, codeplot, idtiers, value_mean, code_analysis, code_analysis_3, group 
  ) %>%
  filter(!code_analysis_3 =="Helotialies_all") %>%
  mutate(code_analysis_3 =  case_when(
    code_analysis_3 == "TairMean" ~ "Mean_air_temp",
    code_analysis_3 == "net_solar_radiation" ~ "Net_solar_radiation",
    code_analysis_3 == "CWS" ~ "CWS",
    code_analysis_3 == "CN_2mm" ~ "CN",
    code_analysis_3 == "Phosphorus" ~ "P",
    code_analysis_3 == "pH - H2O 1:4_2mm" ~ "pH - H2O 1:4_2mm",
    code_analysis_3 == "pH - H2O 1:5_2mm" ~ "pH - H2O 1:5_2mm",
    grepl("Helotiales", group) ~ "fung_relative_abundance",
    TRUE ~ code_analysis_3)) %>%
  filter(grepl("vegetation_Ericaceae_percent|vegetation_nonAMplant_percent|vegetation_species_richness|^CN$|^P$|pH - H2O 1:4_2mm|pH - H2O 1:5_2mm|total_vegetation_cover|vegetation_shannon_diversity|vegetation_nonAMplant_percent|vegetation_Ericaceae_percent", code_analysis_3)| 
           (grepl("Helotiales", group) &
              code_analysis == "fung_relative_abundance")) %>%
  group_by(id_gradient, codeplot, idtiers, code_analysis_3) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(value_mean =  case_when(
    !is.na(value_mean) & code_analysis_3 %in% c("fung_relative_abundance") ~ log10((value_mean + (0.0000113/2))),
    TRUE ~ value_mean)) %>%
  tidyr::pivot_wider(names_from = code_analysis_3, values_from = value_mean) %>%
  mutate_at(vars(matches("fung_relative_abundance")), ~replace(., is.na(.), 0)) %>%
  mutate(pH = case_when(
    (!is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ (`pH - H2O 1:4_2mm` + `pH - H2O 1:5_2mm`)/2,
    (is.na(`pH - H2O 1:4_2mm`) & !is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:5_2mm`,
    (!is.na(`pH - H2O 1:4_2mm`) & is.na(`pH - H2O 1:5_2mm`)) ~ `pH - H2O 1:4_2mm`,
    TRUE ~ NA)) %>%
  select(-`pH - H2O 1:4_2mm`, -`pH - H2O 1:5_2mm`) %>%
  tidyr::drop_na()

setDT(orchamp_agg_all_wide_global)[orchamp_agg_all_wide_meteo, c("Mean_air_temp", "Net_solar_radiation", "CWS") :=
                                     .(Mean_air_temp, Net_solar_radiation, CWS), on=c("id_gradient", "codeplot")]

## conduct the tansformations
orchamp_agg_all_wide_global_trans = orchamp_agg_all_wide_global %>%
  mutate(Mean_air_temp_log = log(Mean_air_temp),
         Mean_air_temp_poly = poly(Mean_air_temp, 2, raw=T)[,2],
         Net_solar_radiation_orderNorm = predict(orderNorm(Net_solar_radiation, n_logit_fit = min(length(Net_solar_radiation), 10000), warn = TRUE), newdata = NULL, inverse = FALSE, warn = TRUE),
         Net_solar_radiation_log = log(Net_solar_radiation),
         Net_solar_radiation_poly = poly(Net_solar_radiation, 2, raw=T)[,2],
         CWS_sqrt = sqrt(CWS),
         CWS_log = log(CWS),
         CWS_poly = poly(CWS, 2, raw=T)[,2],
         CN_yeojohnson = predict(yeojohnson(CN, eps = 0.001, standardize = TRUE), newdata = NULL, inverse = FALSE),
         CN_log = log(CN),
         CN_poly = poly(CN, 2, raw=T)[,2],
         P_Log_b = predict(log_x(P, a = NULL, b = 10, standardize = TRUE, eps = 0.001, warn = TRUE), newdata = NULL, inverse = FALSE),
         P_log = log(P),
         P_poly = poly(P, 2, raw=T)[,2],
         pH_asinh = asinh(pH),
         pH_log = log(pH),
         pH_poly = poly(pH, 2, raw=T)[,2],
         total_vegetation_cover_log = log(total_vegetation_cover),
         total_vegetation_cover_poly = poly(total_vegetation_cover, 2, raw=T)[,2],
         veg_shannon_div_boxcox = predict(boxcox(vegetation_shannon_diversity+0.0000000001, standardize = TRUE), newdata = NULL, inverse = FALSE),
         veg_shannon_div_log = log(vegetation_shannon_diversity+0.0000000001),
         veg_shannon_div_poly = poly(vegetation_shannon_diversity, 2, raw=T)[,2],
         veg_species_richness_boxcox = predict(boxcox(vegetation_species_richness+0.0000000001, standardize = TRUE), newdata = NULL, inverse = FALSE),
         veg_species_richness_log = log(vegetation_species_richness),
         veg_species_richness_poly = poly(vegetation_species_richness, 2, raw=T)[,2],
         veg_Ericaceae_per_orderNorm = predict(orderNorm(vegetation_Ericaceae_percent, n_logit_fit = min(length(vegetation_Ericaceae_percent), 10000), warn = TRUE), newdata = NULL, inverse = FALSE, warn = TRUE),
         veg_Ericaceae_per_log = log(vegetation_Ericaceae_percent+0.0000000001),
         veg_Ericaceae_per_poly = poly(vegetation_Ericaceae_percent, 2, raw=T)[,2],
         veg_nonAMplant_per_sqrt = sqrt(vegetation_nonAMplant_percent),
         veg_nonAMplant_per_log = log(vegetation_nonAMplant_percent+0.0000000001),
         veg_nonAMplant_per_poly = poly(vegetation_nonAMplant_percent, 2, raw=T)[,2]
  )

orchamp_agg_all_wide_global_trans_scaled = orchamp_agg_all_wide_global_trans %>%
  mutate_at(vars(-matches("id_gradient|codeplot|idtiers|Helotiales_other|isolate_1-39|isolate_1-40|AF214576.1|AF486119.1|AJ292199.1|AJ319077.1|AY078144.1|AY394921.1|AY599235.1|AY599236.1|AY762619.1|AY762620.1|Acema1|Bisps1|Cadsp1|DQ148411.1|DQ485204.1|EF446149.1|EF446150.1|EU938674.1|EU938675.1|FM180476.1|FM180477.1|FM207641.1|FN178469.1|FN178471.1|Glomeromycota|Hymvar1|JAHLEZ010000069.1|JF304942.1|JN882306.1|KF836358.1|KT827742.1|KT835080.1|KU382495.1|KU886589.1|LT821494.1|MF943108.1|MH860824.1|MT276008.1|Melbi2|Melva1|OL742649.1|OM238143.1|Oidma1|PRJEB12348|PRJNA407395|Rhesp1|AY176753.1|Phisc1")), ~scale(., center = TRUE, scale = TRUE)[,1])


options(na.action = "na.omit")

## Run the dredge for Helotiales_all
globalmodel = lmer(data=orchamp_agg_all_wide_global_trans_scaled,
                   fung_relative_abundance ~ Mean_air_temp + Mean_air_temp_log + Mean_air_temp_poly +
                     Net_solar_radiation_orderNorm + Net_solar_radiation_log + Net_solar_radiation_poly +
                     CWS_sqrt + CWS_log + CWS_poly +
                     CN_yeojohnson + CN_log + CN_poly +
                     P_Log_b + P_log + P_poly +
                     pH_asinh + pH_log + pH_poly +
                     veg_species_richness_boxcox + veg_species_richness_log + veg_species_richness_poly +
                     veg_Ericaceae_per_orderNorm + veg_Ericaceae_per_log + veg_Ericaceae_per_poly +
                     veg_nonAMplant_per_sqrt + veg_nonAMplant_per_log + veg_nonAMplant_per_poly +
                     (1|codeplot) + (1|id_gradient),  REML = F
)


nullmodel <- MuMIn:::.nullFitRE(globalmodel)
options(na.action = "na.fail")
globalmodel_dredge <- dredge(globalmodel,
                             m.lim = c(0,5),
                             rank = "AIC", 
                             extra =list(R2 = function(x) {
                               r.squaredGLMM(x, null = nullmodel)}),
                             subset =
                               !({Mean_air_temp} && {Mean_air_temp_log} && {Mean_air_temp_poly}) &&
                               !({Mean_air_temp} && {Mean_air_temp_log}) &&
                               !({Mean_air_temp}&& {Mean_air_temp_poly}) &&
                               !({Mean_air_temp_log} && {Mean_air_temp_poly}) &&

                               !({Net_solar_radiation_orderNorm} && {Net_solar_radiation_log} && {Net_solar_radiation_poly}) &&
                               !({Net_solar_radiation_orderNorm} && {Net_solar_radiation_log}) &&
                               !({Net_solar_radiation_orderNorm}&& {Net_solar_radiation_poly}) &&
                               !({Net_solar_radiation_log} && {Net_solar_radiation_poly}) &&

                               !({CWS_sqrt} && {CWS_log} && {CWS_poly}) &&
                               !({CWS_sqrt} && {CWS_log}) &&
                               !({CWS_sqrt}&& {CWS_poly}) &&
                               !({CWS_log} && {CWS_poly}) &&

                               !({CN_yeojohnson} && {CN_log} && {CN_poly}) &&
                               !({CN_yeojohnson} && {CN_log}) &&
                               !({CN_yeojohnson}&& {CN_poly}) &&
                               !({CN_log} && {CN_poly}) &&

                               !({P_Log_b} && {P_log} && {P_poly}) &&
                               !({P_Log_b} && {P_log}) &&
                               !({P_Log_b}&& {P_poly}) &&
                               !({P_log} && {P_poly}) &&

                               !({pH_asinh} && {pH_log} && {pH_poly}) &&
                               !({pH_asinh} && {pH_log}) &&
                               !({pH_asinh}&& {pH_poly}) &&
                               !({pH_log} && {pH_poly}) &&
                               
                               !({veg_species_richness_boxcox} && {veg_species_richness_log} && {veg_species_richness_poly}) &&
                               !({veg_species_richness_boxcox} && {veg_species_richness_log}) &&
                               !({veg_species_richness_boxcox}&& {veg_species_richness_poly}) &&
                               !({veg_species_richness_log} && {veg_species_richness_poly}) &&

                               !({veg_Ericaceae_per_orderNorm} && {veg_Ericaceae_per_log} && {veg_Ericaceae_per_poly}) &&
                               !({veg_Ericaceae_per_orderNorm} && {veg_Ericaceae_per_log}) &&
                               !({veg_Ericaceae_per_orderNorm}&& {veg_Ericaceae_per_poly}) &&
                               !({veg_Ericaceae_per_log} && {veg_Ericaceae_per_poly}) &&

                               !({veg_nonAMplant_per_sqrt} && {veg_nonAMplant_per_log} && {veg_nonAMplant_per_poly}) &&
                               !({veg_nonAMplant_per_sqrt} && {veg_nonAMplant_per_log}) &&
                               !({veg_nonAMplant_per_sqrt}&& {veg_nonAMplant_per_poly}) &&
                               !({veg_nonAMplant_per_log} && {veg_nonAMplant_per_poly})
)

saveRDS(globalmodel_dredge, file = "data/globalmodel_dredge.rds")

# Figure 4 ----------------------------------------------
# import dredge results if necessary
# dredges <- readRDS("data/dredges_20260424.rds")
# globalmodel_dredge <- readRDS("data/globalmodel_dredge_20260424.rds")

# combine dredges and globalmodel_dredge
globalmodel_dredge_list = list(globalmodel_dredge)
names(globalmodel_dredge_list) <- c("Helotiales_all")
dredges <- append(dredges, globalmodel_dredge_list)

dredges_resultats <- purrr::map_df(dredges, ~as.data.frame(.x[1,]), .id="id")

dredges_resultats = dredges_resultats %>%
  merge(x=., y=colors_pca, by.x="id", by.y="code_analysis_3", all.x = TRUE) %>%
  mutate(isolate_nom_Bruyant_2025 = case_when(
    id == "AJ292199.1" ~ "AJ292199.1 Meliniomyces vraolstadiae strain G1", 
    id == "Melbi2" ~ "Melbi2 Meliniomyces bicolor E", 
    id == "DQ485204.1" ~ "DQ485204.1 Cadophora finlandica isolate PRF15", 
    id == "AF486119.1" ~ "AF486119.1 Phialophora finlandia strain CBS 444.86 16S", 
    id == "Melva1" ~ "Melva1 Meliniomyces variabilis T", 
    id == "AY762619.1" ~ "AY762619.1 Meliniomyces variabilis strain UAMH 8861", 
    id == "KT835080.1" ~ "KT835080.1 Pezoloma ericae LivJKC02", 
    id == "AJ319077.1" ~ "AJ319077.1 Hymenoscyphus ericae isolate UAMH 5828", 
    id == "KT827742.1" ~ "KT827742.1 Pezoloma ericae EriJKE06", 
    id == "OM238143.1" ~ "OM238143.1 Hyaloscypha hepaticicola strain UAMH 7357", 
    id == "KF836358.1" ~ "KF836358 Cairneyella variabilis VPRI 42388 T", 
    id == "AY599235.1" ~ "AY599235.1 Grass root mycorrhizal sp. PPO-1", 
    id == "AY599236.1" ~ "AY599236.1 Grass root mycorrhizal sp. PPO-2", 
    id == "Phisc1" ~ "Phisc1 Phialocephala scopiformis", 
    id == "FN178471.1" ~ "FN178471.1 Ascomycota sp. Da MaP I9 Tapesia sp.", 
    id == "FN178469.1" ~ "FN178469.1 Ascomycota sp. Da PoT1 C7", 
    id == "Acema1" ~ "Acema1 Acephala macrosclerotiorum EW76-UTF0540", 
    id == "EU938674.1" ~ "EU938674.1 Phialocephala fortinii strain CBS 179.46", 
    id == "EF446149.1" ~ "EF446149.1 Phialocephala fortinii isolate PFO F", 
    id == "PRJEB12348" ~ "GCA_900073065.1 Phialocephala subalpina UAMH 11012", 
    id == "FM207641.1" ~ "FM207641.1 Ascomycota sp. Da Bac I4", 
    id == "KU886589.1" ~ "KU886589.1 Cadophora sp. Fe06", 
    id == "isolate_1-39" ~ "isolate_1-39", 
    id == "Rhesp1" ~ "Rhesp1 Rhexocercosporidium sp. MPI-PUGE-AT-0058", 
    id == "DQ148411.1" ~ "DQ148411.1 Leptodontidium sp. EF-M", 
    id == "JAHLEZ010000069.1" ~ "GCA_020726505.1 Leptodontidium sp. MPI-SDFR-AT-0119", 
    id == "AF214576.1" ~ "AF214576.1 Leptodontidium orchidicola strain UAMH 8152", 
    id == "OL742649.1" ~ "OL742649.1 Leptodophora orchidicola UAMH 5422", 
    id == "Oidma1" ~ "Oidma1 Oidiodendron maius Zn", 
    id == "KU382495.1" ~ "KU382495.1 Oidiodendron maius OM19", 
    id == "JF304942.1" ~ "JF304942.1 Oidiodendron maius isolate OMA B", 
    id == "AY176753.1" ~ "AY176753.1 Cryptosporiopsis rhizophila strain CBS 109839", 
    id == "Hymvar1" ~ "Hymvar1 Hymenoscyphus varicosporoides PMI453", 
    id == "MF943108.1" ~ "MF943108.1 Hyaloscyphaceae sp. isolate 1M05", 
    id == "PRJNA407395" ~ "GCA_002803225.1 Gamarada debralockiae T6G9 T", 
    id == "FM180476.1" ~ "FM180476.1 Pezeloma ericae AC2 IMI 392567", 
    id == "isolate_1-40" ~ "isolate_1-40", 
    id == "Bisps1" ~ "Bisps1 Calycina sp. as Bisporella", 
    TRUE ~ NA
  ),
  order_phylo = case_when(
    id == "Helotiales_all" ~ 0.25, 
    id == "Helotiales_other" ~ 0.5, 
    id == "Glomeromycota" ~ 0.75, 
    id == "AF486119.1" ~ 1,
    id == "DQ485204.1" ~ 2,
    id == "AY762619.1" ~ 3,
    id == "AJ292199.1" ~ 4,
    id == "AJ319077.1" ~ 5,
    id == "KT835080.1" ~ 6,
    id == "Melbi2" ~ 7,
    id == "OM238143.1" ~ 8,
    id == "KF836358.1" ~ 9,
    id == "Melva1" ~ 10,
    id == "KT827742.1" ~ 11,
    id == "FN178469.1" ~ 12,
    id == "PRJEB12348" ~ 13,
    id == "EU938674.1" ~ 14,
    id == "AY599235.1" ~ 15,
    id == "JAHLEZ010000069.1" ~ 16,
    id == "EF446149.1" ~ 17,
    id == "FM207641.1" ~ 18,
    id == "FN178471.1" ~ 19,
    id == "AY599236.1" ~ 20,
    id == "OL742649.1" ~ 21,
    id == "Acema1" ~ 22,
    id == "KU886589.1" ~ 23,
    id == "isolate_1-39" ~ 24,
    id == "Rhesp1" ~ 25,
    id == "DQ148411.1" ~ 26,
    id == "AF214576.1" ~ 27,
    id == "Phisc1" ~ 28,
    id == "MF943108.1" ~ 29,
    id == "Hymvar1" ~ 30,
    id == "Bisps1" ~ 31,
    id == "AY176753.1" ~ 32,
    id == "Oidma1" ~ 33,
    id == "JF304942.1" ~ 34,
    id == "FM180476.1" ~ 35,
    id == "PRJNA407395" ~ 36,
    id == "KU382495.1" ~ 37,
    id == "isolate_1-40" ~ 38,
    TRUE ~ NA
  ),
  clade = case_when(
    grepl("AJ292199.1|Melbi2|DQ485204.1|AF486119.1|Melva1|AY762619.1|KT835080.1|AJ319077.1|KT827742.1|OM238143.1|KF836358.1", id)  ~ "Hyaloscyphoid clade", 
    grepl("AY599235.1|AY599236.1|Phisc1|FN178471.1|FN178469.1|Acema1|EU938674.1|EF446149.1|PRJEB12348|FM207641.1|KU886589.1|isolate_1-39|Rhesp1|DQ148411.1|JAHLEZ010000069.1|AF214576.1|OL742649.1", id)  ~ "Mollisioid clade", 
    grepl("Oidma1|KU382495.1|JF304942.1|AY176753.1|Hymvar1|MF943108.1|PRJNA407395|FM180476.1|Bisps1", id)  ~ "Pezizelloid clade", 
    grepl("isolate_1-40", id)  ~ "Vandijckellaceae – Calloriaceae clade", 
    grepl("Helotiales_all|Helotiales_other|Glomeromycota", id)  ~ "Misc", 
    TRUE ~ NA
  ),
  lifestyle = case_when(
    id == "Helotiales_all" ~ "Mixed",
    TRUE ~ lifestyle
  ),
  color_lifestyle = case_when(
    id == "Helotiales_all" ~ "black",
    TRUE ~ color_lifestyle
  ),
  lifestyle_JA = case_when(
    id == "Helotiales_all" ~ "Mixed",
    TRUE ~ lifestyle_JA
  ),
  color_lifestyle_JA = case_when(
    id == "Helotiales_all" ~ "black",
    TRUE ~ color_lifestyle_JA
  ))

dredges_summary_list = list()
library(purrr)
for( i in names(dredges)) {
    dredges_summary_list[[i]] <- summary(get.models(dredges[[i]], subset = 1)[[1]])
    }

dredges_summary <- purrr::map_df(dredges_summary_list, ~as.data.frame(.x$coefficients), .id="id") %>%
  tibble::rownames_to_column(., "variable") %>%
  mutate(variable_2 = gsub("_poly|_log|_yeojohnson|_boxcox|_orderNorm|_asinh|_sqrt", "", variable)) %>%
  mutate(variable_2 = gsub("\\...[^.]*$", "", gsub("_poly|_log|_yeojohnson|_boxcox|_orderNorm|_asinh|_sqrt", "", variable)),
         #transformation = stringr::str_extract(variable, "(?<=\\_).*?(?=\\...)")
         transformation = case_when(
           grepl("_baseline_|_full_model_|id_gradient|(Intercept)", variable) ~ NA,
           variable %in% unique(id) ~ NA,
           grepl("_log", variable) ~ "log",
           grepl("_poly", variable) ~ "poly",
           #grepl("_yeojohnson|_orderNorm|_asinh|_boxcox", variable) ~ "gaussian",
           TRUE ~ "gaussian"
         ),
         p_value_asterisks = case_when(
           (`Pr(>|t|)` > 0.05 & `Pr(>|t|)` < 0.1) ~ ".",
           (`Pr(>|t|)` > 0.01 & `Pr(>|t|)` < 0.05) ~ "*",
           (`Pr(>|t|)` > 0.001 & `Pr(>|t|)` < 0.01) ~ "**",
           (`Pr(>|t|)` < 0.001) ~ "***",
           TRUE ~ NA
         )) %>%
  select(id, variable_2, transformation, variable, everything()) %>%
  mutate(isolate_nom_Bruyant_2025 = case_when(
    id == "AJ292199.1" ~ "AJ292199.1 Meliniomyces vraolstadiae strain G1", 
    id == "Melbi2" ~ "Melbi2 Meliniomyces bicolor E", 
    id == "DQ485204.1" ~ "DQ485204.1 Cadophora finlandica isolate PRF15", 
    id == "AF486119.1" ~ "AF486119.1 Phialophora finlandia strain CBS 444.86 16S", 
    id == "Melva1" ~ "Melva1 Meliniomyces variabilis T", 
    id == "AY762619.1" ~ "AY762619.1 Meliniomyces variabilis strain UAMH 8861", 
    id == "KT835080.1" ~ "KT835080.1 Pezoloma ericae LivJKC02", 
    id == "AJ319077.1" ~ "AJ319077.1 Hymenoscyphus ericae isolate UAMH 5828", 
    id == "KT827742.1" ~ "KT827742.1 Pezoloma ericae EriJKE06", 
    id == "OM238143.1" ~ "OM238143.1 Hyaloscypha hepaticicola strain UAMH 7357", 
    id == "KF836358.1" ~ "KF836358 Cairneyella variabilis VPRI 42388 T", 
    id == "AY599235.1" ~ "AY599235.1 Grass root mycorrhizal sp. PPO-1", 
    id == "AY599236.1" ~ "AY599236.1 Grass root mycorrhizal sp. PPO-2", 
    id == "Phisc1" ~ "Phisc1 Phialocephala scopiformis", 
    id == "FN178471.1" ~ "FN178471.1 Ascomycota sp. Da MaP I9 Tapesia sp.", 
    id == "FN178469.1" ~ "FN178469.1 Ascomycota sp. Da PoT1 C7", 
    id == "Acema1" ~ "Acema1 Acephala macrosclerotiorum EW76-UTF0540", 
    id == "EU938674.1" ~ "EU938674.1 Phialocephala fortinii strain CBS 179.46", 
    id == "EF446149.1" ~ "EF446149.1 Phialocephala fortinii isolate PFO F", 
    id == "PRJEB12348" ~ "GCA_900073065.1 Phialocephala subalpina UAMH 11012", 
    id == "FM207641.1" ~ "FM207641.1 Ascomycota sp. Da Bac I4", 
    id == "KU886589.1" ~ "KU886589.1 Cadophora sp. Fe06", 
    id == "isolate_1-39" ~ "isolate_1-39", 
    id == "Rhesp1" ~ "Rhesp1 Rhexocercosporidium sp. MPI-PUGE-AT-0058", 
    id == "DQ148411.1" ~ "DQ148411.1 Leptodontidium sp. EF-M", 
    id == "JAHLEZ010000069.1" ~ "GCA_020726505.1 Leptodontidium sp. MPI-SDFR-AT-0119", 
    id == "AF214576.1" ~ "AF214576.1 Leptodontidium orchidicola strain UAMH 8152", 
    id == "OL742649.1" ~ "OL742649.1 Leptodophora orchidicola UAMH 5422", 
    id == "Oidma1" ~ "Oidma1 Oidiodendron maius Zn", 
    id == "KU382495.1" ~ "KU382495.1 Oidiodendron maius OM19", 
    id == "JF304942.1" ~ "JF304942.1 Oidiodendron maius isolate OMA B", 
    id == "AY176753.1" ~ "AY176753.1 Cryptosporiopsis rhizophila strain CBS 109839", 
    id == "Hymvar1" ~ "Hymvar1 Hymenoscyphus varicosporoides PMI453", 
    id == "MF943108.1" ~ "MF943108.1 Hyaloscyphaceae sp. isolate 1M05", 
    id == "PRJNA407395" ~ "GCA_002803225.1 Gamarada debralockiae T6G9 T", 
    id == "FM180476.1" ~ "FM180476.1 Pezeloma ericae AC2 IMI 392567", 
    id == "isolate_1-40" ~ "isolate_1-40", 
    id == "Bisps1" ~ "Bisps1 Calycina sp. as Bisporella", 
    id == "Helotiales_all" ~ "Helotiales_all", 
    TRUE ~ NA
  ),
  order_phylo = case_when(
    id == "Helotiales_all" ~ 0.25, 
    id == "Helotiales_other" ~ 0.5, 
    id == "Glomeromycota" ~ 0.75, 
    id == "AF486119.1" ~ 1,
    id == "DQ485204.1" ~ 2,
    id == "AY762619.1" ~ 3,
    id == "AJ292199.1" ~ 4,
    id == "AJ319077.1" ~ 5,
    id == "KT835080.1" ~ 6,
    id == "Melbi2" ~ 7,
    id == "OM238143.1" ~ 8,
    id == "KF836358.1" ~ 9,
    id == "Melva1" ~ 10,
    id == "KT827742.1" ~ 11,
    id == "FN178469.1" ~ 12,
    id == "PRJEB12348" ~ 13,
    id == "EU938674.1" ~ 14,
    id == "AY599235.1" ~ 15,
    id == "JAHLEZ010000069.1" ~ 16,
    id == "EF446149.1" ~ 17,
    id == "FM207641.1" ~ 18,
    id == "FN178471.1" ~ 19,
    id == "AY599236.1" ~ 20,
    id == "OL742649.1" ~ 21,
    id == "Acema1" ~ 22,
    id == "KU886589.1" ~ 23,
    id == "isolate_1-39" ~ 24,
    id == "Rhesp1" ~ 25,
    id == "DQ148411.1" ~ 26,
    id == "AF214576.1" ~ 27,
    id == "Phisc1" ~ 28,
    id == "MF943108.1" ~ 29,
    id == "Hymvar1" ~ 30,
    id == "Bisps1" ~ 31,
    id == "AY176753.1" ~ 32,
    id == "Oidma1" ~ 33,
    id == "JF304942.1" ~ 34,
    id == "FM180476.1" ~ 35,
    id == "PRJNA407395" ~ 36,
    id == "KU382495.1" ~ 37,
    id == "isolate_1-40" ~ 38,
    TRUE ~ NA
  ),
  clade = case_when(
    grepl("AJ292199.1|Melbi2|DQ485204.1|AF486119.1|Melva1|AY762619.1|KT835080.1|AJ319077.1|KT827742.1|OM238143.1|KF836358.1", id)  ~ "Hyaloscyphoid clade", 
    grepl("AY599235.1|AY599236.1|Phisc1|FN178471.1|FN178469.1|Acema1|EU938674.1|EF446149.1|PRJEB12348|FM207641.1|KU886589.1|isolate_1-39|Rhesp1|DQ148411.1|JAHLEZ010000069.1|AF214576.1|OL742649.1", id)  ~ "Mollisioid clade", 
    grepl("Oidma1|KU382495.1|JF304942.1|AY176753.1|Hymvar1|MF943108.1|PRJNA407395|FM180476.1|Bisps1", id)  ~ "Pezizelloid clade", 
    grepl("isolate_1-40", id)  ~ "Vandijckellaceae – Calloriaceae clade", 
    grepl("Helotiales_all|Helotiales_other|Glomeromycota", id)  ~ "Misc", 
    TRUE ~ NA
  ))


# get the variable importance for each lineage

#Helotiales_all
# explore the best model identified by the dredge
subset(dredges$Helotiales_all, delta < 1) %>% View()
dredges$Helotiales_all[1,]
Helotiales_all_avg = model.avg(subset(dredges$Helotiales_all, delta < 2))
Helotiales_all_avg$coefficients %>% View()
Helotiales_all_avg$importance
summary_ = summary(Helotiales_all_avg)
summary_
# use the best model to create a LMM
model_Helotiales_all = lmer(data=orchamp_agg_all_wide_global_trans_scaled,
                            fung_relative_abundance ~ CN_log + Mean_air_temp_log + pH_log + veg_Ericaceae_per_log + veg_nonAMplant_per_log + 
                            (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Helotiales_all_exp <- explain(model_Helotiales_all,
                            data = as.data.frame(orchamp_agg_all_wide_global_trans_scaled[, c("id_gradient", "codeplot", "fung_relative_abundance", "CN_log", 
                                                                                "Mean_air_temp_log", "pH_log", "veg_Ericaceae_per_log",
                                                                                "veg_nonAMplant_per_log")]), 
                            y = as.numeric(orchamp_agg_all_wide_global_trans_scaled$fung_relative_abundance),
                            label = "Linear Regression Model") #, is_multiclass = TRUE
plot(model_profile(Helotiales_all_exp))  


#isolate_1-39
subset(dredges$`isolate_1-39`, delta < 1) %>% View()
dredges$`isolate_1-39`[1,]
isolate_1_39_avg = model.avg(subset(dredges$`isolate_1-39`, delta < 2))
isolate_1_39_avg$coefficients %>% View()
isolate_1_39_avg$importance
summary_ = summary(isolate_1_39_avg)
summary_
model_isolate_1_39 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        `isolate_1-39` ~ Mean_air_temp_log + Net_solar_radiation_orderNorm + pH_poly + 
                          veg_Ericaceae_per_log + veg_nonAMplant_per_poly + (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

isolate_1_39_exp <- explain(model_isolate_1_39,
                     data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "isolate_1-39", "Mean_air_temp_log", 
                                                                                "Net_solar_radiation_orderNorm", "pH_poly",  
                                                                                  "veg_Ericaceae_per_log",  "veg_nonAMplant_per_poly")]), 
                     y = as.numeric(orchamp_agg_all_wide_trans_scaled$`isolate_1-39`),
                     label = "Linear Regression Model")
plot(model_profile(isolate_1_39_exp))  


#isolate_1-40
subset(dredges$`isolate_1-40`, delta < 2) %>% View()
dredges$`isolate_1-40`[1,]
isolate_1_40_avg = model.avg(subset(dredges$`isolate_1-40`, delta < 2))
isolate_1_40_avg$coefficients %>% View()
summary(isolate_1_40_avg)

model_isolate_1_40 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        `isolate_1-40` ~ Mean_air_temp_log + pH_log + veg_Ericaceae_per_poly + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

isolate_1_40_exp <- explain(model_isolate_1_40,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot","isolate_1-40", "Mean_air_temp_log", "pH_log", "veg_Ericaceae_per_poly", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$`isolate_1-40`),
                          label = "Linear Regression Model")
plot(model_profile(isolate_1_40_exp))  



#AF214576.1
subset(dredges$AF214576.1, delta < 2) %>% View()
dredges$AF214576.1[1,]
AF214576.1_avg = model.avg(subset(dredges$AF214576.1, delta < 2))
AF214576.1_avg$coefficients %>% View()
summary(AF214576.1_avg)

model_AF214576.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AF214576.1 ~ Net_solar_radiation_poly + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AF214576.1_exp <- explain(model_AF214576.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Net_solar_radiation_poly", "veg_nonAMplant_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AF214576.1),
                          label = "Linear Regression Model")
plot(model_profile(AF214576.1_exp))  


#AF486119.1
subset(dredges$AF486119.1, delta < 2) %>% View()
dredges$AF486119.1[1,]
AF486119.1_avg = model.avg(subset(dredges$AF486119.1, delta < 2))
AF486119.1_avg$coefficients %>% View()
summary(AF486119.1_avg)

model_AF486119.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AF486119.1 ~ P_poly + pH_log + veg_nonAMplant_per_sqrt  +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AF486119.1_exp <- explain(model_AF486119.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AF486119.1", "P_poly", "pH_log", "veg_nonAMplant_per_sqrt")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AF486119.1),
                          label = "Linear Regression Model")
plot(model_profile(AF486119.1_exp))  


#AJ292199.1
subset(dredges$AJ292199.1, delta < 2) %>% View()
dredges$AJ292199.1[1,]
AJ292199.1_avg = model.avg(subset(dredges$AJ292199.1, delta < 2))
AJ292199.1_avg$coefficients %>% View()
summary(AJ292199.1_avg)

model_AJ292199.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AJ292199.1 ~ CN_poly + P_poly +  pH_log + 
                          veg_Ericaceae_per_poly + veg_nonAMplant_per_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AJ292199.1_exp <- explain(model_AJ292199.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AJ292199.1", "CN_poly", "P_poly",  "pH_log", 
                                                                                       "veg_Ericaceae_per_poly", "veg_nonAMplant_per_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AJ292199.1),
                          label = "Linear Regression Model")
plot(model_profile(AJ292199.1_exp))  


#AJ319077.1
subset(dredges$AJ319077.1, delta < 2) %>% View()
dredges$AJ319077.1[1,]
AJ319077.1_avg = model.avg(subset(dredges$AJ319077.1, delta < 2))
AJ319077.1_avg$coefficients %>% View()
summary(AJ319077.1_avg)

model_AJ319077.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AJ319077.1 ~ Mean_air_temp_log + Net_solar_radiation_orderNorm + veg_Ericaceae_per_orderNorm +
                          veg_nonAMplant_per_poly + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AJ319077.1_exp <- explain(model_AJ319077.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AJ319077.1", "Mean_air_temp_log", "Net_solar_radiation_orderNorm", "veg_Ericaceae_per_orderNorm",
                                                                                       "veg_nonAMplant_per_poly", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AJ319077.1),
                          label = "Linear Regression Model")
plot(model_profile(AJ319077.1_exp))  


#AY176753.1
subset(dredges$AY176753.1, delta < 2) %>% View()
dredges$AY176753.1[1,]
AY176753.1_avg = model.avg(subset(dredges$AY176753.1, delta < 2))
AY176753.1_avg$coefficients %>% View()
summary(AY176753.1_avg)

model_AY176753.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AY176753.1 ~ CN_poly + CWS_sqrt + Mean_air_temp_poly +
                        veg_Ericaceae_per_log + veg_species_richness_boxcox +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AY176753.1_exp <- explain(model_AY176753.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AY176753.1", "CN_poly", "CWS_sqrt", "Mean_air_temp_poly",
                                                                                       "veg_Ericaceae_per_log", "veg_species_richness_boxcox")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AY176753.1),
                          label = "Linear Regression Model")
plot(model_profile(AY176753.1_exp))  


#AY599235.1
subset(dredges$AY599235.1, delta < 2) %>% View()
dredges$AY599235.1[1,]
AY599235.1_avg = model.avg(subset(dredges$AY599235.1, delta < 2))
AY599235.1_avg$coefficients %>% View()
summary(AY599235.1_avg)

model_AY599235.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AY599235.1 ~ CWS_log + Mean_air_temp_log + Net_solar_radiation_log +  
                          pH_log + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AY599235.1_exp <- explain(model_AY599235.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AY599235.1", "CWS_log", "Mean_air_temp_log", "Net_solar_radiation_log",  
                                                                                       "pH_log", "veg_nonAMplant_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AY599235.1),
                          label = "Linear Regression Model")
plot(model_profile(AY599235.1_exp))  


#AY599236.1
subset(dredges$AY599236.1, delta < 2) %>% View()
dredges$AY599236.1[1,]
AY599236.1_avg = model.avg(subset(dredges$AY599236.1, delta < 2))
AY599236.1_avg$coefficients %>% View()
summary(AY599236.1_avg)

model_AY599236.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AY599236.1 ~ CN_poly + Mean_air_temp_log +  pH_poly + veg_nonAMplant_per_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AY599236.1_exp <- explain(model_AY599236.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AY599236.1", "CN_poly",
                                                                                     "Mean_air_temp_log",  "pH_poly", "veg_nonAMplant_per_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AY599236.1),
                          label = "Linear Regression Model")
plot(model_profile(AY599236.1_exp))  


#AY762619.1
subset(dredges$AY762619.1, delta < 2) %>% View()
dredges$AY762619.1[1,]
AY762619.1_avg = model.avg(subset(dredges$AY762619.1, delta < 2))
AY762619.1_avg$coefficients %>% View()
summary(AY762619.1_avg)

model_AY762619.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        AY762619.1 ~ CN_poly + Mean_air_temp_log +  pH_asinh +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

AY762619.1_exp <- explain(model_AY762619.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "AY762619.1", "CN_poly", "Mean_air_temp_log", "pH_asinh")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$AY762619.1),
                          label = "Linear Regression Model")
plot(model_profile(AY762619.1_exp))  


#Glomeromycota
subset(dredges$Glomeromycota, delta < 2) %>% View()
dredges$Glomeromycota[1,]
Glomeromycota_avg = model.avg(subset(dredges$Glomeromycota, delta < 2))
Glomeromycota_avg$coefficients %>% View()
summary(Glomeromycota_avg)

model_Glomeromycota = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Glomeromycota ~ CN_poly + pH_log + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Glomeromycota_exp <- explain(model_Glomeromycota,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Glomeromycota", "CN_poly", 
                                                                                     "pH_log", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Glomeromycota),
                          label = "Linear Regression Model")
plot(model_profile(Glomeromycota_exp))  


#Acema1
subset(dredges$Acema1, delta < 2) %>% View()
dredges$Acema1[1,]
Acema1_avg = model.avg(subset(dredges$Acema1, delta < 2))
Acema1_avg$coefficients %>% View()
summary(Acema1_avg)

model_Acema1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Acema1 ~ CN_yeojohnson + Mean_air_temp_log +  pH_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Acema1_exp <- explain(model_Acema1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Acema1", 
                                                                                     "CN_yeojohnson", "Mean_air_temp_log",  "pH_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Acema1),
                          label = "Linear Regression Model")
plot(model_profile(Acema1_exp))  


#Helotiales_other
subset(dredges$Helotiales_other, delta < 2) %>% View()
dredges$Helotiales_other[1,]
Helotiales_other_avg = model.avg(subset(dredges$Helotiales_other, delta < 2))
Helotiales_other_avg$coefficients %>% View()
summary(Helotiales_other_avg)

model_Helotiales_other = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Helotiales_other ~ Mean_air_temp_log + pH_poly + veg_Ericaceae_per_log + veg_nonAMplant_per_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Helotiales_other_exp <- explain(model_Helotiales_other,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Helotiales_other", "Mean_air_temp_log", 
                                                                                     "pH_poly", "veg_Ericaceae_per_log", "veg_nonAMplant_per_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Helotiales_other),
                          label = "Linear Regression Model")
plot(model_profile(Helotiales_other_exp))  


#Bisps1
subset(dredges$Bisps1, delta < 2) %>% View()
dredges$Bisps1[1,]
Bisps1_avg = model.avg(subset(dredges$Bisps1, delta < 2))
Bisps1_avg$coefficients %>% View()
summary(Bisps1_avg)

model_Bisps1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Bisps1 ~ CN_yeojohnson + CWS_poly + Net_solar_radiation_poly + P_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Bisps1_exp <- explain(model_Bisps1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Bisps1", "CN_yeojohnson", 
                                                                                     "CWS_poly", "Net_solar_radiation_poly", "P_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Bisps1),
                          label = "Linear Regression Model")
plot(model_profile(Bisps1_exp))  


#DQ148411.1
subset(dredges$DQ148411.1, delta < 2) %>% View()
dredges$DQ148411.1[1,]
DQ148411.1_avg = model.avg(subset(dredges$DQ148411.1, delta < 2))
DQ148411.1_avg$coefficients %>% View()
summary(DQ148411.1_avg)

model_DQ148411.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        DQ148411.1 ~ Mean_air_temp_log + pH_log + veg_Ericaceae_per_log + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

DQ148411.1_exp <- explain(model_DQ148411.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "DQ148411.1", "Mean_air_temp_log", "pH_log", 
                                                                                     "veg_Ericaceae_per_log", "veg_nonAMplant_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$DQ148411.1),
                          label = "Linear Regression Model")
plot(model_profile(DQ148411.1_exp))  


#DQ485204.1
subset(dredges$DQ485204.1, delta < 2) %>% View()
dredges$DQ485204.1[1,]
DQ485204.1_avg = model.avg(subset(dredges$DQ485204.1, delta < 2))
DQ485204.1_avg$coefficients %>% View()
summary(DQ485204.1_avg)

model_DQ485204.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        DQ485204.1 ~ CN_log + Mean_air_temp_log + pH_poly + veg_Ericaceae_per_log + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

DQ485204.1_exp <- explain(model_DQ485204.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "DQ485204.1", "CN_log", "Mean_air_temp_log", 
                                                                                     "pH_poly", "veg_Ericaceae_per_log", "veg_nonAMplant_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$DQ485204.1),
                          label = "Linear Regression Model")
plot(model_profile(DQ485204.1_exp))  


#EF446149.1
subset(dredges$EF446149.1, delta < 2) %>% View()
dredges$EF446149.1[1,]
EF446149.1_avg = model.avg(subset(dredges$EF446149.1, delta < 2))
EF446149.1_avg$coefficients %>% View()
summary(EF446149.1_avg)

model_EF446149.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        EF446149.1 ~ CN_poly +  CWS_poly + Net_solar_radiation_log + veg_Ericaceae_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

EF446149.1_exp <- explain(model_EF446149.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "EF446149.1", "CN_poly",  
                                                                                     "CWS_poly", "Net_solar_radiation_log", "veg_Ericaceae_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$EF446149.1),
                          label = "Linear Regression Model")
plot(model_profile(EF446149.1_exp))  


#EU938674.1
subset(dredges$EU938674.1, delta < 2) %>% View()
dredges$EU938674.1[1,]
EU938674.1_avg = model.avg(subset(dredges$EU938674.1, delta < 2))
EU938674.1_avg$coefficients %>% View()
summary(EU938674.1_avg)

model_EU938674.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        EU938674.1 ~ CN_log + Mean_air_temp_log +  pH_poly + veg_Ericaceae_per_log + veg_species_richness_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

EU938674.1_exp <- explain(model_EU938674.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "EU938674.1", "CN_log", "Mean_air_temp_log",  
                                                                                     "pH_poly", "veg_Ericaceae_per_log", "veg_species_richness_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$EU938674.1),
                          label = "Linear Regression Model")
plot(model_profile(EU938674.1_exp))  


#FM180476.1
subset(dredges$FM180476.1, delta < 2) %>% View()
dredges$FM180476.1[1,]
FM180476.1_avg = model.avg(subset(dredges$FM180476.1, delta < 2))
FM180476.1_avg$coefficients %>% View()
summary(FM180476.1_avg)

model_FM180476.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        FM180476.1 ~ CN_yeojohnson + pH_poly + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

FM180476.1_exp <- explain(model_FM180476.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "FM180476.1", "CN_yeojohnson", "pH_poly", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$FM180476.1),
                          label = "Linear Regression Model")
plot(model_profile(FM180476.1_exp))  


#FM207641.1
subset(dredges$FM207641.1, delta < 2) %>% View()
dredges$FM207641.1[1,]
FM207641.1_avg = model.avg(subset(dredges$FM207641.1, delta < 2))
FM207641.1_avg$coefficients %>% View()
summary(FM207641.1_avg)

model_FM207641.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        FM207641.1 ~ CN_poly + CWS_log + Mean_air_temp_poly + veg_nonAMplant_per_log + veg_species_richness_boxcox +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

FM207641.1_exp <- explain(model_FM207641.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "CN_poly", 
                                                                                     "CWS_log", "Mean_air_temp_poly", "veg_nonAMplant_per_log", "veg_species_richness_boxcox")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$FM207641.1),
                          label = "Linear Regression Model")
plot(model_profile(FM207641.1_exp))  


#FN178469.1
subset(dredges$FN178469.1, delta < 2) %>% View()
dredges$FN178469.1[1,]
FN178469.1_avg = model.avg(subset(dredges$FN178469.1, delta < 2))
FN178469.1_avg$coefficients %>% View()
summary(FN178469.1_avg)

model_FN178469.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        FN178469.1 ~ CWS_log + Net_solar_radiation_log + pH_log + veg_species_richness_boxcox +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

FN178469.1_exp <- explain(model_FN178469.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "FN178469.1", "CWS_log", 
                                                                                     "Net_solar_radiation_log", "pH_log", "veg_species_richness_boxcox")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$FN178469.1),
                          label = "Linear Regression Model")
plot(model_profile(FN178469.1_exp))  


#FN178471.1
subset(dredges$FN178471.1, delta < 2) %>% View()
dredges$FN178471.1[1,]
FN178471.1_avg = model.avg(subset(dredges$FN178471.1, delta < 2))
FN178471.1_avg$coefficients %>% View()
summary(FN178471.1_avg)

model_FN178471.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        FN178471.1 ~ Mean_air_temp_log + pH_poly + veg_Ericaceae_per_poly + veg_nonAMplant_per_log + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

FN178471.1_exp <- explain(model_FN178471.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "FN178471.1", "Mean_air_temp_log", "pH_poly", 
                                                                                     "veg_Ericaceae_per_poly", "veg_nonAMplant_per_log", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$FN178471.1),
                          label = "Linear Regression Model")
plot(model_profile(FN178471.1_exp))  


#Hymvar1
subset(dredges$Hymvar1, delta < 2) %>% View()
dredges$Hymvar1[1,]
Hymvar1_avg = model.avg(subset(dredges$Hymvar1, delta < 2))
Hymvar1_avg$coefficients %>% View()
summary(Hymvar1_avg)

model_Hymvar1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Hymvar1 ~ CN_yeojohnson + CWS_log + pH_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Hymvar1_exp <- explain(model_Hymvar1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Hymvar1", 
                                                                                     "CN_yeojohnson", "CWS_log", "pH_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Hymvar1),
                          label = "Linear Regression Model")
plot(model_profile(Hymvar1_exp))  


#JAHLEZ010000069.1
subset(dredges$JAHLEZ010000069.1, delta < 2) %>% View()
dredges$JAHLEZ010000069.1[1,]
JAHLEZ010000069.1_avg = model.avg(subset(dredges$JAHLEZ010000069.1, delta < 2))
JAHLEZ010000069.1_avg$coefficients %>% View()
summary(JAHLEZ010000069.1_avg)

model_JAHLEZ010000069.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        JAHLEZ010000069.1 ~ CN_log + Mean_air_temp_log + Net_solar_radiation_orderNorm + pH_log + veg_species_richness_boxcox +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

JAHLEZ010000069.1_exp <- explain(model_JAHLEZ010000069.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "JAHLEZ010000069.1", "CN_log", "Mean_air_temp_log", 
                                                                                     "Net_solar_radiation_orderNorm", "pH_log", "veg_species_richness_boxcox")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$JAHLEZ010000069.1),
                          label = "Linear Regression Model")
plot(model_profile(JAHLEZ010000069.1_exp))  


#JF304942.1
subset(dredges$JF304942.1, delta < 2) %>% View()
dredges$JF304942.1[1,]
JF304942.1_avg = model.avg(subset(dredges$JF304942.1, delta < 2))
JF304942.1_avg$coefficients %>% View()
summary(JF304942.1_avg)

model_JF304942.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        JF304942.1 ~ CN_poly + CWS_poly + Mean_air_temp_log + pH_log + veg_Ericaceae_per_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

JF304942.1_exp <- explain(model_JF304942.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "CN_poly", "CWS_poly", 
                                                                                     "Mean_air_temp_log", "pH_log", "veg_Ericaceae_per_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$JF304942.1),
                          label = "Linear Regression Model")
plot(model_profile(JF304942.1_exp))  


#KF836358.1
subset(dredges$KF836358.1, delta < 2) %>% View()
dredges$KF836358.1[1,]
KF836358.1_avg = model.avg(subset(dredges$KF836358.1, delta < 2))
KF836358.1_avg$coefficients %>% View()
summary(KF836358.1_avg)

model_KF836358.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        KF836358.1 ~ CWS_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

KF836358.1_exp <- explain(model_KF836358.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "CWS_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$KF836358.1),
                          label = "Linear Regression Model")
plot(model_profile(KF836358.1_exp))  


#KT827742.1
subset(dredges$KT827742.1, delta < 2) %>% View()
dredges$KT827742.1[1,]
KT827742.1_avg = model.avg(subset(dredges$KT827742.1, delta < 2))
KT827742.1_avg$coefficients %>% View()
summary(KT827742.1_avg)

model_KT827742.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        KT827742.1 ~ veg_Ericaceae_per_orderNorm +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

KT827742.1_exp <- explain(model_KT827742.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "KT827742.1", "veg_Ericaceae_per_orderNorm")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$KT827742.1),
                          label = "Linear Regression Model")
plot(model_profile(KT827742.1_exp))  


#KT835080.1
subset(dredges$KT835080.1, delta < 2) %>% View()
dredges$KT835080.1[1,]
KT835080.1_avg = model.avg(subset(dredges$KT835080.1, delta < 2))
KT835080.1_avg$coefficients %>% View()
summary(KT835080.1_avg)

model_KT835080.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        KT835080.1 ~ CN_poly + CWS_poly + Mean_air_temp + pH_asinh + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

KT835080.1_exp <- explain(model_KT835080.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "KT835080.1", "CN_poly", "CWS_poly", 
                                                                                     "Mean_air_temp", "pH_asinh", "veg_nonAMplant_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$KT835080.1),
                          label = "Linear Regression Model")
plot(model_profile(KT835080.1_exp))  


#KU382495.1
subset(dredges$KU382495.1, delta < 2) %>% View()
dredges$KU382495.1[1,]
KU382495.1_avg = model.avg(subset(dredges$KU382495.1, delta < 2))
KU382495.1_avg$coefficients %>% View()
summary(KU382495.1_avg)

model_KU382495.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        KU382495.1 ~ CN_poly + pH_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

KU382495.1_exp <- explain(model_KU382495.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "KU382495.1", 
                                                                                     "CN_poly", "pH_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$KU382495.1),
                          label = "Linear Regression Model")
plot(model_profile(KU382495.1_exp))  


#KU886589.1
subset(dredges$KU886589.1, delta < 2) %>% View()
dredges$KU886589.1[1,]
KU886589.1_avg = model.avg(subset(dredges$KU886589.1, delta < 2))
KU886589.1_avg$coefficients %>% View()
summary(KU886589.1_avg)

model_KU886589.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        KU886589.1 ~ pH_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

KU886589.1_exp <- explain(model_KU886589.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "KU886589.1", "pH_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$KU886589.1),
                          label = "Linear Regression Model")
plot(model_profile(KU886589.1_exp))  


#MF943108.1
subset(dredges$MF943108.1, delta < 2) %>% View()
dredges$MF943108.1[1,]
MF943108.1_avg = model.avg(subset(dredges$MF943108.1, delta < 2))
MF943108.1_avg$coefficients %>% View()
summary(MF943108.1_avg)

model_MF943108.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        MF943108.1 ~ CN_poly + Mean_air_temp_poly + Net_solar_radiation_orderNorm + pH_poly + veg_species_richness_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

MF943108.1_exp <- explain(model_MF943108.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "MF943108.1", "CN_poly", "Mean_air_temp_poly", 
                                                                                     "Net_solar_radiation_orderNorm", "pH_poly", "veg_species_richness_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$MF943108.1),
                          label = "Linear Regression Model")
plot(model_profile(MF943108.1_exp))  
View(model_profile(MF943108.1_exp)$agr_profiles)


#Melbi2
subset(dredges$Melbi2, delta < 2) %>% View()
dredges$Melbi2[1,]
Melbi2_avg = model.avg(subset(dredges$Melbi2, delta < 2))
Melbi2_avg$coefficients %>% View()
summary(Melbi2_avg)

model_Melbi2 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Melbi2 ~ P_poly + pH_log + veg_Ericaceae_per_poly + veg_nonAMplant_per_poly + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Melbi2_exp <- explain(model_Melbi2,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Melbi2", "P_poly", "pH_log", 
                                                                                     "veg_Ericaceae_per_poly", "veg_nonAMplant_per_poly", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Melbi2),
                          label = "Linear Regression Model")
plot(model_profile(Melbi2_exp))  


#Melva1
subset(dredges$Melva1, delta < 2) %>% View()
dredges$Melva1[1,]
Melva1_avg = model.avg(subset(dredges$Melva1, delta < 2))
Melva1_avg$coefficients %>% View()
summary(Melva1_avg)

model_Melva1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Melva1 ~ CWS_log + P_poly + pH_log + veg_Ericaceae_per_poly + veg_species_richness_log +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Melva1_exp <- explain(model_Melva1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "CWS_log", "P_poly", "pH_log", 
                                                                                     "veg_Ericaceae_per_poly", "veg_species_richness_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Melva1),
                          label = "Linear Regression Model")
plot(model_profile(Melva1_exp))  


#OL742649.1
subset(dredges$OL742649.1, delta < 2) %>% View()
dredges$OL742649.1[1,]
OL742649.1_avg = model.avg(subset(dredges$OL742649.1, delta < 2))
OL742649.1_avg$coefficients %>% View()
summary(OL742649.1_avg)

model_OL742649.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        OL742649.1 ~ Net_solar_radiation_orderNorm + veg_Ericaceae_per_orderNorm + veg_nonAMplant_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

OL742649.1_exp <- explain(model_OL742649.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "OL742649.1", "Net_solar_radiation_orderNorm", 
                                                                                     "veg_Ericaceae_per_orderNorm", "veg_nonAMplant_per_poly" )]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$OL742649.1),
                          label = "Linear Regression Model")
plot(model_profile(OL742649.1_exp))  


#OM238143.1
subset(dredges$OM238143.1, delta < 2) %>% View()
dredges$OM238143.1[1,]
OM238143.1_avg = model.avg(subset(dredges$OM238143.1, delta < 2))
OM238143.1_avg$coefficients %>% View()
summary(OM238143.1_avg)

model_OM238143.1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        OM238143.1 ~ CN_yeojohnson + Mean_air_temp_log + P_poly + pH_poly + veg_nonAMplant_per_log + 
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

OM238143.1_exp <- explain(model_OM238143.1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "OM238143.1", "CN_yeojohnson", "Mean_air_temp_log", 
                                                                                     "P_poly", "pH_poly", "veg_nonAMplant_per_log")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$OM238143.1),
                          label = "Linear Regression Model")
plot(model_profile(OM238143.1_exp))  


#Oidma1
subset(dredges$Oidma1, delta < 2) %>% View()
dredges$Oidma1[1,]
Oidma1_avg = model.avg(subset(dredges$Oidma1, delta < 2))
Oidma1_avg$coefficients %>% View()
summary(Oidma1_avg)

model_Oidma1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Oidma1 ~ CN_yeojohnson + Net_solar_radiation_orderNorm + pH_poly + veg_Ericaceae_per_orderNorm +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Oidma1_exp <- explain(model_Oidma1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Oidma1", "CN_yeojohnson", 
                                                                                     "Net_solar_radiation_orderNorm", "pH_poly", "veg_Ericaceae_per_orderNorm")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Oidma1),
                          label = "Linear Regression Model")
plot(model_profile(Oidma1_exp))  


#PRJEB12348
subset(dredges$PRJEB12348, delta < 2) %>% View()
dredges$PRJEB12348[1,]
PRJEB12348_avg = model.avg(subset(dredges$PRJEB12348, delta < 2))
PRJEB12348_avg$coefficients %>% View()
summary(PRJEB12348_avg)

model_PRJEB12348 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        PRJEB12348 ~ CN_poly + veg_Ericaceae_per_log + veg_nonAMplant_per_poly + veg_species_richness_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

PRJEB12348_exp <- explain(model_PRJEB12348,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "PRJEB12348", "CN_poly", "veg_Ericaceae_per_log", 
                                                                                     "veg_nonAMplant_per_poly", "veg_species_richness_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$PRJEB12348),
                          label = "Linear Regression Model")
plot(model_profile(PRJEB12348_exp))  


#PRJNA407395
subset(dredges$PRJNA407395, delta < 2) %>% View()
dredges$PRJNA407395[1,]
PRJNA407395_avg = model.avg(subset(dredges$PRJNA407395, delta < 2))
PRJNA407395_avg$coefficients %>% View()
summary(PRJNA407395_avg)

model_PRJNA407395 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        PRJNA407395 ~ CN_yeojohnson + Mean_air_temp_log + Net_solar_radiation_poly + veg_Ericaceae_per_log + veg_species_richness_boxcox +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

PRJNA407395_exp <- explain(model_PRJNA407395,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "PRJNA407395", "CN_yeojohnson", "Mean_air_temp_log", 
                                                                                     "Net_solar_radiation_poly", "veg_Ericaceae_per_log", "veg_species_richness_boxcox")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$PRJNA407395),
                          label = "Linear Regression Model")
plot(model_profile(PRJNA407395_exp))  


#Phisc1
subset(dredges$Phisc1, delta < 2) %>% View()
dredges$Phisc1[1,]
Phisc1_avg = model.avg(subset(dredges$Phisc1, delta < 2))
Phisc1_avg$coefficients %>% View()
summary(Phisc1_avg)

model_Phisc1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Phisc1 ~ P_poly + pH_log + veg_Ericaceae_per_poly +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Phisc1_exp <- explain(model_Phisc1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Phisc1", "P_poly", 
                                                                                     "pH_log", "veg_Ericaceae_per_poly")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Phisc1),
                          label = "Linear Regression Model")
plot(model_profile(Phisc1_exp))  


#Rhesp1
subset(dredges$Rhesp1, delta < 2) %>% View()
dredges$Rhesp1[1,]
Rhesp1_avg = model.avg(subset(dredges$Rhesp1, delta < 2))
Rhesp1_avg$coefficients %>% View()
summary(Rhesp1_avg)

model_Rhesp1 = lmer(data=orchamp_agg_all_wide_trans_scaled,
                        Rhesp1 ~ veg_Ericaceae_per_orderNorm + veg_nonAMplant_per_sqrt +
                          (1|codeplot) + (1|id_gradient), REML = F, na.action=na.fail)

Rhesp1_exp <- explain(model_Rhesp1,
                          data = as.data.frame(orchamp_agg_all_wide_trans_scaled[, c("id_gradient", "codeplot", "Rhesp1", 
                                                                                     "veg_Ericaceae_per_orderNorm", "veg_nonAMplant_per_sqrt")]), 
                          y = as.numeric(orchamp_agg_all_wide_trans_scaled$Rhesp1),
                          label = "Linear Regression Model")
plot(model_profile(Rhesp1_exp))  


### RMSE loss df ####
explainer <- list(Helotiales_all_exp, isolate_1_39_exp, isolate_1_40_exp, AF214576.1_exp, AF486119.1_exp, AJ292199.1_exp, AJ319077.1_exp, AY176753.1_exp, AY599235.1_exp, AY599236.1_exp, AY762619.1_exp, Glomeromycota_exp, Acema1_exp, Helotiales_other_exp, Bisps1_exp, DQ148411.1_exp, DQ485204.1_exp, EF446149.1_exp, EU938674.1_exp, FM180476.1_exp, FM207641.1_exp, FN178469.1_exp, FN178471.1_exp, Hymvar1_exp, JAHLEZ010000069.1_exp, JF304942.1_exp, KF836358.1_exp, KT827742.1_exp, KT835080.1_exp, KU382495.1_exp, KU886589.1_exp, MF943108.1_exp, Melbi2_exp, Melva1_exp, OL742649.1_exp, OM238143.1_exp, Oidma1_exp, PRJEB12348_exp, PRJNA407395_exp, Phisc1_exp, Rhesp1_exp)
names(explainer) <- c("Helotiales_all_exp", "isolate_1_39_exp", "isolate_1_40_exp", "AF214576.1_exp", "AF486119.1_exp", "AJ292199.1_exp", "AJ319077.1_exp", "AY176753.1_exp", "AY599235.1_exp", "AY599236.1_exp", "AY762619.1_exp", "Glomeromycota_exp", "Acema1_exp", "Helotiales_other_exp", "Bisps1_exp", "DQ148411.1_exp", "DQ485204.1_exp", "EF446149.1_exp", "EU938674.1_exp", "FM180476.1_exp", "FM207641.1_exp", "FN178469.1_exp", "FN178471.1_exp", "Hymvar1_exp", "JAHLEZ010000069.1_exp", "JF304942.1_exp", "KF836358.1_exp", "KT827742.1_exp", "KT835080.1_exp", "KU382495.1_exp", "KU886589.1_exp", "MF943108.1_exp", "Melbi2_exp", "Melva1_exp", "OL742649.1_exp", "OM238143.1_exp", "Oidma1_exp", "PRJEB12348_exp", "PRJNA407395_exp", "Phisc1_exp", "Rhesp1_exp")

ftr_imp <- purrr::map(explainer, ~ feature_importance(.x, loss_function = DALEX::loss_root_mean_square, B = 499))


# create dataframe from variable importances
RMSE_list = list()
for( i in names(ftr_imp)) {
  RMSE_list[[i]] <- ftr_imp[[i]] %>% as.data.frame %>% 
    group_by(variable) %>% 
    summarize(dropout_loss_mean = mean(dropout_loss)) %>% 
    ungroup() 
}
RMSE_resultats <- purrr::map_df(RMSE_list, ~as.data.frame(.x), .id="id") %>%
  mutate(id = gsub("_exp|fi_", "", id),
         variable_2 = gsub("_poly|_log|_yeojohnson|_boxcox|_orderNorm|_asinh|_sqrt", "", variable),
         transformation = case_when(
           grepl("_baseline_|_full_model_|id_gradient", variable) ~ NA,
           variable %in% unique(id) ~ NA,
           grepl("_log", variable) ~ "log",
           grepl("_poly", variable) ~ "poly",
           TRUE ~ "gaussian"
         )) %>%
  mutate(  id = case_when(
    id == "isolate_1_39" ~ "isolate_1-39",
    id == "isolate_1_40" ~ "isolate_1-40",
    TRUE ~ id
  )) %>%
  merge(x=., y=colors_pca, by.x="id", by.y="code_analysis_3", all.x = TRUE) %>%
  mutate(isolate_nom_Bruyant_2025 = case_when(
    id == "AJ292199.1" ~ "AJ292199.1 Meliniomyces vraolstadiae strain G1", 
    id == "Melbi2" ~ "Melbi2 Meliniomyces bicolor E", 
    id == "DQ485204.1" ~ "DQ485204.1 Cadophora finlandica isolate PRF15", 
    id == "AF486119.1" ~ "AF486119.1 Phialophora finlandia strain CBS 444.86 16S", 
    id == "Melva1" ~ "Melva1 Meliniomyces variabilis T", 
    id == "AY762619.1" ~ "AY762619.1 Meliniomyces variabilis strain UAMH 8861", 
    id == "KT835080.1" ~ "KT835080.1 Pezoloma ericae LivJKC02", 
    id == "AJ319077.1" ~ "AJ319077.1 Hymenoscyphus ericae isolate UAMH 5828", 
    id == "KT827742.1" ~ "KT827742.1 Pezoloma ericae EriJKE06", 
    id == "OM238143.1" ~ "OM238143.1 Hyaloscypha hepaticicola strain UAMH 7357", 
    id == "KF836358.1" ~ "KF836358 Cairneyella variabilis VPRI 42388 T", 
    id == "AY599235.1" ~ "AY599235.1 Grass root mycorrhizal sp. PPO-1", 
    id == "AY599236.1" ~ "AY599236.1 Grass root mycorrhizal sp. PPO-2", 
    id == "Phisc1" ~ "Phisc1 Phialocephala scopiformis", 
    id == "FN178471.1" ~ "FN178471.1 Ascomycota sp. Da MaP I9 Tapesia sp.", 
    id == "FN178469.1" ~ "FN178469.1 Ascomycota sp. Da PoT1 C7", 
    id == "Acema1" ~ "Acema1 Acephala macrosclerotiorum EW76-UTF0540", 
    id == "EU938674.1" ~ "EU938674.1 Phialocephala fortinii strain CBS 179.46", 
    id == "EF446149.1" ~ "EF446149.1 Phialocephala fortinii isolate PFO F", 
    id == "PRJEB12348" ~ "GCA_900073065.1 Phialocephala subalpina UAMH 11012", 
    id == "FM207641.1" ~ "FM207641.1 Ascomycota sp. Da Bac I4", 
    id == "KU886589.1" ~ "KU886589.1 Cadophora sp. Fe06", 
    id == "isolate_1-39" ~ "isolate_1-39", 
    id == "Rhesp1" ~ "Rhesp1 Rhexocercosporidium sp. MPI-PUGE-AT-0058", 
    id == "DQ148411.1" ~ "DQ148411.1 Leptodontidium sp. EF-M", 
    id == "JAHLEZ010000069.1" ~ "GCA_020726505.1 Leptodontidium sp. MPI-SDFR-AT-0119", 
    id == "AF214576.1" ~ "AF214576.1 Leptodontidium orchidicola strain UAMH 8152", 
    id == "OL742649.1" ~ "OL742649.1 Leptodophora orchidicola UAMH 5422", 
    id == "Oidma1" ~ "Oidma1 Oidiodendron maius Zn", 
    id == "KU382495.1" ~ "KU382495.1 Oidiodendron maius OM19", 
    id == "JF304942.1" ~ "JF304942.1 Oidiodendron maius isolate OMA B", 
    id == "AY176753.1" ~ "AY176753.1 Cryptosporiopsis rhizophila strain CBS 109839", 
    id == "Hymvar1" ~ "Hymvar1 Hymenoscyphus varicosporoides PMI453", 
    id == "MF943108.1" ~ "MF943108.1 Hyaloscyphaceae sp. isolate 1M05", 
    id == "PRJNA407395" ~ "GCA_002803225.1 Gamarada debralockiae T6G9 T", 
    id == "FM180476.1" ~ "FM180476.1 Pezeloma ericae AC2 IMI 392567", 
    id == "isolate_1-40" ~ "isolate_1-40", 
    id == "Bisps1" ~ "Bisps1 Calycina sp. as Bisporella", 
    TRUE ~ NA
  ),
  Name_JA = case_when(
    id == "isolate_1-39" ~ "Helotiales_sp_1-39", 
    id == "isolate_1-40" ~ "T._maxilliforme_1-40", 
    id == "AF214576.1" ~ "L._orchidicola_UAMH8152", 
    id == "AF486119.1" ~ "P._finlandica_CBS444.86", 
    id == "AJ292199.1" ~ "M._vraolstadiae_G1", 
    id == "AJ319077.1" ~ "H._ericae_UAMH5828", 
    id == "AY176753.1" ~ "C._rhizophila_CBS109839", 
    id == "AY599235.1" ~ "Grass_root_mycorrhizal_PPO-1", 
    id == "AY599236.1" ~ "Grass_root_mycorrhizal_PPO-2", 
    id == "AY762619.1" ~ "M._variabilis_UAMH8861", 
    id == "Glomeromycota" ~ "Glomeromycota", 
    id == "Acema1" ~ "A._macrosclerotiorum_Acema1", 
    id == "Helotiales_other" ~ "Helotiales_other", 
    id == "Helotiales_all" ~ "Helotiales_all", 
    id == "Bisps1" ~ "Calycina_sp._Bisps1", 
    id == "DQ148411.1" ~ "Leptodontidium_sp._EF-M", 
    id == "DQ485204.1" ~ "C._finlandica_PRF15", 
    id == "EF446149.1" ~ "P._fortinii_PFO_F", 
    id == "EU938674.1" ~ "P._fortinii_CBS179.46", 
    id == "FM180476.1" ~ "P._ericae_AC2_IMI_392567", 
    id == "FM207641.1" ~ "Ascomycota_sp._Da_Bac_I4", 
    id == "FN178469.1" ~ "Ascomycota_sp._Da_PoT1_C7", 
    id == "FN178471.1" ~ "Tapesia_sp._Da_MaP_I9", 
    id == "Hymvar1" ~ "H._varicosporoides_Hymvar1", 
    id == "JAHLEZ010000069.1" ~ "Leptodontidium_sp._MPI-SDFR-AT-0119", 
    id == "JF304942.1" ~ "O._maius_OMA_B", 
    id == "KF836358.1" ~ "C._variabilis_VPRI_42388_T", 
    id == "KT827742.1" ~ "P._ericae_EriJKE06", 
    id == "KT835080.1" ~ "P._ericae_LivJKC02", 
    id == "KU382495.1" ~ "O._maius_OM19", 
    id == "KU886589.1" ~ "Cadophora_sp._Fe06", 
    id == "MF943108.1" ~ "Hyaloscyphaceae_sp._1M05", 
    id == "Melbi2" ~ "M._bicolor_E_Melbi2", 
    id == "Melva1" ~ "M._variabilis_T_Melva1", 
    id == "OL742649.1" ~ "L._orchidicola_UAMH5422", 
    id == "OM238143.1" ~ "H._hepaticicola_UAMH7357", 
    id == "Oidma1" ~ "O._maius_Zn_Oidma1", 
    id == "PRJEB12348" ~ "P._subalpina_UAMH-11012", 
    id == "PRJNA407395" ~ "G._debralockiae_T6G9_T", 
    id == "Phisc1" ~ "P._scopiformis_Phisc1", 
    id == "Rhesp1" ~ "Rhexocercosporidium_sp_Rhesp1", 
    TRUE ~ NA
  ),
  order_phylo = case_when(
    id == "Helotiales_all" ~ 0.25, 
    id == "Helotiales_other" ~ 0.5, 
    id == "Glomeromycota" ~ 0.75, 
    id == "AF486119.1" ~ 1,
    id == "DQ485204.1" ~ 2,
    id == "AY762619.1" ~ 3,
    id == "AJ292199.1" ~ 4,
    id == "AJ319077.1" ~ 5,
    id == "KT835080.1" ~ 6,
    id == "Melbi2" ~ 7,
    id == "OM238143.1" ~ 8,
    id == "KF836358.1" ~ 9,
    id == "Melva1" ~ 10,
    id == "KT827742.1" ~ 11,
    id == "FN178469.1" ~ 12,
    id == "PRJEB12348" ~ 13,
    id == "EU938674.1" ~ 14,
    id == "AY599235.1" ~ 15,
    id == "JAHLEZ010000069.1" ~ 16,
    id == "EF446149.1" ~ 17,
    id == "FM207641.1" ~ 18,
    id == "FN178471.1" ~ 19,
    id == "AY599236.1" ~ 20,
    id == "OL742649.1" ~ 21,
    id == "Acema1" ~ 22,
    id == "KU886589.1" ~ 23,
    id == "isolate_1-39" ~ 24,
    id == "Rhesp1" ~ 25,
    id == "DQ148411.1" ~ 26,
    id == "AF214576.1" ~ 27,
    id == "Phisc1" ~ 28,
    id == "MF943108.1" ~ 29,
    id == "Hymvar1" ~ 30,
    id == "Bisps1" ~ 31,
    id == "AY176753.1" ~ 32,
    id == "Oidma1" ~ 33,
    id == "JF304942.1" ~ 34,
    id == "FM180476.1" ~ 35,
    id == "PRJNA407395" ~ 36,
    id == "KU382495.1" ~ 37,
    id == "isolate_1-40" ~ 38,
    TRUE ~ NA
  ),
  clade = case_when(
    grepl("AJ292199.1|Melbi2|DQ485204.1|AF486119.1|Melva1|AY762619.1|KT835080.1|AJ319077.1|KT827742.1|OM238143.1|KF836358.1", id)  ~ "Hyaloscyphoid clade", 
    grepl("AY599235.1|AY599236.1|Phisc1|FN178471.1|FN178469.1|Acema1|EU938674.1|EF446149.1|PRJEB12348|FM207641.1|KU886589.1|isolate_1-39|Rhesp1|DQ148411.1|JAHLEZ010000069.1|AF214576.1|OL742649.1", id)  ~ "Mollisioid clade", 
    grepl("Oidma1|KU382495.1|JF304942.1|AY176753.1|Hymvar1|MF943108.1|PRJNA407395|FM180476.1|Bisps1", id)  ~ "Pezizelloid clade", 
    grepl("isolate_1-40", id)  ~ "Vandijckellaceae – Calloriaceae clade", 
    grepl("Helotiales_all|Helotiales_other|Glomeromycota", id)  ~ "Misc", 
    TRUE ~ NA
  ),
  lifestyle = case_when(
    id == "Helotiales_all" ~ "Mixed",
    TRUE ~ lifestyle
  ),
  color_lifestyle = case_when(
    id == "Helotiales_all" ~ "black",
    TRUE ~ color_lifestyle
  ),
  lifestyle_JA = case_when(
    id == "Helotiales_all" ~ "Mixed",
    TRUE ~ lifestyle_JA
  ),
  color_lifestyle_JA = case_when(
    id == "Helotiales_all" ~ "black",
    TRUE ~ color_lifestyle_JA
  ))

# add the variable importances to the dredge results
setDT(RMSE_resultats )[dredges_summary, c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)", "p_value_asterisks") := 
                                .(Estimate, `Std. Error`, df, `t value`, `Pr(>|t|)`, p_value_asterisks), on=c("id", "variable_2", "transformation")]


RMSE_resultats = RMSE_resultats %>%
  mutate(poly_type = case_when(
    id == "Acema1" & variable == "pH_poly" ~ "poly_exp", 
    id == "AF214576.1" & variable == "Net_solar_radiation_poly" ~ "poly_parabolic", 
    id == "AF214576.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "AF486119.1" & variable == "P_poly" ~ "poly_parabolic", 
    id == "AJ292199.1" & variable == "CN_poly" ~ "poly_parabolic", 
    id == "AJ292199.1" & variable == "P_poly" ~ "poly_exp", 
    id == "AJ292199.1" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "AJ319077.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "AY176753.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "AY176753.1" & variable == "Mean_air_temp_poly" ~ "poly_exp", 
    id == "AY599235.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "AY599236.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "AY599236.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "AY762619.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "Bisps1" & variable == "CWS_poly" ~ "poly_parabolic", 
    id == "Bisps1" & variable == "Net_solar_radiation_poly" ~ "poly_parabolic", 
    id == "Bisps1" & variable == "P_poly" ~ "poly_parabolic", 
    id == "DQ148411.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "DQ485204.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "DQ485204.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_exp", 
    id == "EF446149.1" & variable == "CN_poly" ~ "poly_parabolic", 
    id == "EF446149.1" & variable == "CWS_poly" ~ "poly_parabolic", 
    id == "EF446149.1" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "EU938674.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "EU938674.1" & variable == "veg_species_richness_poly" ~ "poly_exp", 
    id == "FM180476.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "FM207641.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "FM207641.1" & variable == "Mean_air_temp_poly" ~ "poly_exp", 
    id == "FN178471.1" & variable == "pH_poly" ~ "poly_parabolic", 
    id == "FN178471.1" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "Glomeromycota" & variable == "CN_poly" ~ "poly_exp", 
    id == "Helotiales_other" & variable == "pH_poly" ~ "poly_exp", 
    id == "Hymvar1" & variable == "pH_poly" ~ "poly_parabolic", 
    id == "isolate_1-39" & variable == "pH_poly" ~ "poly_parabolic", 
    id == "isolate_1-39" & variable == "veg_nonAMplant_per_poly" ~ "poly_exp", 
    id == "isolate_1-40" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "JF304942.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "JF304942.1" & variable == "CWS_poly" ~ "poly_exp", 
    id == "KT835080.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "KT835080.1" & variable == "CWS_poly" ~ "poly_exp", 
    id == "KT835080.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_exp", 
    id == "KU382495.1" & variable == "CN_poly" ~ "poly_parabolic", 
    id == "KU886589.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "Melbi2" & variable == "P_poly" ~ "poly_exp", 
    id == "Melbi2" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "Melbi2" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "Melva1" & variable == "P_poly" ~ "poly_parabolic", 
    id == "Melva1" & variable == "veg_Ericaceae_per_poly" ~ "poly_parabolic", 
    id == "MF943108.1" & variable == "CN_poly" ~ "poly_exp", 
    id == "MF943108.1" & variable == "Mean_air_temp_poly" ~ "poly_exp", 
    id == "MF943108.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "MF943108.1" & variable == "veg_species_richness_poly" ~ "poly_parabolic", 
    id == "Oidma1" & variable == "pH_poly" ~ "poly_exp", 
    id == "OL742649.1" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "OM238143.1" & variable == "P_poly" ~ "poly_parabolic", 
    id == "OM238143.1" & variable == "pH_poly" ~ "poly_exp", 
    id == "Phisc1" & variable == "P_poly" ~ "poly_exp", 
    id == "Phisc1" & variable == "veg_Ericaceae_per_poly" ~ "poly_exp", 
    id == "PRJEB12348" & variable == "CN_poly" ~ "poly_exp", 
    id == "PRJEB12348" & variable == "veg_nonAMplant_per_poly" ~ "poly_parabolic", 
    id == "PRJEB12348" & variable == "veg_species_richness_poly" ~ "poly_parabolic", 
    id == "PRJNA407395" & variable == "Net_solar_radiation_poly" ~ "poly_exp", 
    TRUE ~ NA
  ),
  transformation_2 = case_when(
    transformation == "poly" ~ poly_type,
    TRUE ~ transformation
  ),
  transformation_image = case_when(
    (transformation_2 == "gaussian" & Estimate <0) ~ "linear_negative.png",
    (transformation_2 == "gaussian" & Estimate >0) ~ "linear_positive.png",
    (transformation_2 == "log" & Estimate <0) ~ "log_negative_2.png",
    (transformation_2 == "log" & Estimate >0) ~ "log_positive_2.png",
    (transformation_2 == "poly_exp" & Estimate <0) ~ "poly_exp_negative.png",
    (transformation_2 == "poly_exp" & Estimate >0) ~ "poly_exp_positive.png",
    (transformation_2 == "poly_parabolic" & Estimate <0) ~ "poly_concave.png",
    (transformation_2 == "poly_parabolic" & Estimate >0) ~ "poly_convex.png",
    TRUE ~ NA
  ))


# to calculate the variable importances as a percentage of the total, aka relative loss mean
RMSE_resultats_total = RMSE_resultats %>%
  filter(!(variable %in% c("_baseline_", "_full_model_", "id_gradient","isolate_1-39","isolate_1-40", "fung_relative_abundance")|
             variable %in% unique(.$id))) %>% 
  group_by(id) %>%
  summarize(total_dropout_loss_mean = sum(dropout_loss_mean)) %>%
  ungroup()

require(data.table) # v1.9.5+
setDT(RMSE_resultats)[RMSE_resultats_total, c("total_dropout_loss_mean") := 
                        .(total_dropout_loss_mean), on=c("id")]

RMSE_resultats = RMSE_resultats %>%
  mutate(relative_dropout_loss_mean = case_when(
    !(variable %in% c("_baseline_", "_full_model_", "id_gradient","isolate_1-39","isolate_1-40", "fung_relative_abundance")|variable %in% unique(.$id)) ~ (dropout_loss_mean/total_dropout_loss_mean*100),
    TRUE ~ NA))


# Heatmap 
RMSE_resultats$id_col <- paste0("<span style=\"color: ", RMSE_resultats$color_lifestyle_JA, "\">", RMSE_resultats$Name_JA, "</span>")

# sometimes the images cache needs to be cleaned if the images have been changed
#gridpattern::reset_image_cache()

memoise::forget(gridpattern:::img_read_memoised)

setwd("C:/Users/lem/Desktop/Lauren_SymbiLoss/ORCHAMP/Schreiben/figures/courbes")
setwd("/Users/laurengillespie/Desktop/LEM_Mac/courbes")


RMSE_heatmap <- RMSE_resultats %>%
  filter(!(variable %in% c("_baseline_", "_full_model_", "id_gradient", "codeplot", "isolate_1-39","isolate_1-40", "fung_relative_abundance")|
             variable %in% unique(.$id)|
            id == "Glomeromycota"
           )) %>%
  mutate(variable_2 = case_when(
    variable_2 == "veg_Ericaceae_per" ~ "Percent_Ericaceae",
    variable_2 == "veg_nonAMplant_per" ~ "Percent_non_AM",
    variable_2 == "veg_species_richness" ~ "Plant_sp_richness",
    TRUE ~ variable_2),
  dropout_loss_mean = (dropout_loss_mean*100)) %>% 
  mutate(
           variable_2 = factor(variable_2, levels=c("Mean_air_temp", "Net_solar_radiation", "CWS", "CN", "P", "pH", "Percent_Ericaceae", "Percent_non_AM","Plant_sp_richness")), #"id_gradient"
         clade = factor(clade, levels = c("Misc", "Hyaloscyphoid clade", "Mollisioid clade", "Pezizelloid clade", "Vandijckellaceae – Calloriaceae clade")),
          id_col =  fct_reorder(id_col, order_phylo, .desc = TRUE) ) %>% 
  mutate(relative_dropout_loss_mean = case_when(
    Estimate <0 ~ relative_dropout_loss_mean * -1,
    TRUE ~ relative_dropout_loss_mean
  )) %>%
  ggplot(., aes(variable_2, id_col)) +
  geom_tile(aes(fill = relative_dropout_loss_mean)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-100, 100)) + 

  labs(fill="Variable importance (%)") +
  geom_text(aes(label = p_value_asterisks), size=4, hjust = 0.5,
            nudge_x = 0.25) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(labels= c("Mean air temp", "Net solar radiation", "CWS", "CN", "P", "pH", "% Ericaceae", "% non-AM plants","Plant sp. richness"),
                   expand = c(0, 0))  +
  geom_hline(yintercept = c(1.5, 10.5, 27.5, 38.5), col = 'black', lwd=0.6)+
  theme_bw() +
  theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 40, "points"),
        panel.grid.major = element_line(color = "gray85",linetype="longdash",linewidth=0.1),
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 1.11),
        axis.title.y = element_blank(),
        axis.text.y = ggtext::element_markdown(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.spacing.y = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_rect(fill=NA,color=NA) )  +
  geom_image(aes(image = transformation_image), size = 0.025, nudge_x = -0.2, nudge_y = 0.03)  +
  annotate("text", x = -4.1, y = 39.5, label = "Misc.", angle = 90) +  
  annotate("text", x = -4.1, y = 33, label = "Hyaloscyphoid clade", angle = 90) +  
  annotate("text", x = -4.1, y = 19, label = "Mollisioid clade", angle = 90) +  
  annotate("text", x = -4.1, y = 6.5, label = "Pezizelloid clade", angle = 90) +
  annotate("text", x = -3.4, y = 1, label = "Vandijckellaceae –\nCalloriaceae clade", angle = 0) + 
  coord_cartesian(xlim = c(0.5, 9.5),clip="off")  
  
RMSE_heatmap_image = RMSE_heatmap +
  geom_image(aes(image = transformation_image), size = 0.03, nudge_x = -0.2) 

# R21
color_lifestyle_JA   <-  c(
  "Endophyte" = "#E41A1C",
  "Ericoid_Mycorrhizal" = "#377EB8",
  "Other_or_Unknown" = "#4DAF4A",
  "Mixed" = "black"
)

R21_label = dredges_resultats %>%
  filter(!id == "Glomeromycota") %>%
  mutate(id = factor(id, levels=c("Helotiales_all", "Helotiales_other", "PRJNA407395", "OM238143.1", "Oidma1", "KU382495.1", "KF836358.1", "JF304942.1", "FM180476.1", "Rhesp1", "Phisc1", "JAHLEZ010000069.1", "Hymvar1", "Bisps1", "Melbi2", "DQ485204.1", "Acema1", "Melva1", "KT835080.1", "KT827742.1", "Glomeromycota", "PRJEB12348", "OL742649.1", "MF943108.1", "KU886589.1", "isolate_1-40", "isolate_1-39", "FN178471.1", "FN178469.1", "FM207641.1", "EU938674.1", "EF446149.1", "DQ148411.1", "AY762619.1", "AY599236.1", "AY599235.1", "AY176753.1", "AJ319077.1", "AJ292199.1", "AF486119.1", "AF214576.1"))) %>%
  select(R21, id, lifestyle_JA, color_lifestyle_JA) %>%
  mutate(R21 = as.factor(round(R21, digits = 2)),
         R2 = "R2m") %>%
  mutate(R21 = case_when(
    R21 == "0.4" ~ "0.40",
    R21 == "0.1" ~ "0.10",
    R21 == "0.6" ~ "0.60",
    R21 == "0.2" ~ "0.20",
    R21 == "0.5" ~ "0.50",
    R21 == "0" ~ "0.001",
    TRUE ~ R21
  )) 

R21_plot = dredges_resultats %>%
  select(id, lifestyle_JA, color_lifestyle_JA, R21, order_phylo) %>%
  filter(!id == "Glomeromycota") %>%
  mutate(R2 = "R2m",
         id =  fct_reorder(id, order_phylo, .desc = TRUE)) %>%
  ggplot(., aes(x= R2, y=id, label= R21)) +
  geom_point(color = "white") +
  geom_text(data = R21_label, aes(x = R2, y = id, label = R21, color = lifestyle_JA), 
            size = 3,
            inherit.aes = FALSE) +
  scale_color_manual('Lifestyle', values = color_lifestyle_JA,
                     limits = c("Mixed", "Other_or_Unknown", "Ericoid_Mycorrhizal",  "Endophyte"),
                     labels=c("Mixed", "Other or Unknown",  "Ericoid Mycorrhizal",  "Endophyte")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  +
  guides(color = guide_legend(override.aes = aes(label = "●")))


#R22
R22_label = dredges_resultats %>%
  filter(!id == "Glomeromycota") %>%
  mutate(id = factor(id, levels=c("Helotiales_all" , "Helotiales_other", "PRJNA407395", "OM238143.1", "Oidma1", "KU382495.1", "KF836358.1", "JF304942.1", "FM180476.1", "Rhesp1", "Phisc1", "JAHLEZ010000069.1", "Hymvar1", "Bisps1", "Melbi2", "DQ485204.1", "Acema1", "Melva1", "KT835080.1", "KT827742.1", "Glomeromycota", "PRJEB12348", "OL742649.1", "MF943108.1", "KU886589.1", "isolate_1-40", "isolate_1-39", "FN178471.1", "FN178469.1", "FM207641.1", "EU938674.1", "EF446149.1", "DQ148411.1", "AY762619.1", "AY599236.1", "AY599235.1", "AY176753.1", "AJ319077.1", "AJ292199.1", "AF486119.1", "AF214576.1"))) %>%
  select(R22, id, lifestyle_JA, color_lifestyle_JA) %>%
  mutate(R22 = as.factor(round(R22, digits = 2)),
         R2 = "R2c") %>%
  mutate(R22 = case_when(
    R22 == "0.4" ~ "0.40",
    R22 == "0.1" ~ "0.10",
    R22 == "0.6" ~ "0.60",
    R22 == "0.2" ~ "0.20",
    R22 == "0.7" ~ "0.70",
    R22 == "0" ~ "0.001",
    TRUE ~ R22
  )) 

R22_plot = dredges_resultats %>%
  select(id, lifestyle_JA, color_lifestyle, R22, order_phylo) %>%
  filter(!id == "Glomeromycota") %>%
  mutate(R2 = "R2c",
         id =  fct_reorder(id, order_phylo, .desc = TRUE)) %>%
  ggplot(., aes(x= R2, y=id, label= R22)) +
  geom_point(color = "white") +
  geom_text(data = R22_label, aes(x = R2, y = id, label = R22, color = lifestyle_JA), 
            size = 3,
            inherit.aes = FALSE) +
  scale_color_manual('Lifestyle', values = color_lifestyle_JA,
                     limits = c("Mixed", "Other_or_Unknown", "Ericoid_Mycorrhizal",  "Endophyte"),
                     labels=c("Mixed", "Other or Unknown",  "Ericoid Mycorrhizal",  "Endophyte")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  +
  guides(color = guide_legend(override.aes = aes(label = "●")))

 

# Habitat heatmap

habitat_heatmap_df <- orchamp_metabar_blast_agg_years_2 %>%
  filter(code_analysis == "fung_relative_abundance") %>%
  filter(!code_analysis_4 == "Glomeromycota") %>%
  filter(!is.na(CorineLandCover_simplifie_2)) %>%
  mutate(CorineLandCover_simplifie_2 = case_when(
           is.na(CorineLandCover_simplifie_2) & habitat_type == "Milieu forestier" ~ "Forest",
           is.na(CorineLandCover_simplifie_2) & habitat_type == "Milieu ouvert" ~ "Grassland",
           TRUE ~ CorineLandCover_simplifie_2
         ),
         CorineLandCover_simplifie =  case_when(
           CorineLandCover_simplifie == "Forets de coniferes" ~ "Evergreen forest",
           CorineLandCover_simplifie == "Forets de feuillus" ~ "Deciduous forest",
           CorineLandCover_simplifie == "Forets mixtes" ~ "Mixed forest",
           CorineLandCover_simplifie == "Pelouses et patures" ~ "Grassland and pastures",
           CorineLandCover_simplifie == "Pelouses ecorches" ~ "Stripped grassland",
           CorineLandCover_simplifie == "Landes" ~ "Heathland",
           TRUE ~ CorineLandCover_simplifie)) %>% 
  group_by(CorineLandCover_simplifie, code_analysis_4, code_analysis_3, order_phylo) %>%
  summarize(value_mean = sum(value_mean)) %>%
  ungroup() %>%
  mutate(value_mean = log10((value_mean + (0.0000113/2))),
         order_phylo = case_when(
           code_analysis_3 == "Helotiales_all" ~ 0.25, 
           code_analysis_3 == "Helotiales_other" ~ 0.5, 
           code_analysis_3 == "Glomeromycota" ~ 0.75, 
           code_analysis_3 == "AJ292199.1" ~ 1, 
           code_analysis_3 == "Melbi2" ~ 2, 
           code_analysis_3 == "DQ485204.1" ~ 3, 
           code_analysis_3 == "AF486119.1" ~ 5, 
           code_analysis_3 == "Melva1" ~ 7, 
           code_analysis_3 == "AY762619.1" ~ 8, 
           code_analysis_3 == "KT835080.1" ~ 9, 
           code_analysis_3 == "AJ319077.1" ~ 10, 
           code_analysis_3 == "KT827742.1" ~ 11, 
           code_analysis_3 == "OM238143.1" ~ 13, 
           code_analysis_3 == "KF836358.1" ~ 14, 
           code_analysis_3 == "AY599235.1" ~ 16, 
           code_analysis_3 == "AY599236.1" ~ 17, 
           code_analysis_3 == "Phisc1" ~ 19, 
           code_analysis_3 == "FN178471.1" ~ 20, 
           code_analysis_3 == "FN178469.1" ~ 21, 
           code_analysis_3 == "Acema1" ~ 22, 
           code_analysis_3 == "EU938674.1" ~ 25, 
           code_analysis_3 == "EF446149.1" ~ 27, 
           code_analysis_3 == "PRJEB12348" ~ 28, 
           code_analysis_3 == "FM207641.1" ~ 33, 
           code_analysis_3 == "KU886589.1" ~ 34, 
           code_analysis_3 == "isolate_1-39" ~ 35, 
           code_analysis_3 == "Rhesp1" ~ 40, 
           code_analysis_3 == "DQ148411.1" ~ 41, 
           code_analysis_3 == "JAHLEZ010000069.1" ~ 42, 
           code_analysis_3 == "AF214576.1" ~ 43, 
           code_analysis_3 == "OL742649.1" ~ 44, 
           code_analysis_3 == "Oidma1" ~ 48, 
           code_analysis_3 == "KU382495.1" ~ 51, 
           code_analysis_3 == "JF304942.1" ~ 52, 
           code_analysis_3 == "AY176753.1" ~ 53, 
           code_analysis_3 == "Hymvar1" ~ 59, 
           code_analysis_3 == "MF943108.1" ~ 60, 
           code_analysis_3 == "PRJNA407395" ~ 62, 
           code_analysis_3 == "FM180476.1" ~ 64, 
           code_analysis_3 == "isolate_1-40" ~ 65, 
           code_analysis_3 == "Bisps1" ~ 69, 
           TRUE ~ NA
         )) %>% 
  mutate(code_analysis_4 =  forcats::fct_reorder(code_analysis_4, order_phylo, .desc = TRUE) ,
         CorineLandCover_simplifie = factor(CorineLandCover_simplifie, levels=c("Deciduous forest", "Mixed forest", "Evergreen forest", "Heathland", "Grassland and pastures", "Stripped grassland"))
         ) 

habitat_heatmap = habitat_heatmap_df %>% 
  ggplot(., aes(CorineLandCover_simplifie, code_analysis_4, fill = CorineLandCover_simplifie)) +                           
  geom_tile(aes(alpha = value_mean)) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      breaks=c("Stripped grassland", "Grassland and pastures", "Heathland", "Evergreen forest", "Mixed forest", "Deciduous forest")) +
  scale_alpha_continuous() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))  +
  labs(fill="Rel. abun.\n(log10(x+RAmin/2))") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray85",linetype="longdash",linewidth=0.1),
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 1.11),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,  hjust=1),
    panel.spacing.y = unit(0, "lines")
  )
   

habitat_heatmap_legend = habitat_heatmap_df %>% 
  ggplot(., aes(CorineLandCover_simplifie, code_analysis_4)) +                           
  geom_tile(aes(fill = value_mean)) +
  scale_fill_distiller(palette = "Greys", 
                       direction = 1) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))  +
  labs(fill="Relative abundance\n(RA) per habitat\n(log10(x+RAmin/2))" ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray85",linetype="longdash",linewidth=0.1),
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 1.11),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,  
                               hjust=1),
    panel.spacing.y = unit(0, "lines")
  )


RMSE_layout <- "
  1234567
"

RMSE_heatmap + plot_spacer() + free(R21_plot, type = "label")  + plot_spacer() + free(R22_plot, type = "label") + plot_spacer() + habitat_heatmap_legend +
  plot_layout(design = RMSE_layout) +
  plot_layout(widths = c(5, -0.35, 0.5, -0.47, 0.5, -0.30, 1.0), guides = "collect")


ggsave("Figures/Fig_heatmap_relative_RMSE_all.pdf", width = 30, height = 25, units = "cm", bg = "white")

