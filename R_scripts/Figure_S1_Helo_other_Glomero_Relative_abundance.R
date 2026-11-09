## Figure S1. Helotiales ‘other’ and Gomeromycota relative abundances 

#Relative abundance of (a) Helotiales ‘other’ strains and (b) Glomeromycota strains relative to the entire fungal community per plot 
#by increasing altitude. The plot code number indicates the plot elevation (m).


# Required packages -------------------------------------
library("dplyr")
library("ggplot2")
library("forcats")
library("patchwork")

# Load data ---------------------------------------------

orchamp_metabar_blast_agg_years <- readRDS("data/orchamp_metabar_blast_agg_years_20260327.rds")


# Figure S1a --------------------------------------------

### (a) Helotiales ‘other’ strains ####

Helo_other_stackedbarplot = orchamp_metabar_blast_agg_years_2 %>%
  filter(!is.na(isolate) & code_analysis == "fung_relative_abundance") %>%
  mutate(code_analysis_2 = case_when(
    grepl("_sp$", code_analysis_2) ~ gsub("_sp$", "_unclassified", code_analysis_2),
    grepl("unclassified_", code_analysis_2) ~ paste(gsub("unclassified_", "", code_analysis_2), "_unclassified", sep = ""),
    TRUE ~ code_analysis_2)) %>% 
  select(codeplot,  altitude, code_analysis_2, code_analysis_3, value_mean) %>%  
  filter(grepl("Helotiales_other", code_analysis_3)) %>%
  group_by(codeplot, altitude, code_analysis_2) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(codeplot = fct_reorder(codeplot, altitude)) %>%
  ggplot(., aes(x = codeplot, y = (value_mean*100), fill= code_analysis_2)) +
  geom_bar(position="stack", stat="identity", color = "black") + 
  scale_y_continuous(expand = c(0, 0), lim = c(0, 17),
                     breaks = seq(0.0, 17, by = 5)) +
  guides(fill=guide_legend(title="Helotiales Other lineages", ncol = 2)) +
  labs(title = "(a)") +
  ylab("Helotiales Other relative abundance within\nthe entire soil fungal community (%)") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 13,  face = "bold"),  
        axis.text.x = element_text(size = 13),
        plot.title = element_text(hjust = 0),
        axis.title.y = element_blank()) +
  coord_flip()

# Figure S1b --------------------------------------------

### (b) Glomeromycota strains ####

palette_glomero <- c(
  "Acaulospora_alpina" =  "#D81E00",
  "Acaulospora_brasiliensis" =  "#FDAF91",
  "Acaulospora_cavernata" =  "#54278F",
  "Acaulospora_unclassified" =  "#9996C3",
  "Acaulospora_tortuosa" =  "#08519C",
  "Acaulosporaceae_unclassified" =  "#6BAED6",
  "Ambispora_unclassified" =  "#09306B",
  "Ambisporaceae_unclassified" =  "#797900",
  "Archaeospora_unclassified" =  "#008A8B",
  "Archaeosporaceae_unclassified" =  "#74C475",
  "Archaeosporales_unclassified" =  "#238B45",
  "Claroideoglomeraceae_unclassified" =  "#83F4F7",
  "Claroideoglomus_etunicatum" =  "#CAFF70",
  "Claroideoglomus_unclassified" =  "#02441B",
  "Complexispora_mediterranea" =  "#76DAB3",
  "Diversispora_epigaea" =  "#8AAF9A",
  "Diversispora_sabulosa" =  "#9C669C",
  "Diversispora_unclassified" =  "#FFD966",
  "Dominikia_difficilevidera" =  "#FF56D6",
  "Dominikia_unclassified" =  "#AFABAB",
  "Entrophospora_glacialis" =  "grey85",
  "Epigeocarpum_japonicum" =  "#ffcb3d",
  "Funneliformis_pilosus" =  "#5e7fff",
  "Glomeraceae_unclassified" =  "#7da800",
  "Glomeromycetes_unclassified" =  "#7b009f",
  "Glomeromycota_unclassified" =  "#a99000",
  "Glomus_indicum" =  "#0277bc",
  "Glomus_invermaium" =  "#8c6900",
  "Glomus_unclassified" =  "#ccabff",
  "Glomus_tetrastratosum" =  "#ff83cf",
  "Microdominikia_litorea" =  "#00511c",
  "Microdominikia_unclassified" =  "#ff728c",
  "Paraglomerales_unclassified" =  "#ff6700",
  "Paraglomus_laccatum" =  "#e5ffa1",
  "Paraglomus_unclassified" =  "#860042",
  "Rhizoglomus_unclassified" =  "#bed8ff",
  "Rhizophagus_unclassified" =  "#4DAF4A",
  "Scutellospora_unclassified" =  "#A65628",
  "Septoglomus_mediterraneum" = "#EE3B3B",
  "Complexispora_unclassified" = "#66CD00",
  "Diversisporaceae_unclassified" = "#EEE9BF",
  "Entrophosporaceae_unclassified" = "#CD6839",
  "Funneliformis_unclassified" = "#EEEE00",
  "Gigasporales_unclassified" = "#7A67EE",
  "Glomerales_unclassified" = "#D02090",
  "Pacispora_unclassified" = "#D2B48C",
  "Septoglomus_unclassified" = "#BBFFFF"
  )

Glomero_stackedbarplot = orchamp_metabar_blast_agg_years_2 %>%
  filter(!is.na(isolate)) %>%
  mutate(code_analysis_2 = case_when(
    grepl("_sp$", code_analysis_2) ~ gsub("_sp$", "_unclassified", code_analysis_2),
    grepl("unclassified_", code_analysis_2) ~ paste(gsub("unclassified_", "", code_analysis_2), "_unclassified", sep = ""),
    TRUE ~ code_analysis_2)) %>% 
  select(codeplot,  altitude, code_analysis_2, code_analysis_3, value_mean) %>%  
  filter(grepl("Glomero", code_analysis_3)) %>%
  group_by(codeplot, altitude, code_analysis_2) %>% 
  summarize( value_mean = sum(value_mean, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(codeplot = fct_reorder(codeplot, altitude)) %>%
  ggplot(., aes(x = codeplot, y = (value_mean*100), fill= code_analysis_2)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 65),
                     breaks = seq(0.0, 65, by = 10)) +
  scale_fill_manual(values=palette_glomero) + 
  guides(fill=guide_legend(title="Glomeromycota lineages", ncol = 1)) +
  labs(title = "(b)") +
  ylab("Glomeromycota relative abundance within\nthe entire soil fungal community (%)") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 13,  face = "bold"),  
        axis.text.x = element_text(size = 13),
        plot.title = element_text(hjust = 0),
        axis.title.y = element_blank()) +
  coord_flip()


### Increasing altitude triangle. ####
increasing_altitude_Helo_other_Glomero = data.frame(x =c(1,2,3), y = c(1:174)) %>%
  ggplot(., aes(x = x, y = y)) +
  geom_point(color = "white") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), clip = "off") +
  annotate("segment", x = 2, y = 3, xend = 2, yend = 11, 
           linewidth = 1.25, linejoin = "mitre",
           arrow = arrow(type = "closed", length = unit(0.2, "npc"))) +
  annotate("text", x=2, 
           y=0, label = "Increasing altitude", color = "black", 
           angle = 90, size = 5, fontface = "bold", hjust = 1 ) +
  annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(Inf, -Inf, Inf), fill = "darkgray" )+
  theme_void()+
  theme(plot.margin = margin(0.5,0.15,0.5,0.15, "cm") )


### Combine figures ###
Fung_Glomo_Helo_other <- "
  123
"
increasing_altitude_Helo_other_Glomero + Helo_other_stackedbarplot + Glomero_stackedbarplot +
  plot_layout(design = Fung_Glomo_Helo_other) +
  plot_layout(widths = c(0.6, 5, 5))+
  plot_annotation(theme = theme(plot.title = element_text(size = 18)))


ggsave("FigS1.pdf", width = 58, height = 60, units = "cm")


