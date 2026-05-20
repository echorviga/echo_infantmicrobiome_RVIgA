
#######################################
##### Figure 2: Longitudinal development of the infant gut microbiome.
########################################


######## Set up ############
library(tidyverse); library(dplyr); library(data.table);
library(rstatix); library(ggpubr); library(vegan) ; library(ggplot2); 
library(glue);library(ape);library(grid)

#### ASV and Metadata for Figures 1A and 1B #####

in.full.meta <- read.csv("../Microbiome_RVV/paper/working_code/code_upload/figure2/metadata_fig2_AB.csv")%>%
  mutate(sample_id=sid) %>%
  column_to_rownames(., var="sample_id")

in.full.ra <- read.csv("../Microbiome_RVV/paper/working_code/code_upload/figure2/asv_fig2_AB.csv") %>%
  column_to_rownames(., var="sid") %>%
  as.matrix()


############################

########## Panel A: Alpha diversity ############

shannondiv_all<- in.full.ra%>%
  diversity(., index = "shannon") 

shandiv_df_all <- shannondiv_all %>% 
  enframe() %>%
  rename(sid = name,
         shannon = value) %>%
  merge(x=.,y=in.full.meta, by= "sid")%>%
  mutate(month = factor(month, levels=c(0,1,6,12), labels=c("Birth", "M1", "M6", "M12")))


### Linear mixed effects mode, using continous breastfeeding duration, as feedingg method is not recorded for birth timepoint

model <- lmerTest::lmer(shannon ~ month + (1 | participant_id) + bf_duration_thru12_w , data = shandiv_df_all)

em_results <- emmeans::emmeans(model, pairwise ~ month, adjust = "BH")  

summary(em_results)

### Figure 1A

format_p_value_interpunct <- function(p, digits = 6) {
  formatted_p <- ifelse(p < 0.0001,
                        "<0·0001",
                        gsub("\\.", "\u00B7", sprintf(paste0("%.", digits, "f"), p))) 
  return(formatted_p)
}

comparisons<- as.data.frame(em_results$contrasts)%>%
  filter(contrast %>% row_number() %in% c(1, 5, 6))%>%
  mutate( x_start = case_when(contrast== "Birth - M1" ~ 1,
                              contrast== "M1 - M6" ~ 2,
                              contrast== "M6 - M12" ~ 3),
          x_end = case_when(contrast== "Birth - M1" ~ 2,
                            contrast== "M1 - M6" ~ 3,
                            contrast== "M6 - M12" ~ 4),
          y_position = case_when(contrast== "Birth - M1" ~ 4,
                                 contrast== "M1 - M6" ~ 4.5,
                                 contrast== "M6 - M12" ~ 5),
          label = glue("p {format_p_value_interpunct(p.value)}"))

mc<-list(c("Birth", "M1"), c("M1", "M6"), c("M6", "M12"))

visit.cols = c("Birth"= "#999999",
               "M1"= "#FF7256",
               "M6"= "#56B4E9",
               "M12"= "darkblue")

shandiv_df_all %>%
  ggplot(.,aes(x=month, y=shannon)) +
  geom_boxplot(outlier.shape=NA, aes(color=month))+
  geom_jitter(position=position_jitter(), aes(color=month), size=2) +
  scale_color_manual(values = visit.cols)+ 
  theme_bw()+
  labs(y= " Shannon Diversity")+
  geom_segment(data = comparisons,
               aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
               inherit.aes = FALSE) +
  geom_text(data = comparisons,
            aes(x = (x_start + x_end)/2, y = y_position + 0.2, label = label),
            size = 8,
            inherit.aes = FALSE)+
  theme(text = element_text(size=20), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size=20))

###############################################
########## Panel B: Beta diversity ############
###############################################

dis_in_all<- vegdist(in.full.ra, method = "bray")

####################
set.seed(2026)

participant_block <- in.full.meta[rownames(as.matrix(dis_in_all)), "participant_id"]


in.full.meta <- in.full.meta[rownames(as.matrix(dis_in_all)),]


adonis2(dis_in_all ~ month, 
        data = in.full.meta, 
        permutations = 20000,
        strata = in.full.meta[rownames(as.matrix(dis_in_all)), "participant_id"])

bd <- betadisper(dis_in_all, 
                 group = factor(in.full.meta[rownames(as.matrix(dis_in_all)), "month"]))

permutest(bd, permutations = 50000)

permutest(bd, permutations = 50000, pairwise = TRUE)

TukeyHSD(bd)

####### Birth- M1 ##########
index_01 <- which(colnames(as.matrix(dis_in_all)) %in% in.full.meta$sid[in.full.meta$month%in%c("0","1")])

dist_in01<- as.matrix(dis_in_all)[index_01,index_01] %>%
  as.dist()

betadisper(dist_in01, group= factor(as.matrix(in.full.meta[in.full.meta$month%in%c("0","1"),"month"]))) %>% anova()

adonis2(dist_in01~month, data=in.full.meta[in.full.meta$month%in%c("0","1"),],permutations = 50000, strata = in.full.meta[in.full.meta$month%in%c("0","1"),]$participant_id)


######### M1 - M6 ########
index_16 <- which(colnames(as.matrix(dis_in_all)) %in% in.full.meta$sid[in.full.meta$month%in%c("1","6")])

dist_in16<- as.matrix(dis_in_all)[index_16,index_16] %>%
  as.dist()

betadisper(dist_in16, group= factor(as.matrix(in.full.meta[in.full.meta$month%in%c("1","6"),"month"]))) %>% anova()


adonis2(dist_in16~month, data=in.full.meta[in.full.meta$month%in%c("1","6"),],
        permutations = 50000, strata = in.full.meta[in.full.meta$month%in%c("1","6"),]$participant_id)


######### M6- M12 ##############
index_612 <- which(colnames(as.matrix(dis_in_all)) %in% in.full.meta$sid[in.full.meta$month%in%c("6","12")])
dist_in612<- as.matrix(dis_in_all)[index_612,index_612] %>%
  as.dist()

betadisper(dist_in612, 
           group= factor(as.matrix(in.full.meta[in.full.meta$month%in%c("6","12"),"month"]))) %>% anova()


adonis2(dist_in612~month, data=in.full.meta[in.full.meta$month%in%c("6","12"),],
        permutations = 50000, strata = in.full.meta[in.full.meta$month%in%c("6","12"),]$participant_id)



##########    Figure 2B  ################
PCoAs<- ape::pcoa(dis_in_all)
barplot(PCoAs$values$Relative_eig[1:10])
# Dimension (i.e., Axis 1 (PCOA1))
Axis1.percent <-
  PCoAs$values$Relative_eig[[1]] * 100
# Dimension (i.e., Axis 2 (PCOA2))
Axis2.percent <-
  PCoAs$values$Relative_eig[[2]] * 100
## convert PCoA result into a data frame
PCoAs.data <-
  data.frame(
    sid = rownames(PCoAs$vectors),
    X = PCoAs$vectors[, 1],
    Y = PCoAs$vectors[, 2])
################
PCoAs.met.dat <-
  left_join(PCoAs.data, in.full.meta, by = "sid")
#############
centroids_all<- PCoAs.met.dat %>%
  group_by(month) %>%
  summarize(X=mean(X), Y=mean(Y)) %>%
  mutate(month = factor(month, levels=c(0,1,6,12), labels=c("Birth", "Month 1", "Month 6", "Month 12")))

### vector for colors ####

visit.cols = c("Birth"= "#999999",
               "Month 1"= "#FF7256",
               "Month 6"= "#56B4E9",
               "Month 12"= "darkblue")

light.visit= c("Birth"= "#CCCCCC",
               "Month 1"= "#FF7256",
               "Month 6"= "#A6DFFF",
               "Month 12"= "#6C85CC")

visit.01 = c("Birth"= "#999999",
             "Month 1"= "#E69F00")


PCoAs.met.dat %>%
  mutate(month = factor(month, levels=c(0,1,6,12), labels=c("Birth", "Month 1", "Month 6", "Month 12"))) %>%
  ggplot(., aes(x = X, y = Y)) +
  geom_point(aes(colour = month), size=3)  +
  scale_color_manual(values = visit.cols)+
  xlab(paste("PCoA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
  ylab(paste("PCoA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
  theme_bw() +
  stat_ellipse(aes(color=month), level=0.90)+ 
  theme(text = element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(1, 'lines')) + 
  geom_point(data=centroids_all, size=7, shape= 24,aes(fill=month))+
  scale_fill_manual(values= visit.cols)

###############################################
#### Panel C: Relative Abundance Over Time ####
###############################################

#####################################################################################
######################### Load in Packages ##########################################

library(data.table); library(tidyverse); library(microshades); library(patchwork)

#####################################################################################
######################### Load in Data ##############################################

## ASV for infants included in this analysis ##
in.dat <- read.csv(asv_fig2C.csv)

## Metadata needed for innfants included in this analysis ##
in.meta <- read.csv(metadata_fig2C.csv)

#####################################################################################
#### Removing zero reads samples #####

#remove samples with 0 reads
in.dat <- as.data.frame(in.dat[which(rowSums(in.dat[,-1]) > 0),]) %>%
  column_to_rownames(var="sample_id") %>% t() %>% as.data.frame() %>%
  rownames_to_column(.,var="X.OTU.ID")

#####################################################################################
##### Creating a taxonomy table #######


KINGDOM = 1; PHYLUM = 2; CLASS = 3; ORDER = 4; FAMILY = 5; GENUS = 6; SPECIES = 7

txnames <- in.dat$X.OTU.ID

taxonomy <- str_split_fixed(txnames, "\\;", 7) # create taxonomy table

rownames(taxonomy) <- txnames # the row names will be the full taxon names

colnames(taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # columns are each part of the classification


#####################################################################################
#######                 Agregating at Genus Level                              ######
### Create data frame that lists genus names and use this to aggregate by genus level

asv.genus<- taxonomy %>%
  as.data.frame() %>%
  rownames_to_column(var = "X.OTU.ID") %>%
  dplyr::select(X.OTU.ID, Genus) %>%
  inner_join(in.dat, by = "X.OTU.ID")

asv.genus <- aggregate(asv.genus[,-c(1:2)], by = list(asv.genus$Genus), FUN = sum)
rownames(asv.genus) <- asv.genus[,"Group.1"] # make genus names the row names
asv.genus[,"Group.1"] <- NULL
asv.genus <- t(asv.genus) # take the transpose so that samples are rows and taxa are columns
### Filter out taxa that appear in fewer than 20% of samples and rename unclassified genus
asv.genus <-asv.genus[,apply(asv.genus, MARGIN = 2, FUN = function (x) {sum(x > 0) > .10*length(x)} )]
asv.genus <- asv.genus[which(rowSums(asv.genus) > 115),]
## remove samples with data only for 2 taxa
asv.genus <- asv.genus[which(apply(asv.genus, 1, function(row) sum(row > 0)>2)),]
asv.genus <- data.table(asv.genus, keep.rownames = T)
asv.genus <- asv.genus %>%
  rename( "unclassified"="__") %>%
  select(-c("unclassified","g__uncultured")) %>%
  column_to_rownames(var="rn")


##### Genus aggregated relative abundance, long format dataframe ##########
asv.ra.l <- t(apply(asv.genus ,1, function(x) x/sum(x)))%>% 
  as.data.frame() %>% rownames_to_column(var="sample_id") %>%
  left_join(., in.meta) %>%
  column_to_rownames(var = "sample_id") %>%
  select(-c(sample_type, sample_source)) %>%
  pivot_longer( cols = -c("participant_id", "month"),
                names_to = "Genus", values_to = "rel_abund" )

#####################################################################################
# Phylum and Genus Level Composition
# To visualize the overall microbial landscape while maintaining clarity, taxa with extremely low abundance across all timepoints were filtered out.
# Low abundance is determined after calculating the mean relative abundance of a taxa at each time point #

# Make a vector list of low abundant genus across all timepoints
low.abun.phy <- taxonomy %>% as.data.frame() %>%
  filter(Genus %in% names(asv.genus)) %>%
  select(Genus, Phylum) %>%
  distinct() %>% left_join(., asv.ra.l , by="Genus") %>% 
  group_by(month,Phylum)  %>%
  mutate(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  pivot_wider(names_from = "month", values_from="mean_rel_abund") %>% 
  filter(`1` <1, `6`<1, `12`<1,`0`<1) %>% distinct(Genus) %>% pull(Genus)


### Now remove those low abundance Phylum and recalculate the relative abundance and make long format dataframe
asv.ra.l<- t(apply(asv.genus[, !colnames(asv.genus) %in% low.abun.phy] ,1, function(x) x/sum(x))) %>% 
  as.data.frame() %>% rownames_to_column(var="sample_id")%>%
  left_join(., in.meta) %>%
  column_to_rownames(var = "sample_id") %>%
  select(-c(sample_type, sample_source)) %>%
  pivot_longer( cols = -c("participant_id", "month"),
                names_to = "taxa", values_to = "rel_abund" )

##### Build asv.ra.l with Phylum info #########
asv.ra.l <- taxonomy %>%
  as.data.frame() %>%
  filter(Genus %in% unique(asv.ra.l$taxa)) %>%
  select(Genus, Phylum) %>%
  distinct() %>%
  rename(taxa = Genus) %>%
  right_join(asv.ra.l) %>%
  mutate(
    taxa   = gsub("g__", "", taxa),
    Phylum = gsub("p__", "", Phylum)
  ) %>%
  rename(Abundance = rel_abund, Genus = taxa) %>%
  mutate(combination_variable = paste(Phylum, Genus, sep = "-"))

#### Identify low-abundance genera to collapse into "Other" ########
low.f <- asv.ra.l %>%
  group_by(month, Genus, Phylum) %>%                        # add Phylum here
  summarize(mean_rel_abund = 100 * mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = month, values_from = mean_rel_abund) %>%
  rename(Birth = `0`, M1 = `1`, M6 = `6`, M12 = `12`) %>%
  filter((Phylum == "Actinobacteriota" & Birth < 1.1 & M1 < 1.1 & M6 < 1.1 & M12 < 1.1) |
           (Phylum == "Firmicutes"       & Birth < 5   & M1 < 5   & M6 < 5   & M12 < 5)) %>%
  pull(Genus) %>%
  unique()

#### Build lookup table for Top_Genus / Top_Phylum #########
group_lookup <- asv.ra.l %>%
  mutate(Genus = if_else(Genus %in% low.f, "Other", Genus), group = paste(Phylum, Genus, sep = "-")) %>%
  distinct(group, Genus, Phylum) %>%
  rename(Top_Genus = Genus, Top_Phylum = Phylum)

#####################################################################################
##### Final summarized table ########
## This table is the final dataframe with the relative abundance formatted for alluvial plot
plot.df <- asv.ra.l %>%
  rename(Sample = participant_id) %>%
  mutate(Genus = if_else(Genus %in% low.f, "Other", Genus), group = paste(Phylum, Genus, sep = "-")) %>%
  group_by(month, Sample, group) %>%
  summarize(rel_abund = sum(Abundance), .groups = "drop") %>%
  group_by(month, group) %>%
  summarize(mean_rel_abund = 100 * mean(rel_abund), .groups = "drop") %>%
  mutate(month = factor(month,
                        levels = c("0", "1", "6", "12"),
                        labels = c("Birth", "Month 1", "Month 6", "Month 12"))) %>%
  left_join(group_lookup, by = "group")


tax.colors <- c(
  # Purples (Actinobacteriota)
  "Actinobacteriota-Bifidobacterium" = "#54278f",  # Dark Indigo
  "Actinobacteriota-Collinsella" = "#756bb1",  # Blue Violet
  "Actinobacteriota-Corynebacterium" = "#9e9ac8",  # Orchid
  "Actinobacteriota-Lawsonella" = "#bcbddc",  # Light Pink
  "Actinobacteriota-Varibaculum" = "#dadaeb",  # Soft Pink
  "Actinobacteriota-Other"  = "#f2f0f7",  # Very Light Pink
  
  # Blues (Bacteroidota)
  "Bacteroidota-Bacteroides" = "#08519c",  # Dark Blue
  "Bacteroidota-Parabacteroides" = "#3182bd",  # Royal Blue
  "Bacteroidota-Porphyromonas" = "#6baed6",  # Light Sky Blue
  "Bacteroidota-Prevotella" = "#bdd7e7",  # Powder Blue
  
  # Greens (Proteobacteria)
  "Proteobacteria-Escherichia-Shigella" = "#006d2c",  # Dark Green
  "Proteobacteria-Haemophilus" = "#31a354",  # Forest Green
  "Proteobacteria-Klebsiella" = "#74c476",  # Medium Sea Green
  "Proteobacteria-Sutterella" = "#bae4b3",  # Light Green
  
  #  Reds (Firmicutes)
  "Firmicutes-Anaerococcus" = "#772C4B",  # Dark Brown
  "Firmicutes-Enterococcus" = "#944155",  # Medium Brown
  "Firmicutes-Finegoldia" = "#B2585C",  # Chocolate
  "Firmicutes-Peptoniphilus" = "#CF705E",  # Orange
  "Firmicutes-Staphylococcus" =  "#E38E66",  # Peach
  "Firmicutes-Streptococcus" =   "#F2BC8F",
  "Firmicutes-Other" = "#feedde"
  
)

tax_colors_df <- data.frame(
  group = names(tax.colors),
  color = unname(tax.colors),
  stringsAsFactors = FALSE)


fill_df <- plot.df %>% 
  select(group,Top_Phylum, Top_Genus) %>% 
  unique() %>% 
  left_join(., tax_colors_df)


legends <- lapply(split(plot.df, plot.df$Top_Phylum), function(x) {
  genome <- unique(x$Top_Phylum)
  x <- x %>%
    mutate(Top_Genus = factor(Top_Genus, 
                              levels = c(setdiff(sort(unique(Top_Genus)), "Other"), "Other")))
  patchwork::wrap_elements(
    full = cowplot::get_legend(
      ggplot(x, aes(month, mean_rel_abund, fill = Top_Genus)) + 
        geom_col(color = "black") +
        scale_fill_manual(name = genome, 
                          values = setNames(fill_df$color[fill_df$Top_Phylum == genome], 
                                            fill_df$Top_Genus[fill_df$Top_Phylum == genome])) +
        theme(legend.justification = c(0, 1), 
              legend.spacing.y = unit(0.1, "cm"),
              text=element_text(size=20, face="italic"))))
})


plot.df <- plot.df %>% mutate(group =  factor(group, 
                                              levels = c( "Actinobacteriota-Bifidobacterium",
                                                          "Actinobacteriota-Collinsella",
                                                          "Actinobacteriota-Corynebacterium",
                                                          "Actinobacteriota-Lawsonella",
                                                          "Actinobacteriota-Varibaculum",
                                                          "Actinobacteriota-Other",
                                                          "Bacteroidota-Bacteroides", 
                                                          "Bacteroidota-Parabacteroides",
                                                          "Bacteroidota-Porphyromonas", 
                                                          "Bacteroidota-Prevotella",
                                                          "Firmicutes-Anaerococcus", 
                                                          "Firmicutes-Enterococcus",
                                                          "Firmicutes-Finegoldia", 
                                                          "Firmicutes-Peptoniphilus",
                                                          "Firmicutes-Staphylococcus", 
                                                          "Firmicutes-Streptococcus",
                                                          "Firmicutes-Other", 
                                                          "Proteobacteria-Escherichia-Shigella",
                                                          "Proteobacteria-Haemophilus", 
                                                          "Proteobacteria-Klebsiella",
                                                          "Proteobacteria-Sutterella")))


p1 <- plot.df %>% ggplot(., aes(x=month, y=mean_rel_abund, fill=group,stratum = group, alluvium = group)) +
  ggalluvial::geom_alluvium(aes(fill = group),color = "darkgray", na.rm = TRUE) +
  ggalluvial::geom_stratum(width=0.6)+
  scale_fill_manual(values=tax.colors)+
  scale_x_discrete(labels = c("Birth" = "Birth", "Month 1" = "M1", "Month 6" = "M6", "Month 12" = "M12"))+
  theme_bw() +
  labs(y="Mean Relative Abundance (%)",tag = "C.")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.tag = element_text(size = 25))

# Distribute legends into multiple columns for better spacing
legend_layout <- wrap_plots(
  wrap_plots(legends[1:4], ncol = 1, nrow=4), ncol = 1) +
  theme(legend.key.size = unit(1.5, "cm"),  
        legend.text = element_text(size = 90))

# Adjust final layout
final_plot <- p1 + legend_layout + 
  plot_layout(guides = "collect", widths = c(2, 1))

final_plot
