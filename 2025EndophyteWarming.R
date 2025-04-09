
library(lmerTest)
library(vegan)
library(hillR)
library(ggplot2)
library(matrixStats)
library(multcomp)
library(emmeans)
library(multcompView)
library(egg)
library(dplyr)
library(MASS) 
library(reshape2) 
library(reshape) 
library(MuMIn)
library(lme4)
library(performance)
library(betareg)
library(lme4)   
library(glmmTMB) 
library(DESeq2)
library(forcats)
library(tidytext)
library(grid)
library(cowplot)
library(ggpubr)
library(tidyr)


mytheme <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                   panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                   panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                   #panel.border=element_blank(), #gets rid of square going around the entire graph
                   panel.border = element_rect(colour = "black", fill=NA, size=1),
                   axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                   axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                   axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                   axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                   axis.text.x = element_text(size=10, color="black"), #size of x-axis text
                   axis.text.y = element_text( size=10, color="black"))




#------ read data --------


meta.raw <- read.csv("WarmingEndo_meta.csv", header = TRUE, row.names = 1)
meta.raw$Plot <- as.factor(meta.raw$Plot)
meta <- meta.raw


data.raw <- read.csv("WarmingEndo_OTU.csv", header = TRUE, row.names = 1)


data.raw.f <- t(data.raw[,0:201])

fun.rel <- decostand(data.raw.f, method = "hellinger")

fun.full <- data.frame(code = meta$code, fun.rel)


ddatf<-vegdist(as.matrix(fun.rel), "jaccard")


meta.all <- data.frame(meta.raw, data.raw.f)


meta$Species <- ordered(meta$Species, levels = c("ACLE", "FETH", "POPR", "SCSC", "SONU", "SPCO", "ELCA_Ep", "ELCA_Em", "FEAR_Ep", "FEAR_Em"))

meta$Species <- ordered(meta$Species, levels = c("FEAR_Ep" ,"FEAR_Em", "ELCA_Ep", "ELCA_Em",  "SPCO", "SONU","SCSC", "POPR", "FETH", "ACLE"))

meta$Sample_Site <- ordered(meta$Sample_Site, levels = c("RMBL" ,"Kessler", "KY"))

meta$Species.tissue.plot <- paste(meta$Species, meta$root_leaves, meta$Plot, sep = ".")


#---- Alpha Diversity ---------


q_values <- c(0, 1, 2)

# Loop over each q value to calculate and store alpha diversity in meta dataframe
for (q in q_values) {
  meta[[paste0("alphaq", q)]] <- apply(data.raw.f, 1, hill_taxa, q = q)
}




leaf <- subset(meta, root_leaves == "leaves")

root <- subset(meta, root_leaves == "root")


# Loop to fit models for each alpha diversity value
modelsl <- list()  
for (q in q_values) {
  formula <- as.formula(paste0("alphaq", q, " ~ Treatment + (1 | Species/Plot/Sample_Site)"))
  modelsl[[paste0("modelq", q)]] <- glmmTMB(formula, data = leaf, family = gaussian)
}

modelsr <- list()  # Create an empty list to store the models
for (q in q_values) {
  formula <- as.formula(paste0("alphaq", q, " ~ Treatment + (1 | Species/Plot/Sample_Site)"))
  modelsr[[paste0("modelq", q)]] <- glmmTMB(formula, data = root, family = gaussian)
}

Anova(modelsl$modelq0)
Anova(modelsl$modelq1)
Anova(modelsl$modelq2)

Anova(modelsr$modelq0)
Anova(modelsr$modelq1)
Anova(modelsr$modelq2)

lc <- subset(leaf, Treatment == "C")
lw <- subset(leaf, Treatment == "W")

rc <- subset(root, Treatment == "C")
rw <- subset(root, Treatment == "W")
rc <- subset(root.omit, Treatment == "C")
rw <- subset(root.omit, Treatment == "W")




#---------- Colonization ------------


root.omit <- na.omit(root)


leaf.omit <- na.omit(leaf)



leaf.omit$log.tot <- log(leaf.omit$total +1)

col_values <- c("septate", "aseptate", "vesicle", "arbuscule")

for (c in col_values) {
  root.omit[[c]] <- root.omit[[c]] + 0.001
}

modelsr.col <- list()  
for (c in col_values) {
  formula <- as.formula(paste0(c, " ~ Treatment + (1 | Species/ Plot/Sample_Site)"))
  modelsr.col[[c]] <- glmmTMB(formula, data = root.omit, family = beta_family())
}



model.total <- glmmTMB(log.tot ~ Treatment + (1 | Species/Sample_Site) +(1|Plot/Sample_Site), data = leaf.omit,
                 family = gaussian)




Anova(model.total)
Anova(modelsr.col$septate, type =2)
Anova(modelsr.col$aseptate, type =2)
Anova(modelsr.col$arbuscule, type =2)
Anova(modelsr.col$vesicle, type =2)





#------------Species (Supplemental)-------------
modelsls <- list()  # Create an empty list to store the models
for (q in q_values) {
  formula <- as.formula(paste0("alphaq", q, " ~ Species + (1 | Plot/Sample_Site)"))
  modelsls[[paste0("modelqs", q)]] <- glmmTMB(formula, data = leaf, family = gaussian)
}

modelsrs <- list()  # Create an empty list to store the models
for (q in q_values) {
  formula <- as.formula(paste0("alphaq", q, " ~ Species + (1 | Plot/Sample_Site)"))
  modelsrs[[paste0("modelqs", q)]] <- glmmTMB(formula, data = root, family = gaussian)
}

Anova(glmmTMB(alphaq0 ~ Species +(1|Plot/Sample_Site), data = leaf.omit,
                         family = gaussian))

Anova(modelsls$modelqs0, type =2)
Anova(modelsls$modelqs1, type =2)
Anova(modelsls$modelqs2, type =2)


Anova(modelsrs$modelqs0, type =2)
Anova(modelsrs$modelqs1, type =2)
Anova(modelsrs$modelqs2, type =2)

l0 <- data.frame(cld(emmeans(modelsls$modelqs0, ~ Species), Letters = c("abcde"), sort = FALSE))
l1 <- data.frame(cld(emmeans(modelsls$modelqs1, ~ Species), Letters = c("abcde"), sort = FALSE))
l2 <- data.frame(cld(emmeans(modelsls$modelqs2, ~ Species), Letters = c("abcde"), sort = FALSE))

r0 <- data.frame(cld(emmeans(modelsrs$modelqs0, ~ Species), Letters = c("abcde"), sort = FALSE))
r1 <- data.frame(cld(emmeans(modelsrs$modelqs1, ~ Species), Letters = c("abcde"), sort = FALSE))
r2 <- data.frame(cld(emmeans(modelsrs$modelqs2, ~ Species), Letters = c("abcde"), sort = FALSE))




model.total.sp <- glmmTMB(log.tot ~ Species +(1|Plot/Sample_Site), data = leaf.omit,
                          family = gaussian)

model.root.sp <- glmmTMB(septate ~ Species +(1|Plot/Sample_Site), data = root.omit,
                         family = gaussian)

model.root.sp.arb <- glmmTMB(arbuscule ~ Species +(1|Plot/Sample_Site), data = root.omit,
                         family = gaussian)

model.root.sp.ves <- glmmTMB(vesicle ~ Species +(1|Plot/Sample_Site), data = root.omit,
                             family = gaussian)

model.root.sp.asep <- glmmTMB(aseptate ~ Species +(1|Plot/Sample_Site), data = root.omit,
                             family = gaussian)


Anova(model.root.sp.arb)


rs <- data.frame(cld(emmeans(model.root.sp, ~ Species), Letters = c("abcde"), sort = FALSE))
rv <- data.frame(cld(emmeans(model.root.sp.ves, ~ Species), Letters = c("abcde"), sort = FALSE))
ra <- data.frame(cld(emmeans(model.root.sp.asep, ~ Species), Letters = c("abcde"), sort = FALSE))




#-------- Beta diversity----------

ti.sp.p.site <- paste(meta$root_leaves, meta$Species, meta$Plot, meta$Sample_Site, sep = ".")
p.site <- paste( meta$Plot, meta$Sample_Site, sep = ".")


meta$strata_factor <- interaction(meta$Sample_Site, meta$Plot)

model <- dbrda(ddatf ~ root_leaves * Treatment + Condition(Sp + strata_factor), data = meta)



anova(model, by = "margin", permutations = 999, strata = meta$strata_factor)


RsquareAdj(model)

summary(model)

scoresfun <- scores(model,display="sites")



meta$MDS1 <- scoresfun[,1]
meta$MDS2 <- scoresfun[,2]



#--- Figure 2--------


quartz()

sep.r <- ggplot(root, aes(x = Treatment, y = (septate*100), col = Treatment)) +
  geom_boxplot(size = 1, width = .9, position=position_dodge(), alpha = 0, aes(col= Treatment))+
  geom_point(size = 2,alpha = 0.3, position=position_jitterdodge(), shape = 24, aes(fill = Treatment)) +
  ylab("")+ 
  xlab("")+ 
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  mytheme+ theme(legend.position = "none",  plot.margin = unit(c(0,0,-.25,0), "cm"),     axis.text.x=element_blank())

sep.l <- ggplot(leaf, aes(x = Treatment, y = (total), col = Treatment)) +
  geom_boxplot(size = 1, width = .9, position=position_dodge(), alpha = 0, aes(col= Treatment))+
  geom_point(size = 2, alpha = 0.3,position=position_jitterdodge(), shape = 21, aes(fill = Treatment)) +
  ylab("")+ 
  xlab("")+ 
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  scale_y_continuous(limits = c(0,26))+
  mytheme+ theme(legend.position = "none",   plot.margin = unit(c(0,0,-.25,0), "cm"),    axis.text.x=element_blank())

q1.l <- ggplot(leaf, aes(x = Treatment, y = alphaq1, col = Treatment)) +
  geom_boxplot(size = 1, width = .9, position=position_dodge(), alpha = 0, aes(col= Treatment))+
  geom_point(size = 2, alpha = 0.3,position=position_jitterdodge(), shape = 21, aes(fill = Treatment)) +
  ylab("")+ 
  xlab("")+ 
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  scale_y_continuous(limits = c(0,70))+
  mytheme+ theme(legend.position = "none", plot.margin = unit(c(0,0,-.25,0),"cm"),      axis.text.x=element_blank())

q1.r <- ggplot(root, aes(x = Treatment, y = alphaq1, col = Treatment)) +
  geom_boxplot(size = 1, width = .9, position=position_dodge(), alpha = 0, aes(col= Treatment))+
  geom_point(size = 2, alpha = 0.3,position=position_jitterdodge(), shape = 24, aes(fill = Treatment)) +
  ylab("")+ 
  xlab("")+ 
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_y_continuous(limits = c(0,85))+
  mytheme+ theme(legend.position = "none",   plot.margin = unit(c(0,0,-.25,0), "cm"),     axis.text.x=element_blank())

fungi.ord <- ggplot(meta, aes(x = (MDS1), y = MDS2, col = Treatment, fill = Treatment)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=root_leaves)) +
  ylab("")+ 
  xlab("")+ 
  scale_shape_manual(values=c(21, 24))+
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  mytheme+ theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm"))

warming_total <- grid.arrange(                       # First row with one plot spaning over 2 columns
  arrangeGrob(sep.l, sep.r, ncol = 1),
  arrangeGrob(q1.l, q1.r, ncol = 1),
  fungi.ord,      
  widths = c(1,1,1.5), # Second row with 2 plots in 2 different columns
  nrow = 1)                       # Number of rows



    

#--------- RII Colonization and alpha diversity----------

meta$col <- ifelse(meta$root_leaves == "leaves", meta$total, meta$septate*100) + 0.001

merge.trt <- meta %>%
  filter(Treatment == "C") %>%
  left_join(
    meta %>%
      filter(Treatment == "W") %>%
      select(Species.tissue.plot, 11:23), 
    by = "Species.tissue.plot"
  ) %>%
  na.omit() %>%
  mutate(
    # Calculate the RII variables
    rii.q1 = (alphaq1.y - alphaq1.x) / (alphaq1.y + alphaq1.x),
    rii.col = (col.y  - col.x ) / (col.y  + col.x ),
    lrr.q1 = log(alphaq1.y/alphaq1.x),
    lrr.col = log(col.y/ col.x ) 
  )




model.riiq1 <- glmmTMB(rii.q1 ~ Sample_Site*root_leaves, data = merge.trt,
                       family = gaussian)

model.riicol <- glmmTMB(rii.col ~ Sample_Site*root_leaves, data = merge.trt,
                       family = gaussian)

model.lrrq1 <- glmmTMB(lrr.q1 ~ Sample_Site*root_leaves, data = merge.trt,
                       family = gaussian)

model.lrrcol <- glmmTMB(lrr.col ~ Sample_Site*root_leaves, data = merge.trt,
                        family = gaussian)

Anova(model.riiq1, test = "Chisq")
Anova(model.riicol, test = "Chisq")
Anova(model.lrrq1, test = "Chisq")
Anova(model.lrrcol, test = "Chisq")




#---- RII for beta diversity--------


# Create distance matrix and convert to data frame
ddatf.df <- data.frame(as.matrix(vegdist(fun.rel), "jaccard"))
ddatf.df$name <- row.names(ddatf.df)

# Melt the distance matrix data frame
ddat.m <- melt(ddatf.df, id = "name")

# Merge with metadata and subset relevant columns
merge.name <- data.frame(name = row.names(meta), meta[, c(2, 5:9,16)])
merge.ddatf.df <- merge(ddat.m, merge.name, by = "name", all.x = TRUE)

# Merge again by variable column
merge.ddatf.df1 <- merge(merge.ddatf.df, merge.name, by.x = "variable", by.y = "name", all.x = TRUE)

# Subset where the species and root_leaves match, and there are duplicates
merge.sub <- subset(merge.ddatf.df1, 
                    Species.x == Species.y & 
                      root_leaves.x == root_leaves.y & 
                      duplicated(value) & 
                      value != 0)



# Create a new column for plotID
merge.sub$plotID <- paste(merge.sub$Species.x, merge.sub$root_leaves.x, sep = ".")

#merge.sub$plotID <- paste(merge.sub$Sample_Site.x, merge.sub$Plot.x, merge.sub$root_leaves.x, sep = ".")


# Subset control and treatment data
merge.sub.c <- subset(merge.sub, Treatment.x == "C" & Treatment.y == "C")[, c(3, 18)]
merge.sub.w <- subset(merge.sub, Treatment.x != Treatment.y & Plot.x == Plot.y)[, c(2:10, 18)]

# merge.sub.c <- subset(merge.sub, Treatment.x == "C" & Treatment.y == "C" )[, c(3, 16)]
# merge.sub.w <- subset(merge.sub, Treatment.x != Treatment.y & Species.x == Species.y )[, c(2:9, 16)]

# Aggregate control data by plotID
ag.c <- aggregate(value ~ plotID, merge.sub.c, mean)
colnames(ag.c)[2] <- "mean.c"

# Merge with treatment data
mp <- merge(merge.sub.w, ag.c, by = "plotID")

# Calculate RII and LRR
mp$RII <- (mp$value - mp$mean.c) / (mp$value + mp$mean.c)
mp$LRR <- log(mp$value / mp$mean)

model.riibet <- glmmTMB(RII ~ Sample_Site.x*root_leaves.x, data = mp,
                        family = gaussian)
model.lrrbet <- glmmTMB(LRR ~ Sample_Site.x*root_leaves.x, data = mp,
                        family = gaussian)


Anova(model.riibet, test = "Chisq")
Anova(model.lrrbet, test = "Chisq")

#--------- Figure 4------------

rii <- ggplot(effect.merge, aes(x = root_leaves, y = RII, col = Sample_Site)) +
  geom_hline(yintercept=0, color = "grey", size = 1)+
  geom_boxplot(size = 1, width = .7, position=position_dodge(width = 0.8), alpha = 0, aes(col= Sample_Site))+
  geom_point(size = 3, alpha = 0.3, position=position_jitterdodge(), aes(  shape = root_leaves)) +
  ylab("")+ 
  xlab("")+ 
  scale_y_continuous(breaks=c(-0.2,0,0.2), limits = c(-0.2, 0.2))+
  mytheme+
  theme(legend.position = "none",   
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_blank())

rii.col <- ggplot(effect.merge, aes(x = Sample_Site, y = rii.col, col = Sample_Site)) +
  geom_hline(yintercept=0, color = "grey", size = 1)+
  geom_boxplot(size = 1, width = .7, position=position_dodge(width = 0.8), alpha = 0, aes(col= Sample_Site))+
  geom_point(size = 3, alpha = 0.3, position=position_jitterdodge(), aes(shape = root_leaves)) +
  ylab("")+ 
  xlab("")+ 
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1), limits = c(-1.25, 1.25))+
  mytheme+
  theme(legend.position = "none",   
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_blank())

rii.q1 <- ggplot(effect.merge, aes(x = root_leaves, y = rii.q1, col = Sample_Site)) +
  geom_hline(yintercept=0, color = "grey", size = 1)+
  geom_boxplot(size = 1, width = .7, alpha = 0, col= 'black')+
  geom_point(size = 3, alpha = 0.3, position=position_jitterdodge(), aes(col = Sample_Site, shape = root_leaves)) +
  ylab("")+ 
  xlab("")+ 
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1), limits = c(-1,1))+
  mytheme+
  theme(legend.position = "none",   
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_blank())


rii.merge <- ggarrange(rii.col, rii.q1, rii, ncol = 1)

ggsave(rii.merge, filename = "241205.rii.merge.png")



#-------- Taxa differences between treatments ----------


taxa_in <- data.raw[202:ncol(data.raw)]


taxa_in$species <- taxa_in$taxonomy


ranks <- c("phylum", "class", "order", "family", "genus")
suffixes <- c("c__", "o__", "f__", "g__", "s__")

# Loop through each rank and apply the gsub to extract the relevant part of the species name
for (i in 1:length(ranks)) {
  taxa_in[[ranks[i]]] <- gsub(paste0(";", suffixes[i], ".*$"), "", taxa_in$species)
}




root.sub <- subset(meta.all, root_leaves == "root" )
leaf.sub <- subset(meta.all, root_leaves != "root" )


root.meta <- root.sub[,0:15]
leaf.meta <- leaf.sub[,0:15]

root.fun <- t(root.sub[,16:ncol(root.sub)])
leaf.fun <- t(leaf.sub[,16:ncol(leaf.sub)])

root_ASVs <- data.frame(root.fun)
leaf_ASVs <- data.frame(leaf.fun)



root_ASVs$genus <- taxa_in[rownames(root_ASVs), "genus"]
leaf_ASVs$genus <- taxa_in[rownames(leaf_ASVs), "genus"]

root_ASVs$sum <- 1
leaf_ASVs$sum <- 1


agg.root <- aggregate(. ~ genus, data=root_ASVs, FUN=sum)
agg.leaf <- aggregate(. ~ genus, data=leaf_ASVs, FUN=sum)







# DESeq2 Analysis function (to avoid repetition)
run_DESeq <- function(count_data, metadata, design_formula) {
  dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = metadata, 
                                design = design_formula, tidy = TRUE)
  dds <- DESeq(dds)
  res <- results(dds)
  return(as.data.frame(res))
}

# Run DESeq2 for root and leaf
res.root <- run_DESeq(agg.root[,-ncol(agg.root)], root.meta, ~Treatment)
res.leaf <- run_DESeq(agg.leaf[,-ncol(agg.leaf)], leaf.meta, ~Treatment)


# Print summaries (optional)
summary(res.root)
summary(res.leaf)







#---------Figure 3 --------



# Define a function to process data for both leaf and root

process_data <- function(data, agg_data) {
                         data$means <- rowMeans(decostand(agg_data[,-c(1,ncol(agg_data))], MARGIN = 2, "total")) * 100
                         data$rich <- agg_data$sum
                         data$pres <- (rowSums(decostand(agg_data[,-c(1,ncol(agg_data))], MARGIN = 2, "pa")) / (ncol(agg_data)-2)) * 100
                         data$Genus <- gsub(".*g__", "", rownames(data))
                         data$species <- gsub(".*s__", "", rownames(data))
                         data$Phylum <- gsub(";c__.*$", "", rownames(data))
                         data$Class <- gsub(";o__.*$", "", rownames(data))
                         data$Order <- gsub(";f__.*$", "", rownames(data))
                         data$Family <- gsub(";g__.*$", "", rownames(data))
                         data$Taxon <- gsub(";s__.*$", "", data$Genus)
                         return(data)
}


# Process leaf data
res.leaf <- process_data(res.leaf, agg.leaf)


# Process root data
res.root <- process_data(res.root, agg.root)

res.root$Compartment <- "Root"
res.leaf$Compartment <- "Leaf"

res.all <- bind_rows(res.root, res.leaf)

taxa.in.merge <- taxa_in %>% dplyr::select(Taxon, Lifestyle) %>% distinct()


res.all.merge <- left_join(res.all, taxa.in.merge, by = c("Taxon" = "Taxon"))

res.all.fig <- res.all.merge %>% filter(padj <= 0.05 & Lifestyle != "Unassigned")

#Table S2
res.all.sig <- res.all.merge %>% filter(padj <= 0.05)




full <- res.all.fig %>%
ggplot() +
  geom_point(alpha =0.8, col = "black", aes(y = reorder_within(Taxon, log2FoldChange, Lifestyle), x = log2FoldChange, size = pres, fill = Lifestyle, shape = Compartment)) +
  ylab("")+ 
  xlab("")+ 
  scale_fill_manual(values=c("#88CCEE",  "#DDCC77", "#117733", "#6699CC", "#AA4499")) +
  scale_shape_manual(values=c(21, 24))+
  geom_vline(xintercept=0, color = "grey", size = 1)+
  mytheme +
  scale_y_reordered()+
  facet_wrap(~Lifestyle, scales = "free_y", ncol = 6, strip.position = "top")+
  theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.5, linetype = 2),
        legend.text=element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_blank())


#------------ Supp Figure -------------



res.neg <- subset(res.all, log2FoldChange < 0 & Lifestyle != "Saprotroph")
res.pos <- subset(res.all, log2FoldChange > 0 & Lifestyle != "Saprotroph")

res.phy <- res.all.sig

res.phy$dir <- ifelse(res.phy$log2FoldChange > 0, "Warmed", "Ambient")

res.phy$group <- paste(res.phy$Lifestyle, res.phy$dir)

#ggsave(full, filename ="ref.full.png")


top_15_genera <-   res.phy %>% 
  group_by(Class) %>%
  tally(means) %>%
  top_n(10) %>% 
  select(Class)




#Keep the top 15 but rename everything else starting at the 16th most abundant as Other 
res.phy1 <- res.phy %>% mutate(Class_renamed=if_else(Class %in% top_15_genera$Class, Class, "Other"))

res.phy1$Lifestyle <- ifelse(is.na(res.phy1$Lifestyle), "Unassigned", res.phy1$Lifestyle)

res.phy1$Class_renamed <- gsub("NA", "Unassigned Ascomycota", res.phy1$Class_renamed)

res.phy1$Class_renamed <- fct_rev(fct_reorder(res.phy1$Class_renamed, res.phy$means, sum))

#Put Other last

res.phy1$Lifestyle <- ifelse(is.na(res.phy1$Lifestyle), "Unassigned", res.phy1$Lifestyle)

res.phy1$Class_renamed <- fct_relevel(res.phy1$Class_renamed, "Other", after = Inf) 

quartz()

ri.a <- ggplot(res.phy1, aes(x=dir, y=rich, fill=Class_renamed)) +
  geom_bar(stat="identity") +
  ylab("") +
  xlab("") +
  mytheme+
  theme(legend.position="right",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_manual(values =c( "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#ffcc33", "#AA4499", 
                               "#44AA99","#66ff66","#882255","#888888",  "#999933",  "#661100","#ff355e", "#6699CC",  "black"))+
  facet_grid(Compartment~Lifestyle)

me.a <- ggplot(res.phy1, aes(x=dir, y=means, fill=Class_renamed)) +
  geom_bar(stat="identity") +
  ylab("") +
  xlab("") +
  mytheme+
  theme(legend.position="right",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_manual(values =c( "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#ffcc33", "#AA4499", 
                               "#44AA99","#66ff66","#882255","#888888",  "#999933",  "#661100","#ff355e", "#6699CC",  "black"))+
  facet_grid(Compartment~Lifestyle)

ri <- ggplot(res.phy1, aes(x=dir, y=rich, fill=Class_renamed)) +
  geom_bar(stat="identity") +
  ylab("") +
  xlab("") +
  mytheme+
  theme(legend.position="right",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  scale_fill_manual(values =c( "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#ffcc33", "#AA4499", 
                               "#44AA99","#66ff66","#882255","#888888",  "#999933",  "#661100","#ff355e", "#6699CC",  "black"))

me <- ggplot(res.phy1, aes(x=dir, y=means, fill=Class_renamed)) +
  geom_bar(stat="identity") +
  ylab("") +
  xlab("") +
  mytheme+
  theme(legend.position="right",
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"))+
  scale_fill_manual(values =c( "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#ffcc33", "#AA4499", 
                               "#44AA99","#66ff66","#882255","#888888",  "#999933",  "#661100","#ff355e", "#6699CC",  "black"))

ri.a.nl <- ri.a + theme(legend.position='none', plot.margin = margin(0, 0.25, 0, -.25, "cm"))
me.a.nl <- me.a + theme(legend.position='none', plot.margin = margin(0, 0.25, 0, -.25, "cm"))
ri.nl <- ri + theme(legend.position='none', plot.margin = margin(0, -.5, 0, -.5, "cm"))
me.nl <- me + theme(legend.position='none', plot.margin = margin(0, -.5, 0, -.5, "cm"))



legend <- get_legend(me)

legend <- ggpubr::get_legend(me)

grid.newpage()
le <- as_ggplot(legend)

#with legend
#ggarrange(ri, me, ncol = 1, align = "v")

#No legend
nofacet.fungi <- ggarrange(ri.nl, me.nl, ncol = 1, align = "v")

facet.fungi <- ggarrange(ri.a.nl, me.a.nl, ncol = 1, align = "v")


quartz()

taxa_total <- grid.arrange(                       # First row with one plot spaning over 2 columns
  nofacet.fungi,
  facet.fungi,
  le,      
  widths = c(1.2,2,1), # Second row with 2 plots in 2 different columns
  nrow = 1)  



ggsave(taxa_total, filename = "241209_taxa_total.png")


me.a <- ggplot(res.phy1, aes(x=dir, y=means, fill=Class_renamed)) +
  geom_bar(stat="identity") +
  ylab("") +
  xlab("") +
  mytheme+
  theme(legend.position="right",
        legend.title=element_blank(),
        strip.background = element_blank())+
  scale_fill_manual(values =c( "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#ffcc33", "#AA4499", 
                               "#44AA99","#66ff66","#882255","#888888",  "#999933",  "#661100","#ff355e", "#6699CC",  "black"))+
  facet_grid(Compartment~Lifestyle)



rd <- res.phy1 %>% filter(Class_renamed == c("Dothideomycetes"))
rs <- res.phy1 %>% filter(Class_renamed == c("Sordariomycetes"))

sum(count(rd$Genus)$freq)/111



#----------------Figure 5-----------


bole <- read.csv("WarmingEndo_metabolome.csv")

fun.full <- data.frame(code = meta$code, fun.rel)





fun.meta.bole <- merge(bole, fun.full, by.x = "code", 
                       by.y = "code", all.x = TRUE, all.y = FALSE)


fun.meta.bole <- na.omit(fun.meta.bole)



meta.fm <- fun.meta.bole[,0:10]

bole.fm <- fun.meta.bole[,11:124]
fun.fm <- fun.meta.bole[,125:ncol(fun.meta.bole)]


lome.rel <- scale(bole.fm) - min(scale(bole.fm))



ddatb<-as.matrix(vegdist(lome.rel, "jaccard"))

ddatf<-as.matrix(vegdist(fun.fm, "jaccard"))


colnames(ddatb) <- meta.fm$code
colnames(ddatf) <- meta.fm$code

db1 <- data.frame(code1 = meta.fm$code, Site1 = meta.fm$Sample_Site, Species1= meta.fm$Sp, Root_leaves1= meta.fm$root_leaves, Treatment1 = meta.fm$Treatment, ddatb)
df1 <- data.frame(code = meta.fm$code, Site = meta.fm$Sample_Site, Species= meta.fm$Sp, Treatment = meta.fm$Treatment, ddatf)

b.melt <- melt(db1, ID = c("code1", "Site1", "Species1", "Root_leaves1", "Treatment"))
f.melt <- melt(df1, ID = c("code", "Site", "Species", "Treatment"))

b.melt$fun.dist <- f.melt$value



dist.full <- merge(b.melt, meta.fm, by.x = "variable", 
                   by.y = "code", all.x = TRUE, all.y = FALSE)

dist.full.sub <- subset(dist.full, Root_leaves1 == root_leaves & Treatment1 == Treatment & variable != code1 & Site1 == Sample_Site & duplicated(dist.full$value) == "TRUE" )


model <- lmer(value ~ fun.dist*Root_leaves1*Treatment1 + (1|Sample_Site), data = dist.full.sub)

Anova(model, type =2, by = 'margin')


quartz()
df_rl <- ggplot(dist.full.sub, aes(x = value, y = fun.dist, col = Treatment)) +
  geom_point(size = 3.5, alpha = 0.4, aes(shape = Root_leaves1, fill = Treatment)) +
  stat_smooth(method = "lm", linewidth = .8, alpha = 0.8, col = "black", linetype = 1, se = TRUE)+
  ylab("")+ 
  xlab("")+ 
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  scale_shape_manual(values=c(21, 24))+
  mytheme+
  facet_grid(Root_leaves1~Treatment)+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")




CR <- subset(dist.full.sub, Root_leaves1 == "root" & Treatment1 == "C")
CL <- subset(dist.full.sub, Root_leaves1 == "leaves" & Treatment1 == "C")
WL <- subset(dist.full.sub, Root_leaves1 == "leaves" & Treatment1 == "W")
WR <- subset(dist.full.sub, Root_leaves1 == "root" & Treatment1 == "W")

model.CL <- lm(scale(fun.dist) ~ scale(value), data = CL)
model.CR <- lm(scale(value) ~ scale(fun.dist), data = CR)
model.WL <- lm(scale(value) ~ scale(fun.dist), data = WL)
model.WR <- lm(scale(value) ~ scale(fun.dist), data = WR)

model.CL <- glm(scale(fun.dist) ~ scale(value), data = CL)
model.CR <- glm(scale(value) ~ scale(fun.dist), data = CR)
model.WL <- glm(scale(value) ~ scale(fun.dist), data = WL)
model.WR <- glm(scale(value) ~ scale(fun.dist), data = WR)




summary(model.CL)
summary(model.CR)
summary(model.WL)
summary(model.WR)



#---- Supp Figure----------

lome.rel.all <- scale(bole[,11:124]) - min(scale(bole[,11:124]))

ddatb<-as.matrix(vegdist(lome.rel.all, "jaccard"))

bole$Plot_num <- as.factor(bole$Plot_num)

mod.bole <- dbrda(ddatb ~ Treatment * root_leaves + Condition(Sp:Sample_Site + Plot_num:Sample_Site), data = bole)

anova(mod.bole, by = "margin", type = 3)

plot(mod.bole)



scoresbole <- scores(mod.bole,display="sites")



bole$MDS1 <- scoresbole[,1]
bole$MDS2 <- scoresbole[,2]

quartz()
bole.ord <- ggplot(bole, aes(x = (MDS1), y = MDS2, col = Treatment, fill = Treatment)) +
  geom_point(size = 3, alpha = 0.6, aes(shape=root_leaves)) +
  ylab("")+ 
  xlab("")+ 
  scale_shape_manual(values=c(21, 24))+
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_fill_manual(values=c("darkseagreen", "chocolate1"))+
  mytheme+ theme(legend.position = "none")


RM.leaves.m <- subset(bole,root_leaves == "leaves" &Sample_Site == "RMBL")
RM.root.m <- subset(bole,root_leaves == "root" &Sample_Site == "RMBL")
KE.leaves.m <- subset(bole,root_leaves == "leaves" &Sample_Site == "Kessler")
KE.root.m <- subset(bole,root_leaves == "root" &Sample_Site == "Kessler")



RML.meta.m <- RM.leaves.m[,0:10]
RML.met.rel.t <- t(RM.leaves.m[,11:124])
RML.met<- t(subset(RML.met.rel.t, rowSums(RML.met.rel.t) > 0 & rowVars(RML.met.rel.t) >0))

RMR.meta.m <- RM.root.m[,0:10]
RMR.met.rel.t <- t(RM.root.m[,11:ncol(RM.root.m)])
RMR.met<- t(subset(RMR.met.rel.t, rowSums(RMR.met.rel.t) > 0 & rowVars(RMR.met.rel.t) >0))

KEL.meta.m <- KE.leaves.m[,0:10]
KEL.met.rel.t <- t(KE.leaves.m[,11:124])
KEL.met<- t(subset(KEL.met.rel.t, rowSums(KEL.met.rel.t) > 0 & rowVars(KEL.met.rel.t) >0))

KER.meta.m <- KE.root.m[,0:10]
KER.met.rel.t <- t(KE.root.m[,11:ncol(KE.root.m)])
KER.met<- t(subset(KER.met.rel.t, rowSums(KER.met.rel.t) > 0 & rowVars(KER.met.rel.t) >0))


RML.met.rel <- scale(RML.met) - min(scale(RML.met))
RMR.met.rel <- scale(RMR.met) - min(scale(RMR.met))

KEL.met.rel <- scale(KEL.met)  - min(scale(KEL.met))
KER.met.rel <- scale(KER.met) - min(scale(KER.met))



RML.ddatb<-as.matrix(vegdist(as.matrix(RML.met.rel), "jaccard"))
RMR.ddatb<-as.matrix(vegdist(as.matrix(RMR.met.rel), "jaccard"))

KEL.ddatb<-as.matrix(vegdist(as.matrix(KEL.met.rel), "jaccard"))
KER.ddatb<-as.matrix(vegdist(as.matrix(KER.met.rel), "jaccard"))




RML.model.m <- dbrda(RML.ddatb ~ Treatment*Sp + Condition(Plot), data = RML.meta.m)
RMR.model.m <- dbrda(RMR.ddatb ~ Treatment*Sp + Condition(Plot), data = RMR.meta.m)

KEL.model.m <- dbrda(KEL.ddatb ~ Treatment*Sp + Condition(Plot), data = KEL.meta.m)
KER.model.m <- dbrda(KER.ddatb ~ Treatment*Sp + Condition(Plot), data = KER.meta.m)


RML.model.m <- dbrda(RML.ddatb ~ Treatment+Sp + Condition(Plot), data = RML.meta.m)
RMR.model.m <- dbrda(RMR.ddatb ~ Treatment+Sp + Condition(Plot), data = RMR.meta.m)

KEL.model.m <- dbrda(KEL.ddatb ~ Treatment+Sp + Condition(Plot), data = KEL.meta.m)
KER.model.m <- dbrda(KER.ddatb ~ Treatment+Sp + Condition(Plot), data = KER.meta.m)




anova(RML.model.m, by = "margin", type =3)
anova(RMR.model.m, by = "margin", type =3)

anova(KEL.model.m, by = "margin", type =3)
anova(KER.model.m, by = "margin", type =3)

summary(RML.model.m)
summary(RMR.model.m)

summary(KEL.model.m)
summary(KER.model.m)



RsquareAdj(RML.model.m)
RsquareAdj(RMR.model.m)

RsquareAdj(KEL.model.m)
RsquareAdj(KER.model.m)

RML.scoresmet <- scores(RML.model.m,display="sites")
RMR.scoresmet <- scores(RMR.model.m,display="sites")

KEL.scoresmet <- scores(KEL.model.m,display="sites")
KER.scoresmet <- scores(KER.model.m,display="sites")

RML.meta.m$MDS1 <- RML.scoresmet[,1]
RML.meta.m$MDS2 <- RML.scoresmet[,2]

RMR.meta.m$MDS1 <- RMR.scoresmet[,1]
RMR.meta.m$MDS2 <- RMR.scoresmet[,2]


KEL.meta.m$MDS1 <- KEL.scoresmet[,1]
KEL.meta.m$MDS2 <- KEL.scoresmet[,2]

KER.meta.m$MDS1 <- KER.scoresmet[,1]
KER.meta.m$MDS2 <- KER.scoresmet[,2]


re.meta.all.m <- bind_rows(RML.meta.m, RMR.meta.m, KEL.meta.m, KER.meta.m)

quartz()


re.meta.all.m$Sample_Site <- ordered(re.meta.all.m$Sample_Site, levels = c("RMBL", "Kessler"))

metabo.species <- ggplot(re.meta.all.m, aes(x = MDS1, y = MDS2, col = Treatment)) +
  geom_point(size = 3, aes(shape=Sp)) +
  ylab("")+ 
  xlab("")+ 
  scale_color_manual(values=c("darkseagreen", "chocolate1"))+
  scale_shape_manual(values=c(0,5,6,7,8,9,10))+
  mytheme+
  facet_grid(root_leaves~Sample_Site)+
  theme(strip.text = element_blank(), legend.position = "none")



bole.all.fig <- ggarrange(bole.ord, metabo.species, ncol = 2)


