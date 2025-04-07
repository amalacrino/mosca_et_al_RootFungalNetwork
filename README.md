# Beneficial fungi are the major driver of root fungal microbiome assembly and network robustness

**Saveria Mosca, Edda Francomano, Meriem Miyassa Aci, Nesma Zakaria Mohamed, Leonardo Schena, Antonino Malacrinò**

## Abstract

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject number `PRJNA1080585`.

Our pipeline included:
* nf-core/ampliseq v2.7.1 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* MAFFT [https://academic.oup.com/nar/article/30/14/3059/2904316](https://academic.oup.com/nar/article/30/14/3059/2904316)
* FastTree [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)
* R v4.4.1 [https://www.R-project.org/](https://www.R-project.org/)

# Data processing

Run ampliseq

```bash
nextflow run nf-core/ampliseq -r 2.7.1 -profile singularity \
--input_folder $INDIR \
--FW_primer "CTTGGTCATTTAGAGGAAGTAA" \
--RV_primer "GCTGCGTTCTTCATCGATG" \
--outdir $OUTDIR \
--extension "/*_R{1,2}.fastq.gz" \
--illumina_pe_its \
--dada_ref_tax_custom $taxDB/sh_general_release_dynamic_04.04.2024.fasta \
--skip_dada_addspecies \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--ignore_empty_input_files \
--ignore_failed_trimming \
--ignore_failed_filtering \
--max_cpus 16 \
--max_memory '64.GB'

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```

Build a phyloseq object

```r
library("phyloseq")
library("Biostrings")

creatPSobj <- function(mdata, asvtab, tax, tree, refseqs){
  metadata <- read.table(mdata, sep = '\t', header = T, row.names = 1)
  rownames(metadata) <- gsub('[-]', '.', rownames(metadata))
  data <- read.table(asvtab, sep = '\t', header = T, row.names = 1)
  tax <- read.table(tax, sep = '\t', header = T, row.names = 1)[,c(1:7)]
  ref_tree <- ape::read.tree(tree) 
  refseq <- read.table(refseqs, sep = '\t', header = T)[,c(1,10)]
  seqs.vec <- refseq$sequence 
  names(seqs.vec) <- refseq$ASV_ID
  refseq <- DNAStringSet(seqs.vec)
  ps <- phyloseq(sample_data(metadata),
                 otu_table(data, taxa_are_rows = T),
                 tax_table(as.matrix(tax)),
                 phy_tree(ref_tree),
                 refseq)
  return(ps)
}

ps <- creatPSobj(mdata = "metadata.txt",
                     asvtab = "ASV_table.tsv",
                     tax = "ASV_tax.user.tsv",
                     tree = 'tree.tre',
                     refseqs = "ASV_tax.user.tsv")

rm(list=setdiff(ls(), c("ps")))

save.image(file = 'data.rds')
```

## Data analysis

### Load libraries

```r
library("tidyverse")
library("phyloseq")
library("DESeq2")
library("ggrepel")
library("emmeans")
library("car")
library("lme4")
library("microbiome")
library("ggvenn")
library("decontam")
library("Wrench")
library("RVAideMemoire")
library("picante")
library("decontam")
library("reshape2")
library("ggvenn")
library("RColorBrewer")
library("data.table")
library("psych")
library("microeco")
library("meconetcomp")
library("magrittr")
library("igraph")
library("file2meco")
library("tyRa")
library("minpack.lm")
library("Hmisc")
library("iCAMP")
```

### Load and clean data

```r
load(file = 'data/data.rds')

remove.cont <- function(ps){
  sample_data(ps)$is.neg <- sample_data(ps)$Treatment == "DNA_extraction"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.05)
  cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
  cont.remove <- row.names(cont.remove)
  allTaxa <- taxa_names(ps)
  allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
  ps <- prune_taxa(allTaxa, ps)
  return(ps)
}

ps <- remove.cont(ps)
ps <- phyloseq::subset_samples(ps, Compartment == "Roots")
ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps <- prune_samples(sample_sums(ps) >= 1000, ps)

wrench.norm.ps <- function(ps){
  count_tab <- as.matrix(data.frame(otu_table(ps)))
  group <- paste0(sample_data(ps)$Compartment)
  W <- wrench(count_tab, condition=group)
  norm_factors <- W$nf
  norm_counts <- sweep(count_tab, 2, norm_factors, FUN = '/')
  norm_counts_trim <- data.frame(t(data.frame(norm_counts)))
  #norm_counts_trim[] <- lapply(norm_counts_trim, function(x) DescTools::Winsorize(x, probs = c(0, 0.97), type = 1))
  norm_counts_trim <- data.frame(t(norm_counts_trim))
  norm_counts_trim[norm_counts_trim == 0] <- 1
  norm_counts_trim[norm_counts_trim < 0] <- 1
  norm_counts_trim[is.na(norm_counts_trim)] <- 1
  norm_counts_trim <- log2(norm_counts_trim)
  colnames(norm_counts_trim) <- gsub("\\.", "-", colnames(norm_counts_trim))
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(norm_counts_trim, taxa_are_rows =  TRUE)
  return(ps_norm)
}

ps_n <- wrench.norm.ps(ps)
```

### Multivariate analysis

```r
sampledf <- data.frame(sample_data(ps))
dist.mat <- phyloseq::distance(ps, method = "wunifrac")
dist.mat[dist.mat<0] <- 0
perm <- how(nperm = 999,
            blocks = sampledf$Block)
set.seed(100)
pmv <- adonis2(dist.mat ~ Soil, data = sampledf, permutations = perm)
pmv

pairwise.perm.manova(dist.mat, paste0(sampledf$Soil), nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)

cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = dist.mat, formula = ~ 1)
cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = Soil),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = Soil, shape = Soil), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = "none", legend.justification = c(0.99, 0.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"), labels = c("Agricultural", "Margin", "Uncultivated"), breaks = c("Agricultural", "Margin", "Uncultivated")) +
  scale_shape_manual(name = "Legend", values=c(19,17,15), labels = c("Agricultural", "Margin", "Uncultivated"), breaks = c("Agricultural", "Margin", "Uncultivated")) +
  xlim(-0.4, 0.4) + ylim(-0.3, 0.3)
```

### Diversity

```r
div <- microbiome::alpha(ps)
otus <- as.data.frame(t(otu_table(ps)))

tree <- phy_tree(ps)
div.pd <- pd(otus, tree, include.root = FALSE)

mntd.calc <- function(ps){
  comm <- as.data.frame(t(otu_table(ps)))
  phy <- phy_tree(ps)
  phy.dist <- cophenetic(phy)
  comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
  ks <- sample_data(ps)
  bnti <- cbind(ks, comm.sesmpd)
  return(bnti)
}

mntd.df <- mntd.calc(ps)

div.2 <- cbind(sample_data(ps), div)
div.2 <- cbind(div.2, div.pd, mntd.df$mntd.obs, -mntd.df$mntd.obs.z)
colnames(div.2)[c(32,33)] <- c("mntd.obs", "mntd.obs.z")

model <- lmer(diversity_shannon ~ Soil * (1|Block), data = div.2)
Anova(model)

model <- lmer(PD ~ Soil * (1|Block), data = div.2)
Anova(model)

model <- lmer(mntd.obs ~ Soil * (1|Block), data = div.2)
Anova(model)
ph <- emmeans(model, "Soil")
pairs(ph)

model <- lmer(mntd.obs.z ~ Soil * (1|Block), data = div.2)
Anova(model)
ph <- emmeans(model, "Soil")
pairs(ph)

div.plot.fun <- function(df, var, label, y1, y2){
  plot <- ggplot(div.2, aes(x = Soil, y = get(var), fill = Soil)) +
    theme_bw(base_size = 15) +
    geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
    labs(y = label) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_text(color="black"),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = "none") +
    scale_fill_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"),
                      labels = c("Agricultural", "Margin", "Uncultivated"),
                      breaks = c("Agricultural", "Margin", "Uncultivated")) +
    ylim(y1, y2)
  return(plot)
}

plot1 <- div.plot.fun(div.2, "diversity_shannon", "Shannon diversity", 1.5, 4)
plot2 <- div.plot.fun(div.2, "PD", "Phylogenetic diversity", 0, 70)
plot3 <- div.plot.fun(div.2, "mntd.obs", "MNTD", 0, 1)
plot4 <- div.plot.fun(div.2, "mntd.obs.z", "beta-NTI", -2, 2)
```

### Taxa plot (relative abundance)

```r
glom <- microbiome::aggregate_taxa(ps_n, "Genus")
glom <- microbiome::transform(glom, "compositional")
dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarize(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(Genus %in% filt.gen$Genus))
dat <- dat %>% group_by(Soil, Genus) %>% summarize(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) %>% mutate(Genus = replace(Genus, Genus == "", "unidentified taxa"))
nb.cols <- length(unique(dat$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

taxa_plot <- ggplot(dat, aes(x = as.factor(Soil), y = cs, fill = Genus)) +
  theme_bw(base_size = 14) +
  geom_bar(stat="identity") +
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = mycolors) +
  labs(y = "relative abundance", x="")
```

### Random forest

```r
predictors <- t(otu_table(ps))
dim(predictors)
response <- as.factor(sample_data(ps)$Soil)
rf.data <- data.frame(response, predictors)

set.seed(100)
classify <- randomForest(response~., data = rf.data, ntree = 100)
print(classify)

imp <- importance(classify)
imp <- data.frame(predictors = rownames(imp), imp)
taxt <- tax_table(ps) %>% as.data.frame() %>% tibble::rownames_to_column("ASV") %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom))
imp <- merge(imp, taxt, by.x = "predictors", by.y = "ASV")
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
mlp <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  theme_classic() +
  geom_point(size = 3) +
  scale_x_discrete(labels= imp.20$identification) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.5,
                                          linetype = 2)) +
  xlab(" ") + ylab("Mean Decrease in Gini coefficient") +
  coord_flip()
```

### Network analysis

```r
soil_amp <- phyloseq2meco(ps)

soil_amp_network <- list()
tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Soil == "Agricultural")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$Agricultural <- tmp

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Soil == "Margin")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$Margin <- tmp

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Soil == "Uncultivated")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$Uncultivated <- tmp

soil_amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")

tmp <- cal_network_attr(soil_amp_network)

soil_amp_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

tmp <- robustness$new(soil_amp_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"),
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
tmp$plot(linewidth = 1)

plot.net.fun <- function(rm.strat, measure.var, ylab, xlab, lmTF){
  tmp1 <- tmp$res_table %>% .[.$remove_strategy == rm.strat & .$measure == measure.var, ]
  t1 <- trans_env$new(dataset = NULL, add_data = tmp1)
  t1$dataset$sample_table <- t1$data_env
  net.met.df <- t1$dataset$sample_table
  p <- ggplot(net.met.df, aes(x=remove_ratio, y=value, color = Network)) + 
    geom_point(size = 3, alpha = 0.5) +
    {if(lmTF == T)geom_smooth(method=lm)}+
    {if(lmTF == F)geom_path()} +
    theme_classic() +
    theme(legend.title=element_blank()) +
    xlab(xlab) + 
    ylab(ylab) +
    scale_color_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"),
                      labels = c("Agricultural", "Margin", "Uncultivated"),
                      breaks = c("Agricultural", "Margin", "Uncultivated")) 
  return(p)
}

p1 <- plot.net.fun("edge_rand", "Eigen", "Network connectivity", "Ratio of randomly removed edges", T)
p2 <- plot.net.fun("edge_strong", "Eigen", "Network connectivity", "Ratio of edges removed in decreasing weight", F)
p3 <- plot.net.fun("node_rand", "Eigen", "Network connectivity", "Ratio of randomly removed nodes", T)
p4 <- plot.net.fun("node_degree_high", "Eigen", "Network connectivity", "Ratio of nodes removed in decreasing degree", F)

hdn_agr <- soil_amp_network$Agricultural$res_node_table %>% as_tibble() %>%
  slice_max(order_by = degree, n = round(nrow(.) * 0.25)) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Agricultural")

hdn_mar <- soil_amp_network$Margin$res_node_table %>% as_tibble() %>%
  slice_max(order_by = degree, n = round(nrow(.) * 0.25)) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Margin")

hdn_unc <- soil_amp_network$Uncultivated$res_node_table %>% as_tibble() %>%
  slice_max(order_by = degree, n = round(nrow(.) * 0.25)) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Uncultivated")

hdn.df <- rbind(hdn_agr, hdn_mar, hdn_unc)

nb.cols <- length(unique(hdn.df$identification))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

pxx <- ggplot(data=hdn.df, 
             aes(x = factor(1), fill = identification)) +
        facet_wrap(vars(soil), nrow = 3, scale='free_y') + 
        theme_void() +
        geom_bar(stat = "count") +
        coord_polar(theta='y') +
        scale_fill_manual(values = mycolors) +
        theme(axis.text.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              strip.text = element_blank(),
              legend.position='right',
              legend.title=element_blank())
```

### Network visualization

```r
image.net.fun <- function(group, lab){
  soil_amp_network <- list()
  tmp <- clone(soil_amp)
  tmp$sample_table %<>% subset(Soil == group)
  tmp$tidy_dataset()
  
  
  taxt <- tax_table(ps) %>% as.data.frame() %>%
    mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
    mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom))
  
  totu <- t(tmp$otu_table)
  c_net_calculate(totu, method = "spearman") -> corr
  c_net_build(corr, r_thres = 0.6, p_thres = 0.05, delete_single = T) -> co_net
  co_net1 <- c_net_set(co_net, taxt["identification"])
  c_net_layout(co_net1, method = nicely()) -> coors
  pdf(paste0("figures/", lab, "_network.pdf"))
  c_net_plot(co_net1, coors, labels_num = 0, legend = F, edge.color = "black",
             vertex.size = 5, vertex.color = "lightgray", main = " ", seed = 100)
  dev.off()
}

image.net.fun("Agricultural", "agr")
image.net.fun("Margin", "mar")
image.net.fun("Uncultivated", "unc")
```

### LEfSe

```r
soil_amp <- phyloseq2meco(ps)
t1 <- trans_diff$new(dataset = soil_amp, method = "lefse", group = "Soil", alpha = 0.001, lefse_subgroup = NULL, taxa_level = "ASV")
lefse.df <- t1$res_diff

lefse.df <- lefse.df %>%
  separate(Taxa, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), 
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom))

lfsgr <- ggplot(lefse.df, aes(x = fct_reorder(ASV, Group, .desc = TRUE), y = LDA, fill = Group)) +
  theme_classic() +
  geom_bar(stat = "identity") +
  
  scale_x_discrete(labels= lefse.df$identification) +
  xlab(" ") + ylab("LDA score") +
  scale_fill_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"),
                     labels = c("Agricultural", "Margin", "Uncultivated"),
                     breaks = c("Agricultural", "Margin", "Uncultivated")) +
  coord_flip()
```

### Sloan null model

```r
ks.tmp <- sample_data(ps) %>% as_tibble()
sample_data(ps)$group <- paste0(ks.tmp$Soil,"-",ks.tmp$Block)
list.gm <- unique(sample_data(ps)$group)

sncm.calc <- function(x){
  ks <- sample_data(ps)[["Soil"]] %in% x
  ps <- prune_samples(samples = ks, ps)
  spp.out <- tyRa::fit_sncm(spp = t(otu_table(ps))@.Data, pool=NULL, taxon=data.frame(tax_table(ps)))
  return(spp.out)}

sncm.results.agr <- sncm.calc("Agricultural")
plot_sncm.agr <- plot_sncm_fit(spp_out = sncm.results.agr) +
  theme_bw(base_size = 15) +
  theme(legend.justification=c(0.01,0.99), legend.position=c(0.01,0.99), legend.title=element_blank())
plot_sncm.agr$layers[[6]] <- NULL
plot_sncm.agr$layers[[5]] <- NULL
sncm.results.agr <- sncm.results.agr$predictions %>% 
  filter(fit_class != "As predicted") %>%
  tibble::rownames_to_column("ASV")
sncm.results.agr %>% filter(ASV %in% hdn_agr$name)

sncm.results.mar <- sncm.calc("Margin")
plot_sncm.mar <- plot_sncm_fit(spp_out = sncm.results.mar) +
  theme_bw(base_size = 15) +
  theme(legend.justification=c(0.01,0.99), legend.position=c(0.01,0.99), legend.title=element_blank())
plot_sncm.mar$layers[[6]] <- NULL
plot_sncm.mar$layers[[5]] <- NULL
sncm.results.mar <- sncm.results.mar$predictions %>% 
  filter(fit_class != "As predicted") %>%
  tibble::rownames_to_column("ASV")
sncm.results.mar %>% filter(ASV %in% hdn_mar$name)

sncm.results.unc <- sncm.calc("Uncultivated")
plot_sncm.unc <- plot_sncm_fit(spp_out = sncm.results.unc) +
  theme_bw(base_size = 15) +
  theme(legend.justification=c(0.01,0.99), legend.position=c(0.01,0.99), legend.title=element_blank())
plot_sncm.unc$layers[[6]] <- NULL
plot_sncm.unc$layers[[5]] <- NULL
sncm.results.unc <- sncm.results.unc$predictions %>% 
  filter(fit_class != "As predicted") %>%
  tibble::rownames_to_column("ASV")
sncm.results.unc %>% filter(ASV %in% hdn_unc$name)

ggvenn(list(agr = sncm.results.agr$ASV,
            mar = sncm.results.mar$ASV,
            unc = sncm.results.unc$ASV))

sncm.results.unc %>% filter(ASV %in% sncm.results.agr$ASV) %>%
  filter(ASV %in% sncm.results.mar$ASV)

sncm.results.agr %>% filter(!ASV %in% sncm.results.unc$ASV) %>%
  filter(!ASV %in% sncm.results.mar$ASV) %>% filter(ASV %in% hdn_agr$name)

sncm.results.mar %>% filter(!ASV %in% sncm.results.unc$ASV) %>%
  filter(!ASV %in% sncm.results.agr$ASV) %>% filter(ASV %in% hdn_mar$name)

sncm.results.unc %>% filter(!ASV %in% sncm.results.agr$ASV) %>%
  filter(!ASV %in% sncm.results.mar$ASV)  %>% filter(ASV %in% hdn_unc$name)

ggvenn(list(agr_sncm = sncm.results.agr$ASV,
            agr_hdn = hdn_agr$name))

ggvenn(list(mar_sncm = sncm.results.mar$ASV,
            mar_hdn = hdn_mar$name))

ggvenn(list(unc_sncm = sncm.results.unc$ASV,
            unc_hdn = hdn_unc$name))

com.a <- sncm.results.agr %>% filter(ASV %in% hdn_agr$name) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Agricultural")

com.b <- sncm.results.mar %>% filter(ASV %in% hdn_mar$name) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Margin")

com.c <- sncm.results.unc %>% filter(ASV %in% hdn_unc$name) %>%
  mutate(across(everything(), ~ ifelse(. %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), "", .))) %>%
  mutate(across(everything(), ~ ifelse(str_trim(.) == "", NA, .))) %>%
  mutate(identification = coalesce(Genus, Family, Order, Class, Phylum, Kingdom)) %>%
  mutate(soil = "Uncultivated")

hdn.df <- rbind(com.a, com.b, com.c)

nb.cols <- length(unique(hdn.df$identification))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

pxx <- ggplot(data=hdn.df, 
              aes(x = factor(1), fill = identification)) +
  facet_wrap(vars(soil), scale='free_y') + 
  theme_void() +
  geom_bar(stat = "count") +
  coord_polar(theta='y') +
  scale_fill_manual(values = mycolors) +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.title=element_blank())
```

### iCAMP

```r
preabs <- as.data.frame(t(otu_table(ps)))
phy <- phy_tree(ps)
sss <- icamp.big(preabs, phy, nworker = 8)
df <- sss$detail$processes$CbMPDiCBraya %>%
  select(-sample2) %>%
  group_by(sample1) %>%
  summarise(across(everything(), mean))

ks.tmp <- as.matrix(sample_data(ps)) %>% as.data.frame()
ks.tmp <- ks.tmp %>% rownames_to_column("sample")
df.2 <- merge(df, ks.tmp, by.x = "sample1", by.y = "sample") %>%
  group_by(Soil) %>%
  summarise(across(everything(), mean)) %>%
  select(1,3:7) %>%
  reshape2::melt()

taxa_plot <- ggplot(df.2, aes(x = as.factor(Soil), y = value, fill = variable)) +
  theme_bw(base_size = 22) +
  geom_bar(stat="identity") +
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="transparent"),
        legend.title=element_blank(), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "relative contribution", x="") +
  coord_flip() +
  scale_fill_manual(name = "Legend", values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
                    breaks = c("Heterogeneous.Selection", "Homogeneous.Selection", "Dispersal.Limitation", "Homogenizing.Dispersal", "Drift.and.Others"),
                    labels = c("Heterog. selection", "Homog. selection", "Dispersal limitation", "Homog. dispersal", "Drift + others"))
```
