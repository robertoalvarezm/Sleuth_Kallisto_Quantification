library(sleuth)
library("gridExtra")
library("cowplot")
library("biomaRt")



sample_id <- dir(file.path("01_Raw_Data/kallisto_quant_output/"))
sample_id


kal_dirs <- file.path("01_Raw_Data/kallisto_quant_output/", sample_id)


s2c <- read.table(file.path("01_Raw_Data/kallisto_demo.tsv"),
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t")
s2c


s2c <- dplyr::mutate(s2c, path = kal_dirs)


s2c

library(DT)
DT::datatable(s2c)


marts <- listMarts()
marts


marts <- listMarts(host = "plants.ensembl.org")
marts


plants_mart <- useMart("plants_mart", host = "plants.ensembl.org" )
listDatasets(plants_mart)

plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )

datatable(listAttributes(plants_mart))


plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )
t2g <- getBM(attributes = c("ensembl_transcript_id",
                            "ensembl_gene_id",
                            "description",
                            "external_gene_name"),
             mart = plants_mart)
t2g

ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


s2c$genotype_variation_s <- as.factor(s2c$genotype_variation_s)
s2c$genotype_variation_s <- relevel(s2c$genotype_variation_s, ref = "wild type")


library(cowplot)

No corras esto si tu compu es muy lenta
so <- sleuth_prep(s2c,
full_model = ~genotype_variation_s,
           target_mapping = ttg,
          read_bootstrap_tpm=TRUE,
         extra_bootstrap_summary = TRUE)
#####
saveRDS(so,file="03_Results/so.RDS")
## Mejor corre esta
readRDS(file="03_Results/so.RDS")
so<-readRDS(file="03_Results/so.RDS")
ggplot2::theme_set(theme_cowplot())
plot_pca(so, color_by = 'genotype_variation_s', text_labels = TRUE)


ggplot2::theme_set(theme_cowplot())
g<-plot_pca(so, color_by = 'genotype_variation_s', text_labels = TRUE)
library(plotly)
ggplotly(g)


plot_pca(so, color_by = 'treatment_s', text_labels = TRUE)


plot_loadings(so, pc_input = 1)

plot_bootstrap(so, 'AT2G34420.1', color_by = 'genotype_variation_s')

plot_bootstrap(so, 'AT2G34420.1', color_by = 'treatment_s')

so2 <- sleuth_fit(so, ~genotype_variation_s, 'full')
so2 <- sleuth_fit(so, ~1, 'reduced')
so2 <- sleuth_lrt(so, 'reduced', 'full')
saveRDS(so2,file="../03_Results/so2.RDS")
so2<-readRDS(file="03_Results/so2.RDS")



ull_results <- sleuth_results(so2, 'reduced:full', 'lrt',
                              show_all = FALSE)
library(DT)
datatable(head(full_results))


wald_test <- colnames(design_matrix(so2))[2]
so <- sleuth_wt(so2, wald_test)

sleuth_live(so)