library(AnnotationForge)
library(dplyr)
library(biomaRt)
listEnsembl()
ensembl = useEnsembl(biomart="ensembl", version = 112, dataset="pabelii_gene_ensembl")

ponAbe2 <- getBM(
    attributes = c('entrezgene_id', 'ensembl_gene_id', 'external_gene_name', 'description',
                   'chromosome_name', 'go_id', "go_linkage_type"),
    mart = ensembl
)

pGene <- ponAbe2[,c('entrezgene_id', 'ensembl_gene_id', 'external_gene_name', 'description')]
pChr <- ponAbe2[,c('entrezgene_id', 'chromosome_name')]
pGO <- ponAbe2[,c('entrezgene_id', 'go_id', 'go_linkage_type')]

colnames(pGene) <- c("GID", "ENSEMBL", "SYMBOL", "GENENAME")
colnames(pChr) <- c("GID", "CHROMOSOME")
colnames(pGO) <- c("GID", "GO", "EVIDENCE")

# remove duplicate rows
pGene_unique <- pGene %>%
    distinct(ENSEMBL, .keep_all = TRUE) %>%
    filter(!is.na(GID))

pChr_unique <- pChr %>% filter(!is.na(GID)) %>% distinct()

pGO_unique <- pGO %>%
    distinct() %>%
    filter(GO != "", EVIDENCE != "") %>%
    filter(!is.na(GID))

makeOrgPackage(gene_info=pGene_unique,
               chromosome=pChr_unique,
               go=pGO_unique,
               version="0.1",
               maintainer = "Yao <yao-chung.chen@fu-berlin.de>",
               author = "Yao <yao-chung.chen@fu-berlin.de>",
               outputDir = '.',
               tax_id = '9601',
               genus = 'Pongo',
               species = "abelii",
               goTable = "go")

#R CMD INSTALL ./org.Pabelii.eg.db
