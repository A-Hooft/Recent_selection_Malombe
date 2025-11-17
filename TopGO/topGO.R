#TOP GO

# Install Bioconductor
if(!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")}

# Install topGO package
BiocManager::install("topGO")

# Install package for graph visualization
BiocManager::install("Rgraphviz")

# ape package allows to read annotation files
install.packages("ape")

BiocManager::install("org.Dr.eg.db")    # Zebrafish annotation package from Bioconductor

# Load packages
library(topGO)
library(ape)
library(Rgraphviz)
library(org.Dr.eg.db)


# Function that creates GO data object for further analyses
make_GO_object <- function(description,AllGenes_path, GenesOfInterest_path,
                           ontology, nodeSize) {
  
  # Create gene universe
  all_genes <- read.csv(AllGenes_path, header=FALSE)
  gene_universe <- as.vector(all_genes$V1)
  gene_universe <- gene_universe[!duplicated(gene_universe)]   # remove duplicates
  
  # Read in list of interesting genes to test
  test_genes <- as.vector(read.csv(GenesOfInterest_path,
                                   sep = ',',header=FALSE)[[1]])
  #test_genes <- as.vector(test_genes[[1]])
  
  # Classify genes into test ('1') and no test ('0'). This will be the input data for GOdata object.
  geneClassify <- factor(as.integer(gene_universe %in% test_genes))
  names(geneClassify) <- gene_universe
  str(geneClassify)
  
  # Build GOdata object    
  GOdata <- new("topGOdata", description = description, ontology = ontology,
                allGenes = geneClassify,
                nodeSize = nodeSize,
                annot = annFUN.org, mapping = "org.Dr.eg.db", ID = "symbol")
}

# Path to files

## Here:
### all genes = all genes annotated by snpEff in my GWAS dataset
setwd("~/Documents/Cichlid_Fish")
all_genes_AstCal <- "all_genes_ensemble.txt"
all_genes_zebra <- "zeb_genes_ensembl.txt"

test_genes <- "gene_list_candidate_473_clean.txt"

# Create GO object using function
GO1 <- make_GO_object(description="SNPs passing threshold",
                      all_genes_AstCal, test_genes, ontology="BP", nodeSize=5)

resultfisher_weight01 <- runTest(GO1,
                                 algorithm = "weight01", 
                                 statistic = "fisher", scoreOrder='decreasing')

n_lt0p01 <- sum(score(resultfisher_weight01)<0.05)
sign_table_fisher <- GenTable(GO1, resultfisher_weight01, orderBy = "resultfisher_weight01", topNodes =n_lt0p01)#


write.table(sign_table_fisher, file = "topGO_fisher_473_windows_11_2025.csv", sep = ",")


