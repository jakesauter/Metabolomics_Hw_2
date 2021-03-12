## **Metabolomics Enrichment Analysis**

A report by: **Jake Sauter** and **Nicholas Bartelo**

date: 3/4/2021


## **Problem 1. Enrichment Analysis**

### **1.1 Load Gene Expression Data**


```r
library(tidyverse)
library(magrittr)
library(knitr)
```


```r
data_env <- new.env()
load('data/geod33272.rda', envir = data_env)
ls(data_env)
```

```
[1] "geod33272" "kegg"     
```

`geod33272$x` represents the gene expression data (genes in rows, samples in columns). 


```r
data_env$geod33272$x %>% 
  .[1:5, 1:3] %>% 
  kable()
```



|                | GSM823316_sample_table.txt| GSM823317_sample_table.txt| GSM823318_sample_table.txt|
|:---------------|--------------------------:|--------------------------:|--------------------------:|
|GE_BrightCorner |               18971.290000|               69310.900000|               1.147218e+05|
|DarkCorner      |                   8.751004|                   2.922934|               2.940515e+00|
|DarkCorner      |                  12.294990|                   2.935335|               2.952986e+00|
|DarkCorner      |                   8.807324|                   3.928323|               2.965929e+00|
|DarkCorner      |                   5.561080|                   2.958856|               2.979015e+00|

`geod33272$status` indicates the cell activation status per sample (0 = naive, 1 = activated). 


```r
data_env$geod33272$status
```

```
 [1] 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
```

`geod33272$ID` contains the gene identifiers for all transcripts. 


```r
data_env$geod33272$ID[1:20,] %>% 
  kable()
```



|agilentID       |entrezID |
|:---------------|:--------|
|GE_BrightCorner |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|DarkCorner      |NA       |
|A_24_P66027     |9582     |
|A_32_P77178     |NA       |
|A_23_P212522    |23200    |
|A_24_P934473    |NA       |
|A_24_P9671      |3301     |
|A_32_P29551     |NA       |
|A_24_P801451    |10919    |
|A_32_P30710     |9349     |
|A_32_P89523     |NA       |

We will use `geod33272$ID$entrezID` as gene identifiers.
This mapping is faciliated below. Note that the first few
probes shown are control probes, and thus have a value 
of `NA` for their entrezID. We now perform the entrezID mapping, 
and only keep data from probes that map to an entrezID.


```r
na_entrez_ids <- is.na(data_env$geod33272$ID$entrezID)
  
genes <-
  data_env$geod33272 %>% 
  {.$x[!na_entrez_ids,]} %>% 
  set_colnames(gsub('_sample_table.txt', '', colnames(.))) %>% 
  log2()

# Map names to entrezID
rownames(genes) <-
  data_env$geod33272$ID$entrezID[!na_entrez_ids]
```


```r
genes[1:5, 1:5] %>% 
  kable()
```



|      | GSM823316| GSM823317| GSM823318| GSM823319| GSM823320|
|:-----|---------:|---------:|---------:|---------:|---------:|
|9582  |  8.138524|  7.789053|  8.675485|  8.909267|  7.962474|
|23200 |  9.385482|  9.966857|  8.962000|  9.319318|  8.270992|
|3301  | 12.646104| 12.578232| 12.252538| 13.433171| 13.670922|
|10919 |  7.503905|  8.673080|  7.720187|  8.361439|  7.857598|
|9349  | 15.848321| 15.458346| 15.616941| 15.420573| 15.777043|


### **1.2 Inspecting KEGG Data**

**The object ‘kegg’ contains KEGG pathway annotations. How many pathways does it include? How many genes are annotated with “hsa00020 Citrate cycle (TCA cycle)”?**


```r
type <- typeof(data_env$kegg)
length <- length(data_env$kegg)

cat("kegg variable is type: ", type, " and contains: ", length, " elements")
```

```
kegg variable is type:  list  and contains:  177  elements
```

It appears that each of these elements is a pathway


```r
data_env$kegg %>% 
  names() %>% 
  head()
```

```
[1] "hsa00010 Glycolysis / Gluconeogenesis"            
[2] "hsa00020 Citrate cycle (TCA cycle)"               
[3] "hsa00030 Pentose phosphate pathway"               
[4] "hsa00040 Pentose and glucuronate interconversions"
[5] "hsa00051 Fructose and mannose metabolism"         
[6] "hsa00052 Galactose metabolism"                    
```

It also appears that each "pathway" stores the ID of genes that 
have been determined to be included in the pathway. 


```r
data_env$kegg[[1]][1:10]
```

```
 [1] "10327"  "124"    "125"    "126"    "127"    "128"    "130"    "130589"
 [9] "131"    "160287"
```

Specifically, from the code below, we can determine that the 
"hsa00020 Citrate cycle (TCA cycle)" has **30** annotated genes.


```r
data_env$kegg %>% 
  .[["hsa00020 Citrate cycle (TCA cycle)"]] %>% 
  length()
```

```
[1] 30
```


## **1.3 Preparing a Genes x Pathways Matrix**

**Prepare a genes X pathways matrix: This matrix should contain a 1 if a gene is member of a pathway and 0 otherwise. As mentioned above, use Entrez Gene IDs for mapping.**



```r
pathways <- 
  data_env$kegg %>% 
  names()

genes_by_pathway <- 
  matrix(nrow = nrow(genes), 
         ncol = length(pathways), 
         dimnames = list(
           rownames(genes),
           pathways
         ), 0)

gene_names <- rownames(genes_by_pathway)
for (pathway in pathways) {
  for (pathway_gene in data_env$kegg[[pathway]]) {
    hits <- which(gene_names == pathway_gene)
    genes_by_pathway[hits, pathway] <- 1
  }
}

# More concise but harder to read / understand
# genes_by_pathway <-
#   sapply(kegg, function(x) {gene.ids %in% x}) %>%
#   matrix(ncol = length(kegg))
```


```r
genes_by_pathway[50:60, 5:10] %>% 
  kable()
```



|       | hsa00051 Fructose and mannose metabolism| hsa00052 Galactose metabolism| hsa00053 Ascorbate and aldarate metabolism| hsa00061 Fatty acid biosynthesis| hsa00071 Fatty acid metabolism| hsa00072 Synthesis and degradation of ketone bodies|
|:------|----------------------------------------:|-----------------------------:|------------------------------------------:|--------------------------------:|------------------------------:|---------------------------------------------------:|
|4901   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|8019   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|285521 |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|7253   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|26529  |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|26020  |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|5639   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|6547   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|2706   |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|339779 |                                        0|                             0|                                          0|                                0|                              0|                                                   0|
|94103  |                                        0|                             0|                                          0|                                0|                              0|                                                   0|

**Provide log2-transformed gene expression data (parameter X.mat), the activation status (y.vec), the pathway assignment matrix (C.mat), and set Pi.mat to an appropriate value of your choice. Why did you set Pi.mat to the value you did? Display the top 10 results (check the package documentation).**

**Pi.mat** -- Either an integer, or a matrix or data.frame containing the permutations. See getPImatrix for the acceptable form of a matrix or data.frame. If Pi.mat is an integer, B, then safe will generate B resamples of X.mat.

By mistake we have originally called `safe::safe` with too large of a 
**Pi.mat** value, leading to the warning: 

>  Warning: only 12870 unique resamples exist
          switching to exhaustive permutation
          
which was very interesting. In order to determine how this number was 
calculated, we were able to probe the source code of `safe::safe` to 
see that it was done the following way.


```r
y.vec <- 
  activation_status <- 
    data_env$geod33272$status

n <- ncol(genes)
count <- choose(n, table(y.vec)[1])

cat(count, ' unique permutations exists')
```

```
12870  unique permutations exists
```

We will now use this number of unique permutations as our 
**Pi.mat** value as this number of permutations is computationally 
tractable and will provide us with the maximum statistical power.



```r
library(safe)
library(doParallel)
registerDoParallel(cores=4)
set.seed(497)

safe_res <-
  safe::safe(X.mat = genes,
             y.vec = activation_status,
             C.mat = genes_by_pathway,
             Pi.mat = count,
             parallel = TRUE)

saveRDS(safe_res, 'data/safe_res.rds')
```


Now that we have successfully run `safe`, our top 10 results
are shown below.


```r
library(safe)
safe_res <- readRDS('data/safe_res.rds')
safe::safe.toptable(safe_res, number = 10) %>%
  kable()
```



|GenesetID                                           | Size| Statistic|P.value |Adj.p.value |Description |
|:---------------------------------------------------|----:|---------:|:-------|:-----------|:-----------|
|hsa03440 Homologous recombination                   |   71|   1944898|0.0002  |0.0069      |NA          |
|hsa03030 DNA replication                            |   68|   1894648|0.0002  |0.0069      |NA          |
|hsa03410 Base excision repair                       |   72|   1865724|0.0002  |0.0069      |NA          |
|hsa03430 Mismatch repair                            |   62|   1742655|0.0002  |0.0069      |NA          |
|hsa03420 Nucleotide excision repair                 |   68|   1725735|0.0002  |0.0069      |NA          |
|hsa03450 Non-homologous end-joining                 |   32|    828469|0.0002  |0.0069      |NA          |
|hsa00280 Valine, leucine and isoleucine degradation |   68|   1707931|0.0009  |0.0223      |NA          |
|hsa00190 Oxidative phosphorylation                  |  179|   4196949|0.0010  |0.0223      |NA          |
|hsa00630 Glyoxylate and dicarboxylate metabolism    |   29|    722971|0.0012  |0.0244      |NA          |
|hsa00310 Lysine degradation                         |   80|   1682239|0.0018  |0.0313      |NA          |

### **1.5 Significant Pathway Visualization**

Finally, our `safeplot`:


```r
safeplot(safe_res)
```

```
Shaded limits are the largest local statistics with p > 0.05
   Limits = ( -2.206 , 2.327 )
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


## **Problem 2: Aggregation Based Analysis**

### **2.1 Loading Pre-processed Data**

**Load the attached metabolomics data (Simulated_metabolomics_data_preprocessed.Rdata). The dataset has already been preprocessed by quotient normalization and logged. Metabolites are in columns(!) and samples are in rows of the dat data frame. Annotations of the metabolites can be found in annotations and information on gender is stored in variable gender, with 1 = male and 2 = female.**


```r
data_env <- new.env()
load('./data/Simulated_metabolomics_data_preprocessed.RData', envir = data_env)
ls(data_env)
```

```
[1] "annotations" "dat"         "gender"     
```


### **2.2 Calculate Eigenmetabolite**

Goal: Calculate the eigenmetabolite as a pathway representative for each “Sub_pathway”. 

First we must scale the data (mean of 0, sd of 1). 


```r
df <- data_env$dat %>% 
  as.numeric() %>% 
  {(. - mean(.)) / sd(.)} %>%
  matrix(nrow = nrow(data_env$dat), 
         ncol = ncol(data_env$dat), 
         data = .)

print(data_env$dat[1:3,1:3]) %>% kable()
```

```
         M20675    M34400     M33228
[1,] -0.2175067 0.2667876 -0.2435822
[2,]  0.2669905 0.3000532  0.5189495
[3,]  0.4868916 0.5831026  0.6452525
```



|     M20675|    M34400|     M33228|
|----------:|---------:|----------:|
| -0.2175067| 0.2667876| -0.2435822|
|  0.2669905| 0.3000532|  0.5189495|
|  0.4868916| 0.5831026|  0.6452525|

```r
print(df[1:3,1:3]) %>% kable()
```

```
           [,1]      [,2]       [,3]
[1,] -0.2896416 0.2060601 -0.3163313
[2,]  0.2062678 0.2401092  0.4641616
[3,]  0.4313486 0.5298258  0.5934396
```



|           |          |           |
|----------:|---------:|----------:|
| -0.2896416| 0.2060601| -0.3163313|
|  0.2062678| 0.2401092|  0.4641616|
|  0.4313486| 0.5298258|  0.5934396|

```r
cat('mean: ', mean(df), ' sd: ', sd(df))
```

```
mean:  5.304195e-17  sd:  1
```

We now perform PCA on each Sub_pathway and define the first PC as representative of the respective pathway. 


```r
pathways <- 
  data_env$annotations %>% 
  .$Sub_pathway %>% 
  unique()

eigen_metabolites <- matrix(0, 
                            nrow = nrow(data_env$dat), 
                            ncol = length(pathways),
                            dimnames = list(c(), pathways))

explained_variances <- rep(0, length(pathways))
names(explained_variances) <- pathways


for (pathway in pathways) {
  
  metabolites_in_pathway <- 
    data_env$annotations %>% 
    filter(Sub_pathway == pathway) %>% 
    .$name
  
  if (length(metabolites_in_pathway) >= 2) {
    subpathway_data <- 
      data_env$dat[,metabolites_in_pathway] 

    pathway_rep <- prcomp(t(subpathway_data), 
                             center = TRUE, 
                             scale = TRUE, 
                             rank = 1)
    
    sdev <- pathway_rep$sdev^2
    
    eigen_metabolites[,pathway] <- pathway_rep$rotation
    explained_variances[pathway] <- sdev[1] / sum(sdev)
  }
}
```


**How many variables and samples do you have now?**


```r
dim(eigen_metabolites)
```

```
[1] 906  55
```

We now have **55** dimensions for each of the **906** samples that we
started with. Each of these dimensions represents the first principle 
component of the identified metabolic sub pathways.


```r
sub_pathway_eigen_metabolites <- eigen_metabolites
```


### **2.3 Eigenmetabolite Explained Variance**

**Calculate the explained variance of each eigenmetabolite according to the lecture. Plot the explained variances of all eigenmetabolites, e.g. with a histogram. What is the average explained variance? Why do some pathways have an explained variance of 1?**


```r
hist(explained_variances[explained_variances!=0] * 100, 
  main = 'Explained Variance of PC1 Across Pathways', 
  xlab = '% Explained Variance')

mean_exp_var <- mean(explained_variances)*100
abline(v = mean_exp_var, col = 'red', lty = 2, lwd = 3)
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

We see that the average explained variance is a little below
**43%**, indicated by the red dotted line. I hypothesize
that the average explained variance of some eigenmetabolites
is 1 **because those pathways only contain 1 or 2 metabolites**. 
We show this is actually the case below.


```r
pathways <- names(explained_variances)[explained_variances == 1]
n_metabolites <- sapply(pathways, function(x) {
    data_env$annotations %>% filter(Sub_pathway == x) %>% nrow()
  })

data.frame(pathways, n_metabolites) %>% 
  kable()
```



|                                             |pathways                                     | n_metabolites|
|:--------------------------------------------|:--------------------------------------------|-------------:|
|Monoacylglycerol                             |Monoacylglycerol                             |             2|
|Fatty Acid, Branched                         |Fatty Acid, Branched                         |             2|
|Fatty Acid, Amino                            |Fatty Acid, Amino                            |             2|
|Histidine Metabolism                         |Histidine Metabolism                         |             2|
|Polyamine Metabolism                         |Polyamine Metabolism                         |             2|
|Alanine and Aspartate Metabolism             |Alanine and Aspartate Metabolism             |             2|
|Fatty Acid Metabolism (also BCAA Metabolism) |Fatty Acid Metabolism (also BCAA Metabolism) |             2|
|Phospholipid Metabolism                      |Phospholipid Metabolism                      |             2|
|TCA Cycle                                    |TCA Cycle                                    |             2|
|Creatine Metabolism                          |Creatine Metabolism                          |             2|
|Dipeptide                                    |Dipeptide                                    |             2|


### **2.4**

Perform 2) and 3) for “Super_pathways”. Compare the explained variances plots of Sub- and Super-pathways. What do you observe and why?


```r
pathways <- 
  data_env$annotations %>% 
  .$Super_pathway %>% 
  unique()

eigen_metabolites <- matrix(0, 
                            nrow = nrow(data_env$dat), 
                            ncol = length(pathways),
                            dimnames = list(c(), pathways))

explained_variances <- rep(0, length(pathways))
names(explained_variances) <- pathways


for (pathway in pathways) {
  
  metabolites_in_pathway <- 
    data_env$annotations %>% 
    filter(Super_pathway == pathway) %>% 
    .$name
  
  if (length(metabolites_in_pathway) >= 2) {
    subpathway_data <- 
      data_env$dat[,metabolites_in_pathway] 

    pathway_rep <- prcomp(t(subpathway_data), 
                             center = TRUE, 
                             scale = TRUE, 
                             rank = 1)
    
    sdev <- pathway_rep$sdev^2
    
    eigen_metabolites[,pathway] <- pathway_rep$rotation
    explained_variances[pathway] <- sdev[1] / sum(sdev)
  }
}
```

**How many variables and samples do you have now?**


```r
dim(eigen_metabolites)
```

```
[1] 906   9
```

We now have **9** dimensions for each of the **906** samples that we
started with. Each of these dimensions represents the first principle 
component of the identified metabolic super pathways.



```r
hist(explained_variances[explained_variances!=0] * 100, 
  main = 'Explained Variance of PC1 Across Pathways', 
  xlab = '% Explained Variance')

mean_exp_var <- mean(explained_variances)*100
abline(v = mean_exp_var, col = 'red', lty = 2, lwd = 3)
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


```r
super_pathway_eigen_metabolites <- eigen_metabolites
```

**Comparing explained variance plots** 

We notice that the percentage of explained variance is much lower for super-pathways than for sub-pathways. This is seen by the ~20% decrease of average explained variance between the two plots. In addition, there are sub-pathways which have an average explained variance of 1, discussed previously, which is not the case for the super-pathways, which have a maximum explained variance close to 60%. The reason for this might be because each super-pathway contains many different metabolites, whose variances may not all be explained by the first principal component.


### **2.5 Male vs Female Eigenmetabolite Pathway Analysis**

**Analyze the pathway differences between males and females. Use analysis and visualization techniques you learned in previous lectures.**

**Statistically significant Sub-pathways**


```r
sub_pathways <- colnames(sub_pathway_eigen_metabolites)

p_vals <- vector('double', length(sub_pathways))
names(p_vals) <- sub_pathways

for (sub_pathway in sub_pathways) {
  male   <- sub_pathway_eigen_metabolites[data_env$gender==1, sub_pathway]
  female <- sub_pathway_eigen_metabolites[data_env$gender==2, sub_pathway]
  p_vals[sub_pathway] <- t.test(male, female)$p.value
}

adj_p_vals <- p.adjust(p_vals, 'bonferroni')
hist(adj_p_vals, 
     main='Bonferroni corrected P-values', 
     xlab='P values')
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


**Top 10 Sub-Pathways Male vs Female**


```r
library(knitr)
sorted_p_vals <- sort(adj_p_vals)[1:10]
top_ten <- names(sorted_p_vals)

data.frame(Sub_Pathway = top_ten, 
           P_value = format(sorted_p_vals, digits = 3)) %>% 
  kable()
```



|                                                     |Sub_Pathway                                          |P_value  |
|:----------------------------------------------------|:----------------------------------------------------|:--------|
|Unknown                                              |Unknown                                              |9.28e-50 |
|Leucine, Isoleucine and Valine Metabolism            |Leucine, Isoleucine and Valine Metabolism            |8.30e-43 |
|Steroid                                              |Steroid                                              |2.56e-37 |
|Creatine Metabolism                                  |Creatine Metabolism                                  |4.48e-31 |
|Glutamate Metabolism                                 |Glutamate Metabolism                                 |1.71e-23 |
|Phenylalanine and Tyrosine Metabolism                |Phenylalanine and Tyrosine Metabolism                |8.80e-20 |
|Purine Metabolism, (Hypo)Xanthine/Inosine containing |Purine Metabolism, (Hypo)Xanthine/Inosine containing |1.80e-17 |
|Chemical                                             |Chemical                                             |8.65e-16 |
|Long Chain Fatty Acid                                |Long Chain Fatty Acid                                |1.33e-14 |
|Gamma-glutamyl Amino Acid                            |Gamma-glutamyl Amino Acid                            |4.48e-13 |

**Male vs Female for Top Hit Sub-pathway**

Our top hit other than unknown metabolites is for **Leucine, Isoleucine and Valine Metabolism**. 


```r
top_hit <- names(sort(adj_p_vals)[2])

male <- sub_pathway_eigen_metabolites[data_env$gender==1, top_hit]
female <- sub_pathway_eigen_metabolites[data_env$gender==2, top_hit]

male <- c(male, rep(NA, length(female) - length(male)))

df <- 
  data.frame(Male=male, 
             Female=female) %>%
  mutate(index = row_number()) %>% 
  pivot_longer(-c(index), 
               names_to='gender', values_to="concentration") %>% 
  mutate(gender = factor(gender, levels = c('Male', 'Female')))
  


df %>% 
  ggplot() + 
  geom_boxplot(aes(x=gender, 
                   y=concentration, 
                   fill=gender),
                width = 0.5
               ) +  
  theme(legend.position = "none") + 
  labs(title=paste('Male vs Female Leucine, Isoleucine and Valine Metabolism Sub-Pathway')) + 
  xlab("") + 
  ylab("Normalized Concentration")
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-31-1.png)<!-- -->


**Statistically significant Super-pathways**


```r
super_pathways <- colnames(super_pathway_eigen_metabolites)

p_vals <- vector('double', length(super_pathways))
names(p_vals) <- super_pathways

for (super_pathway in super_pathways) {
  male   <- super_pathway_eigen_metabolites[data_env$gender==1, super_pathway]
  female <- super_pathway_eigen_metabolites[data_env$gender==2, super_pathway]
  p_vals[super_pathway] <- t.test(male, female)$p.value
}

adj_p_vals <- p.adjust(p_vals, 'bonferroni')
hist(adj_p_vals, 
     main='Bonferroni corrected P-values', 
     xlab='P values')
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-32-1.png)<!-- -->


**Top 3 Super-Pathways Male vs Female**


```r
library(knitr)
sorted_p_vals <- sort(adj_p_vals)[1:3]
top_three <- names(sorted_p_vals)

data.frame(Super_Pathway = top_three, 
           P_value = format(sorted_p_vals, digits = 3)) %>% 
  kable()
```



|           |Super_Pathway |P_value  |
|:----------|:-------------|:--------|
|Amino Acid |Amino Acid    |9.77e-69 |
|Unknown    |Unknown       |1.94e-50 |
|Peptide    |Peptide       |7.16e-47 |

**Male vs Female for Top Hit Super-pathway**

Our top hit other than unknown metabolites is for **Amino Acid**. 


```r
top_hit <- names(sort(adj_p_vals)[1])

male <- super_pathway_eigen_metabolites[data_env$gender==1, top_hit]
female <- super_pathway_eigen_metabolites[data_env$gender==2, top_hit]

male <- c(male, rep(NA, length(female) - length(male)))

df <- 
  data.frame(Male=male, 
             Female=female) %>%
  mutate(index = row_number()) %>% 
  pivot_longer(-c(index), 
               names_to='gender', values_to="concentration") %>% 
  mutate(gender = factor(gender, levels = c('Male', 'Female')))
  


df %>% 
  ggplot() + 
  geom_boxplot(aes(x=gender, 
                   y=concentration, 
                   fill=gender),
                width = 0.5
               ) +  
  theme(legend.position = "none") + 
  labs(title=paste('Male vs Female Amino Acid Super-Pathway')) + 
  xlab("") + 
  ylab("Normalized Concentration")
```

![](R/enrichment_analysis_files/figure-html/unnamed-chunk-34-1.png)<!-- -->













