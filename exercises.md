---
title: "Exercise for Lecture 7"
author: "Belen Jimenez-Mena"
date: "13/10/2021"
output: word_document
---

```{r setup, echo=F}
knitr::opts_chunk$set(echo = T)
```

## R Markdown

To complement Lecture 7 ("Elements of Population Genetics for Conservation and Management"), we are going to do a full analysis on a real genomic dataset (i.e. VCF-type), in which we are going to analyse the genetic diversity, levels of differentiation and population structure of an aquatic species. I based these exercises on: <http://https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/m>, so refer to that for more exercises.

For this part you require the use of Computerome cluster (you are already familiar with it), where we will be using the programming language R. 

## Data
For the exercises we will use a genomic dataset of an lobster population, published in a recent article: Single nucleotide polymorphisms reveal a genetic cline across the north‐east Atlantic and enable powerful population assignment in the European lobster by Jenkins et al. (2019). We extracted the dataset and stored it in the cluster. The data can be found in '/home/franb/EXERCICES/LECTURE06/'. We'll show you how to copy it into your own folder so you can follow the practical. We have also uploaded the dataset at DTULearn.

The lobster data that we present here represents variation data from all over the genome generated with a restriction‐site‐associated DNA sequencing. We will apply the tools we have seen throughout the lecture to this dataset to obtain insights for conservation/management in this population.

First, create a folder where you will put the data and the results generated in each of the exercises.
```{r, eval=F}
mkdir LECTURE07
cd LECTURE07
mkdir Results
mkdir Data
```

You can make a copy of the dataset by using the following command in the terminal, from your directory.
```{r, eval=F}
cp /home/LECTURE_07/* LECTURE07/Data/.
```

We also need to load a few modules in the cluster for R to work:
```{r, eval=F}
module load udunits/2.2.26
module load intel/redist/2019_update2
module load intel/perflibs/64
module load lapack/3.8.0
module load gcc/8.2.0
module load proj.4/4.9.3
module load gdal/2.2.3
module load R/4.0.0
```

## Case study
You have genetic data from the European lobster (Homarus gammarus), obtained at different locations around Europe. The company you work for has an interest in conservation of this species, and would like to know more about how is its genetic diversity, are they different populations and how are they different from each other. You are the geneticist in charge of giving back the answers. Let's start!


### R: the programming language for today
As mentioned, we are going to use the programming language R for the analysis and visualization. R is a programming language suited for statistical computing that has been developed by the scientific community and it is widely used for data analysis.

We can call R from the terminal. We do this by simply using
```{r, eval=F}
R
```

If you want to exit R, you just have to type:
```{r, eval=F}
q()
```

and decide whether you want to save or not your session.


Once you are inside R, we need a few packages in order to do the analysis. In R, appart from coding our own functions and programs, we can download packages and use the functions that have been developed and saved there. This is one of the advantages of using R, its active user community that helps us not re-invent the wheel. To install an R package, we have just to type (remember to call R before, aka. the previous command, otherwise, the following commands won't work in the BASH terminal):

## Packages needed for these exercises
There are thousands of helpful R packages for you to use. For the analysis at this lecture, we will be using the following packages. If you work in the cluster, you don't have to install the packages, as they are already installed (you only need to load them using require(), but if you had to install them again, you just uncomment the install.packages() commands). If you work from your laptop, you'd need to install them if you don't have them yet.

```{r, eval=T, message = FALSE}
library(adegenet)
library(poppr)
library(dplyr)
library(hierfstat)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
```

## Loading the dataset
There are many ways to load a dataset in R. For this exercise we will be using a csv file containing SNP (single nucleotide polymorphism) genotypes. You can also work directly with VCF file (and you can use the vcfR package to load that type into R). In our case for today we can load the dataset by:
```{r, eval=T}
# if you need, set your working directory:
setwd("H:/DTU/5. Teaching-Presentations/2021_10_ConservationGenetics/")
# import the dataset
lobster = read.csv("Lobster_SNP_Genotypes.csv")
str(lobster)
```

Convert data.frame from long to wide format. The wide format contains one row for each individual and one column for each locus as well as a column for the ID and site labels.
```{r, eval=T}
lobster_wide = reshape(lobster, idvar = c("ID","Site"), timevar = "Locus", direction = "wide", sep = "")

# Remove "Genotype" from column names
colnames(lobster_wide) = gsub("Genotype", "", colnames(lobster_wide))
```

Subset genotypes and only keep SNP loci used in Jenkins et al. 2019.
```{r, eval=T}
# Subset genotypes
snpgeno = lobster_wide[ , 3:ncol(lobster_wide)]

# Keep only SNP loci used in Jenkins et al. 2019
snps_to_remove = c("25580","32362","41521","53889","65376","8953","21197","15531","22740","28357","33066","51507","53052","53263","21880","22323","22365")
snpgeno = snpgeno[ , !colnames(snpgeno) %in% snps_to_remove]
```


Create vectors of individual and site labels.
```{r, eval=T}
ind = as.character(lobster_wide$ID) # individual ID
site = as.character(lobster_wide$Site) # site ID
```

We will convert the data.frame to a special format to analyse large genomic datasets, a "genlight" and "genind" objects. Check that the genotypes for the first five individuals and loci are as expected.
```{r, eval=T}
lobster_gen = df2genind(snpgeno, ploidy = 2, ind.names = ind, pop = site, sep = "")
lobster_gen$tab[1:5, 1:10]
```

## First look
First, let's look how our data looks like. This is the first step when analysing a genomic dataset. Genomic datasets are quite large and therefore it is a bit difficult ("humanly impossible") to check each data row one by one, by hand. That's why we write scripts and use functions, so we can automate the checks and the analysis. It reduces errors (we are humans!) and assures reproducibility.

By typing the name where we store the dataset (see earlier commands),
```{r, eval=F}
lobster_gen
```

it will show us several lines of information, to summarize all the dataset. Basically we will be able to see a summary of what is inside "lobster_gen".
```{r, eval=T}
lobster_gen
```

EXERCISE: What does all this mean? Please spend a few minutes familiarizing yourself with the information in each entry of the dataset. How many genotypes do we have? What does that mean? How many SNPs have each individual? How much is the % of missing data? Discuss with your partner. 


## Filtering your dataset 
As a second step in a population genetics/genomic analysis, we would need to filter the SNP dataset. There are many ways of filtering a dataset, and it all depends what we are interested in. One thing to remember is that whatever we choose for filtering steps, we would need to report all the steps, so other scientists can replicate our analysis and understand why the results are the way they are. Some of the parameters that genomicists filter their data on is in the % of missing data, or minor allele frequencies. For the sake of simplicity, we will only apply minimum filtering: removing loci with a large amount of missing data, and removing individuals with also higher missing data values.

Calculate the percentage of complete genotypes per loci in the lobster SNP data set.
```{r, eval=T}
locmiss_lobster = propTyped(lobster_gen, by = "loc")
locmiss_lobster[which(locmiss_lobster < 0.80)] # print loci with < 80% complete genotypes


# Barplot
barplot(locmiss_lobster, ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.7)
```

Calculate the percentage of complete genotypes per individual in the lobster SNP data set.
```{r, eval=T}
indmiss_lobster = propTyped(lobster_gen, by = "ind")
indmiss_lobster[ which(indmiss_lobster < 0.80) ] # print individuals with < 80% complete genotypes
```


Remove individuals with > 20% missing genotypes.
```{r, eval=T}
lobster_gen = missingno(lobster_gen, type = "geno", cutoff = 0.20)
```

## Check genotypes are unique
Check all individual genotypes are unique. Duplicated genotypes can result from unintentionally sampling the same individual twice or from sampling clones.

```{r, eval=T}
# Print the number of multilocus genotypes
mlg(lobster_gen)
```

## Identify duplicated genotypes
```{r, eval=T}
dups_lobster = mlg.id(lobster_gen)
for (i in dups_lobster){ # for each element in the list object
  if (length(dups_lobster[i]) > 1){ # if the length is greater than 1
    print(i) # print individuals that are duplicates
  }
}
```

## Remove duplicated genotypes.
```{r, eval=T}
# Create a vector of individuals to remove
lob_dups = c("Laz4","Eye15","Eye16","Eye01","Laz2","Eye08","Gul101","Eye25","Iom02","Hel07","Eye27","Eye05","Eye06","Eye23","Eye22","Eye11","Cro08","Tar1","Eye14","Tar3","Lyn04","Lyn15","Eye07","Eye02","Eye20")
```

```{r, eval=T}
# Create a vector of individual names without the duplicates
lob_Nodups = indNames(lobster_gen)[! indNames(lobster_gen) %in% lob_dups]
```


```{r, eval=T}
# Create a new genind object without the duplicates
lobster_gen = lobster_gen[lob_Nodups, ]
```

## Check loci are still polymorphic after filtering
```{r, eval=T}
isPoly(lobster_gen) %>% summary
```

## Summary statistics
Basic info again:
```{r, eval=T}
lobster_gen
```


EXERCISE: Have a look at the number of alleles per locus. Previously, have a thought about what do you expect to find? 
```{r, eval=T}
table(lobster_gen$loc.fac)
```

How many individuals do you have in each site?
```{r, eval=T}
summary(lobster_gen$pop)
```


## Allelic richness per site across all loci
```{r, eval=T}
allelic.richness(genind2hierfstat(lobster_gen))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)
```

EXERCISE: Are these numbers lower or higher than you expected? 

## Heterozygosity per site

```{r, eval=T}
# Calculate basic stats using hierfstat
basic_lobster = basic.stats(lobster_gen, diploid = TRUE)
```


```{r, eval=T}
# Mean observed heterozygosity per site
Ho_lobster = apply(basic_lobster$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Ho_lobster
```


```{r, eval=T}
# Mean expected heterozygosity per site
He_lobster = apply(basic_lobster$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
He_lobster
```


Visualise heterozygosity per site

```{r, eval=T}
# Create a data.frame of site names, Ho and He and then convert to long format
Het_lobster_df = data.frame(Site = names(Ho_lobster), Ho = Ho_lobster, He = He_lobster) %>%
  melt(id.vars = "Site")
```

```{r, eval=T}
# Custom theme for ggplot2
custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
  )

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])
```

```{r, eval=T}
# Lobster heterozygosity barplot
ggplot(data = Het_lobster_df, aes(x = Site, y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black")+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  ggtitle("European lobster")+
  custom_theme
```


EXERCISE: Do expected and observed Heterozygosity greatly differ, or the differences are small? Are there great differences among sites?

## Fst
Compute pairwise FST (Weir & Cockerham 1984). For this analysis we are going to subselect a few locations so we can speed up the computation:
```{r, eval=T}
# Subset data sets to reduce computation time
lobster_gen_sub = popsub(lobster_gen, sublist = c("Ale","Ber","Brd","Pad","Sar17","Vig"))

# Compute pairwise Fsts
lobster_fst = genet.dist(lobster_gen_sub, method = "WC84")
lobster_fst %>% round(digits = 3)
```

Visualise pairwise FST

```{r, eval=T}
# Convert dist object to data.frame
fst.matrix = as.matrix(lobster_fst)
ind = which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
        )
```

EXERCISE: Based on these values, which populations of these ones we have subselected are more differentiated? Which are more similar?

## Principal Component Analysis (PCA)
We will perform a principal component analysis (PCA) where we will visualize the population structure at an individual level. For the principal component analysis (PCA), we will use some functions that handles genomic objects very "quickly", and also we will continue using the subset of samples for which we estimated Fst previously.


```{r, eval=T}
# Replace missing data with the mean allele frequencies
x = tab(lobster_gen_sub, NA.method = "mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))

```

Visualise PCA results.
```{r, eval=T}
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(lobster_gen_sub)

# Add a column with the site IDs
ind_coords$Site = lobster_gen_sub$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(lobster_gen_sub), "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("Lobster PCA")+
  # custom theme
  ggtheme
```


EXERCISE: How do you see the distribution of populations in the plot: are all populations clustered together in the same area, or some are further appart? How many clusters do you see in this oyster population? Discuss among your peers.

EXERCISE: How this PCA reflect the patterns observed in the analysis of Fst?