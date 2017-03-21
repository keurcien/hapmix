# A PCA-Based Method to detect Introgression
Keurcien Luu, Michael Blum  
March 22, 2017  



## Introduction

### **Example**

Populus Balsamifera            |  Populus Trichocarpa
:-------------------------:|:-------------------------:
<img src="img/Populus_balsamifera.jpg" style="width: 350px;"/>  |  <img src="img/Populus_trichocarpa.jpg" style="width: 350px;"/>


## Introduction

### **Example**

<div class="centered">
  <p>Populus hybrid species</p>
  <img src="img/Populus_hybrid.jpg" style="width: 350px;"/> 
</div>

## Method

*"Indeed, many species show chromosome-scale variation in diversity and 
divergence; species phylogenies can differ along the
genome due to incomplete lineage sorting, adaptive introgression and/or local 
adaptation." (Li, 2016)*

## Method

*"Indeed, many species show chromosome-scale variation in diversity and 
divergence; species phylogenies can differ along the
genome due to incomplete lineage sorting, adaptive introgression and/or local 
adaptation." (Li, 2016)*

## Method

### **Step 1: Normalization**

Let $G$ denote the genotype matrix, and $\tilde{G}$ the normalized genotype matrix
such that for SNP $i$ and individual $j$:

$$\tilde{G}_{ij} = \frac{G_{ij} - f_i}{\sqrt{2 f_i (1 - f_i)}}$$
where $f_i$ represents the minor allele frequency of SNP $i$.

### **Step 2: PCA**

$$\tilde{G} = U \Sigma V^T$$

## Method




## Method


```
## Number of SNPs: 1500
## Number of individuals: 150
```

![](presentation_files/figure-html/unnamed-chunk-1-1.png)<!-- -->



## Method


```r
plot(y, option = "scores", pop = pop)
```

![](presentation_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

## Method




```r
plot(x_local, option = "scores", pop = pop)
```

![](presentation_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

## Simulations

Dataset: Populus

### **Step 1: Ancestral Haplotypes**

Beagle: phase the ancestral haplotypes using the three populations (including the 
hybrid species)

### **Step 1: **




## Results

