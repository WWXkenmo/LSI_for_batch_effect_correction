## Latent semantic indexing for single cell batch effects correction
This repo is for batch effects correction through LSI, which is the fully R implements of RPCI [1]

### running example
we used the example data downloaded from the following links
(a) the matrix file ("bcc_scRNA_counts.txt")
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fbcc%5FscRNA%5Fcounts%2Etxt%
(b) meta_cell file ("bcc_all_metadata.txt")
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813%5Fbcc%5Fall%5Fmetadata%2Etxt%2E
users could first processing the data using the RPCI_example_code, to check the RPCI generate cell embedding vector
```
head(data0@DimReduction$cell.pls)
```
our provided implements could be run through following code
```
embed <- LSI(list(dat6, dat1, dat2, dat3, dat4, dat5, dat7, dat8),var0,18,50,center=TRUE)
```
to make comparison
```
cor(data0@DimReduction$cell.pls[,dims],embed[,dims])
```
I need to note that once the dims chose to > 18, two embedding would be different, the reason is that when doing second round SVD, the inner producted matrix is generated base on the first round truncted singular vector, therefore the latent dimensionality is constraint by the first round SVD (which is, in our case, 18). So if we chose a dimensions > 18 to perform SVD (50), the >=19th singular vector would be random.  

[1] Liu Y, Wang T, Zhou B, et al. Robust integration of multiple single-cell RNA sequencing datasets using a single reference space[J]. Nature biotechnology, 2021: 1-8.
