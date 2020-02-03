# HM450

## TCGA Infinium HM450 microarray

10102 samples from 33 cancer types in TCGA. 700+ tumor adjacent normals.

```
~/zhou_lab/projects/20200202_TCGA_HM450
```
Original data at

```
~/zhou_lab/HFS10T/2019_08_28_HPC_Laird_primary/2016_04_05_TCGA_pancan_renormalization
```

Metadata:

- Cancer type annotation: `zhou_lab/projects/20200202_TCGA_HM450/merged_mapping`


## Blood

#### GSE35069 (Sorted Blood)

Blood sorted into CD4-T, CD8-T, CD19-B, CD56-NK, Granulocytes, Monocytes. Two PBMC and two whole blood samples.
60 samples, 485577 probes, IDATs available.

`/mnt/isilon/zhou_lab/HFS10T/2019_08_28_HPC_Laird_secondary/2016_12_26_TCGA_WGBS/hm450/Reinus_Blood/betas.rda`

["Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility."](https://www.ncbi.nlm.nih.gov/pubmed/22848472)

### GSE36064

Peripheral blood leukocytes from healthy children. 78 samples, 485577 probes. Covariates: sex, ethnicity, agebymonth. Intensity matrix available, no IDATs. 

[Age-associated DNA methylation in pediatric populations.](https://www.ncbi.nlm.nih.gov/pubmed/22300631)

`/mnt/isilon/zhou_lab/HFS10T/2019_08_28_HPC_Laird_secondary/2016_12_26_TCGA_WGBS/hm450/GSE36064/betas.rda`

### GSE40279

`/secondary/projects/shen/projects/2015_04_24_blood/betas.rda`


## Tumors

#### GSE129364

72 FFPE primary colorectal adenoma without recurrence (n=30), primary adenoma with recurrence at the same location (n=19), so-called “matched pair samples” (n=10; comprising the primary adenoma and the recurrent adenoma) and normal mucosa specimens (n=3). IDATs available.

`~/zhou_lab/projects/20191212_GEO_datasets/GSE129364`

[Genome-wide DNA methylation analysis of colorectal adenomas with and without recurrence reveals an association between cytosine-phosphate-guanine methylation and histological subtypes.](https://www.ncbi.nlm.nih.gov/pubmed/31334584)


#### GSE120878

FFPE samples for 89 primary invasive melanoma and 73 nevi. IDATs available.

`~/zhou_lab/projects/20191212_GEO_datasets/GSE120878`

[Identification of a Robust Methylation Classifier for Cutaneous Melanoma Diagnosis.](https://www.ncbi.nlm.nih.gov/pubmed/30529013)

#### GSE108576

FFPE samples for 30 breast cancer brain metastases, 18 lung cancer brain metastases, 37 melanoma brain metastases, and 4 samples with brain metastases from patients with uncertain primary. IDATs available.

`~/zhou_lab/projects/20191212_GEO_datasets/GSE108576`

"DNA methylation analysis of brain metastasis"
Unpublished as of 20200202

#### GSE104293

132 IDH-mutant low grade glioma. 87 FFPE samples and 45 fresh frozen samples from patients subject to radiotherapy (RT) or temozolomide (TMZ). "Two of the CpGs correspond to MGMT and consisted of the 2 probes used in the MGMT-STP27 classifier. The extent of methylation was predictive for PFS in the TMZ, but not the RT arm. LGG IDHmt and 1p and 19q codeleted that overlap largely with oligodendroglioma displayed a higher extent of MGMT methylation than the IDHmt non-codeleted tumors based on several independent LGG datasets."

`~/zhou_lab/projects/20191212_GEO_datasets/GSE104293`

[The DNA methylome of DDR genes and benefit from RT or TMZ in IDH mutant low-grade glioma treated in EORTC 22033](https://www.ncbi.nlm.nih.gov/pubmed/29368212)

#### GSE72251

119 FFPE breast cancer samples from "Cohort 2".

`~/zhou_lab/projects/20191212_GEO_datasets/GSE72251`

[DNA methylation-based immune response signature improves patient diagnosis in multiple cancers.](https://www.ncbi.nlm.nih.gov/pubmed/28714863)