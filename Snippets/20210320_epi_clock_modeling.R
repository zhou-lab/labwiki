
# regressing age on DNA methylation
# elastic net select

meta = read_excel('~/samplesheets//2021/20210226_MouseArray_SampleTableV4_clock.xlsx') %>% dplyr::filter(Tumor_vs_Normal == "Normal")

betas = tbk_data(sprintf('~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210104_mouse_array_data_analysis/tbk_MM285/%s.tbk', meta$IDAT), probes=c('cg36659359_TC11', 'cg36683859_TC21'), max_pval=0.2)
dim(betas)

meta1 = meta %>% mutate(beta = betas[1, match(IDAT, colnames(betas))]) %>% select(IDAT, beta, Mouse_Age_Months, Tissue_Corrected) %>% dplyr::filter(Mouse_Age_Months != "NA") %>% mutate(age = as.numeric(Mouse_Age_Months))

summary(lm(age~beta+Tissue_Corrected, data=meta1))

## initial feature selection based on variance
## sort(apply(mtx, 1, sd, na.rm=T))

levels(as.factor(meta1$Tissue_Corrected))

## glmnet

train(age~Tissue_Corrected, )

## https://daviddalpiaz.github.io/r4sl/elastic-net.html
