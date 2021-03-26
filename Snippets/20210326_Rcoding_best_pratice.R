library(caret)
library(ggplot2)

## sf r:~/samplesheets/2021/20210226_MouseArray_SampleTableV4_clock.xlsx
## sample sheets, should always be at ~/samplesheets
meta = read_excel('~/samplesheets/2021/20210226_MouseArray_SampleTableV4_clock.xlsx') %>% dplyr::filter(Tumor_vs_Normal == "Normal") %>% dplyr::filter(Issue == "NA")
## input data can be sync-ed from HPC
## sf r:~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds
betas <- readRDS("~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210104_mouse_array_data_analysis/20210319_1656_mouse_array_samples.rds")
meta1 = meta %>% dplyr::filter(Mouse_Age_Months != "NA") %>% mutate(age = as.numeric(Mouse_Age_Months)) %>% dplyr::filter(Tissue_Corrected != "NA")

betas1 = betas[rowSums(is.na(betas)) < 10,]

fig = ggplot(df, aes(x, y, colour = label)) +
  geom_point() + scale_colour_manual(values=setNames(df$colour,df$label)) + 
  xlab("Reported Age (Month)") + ylab("Predicted Age (Month)") +
    geom_abline()

## make sure you save the figure to Dropbox
## the following link saves some typing
## ln -s ~/Dropbox/ZhouLab/Lab\ Gallery/zhouw3 ~/gallery
ggsave("~/gallery/20210326_epigenetic_clock_figure.pdf", plot=fig)

## any intermediate files need to sit in the same zhoulab folder
## so that you can sync to HPC easily through
## sf ~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds
saveRDS(fit_test_1, "~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds")


