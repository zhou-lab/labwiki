library(caret)
library(ggplot2)

## set up sf here https://github.com/zhou-lab/labwiki/blob/master/Snippets/20210326_sync_HPC_data.sh

## Note 1: same location for samplesheets
## sf r:~/samplesheets/2021/20210226_MouseArray_SampleTableV4_clock.xlsx
## sample sheets, should always be at ~/samplesheets
meta = read_excel('~/samplesheets/2021/20210226_MouseArray_SampleTableV4_clock.xlsx') %>% dplyr::filter(Tumor_vs_Normal == "Normal") %>% dplyr::filter(Issue == "NA")

## Note 2: input data path should follow the same folder structure so it can be sync-ed to and from HPC
## sf r:~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds
betas <- readRDS("~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210104_mouse_array_data_analysis/20210319_1656_mouse_array_samples.rds")

meta1 = meta %>% dplyr::filter(Mouse_Age_Months != "NA") %>% mutate(age = as.numeric(Mouse_Age_Months)) %>% dplyr::filter(Tissue_Corrected != "NA")
betas1 = betas[rowSums(is.na(betas)) < 10,]
fig = ggplot(df, aes(x, y, colour = label)) +
  geom_point() + scale_colour_manual(values=setNames(df$colour,df$label)) + 
  xlab("Reported Age (Month)") + ylab("Predicted Age (Month)") +
    geom_abline()

## Note 3: figures automatically saved to Dropbox
## the following link saves some typing
## ln -s ~/Dropbox/ZhouLab/Lab\ Gallery/zhouw3 ~/gallery
ggsave("~/gallery/20210326_epigenetic_clock_figure.pdf", plot=fig)

## Note 4: any intermediate files need to sit in the same zhoulab folder
## so that you can sync to HPC easily through
## sf ~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds
saveRDS(fit_test_1, "~/zhoulab/labprojects/20200228_Mouse_Array_Project/20210326_epiclocks/20210325_epigenetic_clock_model_3.rds")


