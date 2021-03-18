cd sesame

## generate document
wzdocument .

wzbuild1_vignette.R
wzbuild2_test.R
wzbuild3_check.R
wzbuild4_bioccheck.R

cd ..
R CMD BUILD sesame
