module load XZ/5.2.4-GCCcore-7.3.0
module load libreadline/7.0-GCCcore-7.3.0
module load bzip2/1.0.6-GCCcore-7.3.0
module load cURL/7.60.0-GCCcore-7.3.0
module load libxml2/2.9.8-GCCcore-7.3.0
module load libpng/1.6.37

wget https://cran.r-project.org/src/base-prerelease/R-devel.tar.gz
cd /mnt/isilon/zhoulab/labsoftware/shared_Renv/versions/4.1.devel/R-devel
./configure --prefix=/mnt/isilon/zhoulab/labsoftware/shared_Renv/versions/4.1.devel -with-recommended-packages=no --without-x --with-cairo --with-libpng --with-libtiff --with-jpeglib --enable-R-shlib --with-pcre1
make
make install

# if you run into Rhdf5lib error, then install the following manually from github
BiocManager::install("grimbough/Rhdf5lib", ref = "cpp-flags")


### run the following in R
Pacs <- c(
    'devtools',
    'BiocManager',
    'reshape2',
    "viridis",
    'tidyverse')

install.packages(pacs, repo = "https://cloud.r-project.org")

options(timeout=1000)

biocPacs <- c(
    'BiocStyle',
    'sesame',
    'BiocCheck',
    'IlluminaHumanMethylation450kmanifest',
    'FlowSorted.CordBloodNorway.450k',
    'ggsci',
    'FlowSorted.Blood.450k',
    'IlluminaHumanMethylation450kanno.ilmn12.hg19'
)

BiocManager::install(biocPacs)
