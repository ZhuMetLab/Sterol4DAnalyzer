# Sterol4DAnalyzer

Sterol4DAnalyzer is an R package for identifying Sterol lipids with 4 dimensional match 
(m/z, RT, CCS and MSMS).

## INSTALL

#### install "devtools" to help installing from github

``` R
install.packages('devtools')
```

#### install depended [SpectraTools](https://github.com/ZhuMetLab/SpectraTools)

``` R
# install dependencies
install.packages('BiocManager')
BiocManager::install(c('BiocParallel', 'ggplot2', 'ggrepel'))
devtools::install_github('ZhuMetLab/SpectraTools')
```

#### install Steorol4DAnalyzer

``` R
# install dependencies
install.packages('BiocManager')
devtools::install_github('ZhuMetLab/Sterol4DAnalyzer')
```

#### run Steorol4DAnalyzer
``` $
library(Steorol4DAnalyzer)

```
