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

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
