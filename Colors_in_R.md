## [grDevices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/palettes.html)

This is the default color palettes defined in R, including `hcl.colors`, `rainbow`, `heat.colors`, `terrain.colors`, `topo.colors`, `cm.colors`. use `colors()` to show the color names R knows about.

```R
heat.colors(4, alpha=1)
```

## [colorspace](https://cran.r-project.org/web/packages/colorspace/vignettes/colorspace.html)

This a very popular package that provides `diverge_hcl`, `diverge_hsl`, `terrain_hcl`, `sequential_hcl`, `rainbow_hcl`, etc.

```R
rainbow_hcl(4)
## "#E495A5" "#ABB065" "#39BEB1" "#ACA4E2“
```

## [pals](https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html)

This is a new collection with some of the most popular palette-generating functions, e.g., `parula`, `turbo`.

Also see [the bivariate color from pals](https://cran.r-project.org/web/packages/pals/vignettes/bivariate_choropleths.html)

## [RColorBrewer](http://applied-r.com/rcolorbrewer-palettes/)

A very popular package, that provides `brewer.pal`.

```R
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()
brewer.pal(8,"Set3")
## [1] "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69"
## [8] "#FCCDE5"
```

[package home](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

## [colorRamps](https://cran.r-project.org/web/packages/colorRamps/index.html)

An old collection which hasn't been much maintained. But still sometimes used.

## [A Cheatsheet](https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf)

## Specify color manually (not recommended)

```R
rgb(r, g, b, maxColorValue=255, alpha=255)
hsv(h, s, v, alpha)
hcl(h, c, l, alpha)
```
