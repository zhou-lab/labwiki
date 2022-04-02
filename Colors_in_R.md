## [pals](https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html) (Recommended!)

This is a new collection with some of the most popular palette-generating functions, e.g., `parula`, `turbo`.

Also see [the bivariate color from pals](https://cran.r-project.org/web/packages/pals/vignettes/bivariate_choropleths.html)

## [grDevices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/palettes.html)

These are the default color palettes defined in R, including `hcl.colors`, `rainbow`, `heat.colors`, `terrain.colors`, `topo.colors`, `cm.colors`. use `colors()` to show the color names R knows about.

```R
heat.colors(4, alpha=1)
```

## [RColorBrewer](http://applied-r.com/rcolorbrewer-palettes/)

A very popular package, that provides the `brewer.pal()` function. However, this package only provides one (though powerful) palette-generating function. The actual palette name is an argument to `brewer.pal`, such as `Set1`, `Set2`, `Paired`, `Dark`, `Accent`, `Spectral`, `RdBu`, `YlGn`, etc.

```R
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()
display.brewer.all(colorblindFriendly=TRUE)
brewer.pal(8,"Set3")
## [1] "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69"
## [8] "#FCCDE5"
```

## [colorspace](https://cran.r-project.org/web/packages/colorspace/vignettes/colorspace.html)

This a very popular package that provides `diverge_hcl`, `diverge_hsl`, `terrain_hcl`, `sequential_hcl`, `rainbow_hcl`, etc.

```R
library(colorspace)
rainbow_hcl(4)
## "#E495A5" "#ABB065" "#39BEB1" "#ACA4E2“
diverge_hcl(7, h = c(246, 40), c = 96, l = c(65, 90))
pal <- choose_palette()
```

[package home](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

## [colorRamps](https://cran.r-project.org/web/packages/colorRamps/index.html)

An old collection which hasn't been much maintained. But still sometimes used. provides `matlab_like`, `matlab_like2`, `ygobb` etc.

## [A Cheatsheet](https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf)

## Specify color manually (not recommended)

```R
rgb(r, g, b, maxColorValue=255, alpha=255)
hsv(h, s, v, alpha)
hcl(h, c, l, alpha)
```
