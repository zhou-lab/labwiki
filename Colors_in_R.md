Most functions either take a palette-generating function or a vector of colors. 

```R
palette_gen_func = grDevices::colorRampPalette(color_vector)
color_vector = palette_gen_func(n_colors)
```

It is common in R to start with a palette-generating function (like `brewer.pal`), generate a color vector and regenerate a palette-generating function.

Also not very commonly, there is a `colorRamp` function. More can be found [here](https://bookdown.org/rdpeng/exdata/plotting-and-color-in-r.html#colorramp)

## [pals](https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html) (Recommended!)

This is a new collection with some of the most popular palette-generating functions, e.g., `parula`, `turbo`.

```R
> brewer.paired(3) # this is the same as brewer.pal(3, "Paired") but allow more than 12 colors
> brewer.set1(3)
> brewer.set2(3)
> brewer.set3(3)
```

Also see [the bivariate color from pals](https://cran.r-project.org/web/packages/pals/vignettes/bivariate_choropleths.html)

## [BuenColors](https://github.com/caleblareau/BuenColors)

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

[Here](https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html) is a good tutorial on how to use it.

```R
library(colorspace)
rainbow_hcl(4)
## "#E495A5" "#ABB065" "#39BEB1" "#ACA4E2“
diverge_hcl(7, h = c(246, 40), c = 96, l = c(65, 90))
pal <- choose_palette()
hcl_palettes("sequential (multi-hue)", n = 7, plot = TRUE)
sequential_hcl
```

[package home](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

## [colorRamps](https://cran.r-project.org/web/packages/colorRamps/index.html)

An old collection which hasn't been much maintained. But still sometimes used. provides `matlab_like`, `matlab_like2`, `ygobb` etc.

## [A Cheatsheet](https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf)

## [ggsci](https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html)

Some additional color palette used in popular journals.

Also see [ggpubr](https://rpkgs.datanovia.com/ggpubr/reference/get_palette.html)

## Specify color manually (not recommended)

```R
rgb(r, g, b, maxColorValue=255, alpha=255)
hsv(h, s, v, alpha)
hcl(h, c, l, alpha)
```
