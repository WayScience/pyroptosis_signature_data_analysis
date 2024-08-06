

suppressPackageStartupMessages(library(ggplot2))

palette <- c('#D81B60',
    '#1E88E5',
    '#FFC107',
    '#004D40',
    '#519C09',
    '#980CAB')

figure_theme <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
    )
)

figure_theme_wide <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
    )
)


suppressPackageStartupMessages(library(ggplot2))

palette <- c('#D81B60',
    '#1E88E5',
    '#FFC107',
    '#004D40',
    '#519C09',
    '#980CAB')

figure_theme <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
    )
)

figure_theme_wide <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
    )
)

# 9 colors
# col 1 - 3 hues
colorgrad1 <- c("grey", "grey", "#585858")
# col 2 - 5 hues color ramp 5 hues
colorgrad2 <- colorRampPalette(c("pink", "darkred"))(11)
# col 3 - 3 hues
colorgrad3 <- colorRampPalette(c("yellow", "brown"))(3)
# col 4 - 3 hues
colorgrad4 <- colorRampPalette(c("lightblue", "darkblue"))(6)
# col 5 - 2 hues
colorgrad5 <- colorRampPalette(c("lightgreen", "darkgreen"))(2)
# col 6 - 3 hues
colorgrad6 <- colorRampPalette(c("purple", "#2e004b"))(3)
# col 7 - 2 hues
colorgrad7 <- colorRampPalette(c("cyan", "darkcyan"))(4)
# col 8 - 2 hues
colorgrad8 <- colorRampPalette(c("#ebb676", "darkorange"))(2)
# col 9 - 3 hues
colorgrad9 <- colorRampPalette(c("magenta", "#833b83"))(3)

colors <- c(
    'Media' = colorgrad1[1],
    'DMSO 0.1%' = colorgrad1[2],

    'Disulfiram 0.1 uM' = colorgrad4[3],
    'Disulfiram 1.0 uM' = colorgrad4[4],
    'Disulfiram 2.5 uM' = colorgrad4[5],

    'Flagellin 0.1 ug/ml'= colorgrad5[1],
    'Flagellin 1.0 ug/ml' = colorgrad5[2],

    'LPS 0.01 ug/ml' = colorgrad2[1],
    'LPS 0.1 ug/ml' = colorgrad2[2],
    'LPS 1.0 ug/ml' = colorgrad2[3],
    'LPS 10.0 ug/ml' = colorgrad2[4],
    'LPS 100.0 ug/ml' = colorgrad2[5],
    'LPS 1.0 ug/ml + Nigericin 1.0 uM' = colorgrad2[6],
    'LPS 1.0 ug/ml + Nigericin 3.0 uM' = colorgrad2[7],
    'LPS 1.0 ug/ml + Nigericin 10.0 uM' = colorgrad2[8],
    'LPS 100.0 ug/ml + Nigericin 1.0 uM' = colorgrad2[9],
    'LPS 100.0 ug/ml + Nigericin 3.0 uM' = colorgrad2[10],
    'LPS 100.0 ug/ml + Nigericin 10.0 uM' = colorgrad2[11],

    'H2O2 100.0 nM'  = colorgrad7[1],
    'H2O2 100.0 uM' = colorgrad7[2],

    'Thapsigargin 1.0 uM' = colorgrad8[1],
    'Thapsigargin 10.0 uM' = colorgrad8[2],

    'Topotecan 5.0 nM' = colorgrad9[1],
    'Topotecan 10.0 nM' = colorgrad9[2],
    'Topotecan 20.0 nM' = colorgrad9[3]

)
# define the shaps of the points
shape_list <- c(
    19, # Media
    19, # DMSO 0.025%
    19, # DMSO 1.0%
    18, # Disulfiram 0.1 uM
    18, # Disulfiram 1.0 uM
    18, # Disulfiram 2.5 uM
    17, # Z-VAD-FMK 30.0 uM
    17 # Z-VAD-FMK 100.0 uM
)

shapes <- c(
    'Media' = shape_list[1],
    'DMSO 0.025%' = shape_list[2],
    'DMSO 1.0%' = shape_list[3],

    'Disulfiram 0.1 uM' = shape_list[4],
    'Disulfiram 1.0 uM' = shape_list[5],
    'Disulfiram 2.5 uM' = shape_list[6],

    'Z-VAD-FMK 30.0 uM' = shape_list[7],
    'Z-VAD-FMK 100.0 uM' = shape_list[8]
)

# fig 2 colors
colors_2 <- c(
    'DMSO 0.1%, Control' = colorgrad1[2],

    'Flagellin 0.1 ug/ml, Control'= colorgrad5[1],
    'Flagellin 1.0 ug/ml, Pyroptosis' = colorgrad5[2],

    'LPS 0.01 ug/ml, Pyroptosis' = colorgrad2[1],
    'LPS 0.1 ug/ml, Pyroptosis' = colorgrad2[2],
    'LPS 1.0 ug/ml, Pyroptosis' = colorgrad2[3],
    'LPS 10.0 ug/ml, Pyroptosis' = colorgrad2[4],
    'LPS 100.0 ug/ml, Pyroptosis' = colorgrad2[5],
    'LPS 1.0 ug/ml + Nigericin 1.0 uM, Pyroptosis' = colorgrad2[9],
    'LPS 1.0 ug/ml + Nigericin 3.0 uM, Pyroptosis' = colorgrad2[10],
    'LPS 1.0 ug/ml + Nigericin 10.0 uM, Pyroptosis' = colorgrad2[11],

    'H2O2 100.0 nM, Control'  = colorgrad7[1],
    'H2O2 100.0 uM, Control' = colorgrad7[2],

    'Thapsigargin 1.0 uM, Apoptosis' = colorgrad8[1],
    'Thapsigargin 10.0 uM, Apoptosis' = colorgrad8[2]
)

shapes_2 <- c(
    'DMSO 0.1%, Control' = 18,

    'Flagellin 0.1 ug/ml, Control'= 18,
    'Flagellin 1.0 ug/ml, Pyroptosis' = 16,

    'LPS 0.01 ug/ml, Pyroptosis' = 16,
    'LPS 0.1 ug/ml, Pyroptosis' = 16,
    'LPS 1.0 ug/ml, Pyroptosis' = 16,
    'LPS 10.0 ug/ml, Pyroptosis' = 16,
    'LPS 100.0 ug/ml, Pyroptosis' = 16,
    'LPS 1.0 ug/ml + Nigericin 1.0 uM, Pyroptosis' = 16,
    'LPS 1.0 ug/ml + Nigericin 3.0 uM, Pyroptosis' = 16,
    'LPS 1.0 ug/ml + Nigericin 10.0 uM, Pyroptosis' = 16,

    'H2O2 100.0 nM, Control'  = 16,
    'H2O2 100.0 uM, Control' = 16,

    'Thapsigargin 1.0 uM, Apoptosis' = 17,
    'Thapsigargin 10.0 uM, Apoptosis' = 17

)

# fig 2 C colors
colors_3 <- c(
    'DMSO 0.1%' = colorgrad1[2],

    'Flagellin 0.1 ug/ml'= colorgrad5[1],
    'Flagellin 1.0 ug/ml' = colorgrad5[2],

    'LPS 0.01 ug/ml' = colorgrad2[1],
    'LPS 0.1 ug/ml' = colorgrad2[2],
    'LPS 1.0 ug/ml' = colorgrad2[3],
    'LPS 10.0 ug/ml' = colorgrad2[4],
    'LPS 100.0 ug/ml' = colorgrad2[5],
    'LPS 1.0 ug/ml + Nigericin 1.0 uM' = colorgrad2[9],
    'LPS 1.0 ug/ml + Nigericin 3.0 uM' = colorgrad2[10],
    'LPS 1.0 ug/ml + Nigericin 10.0 uM' = colorgrad2[11],

    'H2O2 100.0 nM'  = colorgrad7[1],
    'H2O2 100.0 uM' = colorgrad7[2],

    'Thapsigargin 1.0 uM' = colorgrad8[1],
    'Thapsigargin 10.0 uM' = colorgrad8[2]
)

death_curve_colors <- c(
    'Media' = colorgrad1[1],
    'DMSO 0.1%' = colorgrad1[2],
    'LPS 1.0 ug/mL' = colorgrad2[3],
    'LPS 10.0 ug/mL' = colorgrad2[4],
    'Thapsigargin 1.0 uM' = colorgrad8[1],
    'Thapsigargin 10.0 uM' = colorgrad8[2]
)
