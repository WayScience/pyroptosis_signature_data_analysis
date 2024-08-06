suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))
source("../../utils/figure_themes.r")

death_curve_path <- file.path("..","..","..","data","Death_curve_data.csv")
# figure output_path
output_path <- file.path("..","figures")
if (!dir.exists(output_path)){
  dir.create(output_path, recursive = TRUE)
}
oputput_figure_path <- file.path(output_path,"death_curve.png")

death_df <- read.csv(death_curve_path)
head(death_df)

# death percentage
death_df$percentage_dead_cells <- death_df$Death.Cells/(death_df$Death.Cells + death_df$Live.Cells)	*100
death_df$treatment <- paste0(death_df$Compound," ",death_df$Dose," ",death_df$unit, " ",death_df$time, " ", death_df$`unit.1`)
head(death_df)
# aggregate data by treatment
death_df_agg <- death_df %>% group_by(treatment, Compound, unit) %>% summarise(mean_death = mean(percentage_dead_cells), sd_death = sd(percentage_dead_cells), time = mean(time), dose = mean(Dose))
# make dose as factor
death_df_agg$dose <- factor(death_df_agg$dose,levels = c("0","1","10"))
death_df_agg$treatment <- paste0(death_df_agg$Compound," ",death_df_agg$dose, " ",death_df_agg$unit)
# factorize treatment
death_df_agg$treatment <- factor(death_df_agg$treatment,levels = c(
    'Media 0 ',
    'DMSO 0 ',
    'LPS 1 ug/mL',
    'LPS 10 ug/mL',
    'Thapsigargin 1 uM',
    'Thapsigargin 10 uM'
))
head(death_df_agg)

# gsub
death_df_agg$treatment <- gsub("Media 0 ","Media",death_df_agg$treatment)
death_df_agg$treatment <- gsub("DMSO 0 ","DMSO 0.1%",death_df_agg$treatment)
death_df_agg$treatment <- gsub("LPS 1 ug/mL","LPS 1.0 ug/mL",death_df_agg$treatment)
death_df_agg$treatment <- gsub("LPS 10 ug/mL","LPS 10.0 ug/mL",death_df_agg$treatment)
death_df_agg$treatment <- gsub("Thapsigargin 1 uM","Thapsigargin 1.0 uM",death_df_agg$treatment)
death_df_agg$treatment <- gsub("Thapsigargin 10 uM","Thapsigargin 10.0 uM",death_df_agg$treatment)
# factorize treatment
death_df_agg$time <- paste0(death_df_agg$time," Hours")
death_df_agg$treatment <- factor(death_df_agg$treatment,levels = c(
    'Media',
    'DMSO 0.1%',
    'LPS 1.0 ug/mL',
    'LPS 10.0 ug/mL',
    'Thapsigargin 1.0 uM',
    'Thapsigargin 10.0 uM'
))
head(death_df_agg)
unique(death_df_agg$treatment)

# plot the death curve data
width <- 17
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
death_curve_plot <- (
    ggplot(death_df_agg, aes(x = treatment, y = mean_death, fill = treatment))
    + geom_bar(stat = "identity", position = "dodge")
    + geom_errorbar(aes(ymin = mean_death - sd_death, ymax = mean_death + sd_death), width = 0.25, position = position_dodge(0.9))
    + theme_bw()
    + figure_theme
    + ylim(0, 100)
    + labs(
           x = "Treatment",
           y = "Percentage of dead cells")
    + facet_grid(~time)
    + theme(
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        # legend.position = "none"
    )
    + scale_fill_manual(
        name = "Treatment",
        labels = c(
            'Media',
            'DMSO 0.1%',
            'LPS 1.0 ug/mL',
            'LPS 10.0 ug/mL',
            'Thapsigargin 1.0 uM',
            'Thapsigargin 10.0 uM'
        ),
        values = death_curve_colors)
)
# save the plot
ggsave(oputput_figure_path, plot = death_curve_plot, width = width, height = height, units = "in")
death_curve_plot
