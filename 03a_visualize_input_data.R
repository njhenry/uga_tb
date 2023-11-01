## #######################################################################################
##
## VISUALIZE MODEL INPUTS
##
## AUTHOR: Nat Henry
## CREATED: 1 September 2021
## PURPOSE: Visualize input regions, populations, and case notifications
##
## #######################################################################################

load_libs <- c('data.table', 'sf', 'ggplot2', 'glue')
invisible(lapply(load_libs, library, character.only=TRUE))

DATA_VERSION <- '20231024'

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))

config <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/prepped_data', DATA_VERSION, 'config.yaml')
)
viz_dir <- file.path(config$get_dir_path('prepped_data'), 'viz')
dir.create(viz_dir, showWarnings = FALSE)

# Load notifications dataset
notif_data <- config$read("prepped_data", "notif_data")

## Plot bacteriologically-confirmed and clinically-diagnosed cases over time ------------>

tb_labels <- c('Bacteriologically\nconfirmed', 'Clinically\ndiagnosed', 'Extra-pulmonary')
plot_fill_colors <- RColorBrewer::brewer.pal(name="Pastel1", n=3)
names(plot_fill_colors) <- tb_labels

notif_data_long <- rbindlist(list(
  notif_data[, .(ADM1_EN, ADM2_EN, year, lower = 0, upper = PBC, type = tb_labels[1])],
  notif_data[, .(ADM1_EN, ADM2_EN, year, lower = PBC, upper = PBC + PCD, type = tb_labels[2])],
  notif_data[, .(ADM1_EN, ADM2_EN, year, lower = PBC + PCD, upper = PBC + PCD + EPTB, type = tb_labels[3])]
))

adm1_names <- notif_data_long[, sort(unique(ADM1_EN))]
adm1_plots <- vector('list', length = length(adm1_names))
names(adm1_plots) <- adm1_names

for(adm1_name in adm1_names){
  notifs_long_sub <- notif_data_long[ADM1_EN == adm1_name, ]

  n_districts <- notifs_long_sub[, uniqueN(ADM2_EN) ]
  plot_n_cols <- floor(sqrt(n_districts))

  adm1_plots[[adm1_name]] <- ggplot(data = notifs_long_sub) + 
    facet_wrap('ADM2_EN', ncol = plot_n_cols, scales = 'free_y') + 
    geom_ribbon(
      aes(fill = type, x = year, ymin = lower, ymax = upper),
      alpha = .8, color = '#444444'
    ) + 
    geom_vline(xintercept = 2020, color = '#888888', linewidth = 1, linetype = 3) +
    scale_fill_manual(values = plot_fill_colors) +
    labs(
      title = paste('TB case notifications by district, year, and type:', adm1_name), 
      subtitle = 'Across all age groups',
      x = "Year", y = "Number of notified cases", fill = "Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pdf(file.path(viz_dir, 'notifications_by_type.pdf'), height = 10, width = 10)
lapply(adm1_plots, print) |> invisible()
dev.off()
