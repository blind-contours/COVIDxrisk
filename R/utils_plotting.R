
COVIDxRisk_plot <- function(data, type, thresh, y_label, x_label, path){
  risk_plot <- risk_importance %>%
    arrange(Est) %>%
    # First sort by val. This sort the dataframe but NOT the factor levels
    filter(Est > 1.0) %>%
    mutate(name = factor(Variable, levels = Variable)) %>%
    # This trick update the factor levels
    ggplot(aes(x = name, y = Est)) +
    geom_point(size = 4, color = "orange") +
    geom_errorbar(aes(ymin = `Lower_CI`, ymax = `Upper_CI`), width = .1) +
    coord_flip() +
    theme_bw(base_size = 12) +
    ylab("Model Risk Ratio") +
    xlab("County Feature")

  ggsave(here(paste("Figures/", "varimp_", label, ".png", sep = "")), risk_plot, width = 8, height = 6)
}
