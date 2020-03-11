#' Plotting the classical RSF UD
#' @name plotCovariates
#'
#' @param covariates covariates list, that must be a list of "raster-like"
#' elements. It is for the moment supposed that the grid is regular and
#' shared by all elements of covariates
#' @param trajectory_data An optional argument giving the trajectory data to supperpose.
#' @return a ggplot
#' @export
plotCovariates <- function(covariates_list, trajectory_data = NULL){
  J <- length(covariates_list)
  covariates_df <- map_dfr(covariates_list, # To each element of the list, apply
                           Rhabit::rasterToDataFrame, # This function
                           .id = "Covariate") %>% # Keep origin in a Covariate column
    mutate(Covariate = factor(Covariate, labels = paste("Cov.", 1:J)))
  # First, for one covariate
  if(!is.null(trajectory_data)){
    if(is.null(trajectory_data$ID)){
      trajectory_data$ID = rep("1", nrow(trajectory_data))
    }
  }
  plot_covariate <- function(df, cov_name, plot_traj = FALSE) {
    main_plot <-ggplot(data = df, aes(x = x, y = y)) +
      geom_raster(aes(fill = val)) +
      scale_fill_viridis_c(name = cov_name) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = cov_name, x = "X-axis", y = "Y-axis") +
      theme(legend.position = "bottom")
    if(plot_traj){
      return(main_plot +     
               geom_path(data = trajectory_data, aes(group = ID)))
    }
    else{
      return(main_plot)
    }
  }
  covariates_df %>% # Take the data.frame
    group_by(Covariate) %>% # Group by Covariate
    nest() %>% # Make a tibble, the second "column" gathers the data
    mutate(plots = map2(data, Covariate, plot_covariate,
                        plot_traj = !is.null(trajectory_data))) %>% 
    gridExtra::grid.arrange(grobs = .$plots, nrow = 1)
}
