#' Plotting the classical RSF UD
#' @name plotCovariates
#'
#' @param covariates_list covariates list, that must be a list of "raster-like"
#' elements. It is for the moment supposed that the grid is regular and
#' shared by all elements of covariates
#' @param trajectory_data An optional argument giving the trajectory data to supperpose.
#' @importFrom rlang .data
#' @return a ggplot graph
#' @export
plotCovariates <- function(covariates_list, trajectory_data = NULL){
  J <- length(covariates_list)
  covariates_df <- purrr::map_dfr(covariates_list, # To each element of the list, apply
                           Rhabit::rasterToDataFrame, # This function
                           .id = "Covariate") %>% # Keep origin in a Covariate column
    dplyr::mutate(Covariate = factor(.data$Covariate, labels = paste("Cov.", 1:J)))
  # First, for one covariate
  if(!is.null(trajectory_data)){
    if(is.null(trajectory_data$ID)){
      trajectory_data$ID = rep("1", nrow(trajectory_data))
    }
  }
  plot_covariate <- function(df, cov_name, plot_traj = FALSE) {
    main_plot <- ggplot2::ggplot(data = df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$val)) +
      ggplot2::scale_fill_viridis_c(name = .data$cov_name) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(title = cov_name, x = "X-axis", y = "Y-axis") +
      ggplot2::theme(legend.position = "bottom")
    if(plot_traj){
      return(main_plot +     
               ggplot2::geom_path(data = trajectory_data, 
                                  ggplot2::aes(group = .data$ID)))
    }
    else{
      return(main_plot)
    }
  }
  covariates_df %>% # Take the data.frame
    dplyr::group_by(.data$Covariate) %>% # Group by Covariate
    tidyr::nest() %>% # Make a tibble, the second "column" gathers the data
    dplyr::mutate(plots = purrr::map2(.data, .data$Covariate, plot_covariate,
                        plot_traj = !is.null(trajectory_data))) %>% 
    gridExtra::grid.arrange(grobs = .data$plots, nrow = 1)
}
