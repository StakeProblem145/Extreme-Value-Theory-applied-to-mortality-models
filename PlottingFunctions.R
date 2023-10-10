library(tidyverse)
library(ggthemes)


MODEL_PLOT_THEME <- list(
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = "bottom",
    legend.background = element_rect(color = NA)
  )
)


plotLogMortality <-
  function(dataset,
           title = "MISSING TITLE",
           xlim = c(70, 110),
           ylim = c(-4.5, 0.5),
           save = FALSE,
           errorBar = TRUE) {
    plot <- ggplot(dataset$plottingData, aes(x = Age)) +
      geom_line(aes(y = mod), color = "#F8766D", lwd=1) +
      geom_point(aes(y = obs), alpha = 1 / 3) +
      xlim(xlim) +
      ylim(ylim) +
      #geom_col(width=0.8)+
      ggtitle(paste(title)) +
      xlab("Age") +
      ylab("log(m_x)") +
      theme_clean() +
      MODEL_PLOT_THEME
    
    if(errorBar){
      plot <- plot +
        geom_errorbar(aes(ymin = dataset$plottingData$lower, ymax = dataset$plottingData$upper)
                      , width=.3, alpha = 1/2)
    }
    
    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }
    
    plot
  }



plotLogMortalitySplitData <-
  function(dataset,
           title = "MISSING TITLE",
           xlim = c(70, 110),
           ylim = c(-4.5, 0.5),
           hideLegend = TRUE,
           save = FALSE,
           errorBar = TRUE) {
    if(is.null(dataset$plottingData)) {
      W <- dataset
    } else {
      W <- dataset$plottingData
    }
    A <- filter(W, fitted == "Fit")
    B <- filter(W, fitted == "Extrapolation")
    plot <- ggplot() +
      geom_point(aes(x = W$Age, y = W$obs), alpha = 1 / 3) +
      geom_line(aes(x = B$Age, y = B$mod, colour = B$model), linetype = "dashed", lwd=1) +
      geom_line(aes(x = A$Age, y = A$mod, colour = A$model), lwd=1) +
      xlim(xlim) +
      ylim(ylim) +
      #geom_col(width=0.8)+
      ggtitle(paste(title)) +
      xlab("Age") +
      ylab("log(m_x)") +
      labs(color = "Model") +
      theme_clean() +
      MODEL_PLOT_THEME
    
    if(errorBar){
      plot <- plot +
        geom_errorbar(aes(x = A$Age, y = A$mod, ymin = A$lower, ymax = A$upper)
                      , width=.3, alpha = 1/2) +
        geom_errorbar(aes(x = B$Age, y = B$mod, ymin = B$lower, ymax = B$upper)
                      , width=.3, alpha = 1/2)
    }
    
    
    if(hideLegend) {
      plot <- plot +
        theme(legend.position="none")
    }
    
    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }
    
    plot
  }

plotLogMortalitySplitDataWithVLine <-
  function(dataset,
           title = "MISSING TITLE",
           xlim = c(70, 110),
           ylim = c(-4.5, 0.5),
           hideLegend = TRUE,
           save = FALSE,
           errorBar = TRUE,
           vline = c()) {
    if(is.null(dataset$plottingData)) {
      W <- dataset
    } else {
      W <- dataset$plottingData
    }
    A <- filter(W, fitted == "Fit")
    B <- filter(W, fitted == "Extrapolation")
    plot <- ggplot() +
      geom_point(aes(x = W$Age, y = W$obs), alpha = 1 / 3) +
      geom_line(aes(x = B$Age, y = B$mod, colour = B$model), linetype = "dashed", lwd=1) +
      geom_line(aes(x = A$Age, y = A$mod, colour = A$model), lwd=1) +
      geom_vline(xintercept = vline, alpha = 1 / 2, linetype = "dashed") +
      xlim(xlim) +
      ylim(ylim) +
      #geom_col(width=0.8)+
      ggtitle(paste(title)) +
      xlab("Age") +
      ylab("log(m_x)") +
      labs(color = "Model") +
      theme_clean() +
      MODEL_PLOT_THEME

    if(errorBar){
      plot <- plot +
        geom_errorbar(aes(x = A$Age, y = A$mod, ymin = A$lower, ymax = A$upper)
          , width=.3, alpha = 1/2) +
        geom_errorbar(aes(x = B$Age, y = B$mod, ymin = B$lower, ymax = B$upper)
          , width=.3, alpha = 1/2)
    }


    if(hideLegend) {
      plot <- plot +
        theme(legend.position="none")
    }

    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }

    plot
  }



plotLogMortalityExtrapolatedSplitData <-
  function(dataset,
           title = "MISSING TITLE",
           xlim = c(70, 110),
           ylim = c(-4.5, 0.5),
           hideLegend = TRUE,
           save = FALSE,
           errorBar = TRUE) {
    if(is.null(dataset$extrapolatedData)) {
      W <- dataset
    } else {
      W <- dataset$extrapolatedData
    }
    A <- filter(W, fitted == "Fit")
    B <- filter(W, fitted == "Extrapolation")
    plot <- ggplot() +
      geom_point(aes(x = W$Age, y = W$obs), alpha = 1 / 3) +
      geom_line(aes(x = B$Age, y = B$mod, colour = B$model), linetype = "dashed", lwd=1) +
      geom_line(aes(x = A$Age, y = A$mod, colour = A$model), lwd=1) +
      xlim(xlim) +
      ylim(ylim) +
      #geom_col(width=0.8)+
      ggtitle(paste(title)) +
      xlab("Age") +
      ylab("log(m_x)") +
      labs(color = "Model") +
      theme_clean() +
      MODEL_PLOT_THEME

    if(errorBar){
      plot <- plot +
        geom_errorbar(aes(x = A$Age, y = A$mod, ymin = A$lower, ymax = A$upper)
          , width=.3, alpha = 1/2) +
        geom_errorbar(aes(x = B$Age, y = B$mod, ymin = B$lower, ymax = B$upper)
          , width=.3, alpha = 1/2)
    }


    if(hideLegend) {
      plot <- plot +
        theme(legend.position="none")
    }

    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }

    plot
  }

plotResMortality <-
  function(dataset,
           title = "MISSING TITLE",
           xlim = c(70, 110),
           ylim = 1,
           save = FALSE) {
    plot <- ggplot(dataset$residualData, aes(x = Age)) +
      geom_point(aes(y = devRes)) +
      xlim(xlim) +
      ggtitle(paste("Residuals", title)) +
      xlab("Age") +
      ylab("Residuals") +
      theme_clean() +
      MODEL_PLOT_THEME
    if (ylim == 1) {
      plot <- plot +
        ylim(-5, 5) +
          geom_hline(yintercept = -1.96,
                     linetype = "dashed",
                     color = "red") +
          geom_hline(yintercept = 1.96,
                     linetype = "dashed",
                     color = "red")
    }
    if (ylim == 2) {
      plot <- plot +
        ylim(-10, 10) +
          geom_hline(yintercept = -1.96,
                     linetype = "dashed",
                     color = "red") +
          geom_hline(yintercept = 1.96,
                     linetype = "dashed",
                     color = "red") 
          # geom_hline(yintercept = -5,
          #            linetype = "dotted",
          #            color = "red") +
          # geom_hline(yintercept = 5,
          #            linetype = "dotted",
          #            color = "red")
    }
    if (ylim == 3) {
      plot <- plot +
        geom_hline(yintercept = -1.96,
                   linetype = "dashed",
                   color = "red") +
        geom_hline(yintercept = 1.96,
                   linetype = "dashed",
                   color = "red") 
      # geom_hline(yintercept = -5,
      #            linetype = "dotted",
      #            color = "red") +
      # geom_hline(yintercept = 5,
      #            linetype = "dotted",
      #            color = "red")
    }
    
    
    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }
    
    plot
  }



plotParameter <- 
  function(parameterDataFrame,
           parameter,
           title = "MISSING TITLE",
           # xlim = c(1890, 1910),
           ylim = NULL,
           errorBar = TRUE,
           hideLegend = TRUE,
           save = FALSE) {
    W <- filter(parameterDataFrame, Parameter == parameter)
    plot <- ggplot(W, aes(x = factor(Cohort))) +
      geom_point(aes(y = Value)) +
      scale_x_discrete() +
      # xlim(xlim) +
      ggtitle(paste(title)) +
      xlab("Cohort") +
      ylab("Value") +
      theme_clean() +
      MODEL_PLOT_THEME
    
    if(!is.null(ylim)){
      plot <- plot + 
        ylim(ylim)
    }
    
    if(hideLegend) {
      plot <- plot +
        theme(legend.position="none")
    }
    
    if(errorBar){
      plot <- plot +
        geom_errorbar(aes(ymin = Value - SD,
                          ymax = Value + SD)
                      , width=.3)
    }
    
    if(save) {
      title <- gsub(" ", "", title, fixed = TRUE)
      plotname <- paste(title, ".png", sep = "")
      ggsave(plotname, plot = plot, width = 30, height = 20, units = "cm", path = "plots")
    }
    
    plot
  }


