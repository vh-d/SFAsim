library(shiny)
library(ggplot2)
library(data.table)

# Alt 1
shinyApp(
  ui = fluidPage(
    sidebarLayout(
      sidebarPanel = 
        sidebarPanel(
          numericInput("length", "length", 10),
          numericInput("ar_coef", "AR(1) coeff",  0.8, step = 0.1),
          numericInput("ar_const", "AR(1) const", 0.2, step = 0.1),
          numericInput("ar_scale", "scale",       0.05, step = 0.01)
        ), 
      mainPanel = 
        mainPanel(
          plotOutput("plot"),
          plotOutput("plot2")
        )
    )
  ),
  server = function(input, output) {
    ar_values <- 
      reactive(
        sim_ar_gamma(input$length, 
                     ar_coef  = input$ar_coef,
                     ar_const = input$ar_const,
                     ar_scale = input$ar_scale)
      )
    
    output$plot <- 
      renderPlot(
        plot(
          ar_values(), 
          main = "AR(1) gamma process",
          ylab = "value",
          ylim = c(0, 3), 
          type = "l"))
    
    output$plot2 <- 
      renderPlot(
        plot(
          1 - exp(-ar_values()), 
          main = "inefficiency",
          ylab = "value",
          ylim = c(0, 1), 
          type = "l"))
  }
)



# Alt 2 - N simulations
shinyApp(
  ui = fluidPage(
    sidebarLayout(
      sidebarPanel = 
        sidebarPanel(
          selectInput("dist", "Distribution", choices = c("sim_ar_gamma", "sim_ar_lnorm", "sim_ar_norm")),
          numericInput("N", "N", 100),
          numericInput("length", "Length", 10),
          numericInput("ar_coef", "AR(1) coeff",  0.8,  step = 0.1),
          numericInput("ar_const", "AR(1) const", 0.2,  step = 0.1),
          numericInput("ar_scale", "scale",       0.05, step = 0.01)
        ), 
      mainPanel = 
        mainPanel(
          plotOutput("plot"),
          plotOutput("plot2")
        )
    )
  ),
  server = function(input, output) {
    output$plot <- 
      renderPlot(
        ggplot(data = 
                 rbindlist(
                   lapply(1:input$N, 
                          function(X) 
                            data.table(
                              i = X,
                              x = 1:input$length, 
                              y = suppressWarnings(
                                do.call(input$dist,
                                        args = list(
                                          l = input$length, 
                                          ar_coef  = input$ar_coef,
                                          ar_const = input$ar_const,
                                          input$ar_scale)))))),
               aes(x = x, 
                   y = y, 
                   group = i)) +
          geom_line() +
          ylim(c(-1, 5))
      )
  }
)
# rbindlist(lapply(X = 1:100, FUN = function(x) data.table(x = 1:10, y = sim_ar_gamma())))
