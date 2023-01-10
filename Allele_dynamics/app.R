
library("shiny")
library("dplyr")
library("tidyr")
library("ggplot2")
library("plotly")
library("shinycssloaders")

source("./Allele_dynamics/functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Allele dynamics within populations"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          selectInput("sc", "How do the alleles interact?",
                      c("The new allele is recessive", "The new allele is dominanat",
                        "The alleles are codominant", "Over/underdominance of the hetrozygot")),  
          sliderInput("bins",
                        "Proportion of new allele:",
                        min = 1,
                        max = 50,
                        value = 30),
            uiOutput("Old"),
            sliderInput("psize", "Population size:",
                        min =  1000,
                        max = 100000,
                        step = 100,
                        value = 10000),
            sliderInput("t", "How sucessful is the new allele?",
                        min =  -1,
                        max = 1,
                        step = 0.001,
                        value = 0.001),
          uiOutput("ht")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           plotOutput("popPlot"),
           plotOutput("allelePlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  gen <- 3000
  
  output$Old <- renderUI({
    sliderInput("Old", "Proportions of old allele", min = 0,
                max = 1,
                value = (1-(input$bins/100)), step = 0.01)
  })
  
  
    output$ht <- renderUI({ if(input$sc == "Over/underdominance of the hetrozygot"){
      sliderInput("s","How successful is the heterozygot?",
                  min = -1,
                  max = 1,
                  step = 0.001,
                  value = 0.001)} })
  

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
    
    output$popPlot <- renderPlotly({
      a <- input$bins
      sc <- input$sc
      ggplotly(p)
    
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
