#Load libraries
libs <- c("shiny","dplyr","tidyr","ggplot2","plotly",
          "shinycssloaders","RColorBrewer")

lapply(libs, library, character.only = TRUE)

source("functions.R")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Allele dynamics within populations"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          selectInput("sc", "How do the alleles interact?",
                      c("The new allele is recessive", "The new allele is dominanat",
                        "The alleles are codominant", "Overderdominance of the hetrozygot",
                        "Underdominance of heterozygot")),  
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
                        step = 0.01,
                        value = 0.01),
          uiOutput("ht"),
          uiOutput("ud")
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
  
  
    output$ht <- renderUI({ if(input$sc == "Overderdominance of the hetrozygot"){
      sliderInput("s","How successful is the heterozygot?",
                  min = 0.01,
                  max = 1,
                  step = 0.005,
                  value = 0.01)} else {
                    if(input$sc == "Underdominance of heterozygot"){
                      sliderInput("s", "How terrible is the heterozygot?",
                                  min = -1,
                                  max = -0.01,
                                  step = 0.005,
                                  value = -0.01)
                    }
                  } })
  

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
      
      #add if else statement to account for the need for t and s
      a <- input$bins
      sc <- input$sc
      ggplotly(p)
    
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
