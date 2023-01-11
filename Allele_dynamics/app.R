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
                        "The alleles are codominant", "Overderdominance of the heterozygot",
                        "Underdominance of the heterozygot")),  
          
          sliderInput("bins",
                      "Proportion of new allele:",
                      min = 0.1,
                      max = 0.9,
                      value = 0.1,
                      step = 0.01),
            uiOutput("Old"),
            sliderInput("psize", "Population size:",
                        min =  1000,
                        max = 1000000,
                        step = 100,
                        value = 10000),
            sliderInput("t", "How sucessful is the new allele?",
                        min =  -1,
                        max = 1,
                        step = 0.005,
                        value = 0.01),
          uiOutput("ht"),
          uiOutput("ud")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          shinycssloaders::withSpinner(plotlyOutput("allelePlot"), 
                                       size = 2, type = 6),
          shinycssloaders::withSpinner(plotlyOutput("popPlot"),
                                       size = 2, type = 6),
          shinycssloaders::withSpinner(plotlyOutput("pop_sizePlot"),
                                       size = 2, type = 6)
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  genr <- 5000
  
  output$Old <- renderUI({
    sliderInput("Old", "Proportions of old allele", min = 0.1,
                max = 0.9,
                value = (1-(input$bins)), step = 0.01)
  })
  
  
    output$ht <- renderUI({ if(input$sc == "Overderdominance of the heterozygot"){
      sliderInput("s","How successful is the heterozygot?",
                  min = 0.01,
                  max = 1,
                  step = 0.005,
                  value = 0.015)} else {
                    if(input$sc == "Underdominance of the heterozygot"){
                      sliderInput("s", "How terrible is the heterozygot?",
                                  min = -1,
                                  max = -0.01,
                                  step = 0.005,
                                  value = -0.01)
                    }
                  } })
  

    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
    
    output$allelePlot <- renderPlotly({
      
      #add if else statement to account for the need for t and s
      
      if(!grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, ps = input$psize)
        
      }else{if(grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, sd = input$s, ps = input$psize)
        }
        
      }
      Sys.sleep(0.5)
      ggplotly(dl[[4]])
      
    
    })
    
    output$popPlot <- renderPlotly({
      
      #add if else statement to account for the need for t and s
      
      if(!grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, ps = input$psize)
        
      }else{if(grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, sd = input$s, ps = input$psize)
      }
        
      }
      Sys.sleep(0.5)
      ggplotly(dl[[5]])
      
      
    })
    
    output$pop_sizePlot <- renderPlotly({
      
      #add if else statement to account for the need for t and s
      
      if(!grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, ps = input$psize)
        
      }else{if(grepl("heterozygot", input$sc)){
        
        dl <- dnp(qd = input$bins, td = input$t,
                  scd = input$sc, g = genr, sd = input$s, ps = input$psize)
      }
        
      }
      Sys.sleep(0.5)
      ggplotly(dl[[6]])
      
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
