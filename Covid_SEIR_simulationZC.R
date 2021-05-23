#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Covid SEIR interactive Simulation"),
   
   # Sidebar with a slider input for beta 
   sidebarLayout(
      sidebarPanel(
         sliderInput("betas",
                     "Value of Beta:",
                     min = 0.1,
                     max = 0.3,
                     value = 0.22)
      ),
      
      # Show a plot of the plot of simulation
      mainPanel(
         plotOutput("seirPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  library(dplyr)
  library(lubridate)
  library(readr)
  library(deSolve)
  library(rstudioapi)
  directory <- dirname(rstudioapi::getSourceEditorContext()$path) # doesn't work if ran in Terminal
  
  cases <- read.csv(file.path(directory, "data/cases_by_status_and_phu.csv"))
  cases <- cases %>% filter(PHU_NAME == "TORONTO")
  cases_sort <- cases[order(cases$FILE_DATE),]
    
  infectious_period <- 10 # from https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
  latent_period <- 3 #https://link.springer.com/article/10.1007/s10237-020-01332-5
  #recovery_time <- 14
  
  #Data as of 7-06
  active <- 570
  resolved <- 11791
  deaths <- 1070
  N <- 2.3 * 1000000 #population of toronto
  S <- N - active - resolved - deaths # number of susceptible individuals 
  D <- deaths
  
  #calculating beta
  #kirkby recommends using sampling time < time it takes one individual to infect another
  #calculation for transmission rate by https://www.nature.com/articles/s41598-017-09209-x
  
  t <- 7 #sampling interval 
  Inew <- 388 #number of new infections since last sample (from 06-30 to 07-06)
  I <- 570 #number of infected individuals
  
  #parameters
  gamma <- 1/infectious_period
  delta <- 1/latent_period
  
  #Ro <-  beta/gamma
  
  new_case <- read.csv(file.path(directory, "data/COVID19.csv"))
  
  #code for SEIR
  seir <- function(beta, gamma, delta, S, I, resolved, N) {
    parameter_list <- c(beta, gamma, delta)
    W <- S
    X <- I
    Y <- resolved
    Z <- N - W - X - Y - deaths 
    D <- D/N
    
    seir_model = function (current_timepoint, state_values, parameters)
    {
      # create state variables (local variables)
      S = state_values [1]        # susceptibles
      E = state_values [2]        # exposed
      I = state_values [3]        # infectious
      R = state_values [4]        # recovered
      
      with ( 
        as.list (parameters),     # variable names within parameters can be used 
        {
          # compute derivatives
          dS = (-beta * S * I) 
          dE = (beta * S * I) - (delta * E) 
          dI = (delta * E) - (gamma * I) 
          dR = (gamma * I) 
          
          # combine results
          results = c (dS, dE, dI, dR)
          list (results)
        }
      )
    }
    
    initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)
    timepoints = seq (0, 200, by=1)
    output = lsoda(initial_values, timepoints, seir_model, parameter_list)
    return(output)
  }
  
  output$seirPlot <- renderPlot({
    output1 <- reactive(seir(input$betas, gamma, delta, S, I, resolved, N))
    matplot(output1()[,'time'], output1()[,c('S','E','I','R')], 
            col=c('blue','pink','red','green'), type = 'l', lty=2,
            xlab = 'Time', ylab = 'S, E, I, R')
    legend(140, 1, legend = c("Susceptible", "Exposed", "Infected", "Recovered"), 
           col = c("blue", "pink", "red", "green"), lty =1, cex= 0.8)    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

