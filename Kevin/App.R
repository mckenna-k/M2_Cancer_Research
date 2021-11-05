library(shiny)
library(shinydashboard)
library(tidyverse)
library(survival)
library(survminer)

data = read.csv("../data/TCGA-CDR_data.csv")

data = data %>% na_if('#N/A') %>% na_if('[Not Available]') %>% na_if('[Not Applicable]')%>% na_if('[Unknown]')%>% na_if('[Not Evaluated]')
data = data %>% mutate_at(c(2,4,13,16,17,22,26:33),as.character)
data = data %>% mutate_at(c(4,13,16,17,22,26:33),as.numeric)
data_reduc = data[c(1,3:6,9,11,13,16,17,26:33)]
data_arranged = arrange(data_reduc, OS.time)
data_arranged = data_arranged %>% mutate_at(11,as.numeric)

ui <- dashboardPage(
  dashboardHeader(title = 'Cancer trends'),
  dashboardSidebar(
    
  ),
  dashboardBody(fluidRow(
    box(status = 'primary', uiOutput('Cancer'), width = 6),
    box(status = 'primary', uiOutput('Stats'), width = 6),
    box(status = 'primary', uiOutput('Years'), width = 6),
    fluidRow(box(plotOutput('plot'), width = 12,
    )),
    fluidRow(box(plotOutput('plot2'), width = 12))
    
  ))
)

server <- function(input, output){
  
  output$plot <- renderPlot({
    y = Surv(data_arranged$DSS.time[data_arranged$type==input$cancers], event = data_arranged$DSS[data_arranged$type==input$cancers])
    z = Surv(data_arranged$PFI.time[data_arranged$type==input$cancers], event = data_arranged$PFI[data_arranged$type==input$cancers])
    DSS.fit <- survfit(y~1, conf.type = 'plain')
    PFI.fit <- survfit(z~1, conf.type = 'plain')
    check = list(DSS = DSS.fit, PFI = PFI.fit)
    ggsurvplot(check, data = data_arranged, combine = T, conf.int = T, xscale = 365.25, xlim = c(0,input$years*365.25))
      
    
  })
  
  output$plot2 <- renderPlot({
    ggplot(data_arranged)+
    geom_histogram(aes(input$stat))
  })
  output$Cancer <- renderUI({
    selectInput(
      'cancers',
      label = h3('Choix de type de Cancer'), 
      choices = unique(factor(data$type))
    )
  })
  output$Stats <- renderUI({
    selectInput(
      'stat',
      label = h3('Choix de type de Cancer'), 
      choices = colnames(data_arranged)
    )
  })
  output$Years <- renderUI({
    numericInput(
      "years", label = h3("Nombre d'annees"), value = 10)
  })
  
}
shinyApp(ui,server)
