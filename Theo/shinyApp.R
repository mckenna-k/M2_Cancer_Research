library(shiny)
library(shinydashboard)
library(survival)

data = read.csv("data/TCGA-CDR_data.csv")

data <- data[,c(-1,-2)]
quant <- c("age_at_initial_pathologic_diagnosis","initial_pathologic_dx_year","birth_days_to","last_contact_days_to","death_days_to","new_tumor_event_dx_days_to","OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")
for( i in quant){
  data[[i]] <- as.numeric(data[[i]])
}
cat <- colnames(data[,!colnames(data) %in% quant])
data[data=="[Discrepancy]" | data=="[Not Applicable]" | data=="[Not Available]" | data=="[Unknown]" | data=="#N/A" | data=="=#N/D"] <- "NA"
data$ajcc_pathologic_tumor_stage <- as.factor(data$ajcc_pathologic_tumor_stage)
levels(data$ajcc_pathologic_tumor_stage) <- c("I - II","I - II","I - II","I - II","III - IV","I - II","III - IV","III - IV","NA","other","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","I - II","III - IV","III - IV","III - IV","III - IV","III - IV","III - IV","I - II","III - IV","III - IV","III - IV","III - IV","other")

ui <- dashboardPage(
  dashboardHeader(title = 'Cancer trends'),
  dashboardSidebar(
    selectInput("cancerType","Cancer types",choices=unique(data$type),multiple = TRUE),
    selectInput("endType","Endpoints",choices=c("OS"=1,"DSS"=2,"PFI"=3,"DFI"=4),multiple = TRUE),
    selectInput("Covariate","Covariate",choices=cat),
    numericInput("years","Nombre d'annees", value = 10),
    checkboxInput("invert","Invert covariate and endpoints",value=FALSE)
  ),
  dashboardBody(
    fluidRow(
      uiOutput("survPlotUI")
    )
  )
)

server <- function(input, output){

  output$survPlotUI <- renderUI({
    req(input$cancerType)
    plotOutput('survPlot',height = paste0(400+length(input$cancerType)*200,"px"))
  })

  output$survPlot <- renderPlot({
    req(input$cancerType,input$endType,input$Covariate,input$years)
    df <- construct_data_frame(data,input$cancerType,as.numeric(input$endType),data[[input$Covariate]])

    if(input$invert){
      ggplot(df)+
        geom_step(aes(x=time,y=survVal,color=end))+
        geom_step(aes(x=time,y=survValUp,color=end),linetype="dashed")+
        geom_step(aes(x=time,y=survValDo,color=end),linetype="dashed")+
        facet_grid(rows = vars(type),cols=vars(covar))+
        xlab("time (years)")+
        ylab("probability")+
        ggtitle("Courbe de survie en fonction du type de cancer et du type d'évènement")+
        coord_cartesian(ylim=c(0,1),xlim=c(0,input$years))+
        theme_minimal()
    }
    else{
      ggplot(df)+
        geom_step(aes(x=time,y=survVal,color=covar))+
        geom_step(aes(x=time,y=survValUp,color=covar),linetype="dashed")+
        geom_step(aes(x=time,y=survValDo,color=covar),linetype="dashed")+
        facet_grid(rows = vars(type),cols=vars(end))+
        xlab("time (years)")+
        ylab("probability")+
        ggtitle("Courbe de survie en fonction du type de cancer et du type d'évènement")+
        coord_cartesian(ylim=c(0,1),xlim=c(0,input$years))+
        theme_minimal()
    }
  })

}
shinyApp(ui,server)

