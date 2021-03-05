#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Unlist app data
#for(objs in names(envdata)) {
#  assign(objs, envdata[[objs]])
#}

# Define UI for application
ui <- shinyUI(fluidPage(
  
  #Set style for vertical alignment
  tags$head(tags$style(
    HTML('
         #vert { 
            display: flex;
            align-items: center;
            margin-top: 50px;
         }
         '))),
  
   # Application title
   titlePanel(textOutput("title")),
   fluidRow(tags$hr(style="border-color: black;")),
   #Summary data
   fluidRow(column(4, align="left", h3(tags$b("Script options"))), column(8, align="center", h3(tags$b("Quality Control Plots")))),
   fluidRow(column(4, 
                   wellPanel(
                     uiOutput("input"),
                     uiOutput("output"),
                     uiOutput("yaml"),
                     uiOutput("minfeatures"),
                     uiOutput("maxfeatures"),
                     uiOutput("maxpercentmt"),
                     uiOutput("res")
                     ),
                   wellPanel(align="center", tableOutput("sum"))),
            column(8, align="right", id="vert", plotOutput("initQC", width="100%", height="600px"))),
  
   #PCA plots
   fluidRow(tags$hr(style="border-color: black;")),
   fluidRow(column(6, align="center", h3(tags$b("Selection of Principle Components"))), column(6, align="center", h3(tags$b("Plotting Principle Components")))),
   fluidRow(column(6, align="center", id="vert", plotOutput("npcs", width="80%", height="400px")), column(6, id="vert", align="center", plotOutput("pca", width="80%", height="400px"))),
   fluidRow(tags$hr(style="border-color: black;")),
  
    # Clustree plot
    fluidRow(column(12, align="center",  h3(tags$b("Selection of K from Clustree")))),
    fluidRow(column(12, align="center", plotOutput("clust", width="60%", height="800px"))),
    fluidRow(tags$hr(style="border-color: black;")),
  
    #Demultiplex plots
    fluidRow(column(5, align="center", h3(tags$b(textOutput("doubtitle")))), column(7, align="center", h3(tags$b(textOutput("hashtitle"))))),
    fluidRow(column(5, align="right", id="vert", plotOutput("doublets", width="80%", height="700px")), column(7, id="vert", align="right", plotOutput("hashtags", width="80%", height="700px")))
      
   ))

# Define server logic
server <- shinyServer(function(input, output, session) {
   
  #Set text outputs
  output$title <- renderText({paste("QC and Data Processing Report:", sample)})

  inshort <- ifelse(nchar(dirname(opt2$input)) > 50, 
                    substr(dirname(opt2$input), nchar(dirname(opt2$input))-50, nchar(dirname(opt2$input))),
                    dirname(opt2$input)) %>% sub(".*?/", "/", .)
  
  output$input <- renderText({HTML(paste("<b>","Input Directory:", "</b>", "...", inshort))})
  
  outshort <- ifelse(nchar(opt2$output) > 50, 
                     substr(opt2$output, nchar(opt2$output)-50, nchar(opt2$output)),
                     opt2$output) %>% sub(".*?/", "/", .)
  
  output$output <- renderText({HTML(paste("<b>", "Output Directory:", "</b>", "...", outshort))})
  
  #Set QC plots 
  output$initQC <- renderPlot({
  qc2 <- qc2 + theme(legend.position = "none", title = element_text(size=10)) + labs(title = "Percentage mitochondrial features vs number of molecules detected")
  qc3 <- qc3 + theme(legend.position = "none", title = element_text(size=10)) + labs(title = "Number of unique features vs number of molecules detected")
  plot <- grid.arrange(arrangeGrob(qc1.1, qc1.2, qc1.3, ncol=3), arrangeGrob(qc2, qc3, ncol=2), heights=c(2.5/4, 1.5/4), ncol=1)
  plot})
  
  #Set pre/post filtering table
  row.names(data.meta.summ) <- c(
    "<b>Number of cells</b>",
    "<b>Median nCount_RNA</b>",
    "<b>Minimum nCount_RNA</b>",
    "<b>Maximum nCount_RNA</b>",
    "<b>Median nFeature_RNA</b>",
    "<b>Minimum nFeature_RNA</b>",
    "<b>Maximum nFeature_RNA</b>",
    "<b>Median percent.mt</b>",
    "<b>Minimum percent.mt</b>",
    "<b>Maximum percent.mt</b>"
  )
  colnames(data.meta.summ)[1] <- "Post-filtering"
  output$sum <- renderTable(data.meta.summ, spacing = "l", rownames=TRUE, digits=0, hover=TRUE, sanitize.text.function=function(x){x})
  
  #Set npcs plot
  output$npcs <- renderPlot(npcs$plot + theme(text = element_text(size=20)))
  
  #Set PCA and feature plot event
  output$pca <- renderPlot(pca + theme(text = element_text(size=20), legend.position = "none"))
  
  # Plot clustree data
  output$clust <- renderPlot(clust)
  
})
