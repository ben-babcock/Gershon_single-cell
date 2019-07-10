#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(ggplot2)
library(utils)
library(base)
fonts <- list("Liberation Sans")

capgene <- function(gene){
  paste0(toupper(substring(gene, 1, 1)), substring(gene, 2))
}

simpleFeaturePlot <- function(gene = NULL, 
                          color.use = 'red',
                          threshold = 0,
                          plot.intersection = F,
                          pt.size = 0.5,
                          alpha.colors = 0.9,
                          base.transparency = 0.7,
                          annotate = T,
                          to.include = NULL,
                          inc.veh = F,
                          inc.vis = F,
                          data.set = NULL,
                          scale = 900){
  print(scale)
  if(data.set == "P15"){
  fh <- as.character(GTfh.dict[which(GTfh.dict$GENES == capgene(gene)), "FILE"])
  anno1 <-  annotate("text", size = 7*(scale/1000),
             x = c(42, 0, 40, 26, -34, -41), 
             y = c (9, 42, 23, -30, -23, 0), 
             label = c("Astrocytes", "Endothelial Cells", "Vascular\nFibroblasts",
                       "Mature\nNeurons", "Microglia", "Oligodendrocytes"))
  anno2 <-  annotate("text", size = 10*(scale/1000),
             x = c(-5, -15, 20, 26), 
             y = c (-20, 12, 20, -15), 
             label = c("Node A", "Node B", "Node C", "Node D"))
  anno3 <- NULL
  if('inc.veh' %in% to.include){
    inc.veh <- T
  }
  if('inc.vis' %in% to.include){
    inc.vis <- T
  }
  if(inc.veh && inc.vis){
    frame <- rbind(veh, vis)
  }else if(inc.veh){
    frame <- veh
  }else if(inc.vis){
    frame <- vis
  }else{d <- ggplot() + theme_classic()
  return(d)}
  }else if(data.set == "P7"){
    fh <- as.character(WTfh.dict[which(WTfh.dict$GENES == capgene(gene)), "FILE"])
    anno1 <- annotate("text", size = 7*(scale/1000),
                      x = c(40, -42, -42, -42, 25, 37), 
                      y = c (-1, -10, 23, 5, 40, 20), 
                      label = c("Astrocytes", "Endothelial\nCells", "Vascular\nFibroblasts",
                                "Purkinje\nNeurons", "Microglia", "Oligodendrocytes"))
    anno2 <- annotate("text", size = 10*(scale/1000),
                      x = c(-25, 5), 
                      y = c (0, 4), 
                      label = c("Interneuron\nPrecursors", "Granule\nNeuron\nPrecursors"))
    anno3 <- annotate("text", size = 5*(scale/1000),
             x = c(-27, -25, 10, 8), 
             y = c (-10, 10, -16, 25), 
             label = c("Proliferating", "Differentiating", "Proliferating", "Differentiating"))
  }
  
  colorcode <- col2rgb(color.use)
  colorcode <- rgb(colorcode[1,], colorcode[2,], colorcode[3,], max=255, alpha = (255*alpha.colors))
  #Get gene data
  
  data <- readRDS(fh)
  colors <- c(grey(level = 0.7, alpha = base.transparency), colorcode)
  names(colors) <- c('None', paste0(gene))
  #Get Plotting Coordinates
  frame <- cbind(frame,
                 as.data.frame(matrix(data = 0, 
                                      nrow = nrow(frame),
                                      ncol = 3)))
  colnames(frame) <- c('tSNE_1', 'tSNE_2', gene, 'Plot.Status', 'is.null')
  
  positive_cells <- colnames(data[, data[gene, ] > threshold])
  frame[which(rownames(frame) %in% positive_cells), gene] <- '1'
  
  
  #Plot Intersections: Cells that have multiple of genes vector
  
  frame[which(frame[, gene] > 0), 'Plot.Status'] <- paste0(gene)
  frame[which(frame[, 'Plot.Status'] == 0), 'Plot.Status'] <- 'None'
  frame1 <- frame[which(frame[, 'Plot.Status'] == 'None'), ]
  frame2 <- frame[which(frame[, 'Plot.Status'] != 'None'), ]
  if(annotate){
    d <- ggplot(mapping = aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(data = frame1, color = colors[1], size = pt.size) +
      geom_point(data = frame2, aes(color = Plot.Status), size = pt.size, show.legend = F) +
      scale_colour_manual(values = colors) + theme_classic() + anno1 + anno2 + anno3 +
      annotate("text", size = 6, color = color.use, x = -40, y = 40, label = gene)
    return(d)
  }else{
    d <- ggplot(mapping = aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(data = frame1, color = colors[1], size = pt.size) +
      geom_point(data = frame2, aes(color = Plot.Status), size = pt.size, show.legend = F) +
      scale_colour_manual(values = colors) + theme_classic() +
      annotate("text", size = 6, color = color.use, x = -40, y = 40, label = gene)
    return(d)
  }
  
}



#Load functions and data
veh <- readRDS('data/Veh_xy.rds')
vis <- readRDS('data/Vis_xy.rds')
frame <- readRDS('data/WT_xy.rds')
WTfh.dict <- readRDS('data/WT_expr.dict.rds')
GTfh.dict <- readRDS('data/GT_expr.dict.rds')

format.grep <- function(gg){
  grep(substr(as.character(gg), 1, 3), GTfh.dict$GENES, value = T, ignore.case=T)
}
#data <- read.table('~/downloads/WT_data.tsv', header = T, sep = '\t', data.table = F)
#frame <- read.table('~/downloads/WT_coords.tsv', header = T, sep = '\t', data.table = F)


# Define UI for application that draws a histogram
ui <- fluidPage(width = 8, theme = shinythemes::shinytheme("sandstone"),
  
  # Application title
  titlePanel("Gershon Lab Single-Cell Interactive Explorer"),
  #theme = shinythemes::shinytheme("readable"),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(div(class = "panel panel-danger",
                     style = "box-shadow: 5px 5px 15px -5px rgba(0, 0, 0, 0.3);",
                     div(class = "panel-heading", "Gene Controls"), 
                     div(class = "panel-body", textInput("gene", "Gene", value = 'Rbfox3'),
                         sliderInput("threshold", "Expression Threshold",
                                     min = 0.00,
                                     max = 5.00,
                                     step = 0.05,
                                     value = 1))),
                 div(class = "panel panel-danger",
                     style = "box-shadow: 5px 5px 15px -5px rgba(0, 0, 0, 0.3);",
                     div(class = "panel-heading", "Graphics Controls"), 
                     div(class = "panel-body", selectInput("color", "Gene Color:", 
                                                           c('Red' = "#FF0000FF", 'Orange' = "#FF9900FF", 
                                                             'Yellow' = "#FFCC00FF", "Green" = "#79CC3DFF", 
                                                             "Blue" = "#6699FFFF", "Purple" = "#CC33FFFF")),
                         sliderInput("point.size",
                                     "Point Size",
                                     min = 0.00,
                                     max = 2.00,
                                     step = 0.2,
                                     value = 0.5),
                         checkboxInput("annotate", "Annotate Cells", value = T)))
                 
                 ,
                 downloadButton('downloadPlot', "Download as pdf")
    ), 
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(id = "tab", type = "tabs",
                          tabPanel("P7 Cerebellum", value = "P7", plotOutput("WTfeaturePlot", width = '90%')),
                          tabPanel("P15 Tumor", value = "P15", plotOutput("GTfeaturePlot", width = '90%'), 
                                   fixedPanel(top = "120px", right = "100px", draggable = TRUE, 
                                              style = "width: 250px; z-index: 100000;", div(class = "panel panel-danger", 
                                                                                            style = "box-shadow: 5px 5px 15px -5px rgba(0, 0, 0, 0.3);", 
                                                                                            div(class = "panel-heading", "Treatments to Show:"),
                                                                                            checkboxGroupInput(inputId = "treatment", label ="",
                                                                                                               choiceNames = c("Vehicle", "Vismodegib"),
                                                                                                               choiceValues = c('inc.veh', 'inc.vis'),
                                                                                                               selected = c("inc.veh", "inc.vis"), inline = T))))
    )
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  WTplotInput = function(){
    simpleFeaturePlot(gene = capgene(input$gene), color.use = as.character(input$color), 
                  threshold = as.numeric(input$threshold), 
                  pt.size = as.numeric(input$point.size),
                  annotate = as.logical(input$annotate),
                  to.include = input$treatment,
                  data.set = input$tab,
                  scale = session$clientData$output_WTfeaturePlot_width)
  }
  GTplotInput = function(){
    simpleFeaturePlot(gene = capgene(input$gene), color.use = as.character(input$color), 
                      threshold = as.numeric(input$threshold), 
                      pt.size = as.numeric(input$point.size),
                      annotate = as.logical(input$annotate),
                      to.include = input$treatment,
                      data.set = input$tab,
                      scale = session$clientData$output_GTfeaturePlot_width)
  }
  plotInput = function(){
    simpleFeaturePlot(gene = capgene(input$gene), color.use = as.character(input$color), 
                      threshold = as.numeric(input$threshold), 
                      pt.size = as.numeric(input$point.size),
                      annotate = as.logical(input$annotate),
                      to.include = input$treatment,
                      data.set = input$tab)
  }
  output$WTfeaturePlot <- renderPlot({
    validate(
      need(expr = (capgene(input$gene) %in% format.grep(input$gene)), message = c(paste0(input$gene, ' not found!'),
                                                                                  'Did you mean: ',
                                                                                  paste(format.grep(input$gene))))
    )
    WTplotInput()
  }, height = function() {
    session$clientData$output_WTfeaturePlot_width*0.9
  })
  output$GTfeaturePlot <- renderPlot({
    validate(
      need(expr = (capgene(input$gene) %in% format.grep(input$gene)), message = c(paste0(input$gene, ' not found!'),
                                                                                  'Did you mean: ',
                                                                                  paste(format.grep(input$gene))))
    )
    GTplotInput()
  }, height = function() {
    session$clientData$output_GTfeaturePlot_width*0.9
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0(input$gene, '.pdf') },
    content = function(file) {
      ggsave(file, plot = plotInput(), width = 10, height = 10, device = "pdf")
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)




















