library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(shiny)
library(shinythemes)
library(bslib)
library(shinydashboard)



ui <- fluidPage(
  # App title ----
  titlePanel("Second Trimester Single Cell Data"),
  
  
  theme = bs_theme(bootswatch = 'yeti', base_font = "Arial"),
  #shinythemes::themeSelector(),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      textInput(inputId = "gene",
                label = "Human Gene Symbols (separated by spaces, case-sensitive)",
                value = "TBCK"),
      checkboxGroupInput("cells",
                         "Cell Types for Line Plot",
                         choices = c("cerebral cortex endothelial cell",
                                     "blood vessel endothelial cell",
                                     "native cell",
                                     "microglial cell",
                                     "oligodendrocyte precursor cell",
                                     "forebrain radial glia cell",
                                     "Cajal-Retzius cell",
                                     "progenitor cell",
                                     "glutamatergic neuron",
                                     "cerebral cortex GABAergic interneuron"),
                         selected = c("cerebral cortex endothelial cell",
                                      "blood vessel endothelial cell",
                                      "native cell",
                                      "microglial cell",
                                      "oligodendrocyte precursor cell",
                                      "forebrain radial glia cell",
                                      "Cajal-Retzius cell",
                                      "progenitor cell",
                                      "glutamatergic neuron",
                                      "cerebral cortex GABAergic interneuron")),
      radioButtons("scale", "Custom scale plots (Note: Setting plot limits may cut off data.)", c("No", "Yes"), "No"),
      conditionalPanel(
        condition = "input.scale == 'Yes'",
        numericInput(
          "maxPct", "Maximum percent expressing", min = 0, max = 100, 50
        ),
        numericInput(
          "minAvg", "Minimum average scaled expression", 0
        ),
        numericInput(
          "maxAvg", "Maximum average scaled expression", 0
        )),
      actionButton("go", "Plot")
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(style = "position:absolute;right:10px;top:10px;", 
          actionButton("show_about", "About", position = "top right")),
      # Output: Plot ----
      shinydashboard::box(width = NULL, style='overflow-x: scroll;height:400px;overflow-y: scroll;',
          plotOutput("dotPlot")),
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_dp")),
      br(),
      br(),
      shinydashboard::box(width = NULL, style='overflow-x: scroll;height:400px;overflow-y: scroll;',
          plotOutput("pctLinePlot")),
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_plp")),
      br(),
      br(),
      shinydashboard::box(width = NULL, style='overflow-x: scroll;height:400px;overflow-y: scroll;',
          plotOutput("avgLinePlot")),
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_alp")),
      br(),
      br(),
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_data", "Data")),
      br(),
      br()
      
      
      
    )
  )
)



# Define server logic required to draw a histogram ----
server <- function(input, output) {
  d <- read.csv(file = "data/all_genes_second_trimester_sc_data.csv")
  

  rv_h <- eventReactive(input$go,{
    validate(
      need(length(unique(strsplit(input$gene, split = " ")[[1]])) > 0, paste("Please enter at least one gene."))
    )
    ceiling(sqrt(length(unique(strsplit(input$gene, split = " ")[[1]])))) * 300
  })
  rv_w <- eventReactive(input$go,{
    validate(
      need(length(unique(strsplit(input$gene, split = " ")[[1]])) > 0, paste("Please enter at least one gene."))
    )
    (ceiling(sqrt(length(unique(strsplit(input$gene, split = " ")[[1]])))) * 400) + 300
  })
  
  
  rv_dotplot <- eventReactive(input$go,
                              {
                                genes <- strsplit(input$gene, split = " ")[[1]]
                                validate(
                                  need(length(genes) > 0, paste("Please enter at least one gene."))
                                )
                                for (g in genes){
                                  validate(
                                    need(g %in% d$sym, paste("Gene", g, "not found."))
                                  )
                                }
                                
                                ngenes <- length(genes)

                                d.subset <- d %>% filter(sym %in% genes)
                                
                                d.sub <- d.subset %>% select(pct.exp, avg.exp.scaled, features.plot, id, sym)
                                d.sub$id <- substr(d.sub$id, 1, 2)
          
                                d.sub$features.plot <- factor(d.sub$features.plot, levels = rev(c("cerebral cortex endothelial cell",
                                                                                                        "blood vessel endothelial cell",
                                                                                                        "native cell",
                                                                                                        "microglial cell",
                                                                                                        "oligodendrocyte precursor cell",
                                                                                                        "forebrain radial glia cell",
                                                                                                        "Cajal-Retzius cell",
                                                                                                        "progenitor cell",
                                                                                                        "glutamatergic neuron",
                                                                                                        "cerebral cortex GABAergic interneuron")))
                                
                                p <- ggplot(d.sub, aes(y = features.plot, x = substr(id, 1, 2))) +
                                  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
                                  labs(x = "Week", y = "Cell Type", size = "Percent Expressing", color = "Average Expression Scaled") +
                                  theme_minimal(base_size = 16) +
                                  theme(panel.spacing = unit(0.2, "in")) +
                                  theme(plot.title = element_text(hjust = 0.5)) +
                                  facet_wrap(.~sym, ncol = ceiling(sqrt(length(unique(d.subset$sym)))))
                                if (input$scale == "No"){
                                  p + scale_color_gradient(low = "gray", high = "blueviolet")
                                }
                                else{
                                  p + scale_color_gradient(low = "gray", high = "blueviolet", limits =c(input$minAvg, input$maxAvg)) + lims(size = c(0, input$maxPct))
                                }
                              })

    output$dotPlot <- renderPlot({rv_dotplot()},
                               height = function() {
                                 rv_h()
                                 },
                               width = function() {
                                 rv_w()
                               })

  
  rv_pct_lineplot <- eventReactive(input$go,
                                   {
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     validate(
                                       need(length(genes) > 0, paste("Please enter at least one gene."))
                                     )
                                     for (g in genes){
                                       validate(
                                         need(g %in% d$sym, paste("Gene", g, "not found."))
                                         )
                                     }
                                     
                                     validate(
                                       need(length(input$cells) > 0, paste("Please select at least one cell type."))
                                     )
                                     
                                     d.subset <- d %>% filter(sym %in% genes & features.plot %in% input$cells)
                                     
                                     d.sub <- d.subset %>% select(pct.exp, avg.exp.scaled, features.plot, id, sym)
                                     d.sub$id <- substr(d.sub$id, 1, 2)
                                     
                                     
                                     d.sub$features.plot <- factor(d.sub$features.plot, levels = rev(c("cerebral cortex endothelial cell",
                                                                                                             "blood vessel endothelial cell",
                                                                                                             "native cell",
                                                                                                             "microglial cell",
                                                                                                             "oligodendrocyte precursor cell",
                                                                                                             "forebrain radial glia cell",
                                                                                                             "Cajal-Retzius cell",
                                                                                                             "progenitor cell",
                                                                                                             "glutamatergic neuron",
                                                                                                             "cerebral cortex GABAergic interneuron")))
                                     
                                     
                                     p <- ggplot(d.sub, aes(y = pct.exp, x = substr(id, 1, 2), color = features.plot)) +
                                       geom_line(aes(group = features.plot)) +
                                       geom_point() +
                                       labs(x = "Week", color = "Cell Type", y = "Percent Expressing") +
                                       theme_minimal(base_size = 16) +
                                       theme(panel.spacing = unit(0.2, "in")) +
                                       theme(plot.title = element_text(hjust = 0.5)) +
                                       facet_wrap(.~sym, ncol = ceiling(sqrt(length(unique(d.sub$sym)))))

                                     if (input$scale == "No"){
                                       p
                                     }
                                     else{
                                       p + lims(y = c(0, input$maxPct))
                                     }
                                   })
  
  output$pctLinePlot <- renderPlot({rv_pct_lineplot()},
                                   height = function() {
                                     rv_h()
                                   },
                                   width = function() {
                                     rv_w()
                                   })
  
  
  
  
  rv_avg_lineplot <- eventReactive(input$go,
                                   {
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     validate(
                                       need(length(genes) > 0, paste("Please enter at least one gene."))
                                     )
                                     for (g in genes){
                                       validate(
                                         need(g %in% d$sym, paste("Gene", g, "not found."))
                                       )
                                     }
                                     
                                     validate(
                                       need(length(input$cells) > 0, paste("Please select at least one cell type."))
                                     )
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     d.subset <- d %>% filter(sym %in% genes & features.plot %in% input$cells)
                                     
                                     d.sub <- d.subset %>% select(pct.exp, avg.exp.scaled, features.plot, id, sym)
                                     d.sub$id <- substr(d.sub$id, 1, 2)
                                     
                                     d.sub$features.plot <- factor(d.sub$features.plot, levels = rev(c("cerebral cortex endothelial cell",
                                                                                                             "blood vessel endothelial cell",
                                                                                                             "native cell",
                                                                                                             "microglial cell",
                                                                                                             "oligodendrocyte precursor cell",
                                                                                                             "forebrain radial glia cell",
                                                                                                             "Cajal-Retzius cell",
                                                                                                             "progenitor cell",
                                                                                                             "glutamatergic neuron",
                                                                                                             "cerebral cortex GABAergic interneuron")))
                                     
                                     
                                     p <- ggplot(d.sub, aes(y = avg.exp.scaled, x = substr(id, 1, 2), color = features.plot)) +
                                       geom_line(aes(group = features.plot)) +
                                       geom_point() +
                                       labs(x = "Week", color = "Cell Type", y = "Average Expression Scaled") +
                                       theme_minimal(base_size = 16) +
                                       theme(panel.spacing = unit(0.2, "in")) +
                                       theme(plot.title = element_text(hjust = 0.5)) +
                                       facet_wrap(.~sym, ncol = ceiling(sqrt(length(unique(d.sub$sym)))))

                                     if (input$scale == "No"){
                                       p
                                     }
                                     else{
                                       p + lims(y = c(input$minAvg, input$maxAvg))
                                     }
                                   })
  
  output$avgLinePlot <- renderPlot({rv_avg_lineplot()},
                                   height = function() {
                                     rv_h()
                                   },
                                   width = function() {
                                     rv_w()
                                   })
  
  
  output$download_dp = downloadHandler(
    filename = function() { paste0(input$gene, '_dotplot.pdf') },
    content = function(file) {
      pdf(file, width = rv_w() / 80, height = rv_h() / 80)
      print(rv_dotplot())
      dev.off()
    })
  
  output$download_plp = downloadHandler(
    filename = function() { paste0(input$gene, '_percent_line_plot.pdf') },
    content = function(file) {
      pdf(file, width = rv_w() / 80, height = rv_h() / 80)
      print(rv_pct_lineplot())
      dev.off()
    })
  
  output$download_alp = downloadHandler(
    filename = function() { paste0(input$gene, '_average_line_plot.pdf') },
    content = function(file) {
      pdf(file, width = rv_w() / 80, height = rv_h() / 80)
      print(rv_avg_lineplot())
      dev.off()
    })
  
  output$download_data = output$download <- downloadHandler(
    filename = function(){paste0(input$gene, "_single_cell_data.csv")}, 
    content = function(fname){
      write.csv(rv_dotplot()$data, fname, quote = F)
    }
  )
  
  text_about <- HTML("Data from Bhaduri et al. 2021<br>
    Please cite the following when using these plots:<br>
    Bhaduri, A., Sandoval-Espinosa, C., Otero-Garcia, M. et al. An atlas of cortical arealization identifies dynamic molecular signatures. Nature 598, 200–204 (2021). https://doi.org/10.1038/s41586-021-03910-8
    <br>
    Please use gene symbols from Ensembl version 110")
  observeEvent(input$show_about, {showModal(modalDialog(text_about, title = 'About', size = "l"))})
  
  
}





shinyApp(ui = ui, server = server)











