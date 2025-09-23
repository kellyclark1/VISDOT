#load packages needed to generate plots
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(shiny)
library(shinythemes)
library(bslib)
library(shinydashboard)
library(viridis)
library(scales)


#define the user interface
ui <- fluidPage(
  
  # App title
  titlePanel("NeuroTri2-VISDOT"),
  
  #set the theme and font for the ShinyApp
  theme = bs_theme(bootswatch = 'yeti', base_font = "Arial"),

  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Gene names, as HGNC gene symbols
      textInput(inputId = "gene",
                label = "Human HGNC Gene Symbols (separated by spaces, case-sensitive)",
                value = "TBCK"),
      
      # Input: Cell types- all are selected by default and can be unselected
      checkboxGroupInput("cells",
                         "Cell Types",
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
      
      # Input: Custom scaling- choose whether to set the scale of
      #  the plots (Yes) or whether to allow them to auto-scale (No, default)
      radioButtons("scale",
                   "Custom scale plots (Note: Setting plot limits may cut off data.)",
                   c("No", "Yes"), "No"),
      
      #Input: Scaling options- appears only when "Yes" is selected for custom
      #  scaling, allows input of plot limits
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
      
      #Button to click to plot genes
      actionButton("go", "Plot")
    ),
    
    
    
    # Main panel for displaying outputs
    mainPanel(
      
      #About button to show data and citation details
      div(style = "position:absolute;right:10px;top:10px;", 
          actionButton("show_about", "About", position = "top right")),
      
      
      # Output: Dot plot
      shinydashboard::box(width = NULL,
                          #set size of space for plot and allow scrolling
                          #  if plot exceeds that space
                          style='overflow-x: scroll;height:600px;overflow-y: scroll;',
          plotOutput("dotPlot")),
      
      #Button to download dot plot as a pdf
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_dp")),
      br(),
      br(),
      
      #Output: Percent expression line plot
      shinydashboard::box(width = NULL,
                          #set size of space for plot and allow scrolling
                          #  if plot exceeds that space
                          style='overflow-x: scroll;height:600px;overflow-y: scroll;',
          plotOutput("pctLinePlot")),
      
      #Button to download percent expression line plot as a pdf
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_plp")),
      br(),
      br(),
      
      #Output: Average expression line plot
      shinydashboard::box(width = NULL,
                          #set size of space for plot and allow scrolling
                          #  if plot exceeds that space
                          style='overflow-x: scroll;height:600px;overflow-y: scroll;',
          plotOutput("avgLinePlot")),
      
      #Button to download average expression plot as a pdf
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_alp")),
      br(),
      br(),
      
      #Button to download data for selected genes as a csv file
      div(style = "position:absolute;right:10px;", 
          downloadButton("download_data", "Data")),
      br(),
      br()
      
    )
  )
)



# Define the server to generate the outputs
server <- function(input, output) {
  
  # Read in the data file, which is stored in the ShinyApp in the data folder
  d <- read.csv(file = "data/all_genes_second_trimester_sc_data.csv")
  
  # Calculate the appropriate height for plots when the Plot button is clicked
  rv_h <- eventReactive(input$go,{
    
    # If there is no input, print an error message
    validate(
      need(length(unique(strsplit(input$gene, split = " ")[[1]])) > 0,
           paste("Please enter at least one gene."))
    )
    
    # Set the height of a single row of plots (occurs with 1-2 genes)
    if (length(unique(strsplit(input$gene, split = " ")[[1]])) <= 2){
      500
    }else{ # Set the height for multiple rows of plots
      # where unique(strsplit(input$gene, , split = " ")[[1]]) is a vector with
      # the input list of genes, so its length is the number of plots, and 
      # the ceiling of its square root if the number of rows of plots
      # height is 400 pixels per row of plots, plus 100 pixels for labels
    ceiling(sqrt(length(unique(strsplit(input$gene, split = " ")[[1]])))) * 400 + 100
    }
  })
  
  # Calculate the appropriate width for plots when the Plot button is clicked
  rv_w <- eventReactive(input$go,{
    
    # If there is no input, print an error message
    validate(
      need(length(unique(strsplit(input$gene, split = " ")[[1]])) > 0,
           paste("Please enter at least one gene."))
    ) 
    # Set the width by columns of plots
    # where unique(strsplit(input$gene, , split = " ")[[1]]) is a vector with
    # the input list of genes, so its length is the number of plots, and 
    # the ceiling of its square root if the number of columns of plots
    # width is 400 pixels per row of plots, plus 300 pixels for key and labels
    (ceiling(sqrt(length(unique(strsplit(input$gene, split = " ")[[1]])))) * 400) + 300
  })
  
  
  # Define dot plot output
  rv_dotplot <- eventReactive(input$go, # Generate when "Plot" is clicked
                              {
                                # Get the input list of genes as a vector
                                genes <- strsplit(input$gene, split = " ")[[1]]
                                
                                # If there are no genes, print an error message
                                validate(
                                  need(length(genes) > 0,
                                       paste("Please enter at least one gene."))
                                )
                                
                                #If there is a gene input that is not present in
                                #  the data, print an error message
                                for (g in genes){
                                  validate(
                                    need(g %in% d$sym,
                                         paste("Gene", g, "not found."))
                                  )
                                }
                                
                                # Get the number of genes input
                                ngenes <- length(genes)
                                
                                # Print an error message if there are no
                                #  cell types selected
                                validate(
                                  need(length(input$cells) > 0,
                                       paste("Please select at least one cell type."))
                                )
                                
                                # Subset the data to just that of input
                                #  genes and cell types
                                d.subset <- d %>%
                                  dplyr::filter(sym %in% genes &
                                                  features.plot %in% input$cells)
                                
                                # Get just the columns with 
                                d.sub <- d.subset %>%
                                  dplyr::select(pct.exp, # percent expressing,
                                                avg.exp.scaled, # average expression,
                                                features.plot, # cell types,
                                                id, # time points,
                                                sym) # and gene names
                                
                                # From the data for those genes, get:
                                d.sub <- d.subset %>%
                                  dplyr::select(pct.exp, # the percent of cells expressing,
                                                avg.exp.scaled, # the average expression per cell,
                                                features.plot, # the cell types,
                                                id, # the time point
                                                sym) # and the gene names
                                
                                # Redefine the time point as just the number of weeks
                                d.sub$id <- substr(d.sub$id, 1, 2) 
          
                                # Set the order of the cell types
                                d.sub$features.plot <- factor(
                                  d.sub$features.plot,
                                  levels = rev(c("cerebral cortex endothelial cell",
                                                 "blood vessel endothelial cell",
                                                 "native cell",
                                                 "microglial cell",
                                                 "oligodendrocyte precursor cell",
                                                 "forebrain radial glia cell",
                                                 "Cajal-Retzius cell",
                                                 "progenitor cell",
                                                 "glutamatergic neuron",
                                                 "cerebral cortex GABAergic interneuron")))
                                
                                # Make the dot plot with the subset data
                                p <- ggplot(d.sub,
                                            # with y as cell types
                                            aes(y = features.plot,
                                                # and x as time points
                                                x = substr(id, 1, 2))) +
                                  
                                  # Make it a scatter plot with the points
                                  #  sized by percent expressing and colored by
                                  #  average expression
                                  geom_point(aes(size = pct.exp, 
                                                 color = avg.exp.scaled)) +
                                  
                                  # Label axes, size, and color
                                  labs(x = "Week",
                                       y = "Cell Type",
                                       size = "Percent Expressing",
                                       color = "Average Expression\nScaled") +
                                  
                                  # Set plot theme and font sizes
                                  theme_minimal(base_size = 16) +
                                  theme(strip.text = element_text(size = 20),
                                        axis.text = element_text(size = 18),
                                        axis.title = element_text(size = 20),
                                        legend.text = element_text(size = 18),
                                        legend.title = element_text(size = 20)) +
                                  
                                  # Set spacing between plots/genes
                                  theme(panel.spacing = unit(0.2, "in")) +
                                  
                                  # Center plot title
                                  theme(plot.title = element_text(hjust = 0.5)) +
                                  
                                  # Plot each gene as a separate small plot
                                  facet_wrap(.~sym,
                                             # the number of columns is the square root
                                             #  of the number of genes, rounded up to 
                                             #  the nearest whole number
                                             ncol = ceiling(sqrt(length(unique(d.subset$sym))))) +
                                  
                                  #wrap long cell types onto multiple lines
                                  scale_y_discrete(labels = label_wrap(width = 22))
                                
                                
                                # Use the colorblind-friendly viridis palette
                                # Set the color scale for auto plot scaling if 
                                #  custom scaling is not selected
                                if (input$scale == "No"){
                                  p + scale_color_viridis(end = 0.9,
                                                          direction = -1)
                                }
                                # Or limit it by the custom scale limits if
                                #  custom scaling is selected
                                else{
                                  p + scale_color_viridis(end = 0.9,
                                                          direction = -1,
                                                          limits =c(input$minAvg, input$maxAvg)) +
                                    lims(size = c(0, input$maxPct))
                                }
                              })

    # Render the dot plot
    output$dotPlot <- renderPlot({rv_dotplot()},
                               # Calculate the height of the plot with the
                               #  previously defined function
                               height = function() {
                                 rv_h()
                                 },
                               # Calculate the width of the plot with the
                               #  previously defined function
                               width = function() {
                                 rv_w()
                               })

  
   # Generate the percent expressing line plot when "Plot" is clicked
   rv_pct_lineplot <- eventReactive(input$go,
                                   {
                                     # Get the list of input genes as a vector
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     
                                     # Print an error message if there are no genes input
                                     validate(
                                       need(length(genes) > 0,
                                            paste("Please enter at least one gene."))
                                     )
                                     
                                     #Print an error message if there is a gene
                                     #  that is not in the data
                                     for (g in genes){
                                       validate(
                                         need(g %in% d$sym,
                                              paste("Gene", g, "not found."))
                                         )
                                     }
                                     
                                     # Print an error message if there are no
                                     #  cell types selected
                                     validate(
                                       need(length(input$cells) > 0,
                                            paste("Please select at least one cell type."))
                                     )
                                     
                                     # Subset the data to just that of input
                                     #  genes and cell types
                                     d.subset <- d %>%
                                       dplyr::filter(sym %in% genes &
                                                       features.plot %in% input$cells)
                                     
                                     # Get just the columns with 
                                     d.sub <- d.subset %>%
                                       dplyr::select(pct.exp, # percent expressing,
                                                     avg.exp.scaled, # average expression,
                                                     features.plot, # cell types,
                                                     id, # time points,
                                                     sym) # and gene names
                                     
                                     # Redefine the time point as just number of weeks
                                     d.sub$id <- substr(d.sub$id, 1, 2)
                                     
                                     # Set the order of the cell types
                                     d.sub$features.plot <- factor(
                                       d.sub$features.plot,
                                       levels = rev(c("cerebral cortex endothelial cell",
                                                      "blood vessel endothelial cell",
                                                      "native cell",
                                                      "microglial cell",
                                                      "oligodendrocyte precursor cell",
                                                      "forebrain radial glia cell",
                                                      "Cajal-Retzius cell",
                                                      "progenitor cell",
                                                      "glutamatergic neuron",
                                                      "cerebral cortex GABAergic interneuron")))
                                     
                                     # Set the colors using the viridis palette
                                     category_colors <- viridis(n = 10, end = 0.9)
    
                                     # Set names for the colors as cell types to
                                     #  ensure consistent color coding
                                     names(category_colors) <- c("cerebral cortex endothelial cell",
                                                                 "blood vessel endothelial cell",
                                                                 "native cell",
                                                                 "microglial cell",
                                                                 "oligodendrocyte precursor cell",
                                                                 "forebrain radial glia cell",
                                                                 "Cajal-Retzius cell",
                                                                 "progenitor cell",
                                                                 "glutamatergic neuron",
                                                                 "cerebral cortex GABAergic interneuron")
                                     
                                     # Set different line type for each cell type
                                     #  with cell types as names to ensure consistency
                                     category_lines <- c("cerebral cortex endothelial cell" = "solid",
                                                         "blood vessel endothelial cell" = "11",
                                                         "native cell" = "31",
                                                         "microglial cell" = "61",
                                                         "oligodendrocyte precursor cell" = "1131",
                                                         "forebrain radial glia cell" = "1161",
                                                         "Cajal-Retzius cell" = "3161",
                                                         "progenitor cell" = "111131",
                                                         "glutamatergic neuron" = "111161",
                                                         "cerebral cortex GABAergic interneuron" = "113161")
                                     
                                     # Set different shape for each cell type
                                     #  with cell types as names to ensure consistency
                                     category_shapes <- c("cerebral cortex endothelial cell" = 7,
                                                         "blood vessel endothelial cell" = 16,
                                                         "native cell" = 17,
                                                         "microglial cell" = 18,
                                                         "oligodendrocyte precursor cell" = 22,
                                                         "forebrain radial glia cell" = 25,
                                                         "Cajal-Retzius cell" = 3,
                                                         "progenitor cell" = 4,
                                                         "glutamatergic neuron" = 8,
                                                         "cerebral cortex GABAergic interneuron" = 13)
                                       
                                     # Plot percent expressing line plot with
                                     p <- ggplot(d.sub, aes(y = pct.exp, # y as percent expressing,
                                                            x = substr(id, 1, 2), # x as time point,
                                                            color = features.plot, # color as cell type,
                                                            group = features.plot, # grouped by cell type,
                                                            linetype = features.plot, # line as cell type,
                                                            fill = features.plot, # fill as cell type,
                                                            shape = features.plot)) + # and shape as cell type
                                       geom_line(linewidth = 0.8, # Plot as line plot and set line width
                                                 aes(group = features.plot, # group by cell type
                                                     color = features.plot, # color by cell type
                                                     linetype = features.plot)) + #line type by cell type
                                       geom_point(aes(color = features.plot, #Add points
                                                      fill = features.plot, # fill by cell type
                                                      shape = features.plot), # shape as cell type
                                                  size = 3) + # Set point size
                                       
                                       # Set labels fpr axes and key
                                       labs(x = "Week",
                                            color = "Cell Type",
                                            fill = "Cell Type",
                                            linetype = "Cell Type",
                                            shape = "Cell Type",
                                            y = "Percent Expressing") +
                                       
                                       # Set theme and font sizes
                                       theme_minimal(base_size = 18) +
                                       theme(strip.text = element_text(size = 20),
                                             axis.text = element_text(size = 18),
                                             axis.title = element_text(size = 20),
                                             legend.text = element_text(size = 18),
                                             legend.title = element_text(size = 20)) +
                                       
                                       # Set space between genes/plots
                                       theme(panel.spacing = unit(0.2, "in")) +
                                       
                                       # Center plot title
                                       theme(plot.title = element_text(hjust = 0.5)) +
                                       
                                       # Set key size and spacing
                                       theme(legend.key.width = unit(1.5, "cm"),
                                             legend.key.spacing.y = unit(0.2, "cm")) +
                                       
                                       # Plot each gene separately, with the number of
                                       #  columns being the square root of the number
                                       #  of genes, rounded up to the nearest whole number
                                       facet_wrap(.~sym,
                                                  ncol = ceiling(sqrt(length(unique(d.sub$sym))))) +
                                       
                                       # Set colors as previously defined palette
                                       scale_color_manual(
                                         values = category_colors,
                                         breaks = names(category_colors),
                                         
                                         # Wrap long cell type labels
                                         labels = label_wrap(width = 22)
                                       ) +
                                       
                                       # Set fill as previously defined palette
                                       scale_fill_manual(
                                         values = category_colors,
                                         breaks = names(category_colors),
                                         
                                         # Wrap long cell type labels
                                         labels = label_wrap(width = 22) 
                                       ) +
                                       
                                       # Set line type by cell type
                                       scale_linetype_manual(values = category_lines,
                                                             breaks = names(category_lines),
                                                             # Wrap long cell type labels
                                                             labels = label_wrap(width = 22)) +
                                       
                                       # Set shape by cell type
                                       scale_shape_manual(values = category_shapes,
                                                          breaks = names(category_shapes),
                                                          # Wrap long cell type labels
                                                          labels = label_wrap(width = 22))
                                       
                                     # Plot with auto-scaling is custom scale is not checked
                                     if (input$scale == "No"){
                                       p
                                     }
                                     # Or plot using set limit if custom scale is checked
                                     else{
                                       p + lims(y = c(0, input$maxPct))
                                     }
                                   })
   
  # Render percent expressing line plot on ShinyApp
  output$pctLinePlot <- renderPlot({rv_pct_lineplot()},
                                   
                                   # Calculate and set height of plot based on
                                   #  previously defined function
                                   height = function() {
                                     rv_h()
                                   },
                                   
                                   # Calculate and set width of plot based on
                                   #  previously defined function
                                   width = function() {
                                     rv_w()
                                   })
  
  
  
  # Generate average expression line plot when Plot is clicked
  rv_avg_lineplot <- eventReactive(input$go,
                                   {
                                     # Get the list of input genes
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     
                                     # Print an error message if there are no genes entered
                                     validate(
                                       need(length(genes) > 0,
                                            paste("Please enter at least one gene."))
                                     )
                                     
                                     # Print an error message if an entered gene is not in the data
                                     for (g in genes){
                                       validate(
                                         need(g %in% d$sym,
                                              paste("Gene", g, "not found."))
                                       )
                                     }
                                     
                                     # Print an error message if there are no cell types selected
                                     validate(
                                       need(length(input$cells) > 0,
                                            paste("Please select at least one cell type."))
                                     )
                                     
                                     
                                     genes <- strsplit(input$gene, split = " ")[[1]]
                                     
                                     # Get just the data for the input genes and cell types
                                     d.subset <- d %>%
                                       dplyr::filter(sym %in% genes &
                                                       features.plot %in% input$cells)
                                     
                                     # Select from the data:
                                     d.sub <- d.subset %>%
                                       dplyr::select(pct.exp, # percent expressing,
                                                     avg.exp.scaled, # average expression,
                                                     features.plot, # cell types,
                                                     id, # time points,
                                                     sym) # and genes
                                     
                                     # Reassign time points to number of weeks
                                     d.sub$id <- substr(d.sub$id, 1, 2)
                                     
                                     
                                     # Assign order to the cell types
                                     d.sub$features.plot <-
                                       factor(d.sub$features.plot,
                                              levels = rev(c("cerebral cortex endothelial cell",
                                                             "blood vessel endothelial cell",
                                                             "native cell",
                                                             "microglial cell",
                                                             "oligodendrocyte precursor cell",
                                                             "forebrain radial glia cell",
                                                             "Cajal-Retzius cell",
                                                             "progenitor cell",
                                                             "glutamatergic neuron",
                                                             "cerebral cortex GABAergic interneuron")))
                                     
                                     # Get the viridid color palette
                                     category_colors <- viridis(n = 10, end = 0.9)
                                     
                                     # Assign cell types as names for each color for consistency
                                     names(category_colors) <- c("cerebral cortex endothelial cell",
                                                                 "blood vessel endothelial cell",
                                                                 "native cell",
                                                                 "microglial cell",
                                                                 "oligodendrocyte precursor cell",
                                                                 "forebrain radial glia cell",
                                                                 "Cajal-Retzius cell",
                                                                 "progenitor cell",
                                                                 "glutamatergic neuron",
                                                                 "cerebral cortex GABAergic interneuron")
                                     
                                     # Assign a line type to each cell type
                                     category_lines <- c("cerebral cortex endothelial cell" = "solid",
                                                         "blood vessel endothelial cell" = "11",
                                                         "native cell" = "31",
                                                         "microglial cell" = "61",
                                                         "oligodendrocyte precursor cell" = "1131",
                                                         "forebrain radial glia cell" = "1161",
                                                         "Cajal-Retzius cell" = "3161",
                                                         "progenitor cell" = "111131",
                                                         "glutamatergic neuron" = "111161",
                                                         "cerebral cortex GABAergic interneuron" = "113161")
                                     
                                     # Assign a shape to each cell type
                                     category_shapes <- c("cerebral cortex endothelial cell" = 7,
                                                          "blood vessel endothelial cell" = 16,
                                                          "native cell" = 17,
                                                          "microglial cell" = 18,
                                                          "oligodendrocyte precursor cell" = 22,
                                                          "forebrain radial glia cell" = 25,
                                                          "Cajal-Retzius cell" = 3,
                                                          "progenitor cell" = 4,
                                                          "glutamatergic neuron" = 8,
                                                          "cerebral cortex GABAergic interneuron" = 13)
                                     
                                     # Make the plot
                                     p <- ggplot(d.sub,
                                                 # with y as the average expression,
                                                 aes(y = avg.exp.scaled,
                                                     # x as the time point,
                                                     x = substr(id, 1, 2),
                                                     # color by cell type,
                                                     color = features.plot,
                                                     # group by cell type,
                                                     group = features.plot,
                                                     # line style by cell type,
                                                     linetype = features.plot,
                                                     # fill by cell type,
                                                     fill = features.plot,
                                                     # and shape as cell type
                                                     shape = features.plot)) +
                                       
                                       # Make the line plot and set the line width
                                       geom_line(linewidth = 0.8) +
                                       
                                       # Add the points
                                       geom_point(size = 3) +
                                       
                                       # Label axes and key
                                       labs(x = "Week",
                                            color = "Cell Type",
                                            fill = "Cell Type",
                                            linetype = "Cell Type",
                                            shape = "Cell Type",
                                            y = "Average Expression Scaled") +
                                       
                                       # Set font sizes
                                       theme_minimal(base_size = 16) +
                                       theme(strip.text = element_text(size = 20),
                                             axis.text = element_text(size = 18),
                                             axis.title = element_text(size = 20),
                                             legend.text = element_text(size = 18),
                                             legend.title = element_text(size = 20)) +
                                       
                                       # Set spacing between plots/genes
                                       theme(panel.spacing = unit(0.2, "in")) +
                                       
                                       # Center plot title
                                       theme(plot.title = element_text(hjust = 0.5)) +
                                       
                                       # Plot each gene separately, with the
                                       #  number of columns as the square root
                                       #  of the number of genes, rounded up to
                                       #  the nearest whole number
                                       facet_wrap(.~sym,
                                                  ncol = ceiling(sqrt(length(unique(d.sub$sym))))) +
                                       
                                       # Set key sizing and spacing
                                       theme(legend.key.width = unit(1.5, "cm"),
                                             legend.key.spacing.y = unit(0.2, "cm")) +
                                       
                                       #Assign colors and wrap long cell labels
                                       scale_color_manual(
                                         values = category_colors,
                                         breaks = names(category_colors),
                                         labels = label_wrap( width = 22)
                                       ) +
                                       
                                       #Assign fills and wrap long cell labels
                                       scale_fill_manual(
                                         values = category_colors,
                                         breaks = names(category_colors),
                                         labels = label_wrap( width = 22) # Ensure all levels are included in the legend
                                       ) +
                                       
                                       #Assign line types and wrap long cell labels
                                       scale_linetype_manual(values = category_lines,
                                                             breaks = names(category_lines),
                                                             labels = label_wrap( width = 22)) +
                                       
                                       #Assign shapes and wrap long cell labels
                                       scale_shape_manual(values = category_shapes,
                                                          breaks = names(category_shapes),
                                                          labels = label_wrap( width = 22))
                                     
                                     
                                     # Autoscale or scale by custom scale, depending on input
                                     if (input$scale == "No"){
                                       p
                                     }
                                     else{
                                       p + lims(y = c(input$minAvg, input$maxAvg))
                                     }
                                   })
  
  # Render average expression line plot
  output$avgLinePlot <- renderPlot({rv_avg_lineplot()},
                                   
                                   # Set plot height with previously defined function
                                   height = function() {
                                     rv_h()
                                   },
                                   
                                   # Set plot width with previously defined function
                                   width = function() {
                                     rv_w()
                                   })
  
  
  # Download dot plot as pdf with button
  output$download_dp = downloadHandler(
    filename = function() { paste0(sub(" ", "_", input$gene), '_dotplot.pdf') },
    content = function(file) {
      # Set pdf size
      pdf(file, width = rv_w() / 60, height = rv_h() / 60)
      print(rv_dotplot())
      dev.off()
    })
  
  # Download percent line plot as pdf with button
  output$download_plp = downloadHandler(
    filename = function() { paste0(gsub(" ", "_", input$gene), '_percent_line_plot.pdf') },
    content = function(file) {
      # Set pdf size
      pdf(file, width = rv_w() / 60, height = rv_h() / 60)
      print(rv_pct_lineplot())
      dev.off()
    })
  

  
  # Download average line plot as pdf with button
  output$download_alp = downloadHandler(
    filename = function() { paste0(gsub(" ", "_", input$gene), '_average_line_plot.pdf') },
    content = function(file) {
      # Set pdf size
      pdf(file, width = rv_w() / 60, height = rv_h() / 60)
      print(rv_avg_lineplot())
      dev.off()
    })
  

  # Download csv of data for genes of interest with button
  output$download_data = output$download <- downloadHandler(
    filename = function(){paste0(sub(" ", "_", input$gene), "_single_cell_data.csv")}, 
    content = function(fname){
      d = rv_dotplot()$data
      d <- d %>% dplyr::select(gene = sym, cell.type = features.plot, week = id, pct.exp, avg.exp.scaled)
      write.csv(d, fname, quote = F, row.names = F)
    }
  )
  
  
  # Set text for "About" popup
  text_about <- HTML("Data from Bhaduri et al. 2021<br> <br>
    NeuroTri2-VISDOT is a web-based data visualization tool that allows users to enter gene names to generate invidual or side-by-side plots of gene expression in cell types of their choosing in the second trimester developing human brain using the data published by Bhaduri et al.
    <br> <br>
    Please cite the following when using NeuroTri2-VISDOT:
    <br> <br>
    Bhaduri, A., Sandoval-Espinosa, C., Otero-Garcia, M. et al. An atlas of cortical arealization identifies dynamic molecular signatures. Nature 598, 200–204 (2021). https://doi.org/10.1038/s41586-021-03910-8
    <br> <br>
    Clark, K.J., Lubin, E.E., Gonzalez, E.M., Sangree, A.K., Layo-Carris, D.E., Durham, E.L., Ahrens-Nicklas, R.C., Nomakuchi, T.T., Bhoj, E.J. NeuroTri2-VISDOT: An open-access tool to harness the power of second trimester human single cell data to inform models of Mendelian neurodevelopmental disorders. BioRxiv (2024). https://doi.org/10.1101/2024.02.01.578438.
    <br> <br>
    Please use HGNC gene symbols consistent with Ensembl version 110")
  observeEvent(input$show_about, {showModal(modalDialog(text_about, title = 'About', size = "l"))})
  
  
}




# Run the app
shinyApp(ui = ui, server = server)











