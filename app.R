# loading all the required libraries
library(shiny)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DT)
library(corrplot)
library(shinyWidgets)
library(shinydashboard)
library(shinyBS)
library(tippy)
library(htmltools)
library(shinycustomloader)
library(shinyjs)
library(scatterD3)
library(viridis)

# setting the working directory and loading all the required datasets
#setwd("/srv/shiny-server/RShiny/")
#setwd("/Users/ajitknarendra/Desktop/PEM/RShiny/")
load("input_files_2.RDa")
load("average_expression_matrices_2.RDa")
load("percent_exprs_matrices_2.RDa")
load("avg_cellty_all_age_matrices_2.RDa")

# arranging the cell types in a specific order so that they appear in this order in the drop-down 
# menus for cell type, as well as the dot plot
cellt_order <- c("alveolar.type.1.cells",
                 "alveolar.type.2.cells",
                 "AT1/AT2.like.cells",
                 "AT2/club.like.cells",
                 "basal.cells",
                 "club.cells",
                 "ciliated.cells",
                 "goblet.cells",
                 "PNECs",
                 "matrix.fibroblast.1",
                 "matrix.fibroblast.2",
                 "myofibroblasts",
                 "pericytes",
                 "airway.smooth.muscle",
                 "vascular.smooth.muscle",
                 "chondrocytes",
                 "arteries",
                 "veins",
                 "Cap1",
                 "Cap2",
                 "lymphatics",
                 "bronchial.vessel",
                 "Alveolar.macrophages",
                 "Interstitial.macrophages",
                 "monocytes",
                 "mast.cells",
                 "dendritic.cells",
                 "T.cells",
                 "B.cells",
                 "NK.cells",
                 "enucleated.erythrocytes")
cellt_order_df <- data.frame("CellTypes" = cellt_order)
all_cell_types <- left_join(cellt_order_df, all_cell_types)
epithelial_cells <- inner_join(cellt_order_df, epithelial_cells)
mesenchymal_cells <- inner_join(cellt_order_df, mesenchymal_cells)
endothelial_cells <- inner_join(cellt_order_df, endothelial_cells)
hematopoeitic_cells <- inner_join(cellt_order_df, hematopoeitic_cells)

# creating a dashboardSidebar, adding menu items and subitems to it
sidebar <- dashboardSidebar(
  width = 250,
  tags$head(tags$style(HTML('.shiny-server-account { display: none; }'))),
  
  sidebarMenu(id = "tabs",       
              menuItem("RNA-seq", icon = icon("dna"), startExpanded = TRUE,
                       menuSubItem("Dot Plot", tabName = "dot_plot", icon = icon("braille"),
                                   selected = TRUE),
                       menuSubItem("Heatmap", tabName = "heatmap", icon = icon("fire")),
                       menuSubItem("Gene Expression Data", tabName = "exprs_data", icon = icon("table"))
              ),
              menuItem("About", icon = icon("info-circle"), tabName = "main_page")
  )
)

# creating the dot plot page along with all the contents
# there are 3 tab items in this page - dot plot by cell type, by age and by sample
# dot plot by cell type
dotplot_pg_celltype <- tabPanel("By Cell Type", style = "background-color: #fff; border-color: #D3D3D3;",
                                wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                          fluidRow(
                                            column(3,
                                                   selectizeInput('color_cellt_dotp',
                                                                  label = "Select Color : ",
                                                                  choices = list(
                                                                    "Red",
                                                                    "Blue",
                                                                    "Green",
                                                                    "Orange",
                                                                    "Purple",
                                                                    "Grey",
                                                                    "Purple - Orange",
                                                                    "Orange - Red",
                                                                    "Purple - Red",
                                                                    "Red - Yellow - Blue",
                                                                    "Yellow - Orange - Red",
                                                                    "Viridis",
                                                                    "Spectral",
                                                                    "Inferno",
                                                                    "Plasma",
                                                                    "Warm",
                                                                    "Rainbow"),
                                                                  multiple = FALSE,
                                                                  selected = "Viridis"
                                                   )
                                            ),
                                            column(2),
                                            column(6,
                                                   setSliderColor(color = c("#D1B2D1", "#D1B2D1", "#D1B2D1", "#D1B2D1",
                                                                            "#D1B2D1", "#D1B2D1"),
                                                                  sliderId = c(1,2,3,4,5,6)),
                                                   sliderInput('leftmarg_dotp_cellt',
                                                               label = "Adjust the Left Margin : ",
                                                               min = 1, max = 250,  
                                                               value = 125)
                                            ),
                                            column(12,
                                                   box(title = "Dot Plot for Cell Type", solidHeader = TRUE, 
                                                       collapsible = TRUE, 
                                                       status = "primary", 
                                                       collapsed = FALSE, br(), width = "100%",
                                                       scatterD3Output('dotplot_plot_cellt', width = "100%"),
                                                       tags$style(HTML('#dotplot_plot_cellt .scatterD3 .x.axis text {
                                                                       fill: #000;
                                                                       font-style: italic;
                                                                       }'))
                                                   )
                                            )
                                          )
                                ),
                           
                           wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                     box(title = "Corresponding Data for Cell Type", 
                                         solidHeader = TRUE, collapsible = TRUE, 
                                         status = "primary", 
                                         collapsed = TRUE, br(), width = "100%",
                                         DT::dataTableOutput('dotplot_data_cellt', width = "100%")
                                     ),
                                     br(),
                                     conditionalPanel("output.dotplot_data_cellt", 
                                                      downloadButton("downloadDotPDCt", "Download Corresponding Data", 
                                                                     class = "d_dotbutt")),
                                     tags$head(tags$style(".d_dotbutt{background-color: #D1B2D1;} 
                                                          .d_dotbutt{border-color: #D3D3D3;}")),
                                     br()
                                     )
)

# dot plot by age 
dotplot_pg_age <- tabPanel("By Age", style = "background-color: #fff; border-color: #D3D3D3;",
                           wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                     fluidRow(
                                       column(3,
                                              selectizeInput('age_dotp_list',
                                                             label = "Select Age : ",
                                                             choices = list(
                                                               "All Ages",
                                                               "Age" = age_list
                                                             ),
                                                             multiple = TRUE,
                                                             selected = "All Ages") 
                                       ),
                                       column(3,
                                              selectizeInput('symbol_age_dotp',
                                                             label = "Select Symbol : ",
                                                             choices = list(
                                                               "Gene",
                                                               "Age",
                                                               "None"
                                                             ),
                                                             multiple = FALSE,
                                                             selected = "None")
                                       ),
                                       column(3,
                                              selectizeInput('color_age_dotp',
                                                             label = "Select Color : ",
                                                             choices = list(
                                                               "Red",
                                                               "Blue",
                                                               "Green",
                                                               "Orange",
                                                               "Purple",
                                                               "Grey",
                                                               "Purple - Orange",
                                                               "Orange - Red",
                                                               "Purple - Red",
                                                               "Red - Yellow - Blue",
                                                               "Yellow - Orange - Red",
                                                               "Viridis",
                                                               "Spectral",
                                                               "Inferno",
                                                               "Plasma",
                                                               "Warm",
                                                               "Rainbow"),
                                                             multiple = FALSE,
                                                             selected = "Viridis"
                                              )
                                       ),
                                       column(3,
                                              sliderInput('leftmarg_dotp_age',
                                                          label = "Adjust the Left Margin : ",
                                                          min = 1, max = 250,  
                                                          value = 125)
                                       ),
                                       column(12,
                                              box(title = "Dot Plot for Age", solidHeader = TRUE, 
                                                  collapsible = TRUE, 
                                                  status = "primary", 
                                                  collapsed = FALSE, br(), width = "100%", 
                                                  scatterD3Output('dotplot_plot_age', width = "100%"),
                                                  tags$style(HTML('#dotplot_plot_age .scatterD3 .x.axis text {
                                                                  fill: #000;
                                                                  font-style: italic;
                                                                  }'))
                                              )
                                       )
                                     )
                           ),
                           
                           wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                     box(title = "Corresponding Data for Age", 
                                         solidHeader = TRUE, collapsible = TRUE, 
                                         status = "primary", 
                                         collapsed = TRUE, br(), width = "100%",
                                         DT::dataTableOutput('dotplot_data_age', width = "100%")
                                     ),
                                     br(),
                                     conditionalPanel("output.dotplot_data_age", 
                                                      downloadButton("downloadDotPDAge", "Download Corresponding Data", 
                                                                     class = "d_dotbutt")),
                                     br()
                           )
)

# dot plot by sample 
dotplot_pg_sample <- tabPanel("By Sample",
                              style = "background-color: #fff; border-color: #D3D3D3;",
                              wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                        fluidRow(
                                          column(3,
                                                 uiOutput("gene_samp"),
                                                 uiOutput("cellty_samp")
                                          ),
                                          column(3,
                                                 selectizeInput('sample_exprs_type',
                                                                label = "Select Expression : ",
                                                                choices = list(
                                                                  "Average",
                                                                  "Percent"
                                                                ),
                                                                multiple = FALSE,
                                                                selected = "Average"
                                                 ),
                                                 selectizeInput('color_samp_dotp',
                                                                label = "Select Color : ",
                                                                choices = list(
                                                                  "Red", 
                                                                  "Blue", 
                                                                  "Green", 
                                                                  "Orange", 
                                                                  "Yellow", 
                                                                  "Maroon", 
                                                                  "Purple", 
                                                                  "Violet", 
                                                                  "Grey", 
                                                                  "Black", 
                                                                  "Color by Age"),
                                                                multiple = FALSE,
                                                                selected = "Blue"
                                                 )
                                          ),
                                          column(3,
                                                 HTML(paste("<b>","Static Scatter Plot : ", "</b>")),
                                                 switchInput(
                                                   inputId = "sc_plot_samp",
                                                   label = "<i class=\"fas fa-chart-bar\"></i>",
                                                   onLabel = "Yes",
                                                   offLabel = "No",
                                                   width = "100%"
                                                 ),
                                                 sliderInput('point_size_sample',
                                                             label = "Adjust the Point Size: ",
                                                             min = 1, max = 500,  
                                                             value = 100)
                                          ),
                                          column(3,
                                                 sliderInput('leftmarg_dotp_sample',
                                                             label = "Adjust the Left Margin: ",
                                                             min = 1, max = 100,  
                                                             value = 30)
                                          )
                                        ),
                                        fluidRow(
                                          column(12,
                                                 box(title = "Graph for Sample grouped by Age", 
                                                     solidHeader = TRUE, collapsible = TRUE, 
                                                     status = "primary", 
                                                     collapsed = FALSE, br(), width = "100%",
                                                     uiOutput('dotplot_plot_sample'),
                                                     tags$head(tags$style(".d_dotbutt{background-color: #D1B2D1;} 
                                                                  .d_dotbutt{border-color: #D3D3D3;}"))
                                                 )
                                          )
                                        )
                              ),
                              
                              wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                        box(title = "Corresponding Data for Sample grouped by Age", 
                                            solidHeader = TRUE, collapsible = TRUE, 
                                            status = "primary", 
                                            collapsed = TRUE, br(), width = "100%",
                                            DT::dataTableOutput('dotplot_data_sample', width = "100%")
                                        ),
                                        br(),
                                        conditionalPanel("output.dotplot_data_sample", 
                                                         downloadButton("downloadDotPDSamp", "Download Corresponding Data", 
                                                                        class = "d_dotbutt")),
                                        br()
                              )
)

# combining all the three tabs in the dot plot page                                      
dotplot_page <- tabItem(
  tabName = "dot_plot",
  fluidRow(
    column(3, 
           wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                     wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                               selectizeInput('dp_gene', 
                                              label = HTML("Select Gene <br/> (One or More) : "),
                                              choices = NULL, 
                                              options = list(create = TRUE),
                                              multiple = TRUE, 
                                              selected = NULL
                               ),
                               selectizeInput('dp_cellty',
                                              label = HTML("Select Cell Type <br/> (One or More) : "),
                                              choices = NULL,
                                              options = list(create = TRUE), 
                                              multiple = TRUE,
                                              selected = NULL
                               ),
                               actionButton("reset_def_dp", "Reset to Default", class = "d_dotbutt",
                                            width = "100%", style='padding:4px; font-size:80%')
                               ),
                     wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                               h4(tags$b("Gene Links")),
                               selectizeInput('gene_link_n', 
                                              label = "Select Gene : ",
                                              choices = NULL,
                                              multiple = FALSE, 
                                              selected = NULL
                               ),
                               h5(tags$b("Click here to see the Input Gene in : ")),
                               actionButton("umap_link", "Single Cell UMAP", 
                                            class = "d_dotbutt",
                                            width = "100%", style='padding:4px; font-size:80%'),
                               br(),
                               br(),
                               actionButton("lungep_link", "Region Overlap", 
                                            class = "d_dotbutt", width = "100%",
                                            style='padding:4px; font-size:80%'),
                               br(),
                               br(),
                               actionButton("genecard_link", "GeneCards", 
                                            class = "d_dotbutt",
                                            width = "100%",
                                            style='padding:4px; font-size:80%')
                               )
                               )
                     ),
    column(9,
           tabsetPanel(
             id = "dot_plot_types",
             dotplot_pg_celltype,
             dotplot_pg_age,
             dotplot_pg_sample
             )
             )
  )         
)

# creating the heatmap page along with all the contents
# there are 2 tab items in this page - average expression and percent expression
# function for average/percent expression tab for heatmap
heatmap_tab_func <- function(tabpan_title, cellt_in1, gene_in1, reset_in, top_most_label,
                             cellt_in2, gene_in2, color_in, scale_in, clust_col_in,
                             clust_row_in, plot_out, down_plot_out, plot_data, down_plot_data,
                             heatmap_down_cond, heatmap_data_down_cond) {
  heatmap_tab <- tabPanel(title = tabpan_title, 
                          style = "background-color: #fff; border-color: #D3D3D3;",
                          fluidRow(
                            column(3, 
                                   wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                             wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                                       selectizeInput(cellt_in1,
                                                                      label = "Select Cell Type(s) : ",
                                                                      choices = NULL,
                                                                      options = list(create = TRUE), 
                                                                      multiple = TRUE,
                                                                      selected = NULL
                                                       ),
                                                       selectizeInput(gene_in1, 
                                                                      label = "Select Gene(s) : ",
                                                                      choices = NULL, 
                                                                      options = list(create = TRUE),
                                                                      multiple = TRUE, 
                                                                      selected = NULL
                                                       )
                                             ),
                                             h5(tags$i("Click here to use the Heatmap option given below :")),
                                             actionButton(reset_in, "Reset", 
                                                          class = "d_dotbutt"),
                                             br(),
                                             br(),
                                             wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                                       h5(tags$b(top_most_label
                                                       )),
                                                       selectizeInput(cellt_in2,
                                                                      label = "Select Cell Type : ",
                                                                      choices = NULL,
                                                                      options = list(create = TRUE),
                                                                      selected = NULL,
                                                                      multiple = FALSE
                                                       ),
                                                       sliderInput(gene_in2,
                                                                   label = "Select Top No. of Genes to Display : ",
                                                                   min = 1, max = 500,
                                                                   value = 25)
                                                       )
                                             )
                          ),
                          column(9,
                                 wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                           fluidRow(
                                             column(3,
                                                    selectizeInput(color_in,
                                                                   label = "Select Color : ",
                                                                   choices = list( 
                                                                     "Red",
                                                                     "Blue",
                                                                     "Green",
                                                                     "Grey", 
                                                                     "Orange",
                                                                     "Purple",
                                                                     "Spectral",
                                                                     "Viridis",
                                                                     "Magma" ,
                                                                     "Heat colors"),
                                                                   multiple = FALSE,
                                                                   selected = NULL
                                                    )
                                             ),
                                             column(1),
                                             column(3,
                                                    selectizeInput(scale_in,
                                                                   label = "Scale : ",
                                                                   choices = list("column", "row", "none"),
                                                             #      options = list(create = TRUE),
                                                                   multiple = FALSE,
                                                                   selected = "none"
                                                    )
                                             ),
                                             column(1),
                                             column(4,
                                                    h5(tags$b("Cluster the Heatmap by :")),
                                                    tags$style(".pretty.p-default input:checked~.state 
                                                               label:after {background-color: #D1B2D1 !important;
                                                               font-family: Arial !important;}"),
                                                    fluidRow(
                                                      column(4,
                                                             prettyCheckbox(
                                                               inputId = clust_col_in,
                                                               label = "column", 
                                                               value = FALSE,
                                                               shape = "curve",
                                                               status = "default"
                                                             )
                                                      ),
                                                      column(3,
                                                             prettyCheckbox(
                                                               inputId = clust_row_in,
                                                               label = "row", 
                                                               value = FALSE,
                                                               shape = "curve",
                                                               status = "default"
                                                             )   
                                                             
                                                      )
                                                    )
                                                    )      
                                           ),
                                           box(title = "Heatmap", solidHeader = TRUE, collapsible = TRUE, 
                                               status = "primary", 
                                               collapsed = FALSE, br(), width = "100%",
                                               plotOutput(plot_out, width = "100%")
                                           ),
                                           br(),
                                           conditionalPanel(heatmap_down_cond, 
                                                            downloadButton(down_plot_out, "Download Heatmap", 
                                                                           class = "d_dotbutt")),
                                           br()
                          ),
                          
                          wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                    box(title = "Heatmap Corresponding Data", 
                                        solidHeader = TRUE, collapsible = TRUE, 
                                        status = "primary", 
                                        collapsed = TRUE, br(), width = "100%",
                                        DT::dataTableOutput(plot_data, width = "100%")
                                    ),
                                    br(),
                                    conditionalPanel(heatmap_data_down_cond, 
                                                     downloadButton(down_plot_data, "Download Corresponding Data", 
                                                                    class = "d_dotbutt")),
                                    br()
                          )
  )
  )
  )
  return(heatmap_tab)
}

# average expression tab for heatmap
heatmap_avg_exprs_tab <- heatmap_tab_func("Average Expression", 'cell_type_hp_avg1',
                                          'genes_hp_avg1', "reset_heatmap_avg",
                                          "See Top Most Avg Gene Expression for a Particular Cell Type 
                                          and corresponding Avg Expression of those Genes in the remaining 
                                          Cell Types", 'cell_type_hp_avg2', 'genes_hp_avg2',
                                          'color_hp_avg', 'scale_hp_avg', 'clust_col_hp_avg',
                                          'clust_row_hp_avg', 'heatmap_plot_avg', "downloadPlotavg",
                                          'heatmap_data_avg', "downloadHmDavg", "output.heatmap_plot_avg", 
                                          "output.heatmap_data_avg")
                                  
# percent expression tab for heatmap
heatmap_perc_exprs_tab <- heatmap_tab_func("Percent Expression", 'cell_type_hp2',
                                           'genes_hp2', "reset_heatmap",
                                           "See Top Most % Gene Expression for a Particular Cell Type 
                                            and corresponding % Expression of those Genes in the remaining 
                                           Cell Types", 'cell_type_hp1', 'genes_hp1',
                                           'color_hp', 'scale_hp', 'clust_col_hp',
                                           'clust_row_hp', 'heatmap_plot', "downloadPlot",
                                           'heatmap_data', "downloadHmD", "output.heatmap_plot", 
                                           "output.heatmap_data")

# combining the two heatmap tabs in the heatmap page
heatmap_page <- tabItem(style = "font-family: Arial !important;",
                        tabName = "heatmap",
                        tabsetPanel(id = "heatmap_types",
                                    heatmap_avg_exprs_tab,
                                    heatmap_perc_exprs_tab
                        )
)

# creating the gene expression data page along with all the contents
# there are multiple tab items in this page, depending on the sample number
# each sample then has 2 tabs - one for average gene expression, and the other for percent gene expression
# creating a function for the average/percent gene expression tab for each sample
gene_exprs_tab_func <- function(tabpantitle, boxtitle, dtout_name, downbut_name_dt) {
  gene_exprs_tp <- tabPanel(title = tabpantitle,
           style = "background-color: #fff;
           border-color: #D3D3D3;",
           wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                     box(title = boxtitle, 
                         solidHeader = TRUE, collapsible = TRUE, 
                         status = "primary", 
                         collapsed = FALSE, br(), width = "100%",
                         DT::dataTableOutput(dtout_name, width = "100%"),
                         br(),
                         downloadButton(downbut_name_dt,"Download the (selected) data",
                                        class = "d_dotbutt")
                     )
           )
  )
  return(gene_exprs_tp)
}

# adding this function inside another function created for each sample tab
sample_tabpan_func <- function(samptab_title, tabpan_id,
                               tabpantitle_a, boxtitle_a, dtout_name_a, downbut_name_dt_a,
                               tabpantitle_p, boxtitle_p, dtout_name_p, downbut_name_dt_p) {
  sample_tabpan <- tabPanel(title = samptab_title,
                            style = "background-color: #fff; border-color: #D3D3D3;",
                            fluidRow(
                              column(12,
                                     tabsetPanel(id = tabpan_id,
                                                 gene_exprs_tab_func(tabpantitle_a, boxtitle_a, dtout_name_a, downbut_name_dt_a), 
                                                 gene_exprs_tab_func(tabpantitle_p, boxtitle_p, dtout_name_p, downbut_name_dt_p)
                                     )
                              )
                            )
  )
  return(sample_tabpan)
}

# combining everything in the gene expression data page
gene_exprs_data_page <- tabItem(
  tabName = "exprs_data",
  tabsetPanel(id = "expressiondata_types",
              sample_tabpan_func("Sample : D032", "d032_exp",
                                 "Average Expression", "Average Expression : Sample D032", 
                                 'd032_average', 'd032a_down',
                                 "Percent Expression", "Percent Expression : Sample D032",
                                 'd032_percent', 'd032p_down'),
              sample_tabpan_func("Sample : D046", "d046_exp",
                                 "Average Expression", "Average Expression : Sample D046",
                                 'd046_average', 'd046a_down',
                                 "Percent Expression", "Percent Expression : Sample D046",
                                 'd046_percent', 'd046p_down'),
              sample_tabpan_func("Sample : D062", "d062_exp",
                                 "Average Expression", "Average Expression : Sample D062",
                                 'd062_average', 'd062a_down',
                                 "Percent Expression", "Percent Expression : Sample D062",
                                 'd062_percent', 'd062p_down'),
              sample_tabpan_func("Sample : D088", "d088_exp",
                                 "Average Expression", "Average Expression : Sample D088",
                                 'd088_average', 'd088a_down',
                                 "Percent Expression", "Percent Expression : Sample D088",
                                 'd088_percent', 'd088p_down'),
              sample_tabpan_func("Sample : D122", "d122_exp",
                                 "Average Expression", "Average Expression : Sample D122",
                                 'd122_average', 'd122a_down',
                                 "Percent Expression", "Percent Expression : Sample D122",
                                 'd122_percent', 'd122p_down'),
              sample_tabpan_func("Sample : D139", "d139_exp",
                                 "Average Expression", "Average Expression : Sample D139",
                                 'd139_average', 'd139a_down',
                                 "Percent Expression", "Percent Expression : Sample D139",
                                 'd139_percent', 'd139p_down'),
              sample_tabpan_func("Sample : D150", "d150_exp",
                                 "Average Expression", "Average Expression : Sample D150",
                                 'd150_average', 'd150a_down',
                                 "Percent Expression", "Percent Expression : Sample D150",
                                 'd150_percent', 'd150p_down'),
              sample_tabpan_func("Sample : D175", "d175_exp",
                                 "Average Expression", "Average Expression : Sample D175",
                                 'd175_average', 'd175a_down',
                                 "Percent Expression", "Percent Expression : Sample D175",
                                 'd175_percent', 'd175p_down'),
              sample_tabpan_func("Sample : D231", "d231_exp",
                                 "Average Expression", "Average Expression  : Sample D231",
                                 'd231_average', 'd231a_down',
                                 "Percent Expression", "Percent Expression : Sample D231",
                                 'd231_percent', 'd231p_down')
  )
)
  
# the 'About' section page
about_section_page <- tabItem(tabName = "main_page",
                              fluidRow(style = "font-family: Arial !important;",
                                       column(12,
                                              wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                                        box(title = "About the Experiment", solidHeader = TRUE, collapsible = TRUE, 
                                                            status = "primary",
                                                            collapsed = FALSE,
                                                            br(), width = "100%",
                                                            c("Respiratory failure associated with COVID-19 has placed focus on 
                                                              the lung. Here, we present  single-nucleus  accessible  chromatin  
                                                              profiles  of  90,980  nuclei  and  matched single-nucleus  
                                                              transcriptomes  of  46,500  nuclei  in  healthy  lung  donors  of
                                                              ~30  weeks gestation,  ~3  years  and  ~30  years.  We  mapped 
                                                              398,395  candidate  cis-regulatory elements (cCREs) in lung cell 
                                                              types and linked distal cCREs to putative target genes and
                                                              leveraged this dataset to investigate loci associated with COVID-19.
                                                              At the SARS-CoV-2 host  entry  gene TMPRSS2,  we  identified  distal
                                                              cCREs  with age-increased  activity  in alveolar  type  2  cells 
                                                              which  had  immune  regulatory  signatures and  harbored  variants
                                                              associated with respiratory traits.  At the novel 3p21.31 COVID-19 
                                                              risk locus, a candidate variant overlapped a distal cCRE linked to 
                                                              SLC6A20, a gene expressed in alveolar cells which has known 
                                                              functional association with ACE2. Our findings provide insight into
                                                              the gene regulatory logic of lung cell types including likely 
                                                              drivers of the increased expression of SARS-CoV-2 host genes in 
                                                              adults, and establish a resource (lungepigenome.org) for 
                                                              interpreting lung disease risk."),
                                                            br(),
                                                            br(),
                                                            "Link to the paper : ", tags$a(href="https://elifesciences.org/articles/62522", 
                                                                                           "https://elifesciences.org/articles/62522"),
                                                            br(),
                                                            "Link to all the datasets : ", 
                                                            tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161383",
                                                                   "GEO: GSE161383")
                                                            )
                                                        ),
                                              wellPanel(style = "background-color: #fff; border-color: #D3D3D3;",
                                                        box(title = "Related Links", solidHeader = TRUE, collapsible = TRUE, status = "primary",
                                                            collapsed = TRUE,
                                                            br(), width = "100%",
                                                            "Xin Sun lab homepage : ", tags$a(href="http://xinsunlab.org", "http://xinsunlab.org"),
                                                            br(), 
                                                            "Kyle Gaulton lab : ", tags$a(href="http://www.gaultonlab.org", 
                                                                                          "http://www.gaultonlab.org"),
                                                            br(), 
                                                            "UCSD Center for Epigenomics : ", 
                                                            tags$a(href="https://medschool.ucsd.edu/som/cmm/research/epigenomics/pages/default.aspx", 
                                                                   "https://medschool.ucsd.edu/som/cmm/research/epigenomics/pages/default.aspx"),
                                                            br(), 
                                                            "LungMAP : ", tags$a(href="https://lungmap.net", "https://lungmap.net"),
                                                            br(), 
                                                            "Lungepigenome Website : ", tags$a(href="https://www.lungepigenome.org/", 
                                                                                               "https://www.lungepigenome.org/"),
                                                            br(), br(),
                                                            img(src="lungepigenome.png", width = "50%")
                                                        )
                                              )
                                       )
                              )
)

# combining all the pages in the 'body' of the app
body <- dashboardBody(style = "font-family: Arial !important;",
                      useShinyjs(),
                      tags$head(
                        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
                        tags$style(HTML('.logo {
                                        background-color: #D1B2D1 !important;
                                        font-family: Arial !important;
                                        }
                                        .navbar {
                                        background-color: #D1B2D1 !important;
                                        }
                                        body {
                                        font-family: Arial !important;
                                        }
                                        .box.box-solid.box-primary>.box-header {
                                        color:#000000;
                                        background:#D1B2D1;
                                        }
                                        .box.box-solid.box-primary{
                                        border-bottom-color:#D3D3D3;
                                        border-left-color:#D3D3D3;
                                        border-right-color:#D3D3D3;
                                        border-top-color:#D3D3D3;
                                        }
                                        .box-header .box-title, .box-header>.fa, 
                                        .box-header>.glyphicon, .box-header>.ion {
                                        font-family: "Arial", sans-serif !important;
                                        }
                                        .panel.panel-primary{
                                        border-color:#D3D3D3;
                                        }
                                        .panel.panel-success{
                                        border-color:#D3D3D3;
                                        }
                                        .panel-group>.panel.panel-success>.panel-heading {
                                        color:#000000;
                                        background:#D1B2D1;
                                        }
                                        .panel-group>.panel.panel-primary>.panel-heading {
                                        color:#000000;
                                        background:white;
                                        }
                                        .scatterD3 .gear-menu path {
                                        fill:rgb(177 1 193) !important;
                                        opacity: 1;
                                        }
                                        .scatterD3 .caption-icon path {
                                        fill:rgb(177 1 193) !important;
                                        opacity: 1;
                                        }
                                        .scatterD3 .axis text {
                                        font-family: Arial;
                                        }
                                        .bootstrap-switch .bootstrap-switch-handle-off.bootstrap-switch-primary, 
                                        .bootstrap-switch .bootstrap-switch-handle-on.bootstrap-switch-primary {
                                        color: #fff;
                                        background: #D1B2D1;
                                        }
                                        '))
    ),
  
  tabItems(
    dotplot_page,
    heatmap_page,
    gene_exprs_data_page,
    about_section_page
  )
)

# creating the ui side
ui <- function(request) {
  dashboardPage(
    skin = "black",
    dashboardHeader(title = "Gene Expression Profiling", titleWidth = 355),
    sidebar,
    body
  )
}

# creating the server side
server <- function(input, output, session) {
  
  #################################### ALL INPUT DATA ##########################################
  
  # creating functions for updating SelectizeInput, with and without the 'selected' option
  # without 'selected' option
  updateinp_func <- function(inputid, choices) {
    update_sel_in <- updateSelectizeInput(session = session, 
                                          inputId = inputid, 
                                          choices = choices, 
                                          server = TRUE)
    return(update_sel_in)
  }
  
  # with 'selected' option
  updateinp_func_w_select <- function(inputid, choices, selected_choice) {
    update_sel_in <- updateSelectizeInput(session = session, 
                                          inputId = inputid, 
                                          choices = choices, 
                                          selected = selected_choice,
                                          server = TRUE)
    return(update_sel_in)
  }
  
  # updating the cell type and gene inputs in dot plot
  updateinp_func_w_select('dp_gene', gene_names_old, "TMPRSS2")
  updateinp_func_w_select('dp_cellty', list("All Cell Types", 
                                            "Epithelial Cells" = c(epithelial_cells$CellTypes,
                                                                   "All Epithelial Cells" = 
                                                                     "Epithelial"),
                                            "Mesenchymal Cells" = c(mesenchymal_cells$CellTypes,
                                                                    "All Mesenchymal Cells" = "Mesenchymal"),
                                            "Endothelial Cells" = c(endothelial_cells$CellTypes,
                                                                    "All Endothelial Cells" = "Endothelial"),
                                            "Hematopoeitic Cells" = c(hematopoeitic_cells$CellTypes,
                                                                      "All Hematopoeitic Cells" = "Hematopoeitic")),
                          "All Cell Types")
  updateinp_func_w_select('gene_link_n', gene_names_old, "TMPRSS2")
  
  # rendering the UI Output for the gene selection for dot plot by sample tab
  output$gene_samp <- renderUI({
    req(input$dp_gene)
    selectizeInput(inputId = 'gene_samp_inp', 
                   label = "Select One Gene : ",
                   choices = input$dp_gene,
                   multiple = FALSE, 
                   selected = NULL
    )
  })
  
  # rendering the UI Output for the cell type selection for dot plot by sample tab
  output$cellty_samp <- renderUI({
    req(input$dp_cellty)
    selectizeInput(inputId = 'cellty_samp_inp', 
                   label = "Select One Cell Type : ",
                   choices = 
                     if (input$dp_cellty == "All Cell Types") {
                       cellt_order 
                     } else {
                       cellt_list1 <- data.frame()
                       cellt_list <- data.frame()
                       for (dpcellty in input$dp_cellty) {
                         if (dpcellty %in% c("Epithelial", "Mesenchymal", "Endothelial", "Hematopoeitic")) {
                           cellt_list <- all_cell_types %>% filter(Cells %in% dpcellty)
                           cellt_list1 <- rbind.data.frame(cellt_list1, cellt_list)
                         }
                         if (dpcellty %in% celltypes_old) {
                           cellt_list <- all_cell_types %>% filter(CellTypes %in% dpcellty)
                           cellt_list1 <- rbind.data.frame(cellt_list1, cellt_list)
                         }
                       }
                       cellt_list1$CellTypes 
                     },
                   multiple = FALSE, 
                   selected = NULL
    )
  })
  
  # updating the cell type and gene inputs in heatmap (average expression)
  updateinp_func('cell_type_hp_avg2', list("Cell Type" = cellt_order)) 
  updateinp_func_w_select('cell_type_hp_avg1', list("All Cell Types", "Cell Type" = cellt_order),
                          "All Cell Types")
  updateinp_func_w_select('genes_hp_avg1', gene_names_old, c("TMPRSS2", "ACE2"))
  
  # updating the cell type and gene inputs in heatmap (percent expression)
  updateinp_func('cell_type_hp1', list("Cell Type" = cellt_order))
  updateinp_func_w_select('cell_type_hp2', list("All Cell Types", "Cell Type" = cellt_order),
                          "All Cell Types")
  updateinp_func_w_select('genes_hp2', gene_names_old, c("TMPRSS2", "ACE2"))
  
  ############################### DOT PLOT PAGE ##############################################
  
  # creating dot plot by cell type, by age and by sample
  # first creating a function to encapsulate the 'color' option in the dot plot by cell type and
  # by age, in a reactive environment
  color_dotplot_func <- function(input_col) {
    color_dp_reactive <- reactive({
      switch(input_col,
             "Red" = "interpolateReds",
             "Blue" = "interpolateBlues",
             "Green" = "interpolateGreens",
             "Orange" = "interpolateOranges",
             "Purple" = "interpolatePurples",
             "Grey" = "interpolateGreys",
             "Purple - Orange" = "interpolatePuOr",
             "Orange - Red" = "interpolateOrRd",
             "Purple - Red" = "interpolatePuRd",
             "Red - Yellow - Blue" = "interpolateRdYlBu",
             "Yellow - Orange - Red" = "interpolateYlOrRd",
             "Viridis" = "interpolateViridis",
             "Spectral" = "interpolateSpectral",
             "Inferno" = "interpolateInferno",
             "Plasma" = "interpolatePlasma",
             "Warm" = "interpolateWarm",
             "Rainbow" = "interpolateRainbow"
      )
    })
    return(color_dp_reactive)
  }
  
  # dot plot by cell type and by age 'color' options in reactive environments using the above function
  color_dotplot_cell <- color_dotplot_func(input$color_cellt_dotp)
  color_dotplot_age <- color_dotplot_func(input$color_age_dotp)
  
  # arranging the celltypes in an order to be used in the same order in the dot plot (scatterD3)
  cellt_order1 <- data.frame("CellTypes" = c("alveolar.type.1.cells",
                                             "alveolar.type.2.cells",
                                             "AT1/AT2.like.cells",
                                             "AT2/club.like.cells",
                                             "basal.cells",
                                             "club.cells",
                                             "ciliated.cells",
                                             "goblet.cells",
                                             "PNECs",
                                             "matrix.fibroblast.1",
                                             "matrix.fibroblast.2",
                                             "myofibroblasts",
                                             "pericytes",
                                             "airway.smooth.muscle",
                                             "vascular.smooth.muscle",
                                             "chondrocytes",
                                             "arteries",
                                             "veins",
                                             "Cap1",
                                             "Cap2",
                                             "lymphatics",
                                             "bronchial.vessel",
                                             "Alveolar.macrophages",
                                             "Interstitial.macrophages",
                                             "monocytes",
                                             "mast.cells",
                                             "dendritic.cells",
                                             "T.cells",
                                             "B.cells",
                                             "NK.cells",
                                             "enucleated.erythrocytes"),
                             "CellTy_Num" = c(1:31))
  
  # function for dot plot by cell type (the starting part of which is common for the dot plot and
  # the corresponding data for the same)
  dotplot_cellt_func <- function(gene_input, celltype_input) {
    
    perc_cellt <- as.data.frame(t(exprs_mat[gene_input,]))
    perc_cellt$CellTypes <- rownames(perc_cellt)
    avg_cellt <- as.data.frame(t(exprs_mat_avg[gene_input,]))
    avg_cellt$CellTypes <- rownames(avg_cellt)
    avg_cellt <- tidyr::gather(avg_cellt[,c("CellTypes", gene_input)], 
                               key = "Gene", value = "Avg_Expression",
                               -"CellTypes")
    avg_cellt$comb_sample <- paste0(avg_cellt$CellTypes, "_", avg_cellt$Gene)
    perc_cellt <- tidyr::gather(perc_cellt[,c("CellTypes", gene_input)], 
                                key = "Gene", value = "Percent_Expression",
                                -"CellTypes")
    perc_cellt$comb_sample <- paste0(perc_cellt$CellTypes, "_", perc_cellt$Gene)
    perc_cellt <- perc_cellt %>% select(-c(CellTypes, Gene))
    cellt_exprs_comb <- left_join(avg_cellt, perc_cellt, by = "comb_sample")
    cellt_exprs_comb <- left_join(cellt_exprs_comb, all_cell_types, by = "CellTypes")
    cellt_exprs_comb <- left_join(cellt_order1, cellt_exprs_comb)
    cellt_exprs_comb <- cellt_exprs_comb %>% arrange(desc(CellTy_Num))
    
    if (celltype_input == "All Cell Types") {
      cellt_exprs_comb <- cellt_exprs_comb 
      cellt_exprs_comb$CellTypes <- factor(cellt_exprs_comb$CellTypes, 
                                           levels = unique(cellt_exprs_comb$CellTypes))
    } else {
      cellt_exp_df1 <- data.frame()
      cellt_exp_df <- data.frame()
      for (dpcellty in celltype_input) {
        if (dpcellty %in% c("Epithelial", "Mesenchymal", "Endothelial", "Hematopoeitic")) {
          cellt_exp_df <- cellt_exprs_comb %>% filter(Cells %in% dpcellty)
          cellt_exp_df1 <- rbind.data.frame(cellt_exp_df1, cellt_exp_df)
          cellt_exp_df1$CellTypes <- factor(cellt_exp_df1$CellTypes, 
                                            levels = unique(cellt_exp_df1$CellTypes))
        }
        if (dpcellty %in% celltypes_old) {
          cellt_exp_df <- cellt_exprs_comb %>% filter(CellTypes %in% dpcellty)
          cellt_exp_df1 <- rbind.data.frame(cellt_exp_df1, cellt_exp_df)
        }
      }
      cellt_exprs_comb <- cellt_exp_df1 
    }
    return(cellt_exprs_comb)
    
  }
  
  # dot plot by cell type
  output$dotplot_plot_cellt <- renderScatterD3(
    {
      req(input$dp_gene)
      req(input$dp_cellty)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
       
        cellt_exprs_comb <- dotplot_cellt_func(input$dp_gene, input$dp_cellty)
        scatterD3(x = cellt_exprs_comb$Gene, y = cellt_exprs_comb$CellTypes, 
                  col_var = cellt_exprs_comb$Avg_Expression,
                  size_var = cellt_exprs_comb$Percent_Expression,
                  xlab = "Gene", ylab = "Cell Types", col_lab = "Avg Expression",
                  size_lab = "% Expression",
                  tooltip_text = paste0("Avg Expression: ", 
                                        round(cellt_exprs_comb$Avg_Expression,3), "</strong><br /> ",
                                        "% Expression: ", 
                                        round(cellt_exprs_comb$Percent_Expression,3), "<br />",
                                        "Cell Type: ", cellt_exprs_comb$CellTypes),
                  caption = list(title = paste0("Gene : ", input$dp_gene), 
                                 subtitle = paste0("Cell Type : ", input$dp_cellty)),
                  left_margin = input$leftmarg_dotp_cellt,
                  col_continuous = TRUE,
                  colors = color_dotplot_cell()) 
    }
  )
  
  # dot plot by cell type data
  output$dotplot_data_cellt <- DT::renderDataTable(
    options = list(scrollX = TRUE),{
      req(input$dp_gene)
      req(input$dp_cellty)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      cellt_exprs_comb <- dotplot_cellt_func(input$dp_gene, input$dp_cellty)
      
      cellt_exprs_comb <- cellt_exprs_comb %>% select(-comb_sample)
      cellt_exprs_comb
    }
  )
  
  # download the dot plot by cell type data
  output$downloadDotPDCt <- downloadHandler(
    filename = function() {
      paste("expression_per_celltype_data.txt")
    },
    content = function(file) {
      req(input$dp_gene)
      req(input$dp_cellty)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      cellt_exprs_comb <- dotplot_cellt_func(input$dp_gene, input$dp_cellty)
      
      cellt_exprs_comb <- cellt_exprs_comb %>% select(-comb_sample)
      write.table(cellt_exprs_comb, file, quote = F)
    }
  )
  
  # function for dot plot by age (the starting part of which is common for the dot plot and
  # the corresponding data for the same)
  dotplot_age_func <- function(gene_input, celltype_input) {
    
    avg_age_df <- tidyr::gather(avg_exprs_per_age[,c("CellTypes", gene_input, "Age")], 
                                key = "Gene", value = "Avg_Expression",
                                -c("CellTypes", "Age"))
    avg_age_df$comb_samp <- paste0(avg_age_df$CellTypes, "_", avg_age_df$Age, "_",
                                   avg_age_df$Gene)
    percent_age_df <- tidyr::gather(percent_exprs_per_age[,c("CellTypes", gene_input,
                                                             "Age")], 
                                    key = "Gene", value = "Percent_Expression",
                                    -c("CellTypes", "Age"))
    percent_age_df$comb_samp <- paste0(percent_age_df$CellTypes, "_", percent_age_df$Age, "_",
                                       percent_age_df$Gene)
    percent_age_df <- percent_age_df %>% select(-c(CellTypes, Age, Gene))
    age_comb_df <- left_join(avg_age_df, percent_age_df, by = "comb_samp")
    age_comb_df$Age1 <- age_comb_df$Age
    age_comb_df$Age1[age_comb_df$Age1 == "30wk"] <- "w30"
    age_comb_df$Age1[age_comb_df$Age1 == "3yo"] <- "y3"
    age_comb_df$Age1[age_comb_df$Age1 == "30yo"] <- "y30"
    age_comb_df$Age_Gene <- paste0(age_comb_df$Gene, "_",  age_comb_df$Age1)
    age_comb_df <- left_join(age_comb_df, all_cell_types, by = "CellTypes")
    age_comb_df <- left_join(cellt_order1, age_comb_df)
    age_comb_df <- age_comb_df %>% arrange(desc(CellTy_Num))
    
    if (celltype_input == "All Cell Types") {
      age_comb_df <- age_comb_df 
      age_comb_df$CellTypes <- factor(age_comb_df$CellTypes, 
                                      levels = unique(age_comb_df$CellTypes))
    } else {
      age_exp_df1 <- data.frame()
      age_exp_df <- data.frame()
      for (dpcellty in celltype_input) {
        if (dpcellty %in% c("Epithelial", "Mesenchymal", "Endothelial", "Hematopoeitic")) 
        {
          age_exp_df <- age_comb_df %>% filter(Cells %in% dpcellty)
          age_exp_df1 <- rbind.data.frame(age_exp_df1, age_exp_df)
          age_exp_df1$CellTypes <- factor(age_exp_df1$CellTypes, 
                                          levels = unique(age_exp_df1$CellTypes))
        }
        if (dpcellty %in% celltypes_old) {
          age_exp_df <- age_comb_df %>% filter(CellTypes %in% dpcellty)
          age_exp_df1 <- rbind.data.frame(age_exp_df1, age_exp_df)
        }
      }
      age_comb_df <- age_exp_df1 
    }
    return(age_comb_df)
    
  }
  
  # dot plot by age
  output$dotplot_plot_age <- renderScatterD3(
    {
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$age_dotp_list)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      age_comb_df <- dotplot_age_func(input$dp_gene, input$dp_cellty)
      
      if (input$age_dotp_list == "All Ages") {
        age_comb_df <- age_comb_df 
      } else if (input$age_dotp_list != "All Ages") {
        age_comb_df <- age_comb_df %>% filter(Age %in% input$age_dotp_list)
      }
      
      symb_age_dp <- reactive({
        switch(input$symbol_age_dotp,
               "Gene" = age_comb_df$Gene,
               "Age" = age_comb_df$Age,
               "None" = NULL
        )
      })
      
      scatterD3(x = age_comb_df$Age_Gene, 
                y = age_comb_df$CellTypes, 
                col_var = age_comb_df$Avg_Expression, 
                size_var = age_comb_df$Percent_Expression,
                xlab = "Age & Gene", ylab = "Cell Type",
                col_lab = "Avg Expression", size_lab = "% Expression",
                left_margin = input$leftmarg_dotp_age,
                colors = color_dotplot_age(),
                symbol_var = symb_age_dp(),
                symbol_lab = input$symbol_age_dotp,
                col_continuous = TRUE,
                tooltip_text = paste0("Avg Expression: ", 
                                      round(age_comb_df$Avg_Expression,3), "</strong><br /> ",
                                      "% Expression: ", 
                                      round(age_comb_df$Percent_Expression,3), "</strong><br /> ",
                                      "Cell Type: ", age_comb_df$CellTypes, "</strong><br /> ",
                                      "Gene: ", age_comb_df$Gene, "<br />",
                                      "Age: ", age_comb_df$Age),
                caption = list(title = paste0("Gene : ", input$dp_gene), 
                               subtitle = paste0("Cell Type : ", input$dp_cellty,
                                                 " ; ",
                                                 "Age : ", input$age_dotp_list)))
      
    }
  )
  
  # dot plot by age data
  output$dotplot_data_age <- DT::renderDataTable(
    options = list(scrollX = TRUE),{
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$age_dotp_list)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      age_comb_df <- dotplot_age_func(input$dp_gene, input$dp_cellty)
      
      if (input$age_dotp_list == "All Ages") {
        age_comb_df <- age_comb_df 
      } else if (input$age_dotp_list != "All Ages") {
        age_comb_df <- age_comb_df %>% filter(Age %in% input$age_dotp_list)
      }
      
      age_comb_df <- age_comb_df %>% select(-c(Age_Gene, Age1, comb_samp))
      age_comb_df
      
    }
  )
  
  # download the dot plot by age data
  output$downloadDotPDAge <- downloadHandler(
    filename = function() {
      paste("expression_per_age_data.txt")
    },
    content = function(file) {
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$age_dotp_list)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      age_comb_df <- dotplot_age_func(input$dp_gene, input$dp_cellty)
      
      if (input$age_dotp_list == "All Ages") {
        age_comb_df <- age_comb_df 
      } else if (input$age_dotp_list != "All Ages") {
        age_comb_df <- age_comb_df %>% filter(Age %in% input$age_dotp_list)
      }
      
      age_comb_df <- age_comb_df %>% select(-c(Age_Gene, Age1, comb_samp))
      write.table(age_comb_df, file, quote = F)
      
    }
  )
  
  # for dot plot by sample, since the user has the ability to select either a static scatter plot
  # or an interactive dot plot to be plotted in the same space, rendering UI first
  output$dotplot_plot_sample <- renderUI({
    if (input$sc_plot_samp == FALSE) {
      scatterD3Output('scd3_sample_dp', width = "100%")
    } else if (input$sc_plot_samp == TRUE) {
      tagList(plotOutput('ggplot_sample_dp', width = "100%", height = 600),
              br(),
              downloadButton("downloadggplotsamp", "Download Scatter Plot", 
                             class = "d_dotbutt"))
    }
  })
  
  # function to wrap the earlier part of both the interactive and static dot plot 
  # for average expression
  dotplot_data_in_avg_func <- function(gene_in, cellty_in) {
    samp_scd3_data <- tidyr::gather(exprs_samp_avg[,c("CellTypes", "Age", "Sample",
                                                      gene_in)], 
                                    key = "Gene", value = "Avg_Expression",
                                    -c("CellTypes", "Age", "Sample"))
    samp_scd3_data <- samp_scd3_data %>% filter(CellTypes == cellty_in)
    samp_scd3_data$Age[samp_scd3_data$Age == "30wk"] <- "weeks_30"
    samp_scd3_data$Age[samp_scd3_data$Age == "3yo"] <- "years_3"
    samp_scd3_data$Age[samp_scd3_data$Age == "30yo"] <- "years_30"
    return(samp_scd3_data)
  }
  
  # for percent expression
  dotplot_data_in_perc_func <- function(gene_in, cellty_in) {
    samp_scd3_data_p <- tidyr::gather(exprs_samp_final[,c("CellTypes", "Age", "Sample",
                                                          gene_in)], 
                                      key = "Gene", value = "Percent_Expression",
                                      -c("CellTypes", "Age", "Sample"))
    samp_scd3_data_p <- samp_scd3_data_p %>% filter(CellTypes == cellty_in)
    samp_scd3_data_p$Age[samp_scd3_data_p$Age == "30wk"] <- "weeks_30"
    samp_scd3_data_p$Age[samp_scd3_data_p$Age == "3yo"] <- "years_3"
    samp_scd3_data_p$Age[samp_scd3_data_p$Age == "30yo"] <- "years_30"
    return(samp_scd3_data_p)
  }
  
  # dot plot by sample (interactive)
  output$scd3_sample_dp <- renderScatterD3(
    {
      
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$gene_samp_inp)
      req(input$cellty_samp_inp)
      req(input$sample_exprs_type)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      samp_scd3_data <- dotplot_data_in_avg_func(input$gene_samp_inp, input$cellty_samp_inp)
      samp_scd3_data_p <- dotplot_data_in_perc_func(input$gene_samp_inp, input$cellty_samp_inp)
      
      if (input$sample_exprs_type == "Average") {
        if (input$color_samp_dotp == "Color by Age") {
          scatterD3(x = samp_scd3_data$Age, y = samp_scd3_data$Avg_Expression, 
                    xlab = "Age", 
                    ylab = paste0("Average Expression of ", input$gene_samp_inp,
                                  " in ", input$cellty_samp_inp),
                    left_margin = input$leftmarg_dotp_sample,
                    col_var = samp_scd3_data$Age,
                    col_lab = "Age",
                    point_size = input$point_size_sample,
                    tooltip_text = paste0("Avg Expression: ", 
                                          round(samp_scd3_data$Avg_Expression,
                                                digits = 3), "</strong><br /> ",
                                          "Sample: ", samp_scd3_data$Sample),
                    caption = list(title = paste0("Gene : ", input$gene_samp_inp),
                                   subtitle = paste0("Cell Type : ", input$cellty_samp_inp))
          )
        } else if (input$color_samp_dotp != "Color by Age") {
          scatterD3(x = samp_scd3_data$Age, y = samp_scd3_data$Avg_Expression, 
                    xlab = "Age", 
                    ylab = paste0("Average Expression of ", input$gene_samp_inp,
                                  " in ", input$cellty_samp_inp),
                    left_margin = input$leftmarg_dotp_sample,
                    colors = input$color_samp_dotp,
                    point_size = input$point_size_sample,
                    tooltip_text = paste0("Avg Expression: ", 
                                          round(samp_scd3_data$Avg_Expression,
                                                digits = 3), "</strong><br /> ",
                                          "Sample: ", samp_scd3_data$Sample),
                    caption = list(title = paste0("Gene : ", input$gene_samp_inp),
                                   subtitle = paste0("Cell Type : ", input$cellty_samp_inp))
          )
        }
      } else if (input$sample_exprs_type == "Percent") {
        if (input$color_samp_dotp == "Color by Age") {
          scatterD3(x = samp_scd3_data_p$Age, y = samp_scd3_data_p$Percent_Expression, 
                    xlab = "Age", 
                    ylab = paste0("% Expression of ", input$gene_samp_inp,
                                  " in ", input$cellty_samp_inp),
                    left_margin = input$leftmarg_dotp_sample,
                    col_var = samp_scd3_data_p$Age,
                    col_lab = "Age",
                    point_size = input$point_size_sample,
                    tooltip_text = paste0("% Expression: ", 
                                          round(samp_scd3_data_p$Percent_Expression,
                                                digits = 3), "</strong><br /> ",
                                          "Sample: ", samp_scd3_data_p$Sample),
                    caption = list(title = paste0("Gene : ", input$gene_samp_inp),
                                   subtitle = paste0("Cell Type : ", input$cellty_samp_inp))
          )
        } else if (input$color_samp_dotp != "Color by Age") {
          scatterD3(x = samp_scd3_data_p$Age, y = samp_scd3_data_p$Percent_Expression, 
                    xlab = "Age", 
                    ylab = paste0("% Expression of ", input$gene_samp_inp,
                                  " in ", input$cellty_samp_inp),
                    left_margin = input$leftmarg_dotp_sample,
                    colors = input$color_samp_dotp,
                    point_size = input$point_size_sample,
                    tooltip_text = paste0("% Expression: ", 
                                          round(samp_scd3_data_p$Percent_Expression,
                                                digits = 3), "</strong><br /> ",
                                          "Sample: ", samp_scd3_data_p$Sample),
                    caption = list(title = paste0("Gene : ", input$gene_samp_inp),
                                   subtitle = paste0("Cell Type : ", input$cellty_samp_inp))
          )
        }
      }
      
    })
  
  # wrapping dot plot by sample (static) in a function since it will be used later while 
  # saving the plot as well
  dotplot_samp_stat_func <- function(gene_in, cellty_in, samp_exprs_ty, color_in,
                                     point_size_in, left_marg_in) {
    
    samp_scd3_data <- dotplot_data_in_avg_func(gene_in, cellty_in)
    samp_scd3_data_p <- dotplot_data_in_perc_func(gene_in, cellty_in)
    
    samp_scd3_data_mean <- aggregate(samp_scd3_data[, "Avg_Expression"], 
                                     list(samp_scd3_data$Age), mean)
    colnames(samp_scd3_data_mean) <- c("Age", "Mean")
    samp_scd3_data_sd <- aggregate(samp_scd3_data[, "Avg_Expression"], 
                                   list(samp_scd3_data$Age), sd)
    colnames(samp_scd3_data_sd) <- c("Age", "Std_Dev")
    mean_and_sd <- left_join(samp_scd3_data_sd, samp_scd3_data_mean)
    samp_scd3_data <- left_join(samp_scd3_data, mean_and_sd)
    
    samp_scd3_data_p_mean <- aggregate(samp_scd3_data_p[, "Percent_Expression"], 
                                       list(samp_scd3_data_p$Age), mean)
    colnames(samp_scd3_data_p_mean) <- c("Age", "Mean")
    samp_scd3_data_p_sd <- aggregate(samp_scd3_data_p[, "Percent_Expression"], 
                                     list(samp_scd3_data_p$Age), sd)
    colnames(samp_scd3_data_p_sd) <- c("Age", "Std_Dev")
    mean_and_sd_p <- left_join(samp_scd3_data_p_sd, samp_scd3_data_p_mean)
    samp_scd3_data_p <- left_join(samp_scd3_data_p, mean_and_sd_p)
    
    if (samp_exprs_ty == "Average") {
      if (color_in == "Color by Age") {
        dp_samp <- ggplot(data = samp_scd3_data, aes(y = Avg_Expression, x = Age)) + 
          geom_point(aes(color = Age), size = point_size_in/50) +
          stat_summary(aes(color = Age), fun = mean, fun.min = mean, fun.max = mean,
                       geom = "crossbar", width = 0.3, fatten = 1.5) +
          geom_errorbar(aes(ymin=Mean-Std_Dev, 
                            ymax=Mean+Std_Dev, color = Age), width=.2,
                        position=position_dodge(.9)) + theme_bw() +
          theme(plot.margin = unit(c(5.5,5.5,5.5,(left_marg_in/2)), "pt"),
                plot.title = element_text(hjust = 0.5, face = "italic"),
                axis.text.x = element_text(face = "bold"),
                axis.text.y = element_text(face = "bold"),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12)) +
          ylab(paste0(samp_exprs_ty, " of ", cellty_in)) + 
          ggtitle(gene_in)
      } else if (color_in != "Color by Age") {
        dp_samp <- ggplot(data = samp_scd3_data, aes(y = Avg_Expression, x = Age)) + 
          geom_point(color = color_in, size = point_size_in/50) +
          stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                       geom = "crossbar", width = 0.3, fatten = 1.5) +
          geom_errorbar(aes(ymin=Mean-Std_Dev, 
                            ymax=Mean+Std_Dev), width=.2,
                        position=position_dodge(.9)) + theme_bw() +
          theme(plot.margin = unit(c(5.5,5.5,5.5,(left_marg_in/2)), "pt"),
                plot.title = element_text(hjust = 0.5, face = "italic"),
                axis.text.x = element_text(face = "bold"),
                axis.text.y = element_text(face = "bold"),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12)) +
          ylab(paste0(samp_exprs_ty, " of ", cellty_in)) + 
          ggtitle(gene_in)
      }
    } else if (samp_exprs_ty == "Percent") {
      if (color_in == "Color by Age") {
        dp_samp <- ggplot(data = samp_scd3_data_p, aes(y = Percent_Expression, x = Age)) + 
          geom_point(aes(color = Age), size = point_size_in/50) +
          stat_summary(aes(color = Age), fun = mean, fun.min = mean, fun.max = mean,
                       geom = "crossbar", width = 0.3, fatten = 1.5) +
          geom_errorbar(aes(ymin=Mean-Std_Dev, 
                            ymax=Mean+Std_Dev, color = Age), width=.2,
                        position=position_dodge(.9)) + theme_bw() +
          theme(plot.margin = unit(c(5.5,5.5,5.5,(left_marg_in/2)), "pt"),
                plot.title = element_text(hjust = 0.5, face = "italic"),
                axis.text.x = element_text(face = "bold"),
                axis.text.y = element_text(face = "bold"),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12)) +
          ylab(paste0(samp_exprs_ty, " of ", cellty_in)) + 
          ggtitle(gene_in)
      } else if (color_in != "Color by Age") {
        dp_samp <- ggplot(data = samp_scd3_data_p, aes(y = Percent_Expression, x = Age)) + 
          geom_point(color = color_in, size = point_size_in/50) +
          stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                       geom = "crossbar", width = 0.3, fatten = 1.5) +
          geom_errorbar(aes(ymin=Mean-Std_Dev, 
                            ymax=Mean+Std_Dev), width=.2,
                        position=position_dodge(.9)) + theme_bw() +
          theme(plot.margin = unit(c(5.5,5.5,5.5,(left_marg_in/2)), "pt"),
                plot.title = element_text(hjust = 0.5, face = "italic"),
                axis.text.x = element_text(face = "bold"),
                axis.text.y = element_text(face = "bold"),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12)) +
          ylab(paste0(samp_exprs_ty, " of ", cellty_in)) + 
          ggtitle(gene_in)
      }
    }
    
    return(dp_samp)
    
  }
  
  # dot plot by sample (static scatter plot)
  output$ggplot_sample_dp <- renderPlot({
    
    req(input$dp_gene)
    req(input$dp_cellty)
    req(input$gene_samp_inp)
    req(input$cellty_samp_inp)
    req(input$sample_exprs_type)
    
    for (i in input$dp_gene) {
      if (!(i %in% gene_names_old)) {
        validate("Data on at least one of the input genes are not found in this dataset")
      }
    }
    
    for (i in input$dp_cellty) {
      if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                     "All Hematopoeitic Cells", "All Endothelial Cells",
                     "All Mesenchymal Cells", "All Epithelial Cells"))) {
        validate("Data on at least one of the input cell types are not found in this dataset")
      }
    }
    
    dotplot_samp_stat_func(input$gene_samp_inp, input$cellty_samp_inp, 
                                   input$sample_exprs_type, input$color_samp_dotp,
                                   input$point_size_sample, input$leftmarg_dotp_sample)
    
  })
  
  # download the dot plot by sample (static scatter plot)
  output$downloadggplotsamp <- downloadHandler(
    filename = function() { paste('scatterplot.png') },
    content = function(file) {
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$gene_samp_inp)
      req(input$cellty_samp_inp)
      req(input$sample_exprs_type)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      scp <- dotplot_samp_stat_func(input$gene_samp_inp, input$cellty_samp_inp, 
                             input$sample_exprs_type, input$color_samp_dotp,
                             input$point_size_sample, input$leftmarg_dotp_sample)
      
      ggsave(file, plot = scp)
    }
  )
  
  # function for displaying the dot plot by sample data and downloading it
  dotplot_sample_func <- function(gene_inp, cellt_inp, samp_type_inp) {
    
    samp_scd3_data <- tidyr::gather(exprs_samp_avg[,c("CellTypes", "Age", "Sample",
                                                      gene_inp)], 
                                    key = "Gene", value = "Avg_Expression",
                                    -c("CellTypes", "Age", "Sample"))
    samp_scd3_data <- samp_scd3_data %>% filter(CellTypes == cellt_inp)
    samp_scd3_data_p <- tidyr::gather(exprs_samp_final[,c("CellTypes", "Age", "Sample",
                                                          gene_inp)], 
                                      key = "Gene", value = "Percent_Expression",
                                      -c("CellTypes", "Age", "Sample"))
    samp_scd3_data_p <- samp_scd3_data_p %>% filter(CellTypes == cellt_inp)
    if (samp_type_inp == "Average") {
      samp_dat <- samp_scd3_data
    } else if (samp_type_inp == "Percent") {
      samp_dat <- samp_scd3_data_p
    }
    return(samp_dat)
  }
  
  # dot plot by sample data
  output$dotplot_data_sample <- DT::renderDataTable(
    options = list(scrollX = TRUE),{
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$gene_samp_inp)
      req(input$cellty_samp_inp)
      req(input$sample_exprs_type)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      dotplot_sample_func(input$gene_samp_inp, input$cellty_samp_inp, input$sample_exprs_type)
      
    }
  )
  
  # download the dot plot by sample data
  output$downloadDotPDSamp <- downloadHandler(
    filename = function() {
      paste("expression_per_sample_data.txt")
    },
    content = function(file) {
      req(input$dp_gene)
      req(input$dp_cellty)
      req(input$gene_samp_inp)
      req(input$cellty_samp_inp)
      req(input$sample_exprs_type)
      
      for (i in input$dp_gene) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$dp_cellty) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      write.table(dotplot_sample_func(input$gene_samp_inp, input$cellty_samp_inp, 
                                      input$sample_exprs_type), file, quote = F)
    }
  )
  
  # reset to default for dot plot
  observeEvent(input$reset_def_dp, {
    updateSelectizeInput(session = session, 
                         inputId = 'dp_gene', 
                         choices = gene_names_old, 
                         selected = "TMPRSS2",
                         server = TRUE)
    updateSelectizeInput(session = session, 
                         inputId = 'dp_cellty', 
                         choices =list("All Cell Types", 
                                       "Epithelial Cells" = c(epithelial_cells$CellTypes,
                                                              "All Epithelial Cells" = "Epithelial"),
                                       "Mesenchymal Cells" = c(mesenchymal_cells$CellTypes,
                                                               "All Mesenchymal Cells" = "Mesenchymal"),
                                       "Endothelial Cells" = c(endothelial_cells$CellTypes,
                                                               "All Endothelial Cells" = "Endothelial"),
                                       "Hematopoeitic Cells" = c(hematopoeitic_cells$CellTypes,
                                                                 "All Hematopoeitic Cells" = "Hematopoeitic")), 
                         selected = "All Cell Types",
                         server = TRUE)
  })
  
  # UMAP link, Region Overlap link and GeneCards Website link
  onclick("umap_link", runjs(paste0("window.open('https://www.genome-browser.lungepigenome.org/cellBrowser/?ds=lungsnrna&gene=",
                                    input$gene_link_n, "')")))
  onclick("lungep_link", runjs(paste0("window.open('https://www.lungepigenome.org/region-search/?annotation=HGNC%3A11876&region=",
                                      input$gene_link_n,
                                      "+(homo+sapiens)&genome=GRCh37", "')")))
  onclick("genecard_link", runjs(paste0("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                                        input$gene_link_n, "')")))
  
  # bookmarking the page so that the selected cell type and gene can be seen in the UMAP 
  # (the linkout in UMAP of lungepigenome that goes to the dotplot page of this RShiny App)
  setBookmarkExclude(c("lungep_link", "heatmap_types", "d122_exp", "d231_exp", "method_corr",
                       "reset_heatmap", "d046_exp", "d150_exp", "clear_gene_corr", "color_corrp",
                       "corr_gene", "genes_hp1", "point_size_sample", "clust_row_hp", "d062_exp",
                       "sidebarCollapsed", "d175_exp", "d139_exp", "umap_link", "cell_type_hp1",
                       "cell_type_hp2", "corr_cell", "sc_plot_samp", "correlationplot_types",
                       "genecard_link", "gene_link_n", "sidebarItemExpanded", "genes_hp2",
                       "age_dotp_list", "symbol_age_dotp", "color_age_dotp", "sample_exprs_type",
                       "d088_exp", "type_corr", "tabs", "color_cellt_dotp", "color_samp_dotp",
                       "scale_hp", "leftmarg_dotp_sample", "d032_exp", "expressiondata_types",
                       "dp_cellty", "clear_cell_corr", "dot_plot_types", "leftmarg_dotp_cellt",
                       "leftmarg_dotp_age", "reset_def_dp", "order_corr", "color_hp", 
                       "clust_col_hp"))
  
  observe({
    reactiveValuesToList(input)
    session$doBookmark()
  })
  
  onBookmarked(function(url) {
    updateQueryString(url)
  })
  
  onRestored(function(state) {
    updateSelectizeInput(session, "dp_gene", 
                         selected=state$input$dp_gene, 
                         choices=c(gene_names_old, state$input$dp_gene), 
                         server=TRUE)
  })
  
  #################################### HEATMAP PAGE ############################################
  
  # creating heatmap for average and percent expression data
  # wrapping the color input (reactive) in a function (since it is the same for both the data)
  color_heatmap_in_func <- function(color_inp_exprs_ty) {
    color_in <- reactive({
      switch(color_inp_exprs_ty,
             "Red" = colorRampPalette(brewer.pal(9, "Reds"))(100),
             "Blue" = colorRampPalette(brewer.pal(9, "Blues"))(100),
             "Green" = colorRampPalette(brewer.pal(9, "Greens"))(100),
             "Grey" = colorRampPalette(brewer.pal(9, "Greys"))(100), 
             "Orange" = colorRampPalette(brewer.pal(9, "Oranges"))(100),
             "Purple" = colorRampPalette(brewer.pal(9, "Purples"))(100),
             "Spectral" = colorRampPalette(brewer.pal(9, "Spectral"))(100),
             "Viridis" = viridis(100),
             "Magma" = magma(100),
             "Heat colors" = heat.colors(100)
      )
    })
    return(color_in)
  }
  
  # color input for average and percent expression heatmap
  color_heatmap <- color_heatmap_in_func(input$color_hp)
  color_heatmap_avg <- color_heatmap_in_func(input$color_hp_avg)
  
  # wrapping the heatmap of average and percent expression in a function (since it is quite similar
  # for the two types of data, and will also be used for downloading the heatmap)
  heatmap_output_func <- function(exprs_data_type, cellty_in_op1, gene_in_op1,
                                  clust_row_in, color_in, scale_in, clust_col_in,
                                  cellty_in_op2, gene_in_op2) {
    
    exprs_mat_tr <- as.data.frame(t(exprs_data_type))
    
    if (length(cellty_in_op1) > 0 | length(gene_in_op1) > 0) {
      req(cellty_in_op1)
      req(gene_in_op1)
      if (cellty_in_op1 == "All Cell Types") {
        if ((clust_row_in == TRUE) & (length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes to Cluster")
        } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Genes to Scale")
        } else {
          hp_out <- pheatmap(exprs_data_type[gene_in_op1,], color = color_in, 
                   scale = scale_in, cluster_cols = clust_col_in, 
                   cluster_rows = clust_row_in)
        }
      } else {
        if ((clust_row_in == TRUE) & (length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes to Cluster")
        } else if ((clust_col_in == TRUE) & (length(cellty_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Cell Types to Cluster")
        } else if ((clust_col_in == TRUE & length(cellty_in_op1) < 2) & 
                   (clust_row_in == TRUE & length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes and Cell Types to Cluster")
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == FALSE & length(gene_in_op1) == 1) {
          hp_out <- validate("Please select Two or More Genes or Cell Types")
        } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Genes to Scale")
        } else if ((scale_in == "column") & length(cellty_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Cell Types to Scale")
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == FALSE & length(gene_in_op1) > 1) {
          hp_out <- pheatmap(exprs_mat_tr[cellty_in_op1,gene_in_op1], color = color_in,
                   scale = scale_in, cluster_cols = clust_col_in, 
                   cluster_rows = clust_row_in)
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == TRUE & length(gene_in_op1) > 1) {
          hp_out <- pheatmap(exprs_mat_tr[cellty_in_op1,gene_in_op1], color = color_in,
                             scale = scale_in, cluster_cols = clust_row_in, 
                             cluster_rows = clust_col_in)
        } else {
          hp_out <- pheatmap(exprs_data_type[gene_in_op1,cellty_in_op1], color = color_in,
                   scale = scale_in, cluster_cols = clust_col_in, 
                   cluster_rows = clust_row_in)
        }
      }
    } else if (length(cellty_in_op1) == 0 & length(gene_in_op1) == 0) {
      req(cellty_in_op2)
      req(gene_in_op2)
      exprs_mat1 <- exprs_data_type %>% arrange(desc(exprs_data_type[,cellty_in_op2]))
      exprs_mat1 <- exprs_mat1[1:gene_in_op2,]
      if ((clust_row_in == TRUE) & (gene_in_op2 < 2)) {
        hp_out <- validate("Please select Two or More Genes to Cluster")
      } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
        hp_out <- validate("Please select Two or More Genes to Scale")
      } else {
        hp_out <- pheatmap(exprs_mat1, color = color_in, scale = scale_in, 
                 cluster_cols = clust_col_in, cluster_rows = clust_row_in)
      }
    }
    
    return(hp_out)
    
  }
  
  # average expression heatmap
  output$heatmap_plot_avg <- renderPlot({
    
    for (i in input$genes_hp_avg1) {
      if (!(i %in% gene_names_old)) {
        validate("Data on at least one of the input genes are not found in this dataset")
      }
    }
    
    for (i in input$cell_type_hp_avg1) {
      if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                     "All Hematopoeitic Cells", "All Endothelial Cells",
                     "All Mesenchymal Cells", "All Epithelial Cells"))) {
        validate("Data on at least one of the input cell types are not found in this dataset")
      }
    }
    
    heatmap_output_func(exprs_mat_avg, input$cell_type_hp_avg1, input$genes_hp_avg1,
                        input$clust_row_hp_avg, color_heatmap_avg(), input$scale_hp_avg, 
                        input$clust_col_hp_avg, input$cell_type_hp_avg2, input$genes_hp_avg2)
    
  })
  
  # percent expression heatmap
  output$heatmap_plot <- renderPlot({
    
    for (i in input$genes_hp2) {
      if (!(i %in% gene_names_old)) {
        validate("Data on at least one of the input genes are not found in this dataset")
      }
    }
    
    for (i in input$cell_type_hp2) {
      if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                     "All Hematopoeitic Cells", "All Endothelial Cells",
                     "All Mesenchymal Cells", "All Epithelial Cells"))) {
        validate("Data on at least one of the input cell types are not found in this dataset")
      }
    }
    
    heatmap_output_func(exprs_mat, input$cell_type_hp2, input$genes_hp2,
                        input$clust_row_hp, color_heatmap(), input$scale_hp, input$clust_col_hp,
                        input$cell_type_hp1, input$genes_hp1)
    
  })
  
  # average expression heatmap download
  output$downloadPlotavg <- downloadHandler(
    filename = function() { paste('heatmap_avg.png') },
    content = function(file) {
      png(file)
      
      for (i in input$genes_hp_avg1) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp_avg1) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      heatmap_output_func(exprs_mat_avg, input$cell_type_hp_avg1, input$genes_hp_avg1,
                          input$clust_row_hp_avg, color_heatmap_avg(), input$scale_hp_avg, 
                          input$clust_col_hp_avg, input$cell_type_hp_avg2, input$genes_hp_avg2)
      
      dev.off()
    }
  )
  
  # percent expression heatmap download
  output$downloadPlot <- downloadHandler(
    filename = function() { paste('heatmap.png') },
    content = function(file) {
      png(file)
      
      for (i in input$genes_hp2) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp2) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      heatmap_output_func(exprs_mat, input$cell_type_hp2, input$genes_hp2,
                          input$clust_row_hp, color_heatmap(), input$scale_hp, input$clust_col_hp,
                          input$cell_type_hp1, input$genes_hp1)
      
      dev.off()
    }
  )
  
  # resetting heatmap inputs so that the second type of heatmap can be displayed (in both average
  # and percent expression)
  # average expression
  observeEvent(input$reset_heatmap_avg, {
    shinyjs::reset("cell_type_hp_avg1") 
    shinyjs::reset("genes_hp_avg1") 
  })
  
  # percent expression
  observeEvent(input$reset_heatmap, {
    shinyjs::reset("cell_type_hp2") 
    shinyjs::reset("genes_hp2") 
  })
  
  # function to wrap the average and percent heatmap expression data table (will also be useful
  # while downloading the data)
  heatmap_data_out_func <- function(exprs_data_type, cellty_in_op1, gene_in_op1,
                                         clust_row_in, scale_in, clust_col_in,
                                         cellty_in_op2, gene_in_op2) {
    
    exprs_mat_tr <- as.data.frame(t(exprs_data_type))
    
    if (length(cellty_in_op1) > 0 | length(gene_in_op1) > 0) {
      req(cellty_in_op1)
      req(gene_in_op1)
      if (cellty_in_op1 == "All Cell Types") {
        if ((clust_row_in == TRUE) & (length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes to Cluster")
        } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Genes to Scale")
        } else {
          hp_out <- exprs_data_type[gene_in_op1,]
        }
      } else {
        if ((clust_row_in == TRUE) & (length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes to Cluster")
        } else if ((clust_col_in == TRUE) & (length(cellty_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Cell Types to Cluster")
        } else if ((clust_col_in == TRUE & length(cellty_in_op1) < 2) & 
                   (clust_row_in == TRUE & length(gene_in_op1) < 2)) {
          hp_out <- validate("Please select Two or More Genes and Cell Types to Cluster")
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == FALSE & length(gene_in_op1) == 1) {
          hp_out <- validate("Please select Two or More Genes or Cell Types")
        } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Genes to Scale")
        } else if ((scale_in == "column") & length(cellty_in_op1) < 2) {
          hp_out <- validate("Please select Two or More Cell Types to Scale")
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == FALSE & length(gene_in_op1) > 1) {
          hp_out <- exprs_mat_tr[cellty_in_op1,gene_in_op1]
        } else if (clust_col_in == FALSE & length(cellty_in_op1) == 1 & 
                   clust_row_in == TRUE & length(gene_in_op1) > 1) {
          hp_out <- exprs_mat_tr[cellty_in_op1,gene_in_op1]
        } else {
          hp_out <- exprs_data_type[gene_in_op1,cellty_in_op1]
        }
      }
    } else if (length(cellty_in_op1) == 0 & length(gene_in_op1) == 0) {
      req(cellty_in_op2)
      req(gene_in_op2)
      exprs_mat1 <- exprs_data_type %>% arrange(desc(exprs_data_type[,cellty_in_op2]))
      exprs_mat1 <- exprs_mat1[1:gene_in_op2,]
      if ((clust_row_in == TRUE) & (gene_in_op2 < 2)) {
        hp_out <- validate("Please select Two or More Genes to Cluster")
      } else if ((scale_in == "column") & length(gene_in_op1) < 2) {
        hp_out <- validate("Please select Two or More Genes to Scale")
      } else {
        hp_out <- exprs_mat1
      }
    }
    
    return(hp_out)
    
  }
  
  # average expression heatmap corresponding data table
  output$heatmap_data_avg <- DT::renderDataTable(
    options = list(scrollX = TRUE),
    {
      for (i in input$genes_hp_avg1) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp_avg1) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      heatmap_data_out_func(exprs_mat_avg, input$cell_type_hp_avg1, input$genes_hp_avg1,
                            input$clust_row_hp_avg, input$scale_hp_avg, 
                            input$clust_col_hp_avg, input$cell_type_hp_avg2, input$genes_hp_avg2)
      
    }
  )
  
  # percent expression heatmap corresponding data table
  output$heatmap_data <- DT::renderDataTable(
    options = list(scrollX = TRUE),
    {
      for (i in input$genes_hp2) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp2) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      heatmap_data_out_func(exprs_mat, input$cell_type_hp2, input$genes_hp2,
                            input$clust_row_hp, input$scale_hp, input$clust_col_hp,
                            input$cell_type_hp1, input$genes_hp1)
      
    }
  )
  
  # average expression heatmap corresponding data download
  output$downloadHmDavg <- downloadHandler(
    filename = function() {
      paste("heatmap_data_avg.txt")
    },
    content = function(file) {
      
      for (i in input$genes_hp_avg1) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp_avg1) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      write.table(heatmap_data_out_func(exprs_mat_avg, input$cell_type_hp_avg1, input$genes_hp_avg1,
                                        input$clust_row_hp_avg, input$scale_hp_avg, 
                                        input$clust_col_hp_avg, input$cell_type_hp_avg2, 
                                        input$genes_hp_avg2),
                  file, quote = F)
    }
  )
  
  # percent expression heatmap corresponding data download
  output$downloadHmD <- downloadHandler(
    filename = function() {
      paste("heatmap_data.txt")
    },
    content = function(file) {
      for (i in input$genes_hp2) {
        if (!(i %in% gene_names_old)) {
          validate("Data on at least one of the input genes are not found in this dataset")
        }
      }
      
      for (i in input$cell_type_hp2) {
        if (!(i %in% c(cells_list, celltypes_old, "All Cell Types",
                       "All Hematopoeitic Cells", "All Endothelial Cells",
                       "All Mesenchymal Cells", "All Epithelial Cells"))) {
          validate("Data on at least one of the input cell types are not found in this dataset")
        }
      }
      
      write.table(heatmap_data_out_func(exprs_mat, input$cell_type_hp2, input$genes_hp2,
                                        input$clust_row_hp, input$scale_hp, input$clust_col_hp,
                                        input$cell_type_hp1, input$genes_hp1),
                  file, quote = F)
    }
  )
  
  ######################## GENE EXPRESSION DATA TABLES #######################################
  
  # function for displaying the average and percent gene expression data table by sample
  gene_exprs_data_tab_func <- function(data_tab_exprs) {
    output_name <- DT::renderDataTable(
      options = list(scrollX = TRUE),
      data_tab_exprs
    )
    return(output_name)
  }
  
  # data table output per sample : average expression
  output$d032_average <- gene_exprs_data_tab_func(d032)
  output$d046_average <- gene_exprs_data_tab_func(d046)
  output$d062_average <- gene_exprs_data_tab_func(d062)
  output$d088_average <- gene_exprs_data_tab_func(d088)
  output$d122_average <- gene_exprs_data_tab_func(d122)
  output$d139_average <- gene_exprs_data_tab_func(d139)
  output$d150_average <- gene_exprs_data_tab_func(d150)
  output$d175_average <- gene_exprs_data_tab_func(d175)
  output$d231_average <- gene_exprs_data_tab_func(d231)
  
  # data table output per sample : percent expression
  output$d032_percent <- gene_exprs_data_tab_func(d032_p)
  output$d046_percent <- gene_exprs_data_tab_func(d046_p)
  output$d062_percent <- gene_exprs_data_tab_func(d062_p)
  output$d088_percent <- gene_exprs_data_tab_func(d088_p)
  output$d122_percent <- gene_exprs_data_tab_func(d122_p)
  output$d139_percent <- gene_exprs_data_tab_func(d139_p)
  output$d150_percent <- gene_exprs_data_tab_func(d150_p)
  output$d175_percent <- gene_exprs_data_tab_func(d175_p)
  output$d231_percent <- gene_exprs_data_tab_func(d231_p)
  
  # function for downloading the entire or selected part of the data table
  download_data_func <- function(file_name, data_table, data_tab_id) {
    out_down <- downloadHandler(
      filename = function() {
        paste(file_name, '.csv', sep = '')
      },
      content = function(file){
        write.csv(data_table[input[[data_tab_id]], ],file)
      }
    )
    
    return(out_down)
  }
  
  # download data table output per sample : average expression
  output$d032a_down <- download_data_func('sample_d032_avg_exprs', d032, "d032_average_rows_all")
  output$d046a_down <- download_data_func('sample_d046_avg_exprs', d046, "d046_average_rows_all")
  output$d062a_down <- download_data_func('sample_d062_avg_exprs', d062, "d062_average_rows_all")
  output$d088a_down <- download_data_func('sample_d088_avg_exprs', d088, "d088_average_rows_all")
  output$d122a_down <- download_data_func('sample_d122_avg_exprs', d122, "d122_average_rows_all")
  output$d139a_down <- download_data_func('sample_d139_avg_exprs', d139, "d139_average_rows_all")
  output$d150a_down <- download_data_func('sample_d150_avg_exprs', d150, "d150_average_rows_all")
  output$d175a_down <- download_data_func('sample_d175_avg_exprs', d175, "d175_average_rows_all")
  output$d231a_down <- download_data_func('sample_d231_avg_exprs', d231, "d231_average_rows_all")
  
  # download data table output per sample : percent expression
  output$d032p_down <- download_data_func('sample_d032_percent_exprs', d032_p, "d032_percent_rows_all") 
  output$d046p_down <- download_data_func('sample_d046_percent_exprs', d046_p, "d046_percent_rows_all")
  output$d062p_down <- download_data_func('sample_d062_percent_exprs', d062_p, "d062_percent_rows_all")
  output$d088p_down <- download_data_func('sample_d088_percent_exprs', d088_p, "d088_percent_rows_all")
  output$d122p_down <- download_data_func('sample_d122_percent_exprs', d122_p, "d122_percent_rows_all")
  output$d139p_down <- download_data_func('sample_d139_percent_exprs', d139_p, "d139_percent_rows_all")
  output$d150p_down <- download_data_func('sample_d150_percent_exprs', d150_p, "d150_percent_rows_all")
  output$d175p_down <- download_data_func('sample_d175_percent_exprs', d175_p, "d175_percent_rows_all")
  output$d231p_down <- download_data_func('sample_d231_percent_exprs', d231_p, "d231_percent_rows_all")
  
}

enableBookmarking(store = "url")
shinyApp(ui = ui, server = server)
