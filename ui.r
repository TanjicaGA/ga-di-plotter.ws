suppressPackageStartupMessages({
    library(shiny)
    library(ga.data)
  #  library(bettertrace)
    library(shinyjs)
})

shinyUI(fluidPage(

    useShinyjs(),
    ## Application title
    ## headerPanel("Plate Analyzer"),

    tags$head(
        tags$style("#diPlot {height:60vh !important;}"),
        tags$link( rel="stylesheet", type="text/css", href="extra-styles.css" )
    ),

    tabsetPanel(

        tabPanel(

            "DI-Plot",

            wellPanel(
                fluidRow(

                    column(
                        4,
                        fileInput("bc_file", label=h4("Plate Analyzer v4.0.0"),
                                  accept="text/csv" )
                    ),
                    column(
                        4,
                        selectInput(
                            "kitlot",
                            label=h4("Kitlot"),
                            choices=local({
                                ab <- available.batches()
                                o <- order( grepl("^[RLK]",ab), ab, decreasing=TRUE )
                                as.list(ab[o])
                            }),
                            selected = 1
                        )
                    ),
                    shinyjs::hidden(column(
                                 4,
                                 h4("DI + Bact table"),
                                 downloadButton("downloadResults","Download"),
                                 id="download-button-column"
                             ))
                )
            ),

            fluidRow(

                ## COL 1-9
                column(12,
                       plotOutput("diPlot", height="500px"),
                       ## Controls:
                       wellPanel(
                           fluidRow(
                               column(2,checkboxInput("plot_names", label="Use names", value=FALSE)),
                               column(2,checkboxInput("show_normals", label="BG Normal", value=FALSE)),
                               column(1,checkboxInput("add_di_scale", label="DI scale", value=TRUE)),
                               column(1,checkboxInput("add_threshold", label="DI border", value=TRUE)),
                               column(1,checkboxInput("unit_scale", label="Unit scale", value=TRUE)),
                               column(1,checkboxInput("log_scale", label="Log scale", value=TRUE)),
                               column(1,checkboxInput("qcc30_filter", label="QCC30 QC", value=TRUE))
                           )
                       )
                       )
            ),

            fluidRow(

                column( 2, tableOutput("QCC33table") ),
                column( 2, tableOutput("QCC23table") ),
                column( 2, tableOutput("QCC29table") ),
                column( 2, tableOutput("QCC30table") ),
                column( 2, tableOutput("Othertable") )

            ),

            fluidRow(
                column( 12, tableOutput("BacteriaTable") )
            )

        ),

        tabPanel(
            "QC Overview",

            wellPanel(
                fluidRow(
                    column( 4, div(tableOutput("qcTableQCC30"),class="qcTable" ))
                ),
                fluidRow(
                    column( 4, tableOutput("qcTableQCC29") )
                ),
                fluidRow(
                    column( 4, tableOutput("qcTableQCC23") )
                ),
                fluidRow(
                    column( 4, tableOutput("qcTableQCC33") )
                )
            )


        ),

        tabPanel(

            "Bacteria Abundance",

            wellPanel(

                fluidRow(
                    column(10, tableOutput("ddQcTables")),
                    column(2, uiOutput("probeButton")),
		              
                ),
                fluidRow(
                    column(1,checkboxInput("akkermansia_pos", label="change_Akkermansia_limit", value=FALSE))
                ),

                fluidRow(
                    div(
                        id="dd-qc-container",
                        uiOutput(
                            outputId="dd_qc_plots"
                        )
                    )
                )
            )
        ),

        tabPanel(

            "Probe Annotations",

            fluidRow(
                column(
                    12,
                    dataTableOutput("ProbeAnnotations")
                )
            )

        ),
        tabPanel(

          "DIV plot",
         
         

                fluidRow(
                  column(12,
                         plotOutput("DIvsDIVplot", height="500px")
                  )
                )
            

       ),
       tabPanel(
         
         "Trending",
         fluidRow(
           column(12,
                  plotOutput("trending", height="2000px")
           )
         )
         
         
       )
    )

    ## tags$script(HTML("console.log(output)"))


))
