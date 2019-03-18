library(shiny)
library(ga.data)
library(ga.biocode)
library(devtools)
library(ga.gamap)
library(ga.utils)
library(ga.gamapqc)
library(foreach)

library(bettertrace)

testfile <- "~/GA/Experiments/Q2-011/Data/01-Raw/Q2-011 L1802_II LX1734 KrG.csv"

run.gamap.from.plate.data <- function(x, input, stop.at, ... ) {

    args <- list(
        x=x,
        batch = input$kitlot,
        start.from = "file",
        qc.check.qcc30 = input$qcc30_filter,
        ...
    )

    if(!missing(stop.at))
        args$stop.at <- stop.at

    do.call( gamap, args )

}

## Define server logic required to draw a histogram
shinyServer(function(input, output) {

    plateData <- reactive({

        bc.file <- input$bc_file
        ## bc.file$datapath <- testfile

        if( !is.null(bc.file$datapath) ) {
            gamap(
                x       = bc.file$datapath,
                stop.at = "file"
            )
        }
        else
            NULL

    })

    di <- reactive({

        pd <- plateData()

        if( !is.null(pd) ) {
            run.gamap.from.plate.data( x=pd, input, stop.at="dysbiosis")
        }

    })

    bt <- reactive({

        pd <- plateData()

        if( !is.null(pd) ) {
            b <- gamap.probe.levels(
                x       = pd,
                start.from = "file",
                batch   = input$kitlot,
                qc.check.qcc30 = input$qcc30_filter
            )
            bl <- bacteria.limits()
            colnames(b) <- bl$Bacteria

            sn <- rownames(b)

            b[ !grepl("QCC(30|29)",sn), ]

        }
        else {
            NULL
        }

    })

    din <- reactive({

        pd <- plateData()

        if( !is.null( pd )  ) {
            run.gamap.from.plate.data( x=pd, input)
        }

    })

    qcc <- function(qn) {

        pd <- plateData()

        if(!is.null(pd)) {

            dd <- data.frame(
                QCC    = pd$Sample,
                DI        = din(),
                row.names = NULL
            )
            nn <- dd$QCC

            i <- grepl( paste0("^",qn), nn )

            d <- dd[i,]
            names(d)[1] <- qn

            return( d )

        }
    }

    qc.check <- function( qn, probe, value, platform ) {

        qcset <- gamap.qc.ranges(platform)

        if( qc == "QCC30" ) {

        }

    }

    qctable <- function(qn) {

        pd <- plateData()

        if( !is.null(pd) ) {

            platform <- pd$Platform[1]

            pd2 <- pd[ grep( paste0("^",qn), pd$Sample ),  ]

            i <- grep( "BLANK|UNI05|HYC01", colnames(pd2), value=TRUE )
            pd3 <- pd2[, c("Sample","Row","Col","Well",i)]
            xr <- probe.data( pd2 )

            if( platform == "Lx200" ) {
                pd3 <- pd3[, colnames(pd3) %!~% "BLANK[12]" ]
                xr <- xr[, colnames(xr) %!in% lx200.missing.probes() ]
            }

            pd3$Total <- rowSums(xr, na.rm=TRUE)

            o <- order(
                pd3$Sample %!~% "^QCC30",
                pd3$Sample %!~% "^QCC29",
                pd3$Sample %!~% "^QCC23",
                pd3$Sample %!~% "^QCC33"
            )

            pd4 <- pd3[o,]

            for( v in c("Row","Col","Well") ) {
                pd4[,v] <- paste(pd4[,v])
            }

            cn <- colnames(pd4)

            pd5 <- foreach( i=1:ncol(pd4), .combine = cbind.data.frame ) %do% {
                v <- pd4[,i]
                d <- data.frame(v=v)
                colnames(d) <- cn[i]

                if( cn[i] == "UNI05" ) {
                    d <- cbind.data.frame( d, "" )
                    colnames(d)[2] <- "?"
                }


                d

            }

            colnames(pd5)[ colnames(pd5) == "?" ] <- ""
            pd5

        }

    }

    output$qcTableQCC30 <- renderTable(qctable("QCC30"), digits=1)
    output$qcTableQCC29 <- renderTable(qctable("QCC29"), digits=1)
    output$qcTableQCC23 <- renderTable(qctable("QCC23"), digits=1)
    output$qcTableQCC33 <- renderTable(qctable("QCC33"), digits=1)

    output$QCC23table <- renderTable(qcc("QCC23"))
    output$QCC33table <- renderTable(qcc("QCC33"))
    output$QCC29table <- renderTable(qcc("QCC29"))
    output$QCC30table <- renderTable(qcc("QCC30"))

    output$BacteriaTable <- renderTable( bt(), digits=0, rownames=TRUE )

    i.qcc <- reactive( grepl( "^QCC", names(din()) ) )

    for( n in "Other" ) {
        output[[paste0(n,"table")]] <- renderTable({

            dd <- data.frame(
                Sample    = names(din()),
                DI        = din(),
                row.names = NULL
            )

            dd[!i.qcc(),]

        })
    }

    output$diPlot <- renderPlot({

        bc.file <- input$bc_file
        ## bc.file$datapath <- testfile
        ## bc.file$name <- testfile

        if( !is.null(di()) ) {

            t2 <- attr( di(), "T2" )
            qres <- attr( di(), "Qres" )

            names(t2) <- names(qres) <- names(di())

            di <- di()

            args <-
                list(
                    t2, qres,

                    main = basename( bc.file$name ),

                    unit.scale = input$unit_scale,
                    plot.names = input$plot_names,
                    log.scale = input$log_scale,
                    add.di.scale = input$add_di_scale,
                    add.threshold = input$add_threshold,
                    add.normal.backdrop = input$show_normals
                )

            do.call( plot_di, args )

        }

    })

    output$downloadResults <- downloadHandler(

        filename = function() {
        sub( "\\.csv$", "-results.csv", input$bc_file, ignore.case=TRUE )
    },

    content = function(file) {
        d <- din()
        d <- d[ !grepl("QCC(30|29)",names(d)) ]
        b <- bt()
        bl <- bacteria.limits()
        colnames(b) <- probe.numbers( bl$Probe )
        if( !is.null(d) && !is.null(b) ) {
            dd <- cbind.data.frame( Sample=names(d), DI=d, b )
            write.csv( dd, file, row.names=FALSE )
        }
    }

    )

})
