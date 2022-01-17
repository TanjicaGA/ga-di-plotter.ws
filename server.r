options( stringsAsFactors=FALSE )

suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(devtools)
    library(ga.data)
    library(ga.gamap)
    library(ga.utils)
    library(ga.gamapqc)
    library(foreach)
    library(ggplot2)
    library(purrr)
    library(stringr)
    library(ga.software)
    library(ga.software.dd)
    library(dplyr)
    library(ga.diversityindex)
   # library(bettertrace)
    ## load_all("~/git/R-packages/ga.software/")
    ## load_all("~/git/R-packages/ga.software.dd/")
})

testfile <- file.path( Sys.getenv("BIOINFORMATICS"), "Experiments/AID/AID20/AID200625/Data/PlateScans/AID200625 L2001 LX1819 LAR original.csv" )
testfile_R <- file.path( Sys.getenv("BIOINFORMATICS"), "Production/Kitlot/R2103/Data/PlateScans/Q2-015 R2103 LX1819 LO.csv" )
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


run.div.from.plate.data<-function(x,input,...)
{
  
  args<- list(
    x=x,
    batch=input$kitlot,
    qc.check.qcc30 = input$qcc30_filter,
    ...
    
  )
  do.call(divind, args)
}


translate.probes <- function( probes, mode=c("probe","phylum","bacteria") ) {
    mode <- match.arg( mode )
    prn <- paste( probe.numbers( probes ) )
    switch(
        mode,
        probe = probe.codes(probes),
        phylum = prn,
        bacteria = unlist(report.bacteria.names())[prn]
    )
}

DEBUG <- FALSE


## Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

    probeAnnotations <- c(
        probe="Probe Codes",
        phylum="Phylum Numbers",
        bacteria="Bacteria Names"
    )

    input_file <- function() {
        b <- input$bc_file
        if( DEBUG ) {
            b$datapath <- testfile
            b$name <- basename(testfile)
        }
        b
    }
    

    ## FUNCTIONS
    

    qcc <- function(qn) {

        pd <- req(plateData())

        dd <- data.frame(
            QCC    = pd$Sample,
            DI        = din(),
            div    = div(),
            row.names = NULL
        )
        nn <- dd$QCC

        i <- grepl( paste0("^",qn), nn )

        d <- dd[i,]
        names(d)[1] <- qn

        return( d )
    }
    qctable <- function(qn) {

        pd <- req(plateData())
        pd$DI <- sprintf( "%.2f", req(din()) )
        pd$div <- sprintf( "%.2f", req(div()) )
        platform <- pd$Platform[1]

        pd2 <- pd[ grep( paste0("^",qn), pd$Sample ),  ]

        i <- grep( "BLANK|UNI05|HYC01|DI|div", colnames(pd2), value=TRUE )
        pd3 <- pd2[, c("Sample","Row","Col","Well",i)]

        j <- grepl("BLANK", colnames(pd3))
        if( all(is.na(pd3[,j]))) {
            pd3 <- pd3[, !j]
        }

        xr <- probe.data( pd2 )

        if( platform == "Lx200" ) {
            pd3 <- pd3[, colnames(pd3) %!~% "BLANK[12]" ]
            xr <- xr[, colnames(xr) %!in% lx200.missing.probes() ]
        }

        pd3$Total <- sprintf( "%.0f", rowSums(xr, na.rm=TRUE) )

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
    currentProbeAnnotation <- function() {
        j <- ga.utils::coalesce( input$probe_labels %% 3 + 1, 1 )
        names(probeAnnotations)[ j ]
    }
    ddQcTables <- function() {
        ##ADDED NEW FUNCTION FOR QC RANGES IN R KIT
        
        if(grepl("^R",input$kitlot)){
                
             
               pd<-req(plateData())
	       pd$Platform=rep("lx200.RUOII",length(pd$Platform))	
          ##   	observeEvent(input$bc_file, {print(paste0("R names: ", pd))})      
	     ## cat(file=stderr(),is.recursive(pd))

         }
        else{    
        
        pd <- req(plateData())
	##observeEvent(input$bc_file, {print(paste0("L names: ", colnames(pd)))})
	}
        set.dd.qc.ranges()
        qc <- abundancy.table.qc( pd, start.from="file", batch=input$kitlot, report.per.sample=FALSE, kitlots=input$kitlot, variant="aa",qc.check.qcc30=input$qc_all_filter )
     ##	observeEvent(input$kitlot, {print(paste0("pd is : ", pd))})
        clear.dd.qc.ranges()
        
        
        qc.data <- attr( qc, "qc.data" )[["1"]]

        ## say( jsonlite::toJSON( qc.data, pretty=TRUE, auto_unbox=TRUE ))

        sum.probes <- function( .x, limit, cumsum=TRUE ) {
            l <- length( .x$probes[[limit]] )
            if( cumsum ) {
                n <- as.numeric( sub("\\D","", limit) )
                if( n < 3 )
                    l <- l + sum.probes( .x, paste0("±",n+1), cumsum=cumsum )
            }
            return( l )
        }
        list.probes <- function( .x, limit, cumsum=TRUE ) {
            p <- unique( translate.probes( .x$probes[[limit]], currentProbeAnnotation() ) )
            if( cumsum ) {
                n <- as.numeric( sub("\\D","", limit) )
                if( n < 3 )
                    p <- unique( append(p, list.probes( .x, paste0("±",n+1), cumsum=cumsum ) ))
            }
            return( p )
        }

        set.criteria <- function( set ) {
            l <- bacteria.table.qc.parameters(input$kitlot)[[set]]$limits
            paste( imap_chr( l, ~{ paste0( paste0("[ ±",.y), " <= " , .x, " ]" ) } ), collapse=" & " )
        }

        q <- imap_dfr( qc.data[c("QCC23","QCC33")], ~{
            print( .x )
            cbind.data.frame(
                Set = .y,
                Sample = imap_chr( .x, ~ .x$sample ),
                
                "±1" = imap_int( .x, ~ sum.probes(.x,"±1") ),
                "±2" = imap_int( .x, ~ sum.probes(.x,"±2") ),
                "±3" = imap_int( .x, ~ sum.probes(.x,"±3") ),
                "±1 List of probes" = imap_chr( .x, ~ paste(list.probes(.x,"±1"),collapse=", ") ),
                "±2 List of probes" = imap_chr( .x, ~ paste(list.probes(.x,"±2"),collapse=", ") ),
                "±3 List of probes" = imap_chr( .x, ~ paste(list.probes(.x,"±3"),collapse=", ") ),
                Result = c( "Fail", "Pass" )[ 1+imap_lgl( .x, ~ .x$qc ) ],
                Criteria = set.criteria( .y )
            )
        })

        q$Set[ duplicated(q$Set) ] <- ""
        q

    }
    

    plateData <- reactive({

        bc.file <- req(input_file())
        
        x_run=gamap(  
            x       = bc.file$datapath,
            stop.at = "file"
        )
       x_run
    })
    

    di <- reactive({
        pd <- req( plateData() )
        run.gamap.from.plate.data( x=pd, input, stop.at="dysbiosis")
    })
    bt <- reactive({
        pd <- req( plateData() )

        b <- gamap.probe.levels.ga(
            x              = pd,
            start.from     = "file",
            batch          = input$kitlot,
            qc.check.qcc30 = input$qcc30_filter
        )
        bl <- bacteria.limits( dont.warn.missing.revision = TRUE, revision="rev5" )
        colnames(b) <- bl$Bacteria

        sn <- rownames(b)

        b[ !grepl("QCC(30|29)",sn), ]


    })
    din <- reactive({
        pd <- req( plateData() )
        run.gamap.from.plate.data( x=pd, input)
    })
    
    div<-reactive({
      bc.file <- req(input_file())
      run.div.from.plate.data(x=bc.file$datapath,input)
    })

    ## Page 1: DI table, BT and DI-Plot:

    output$BacteriaTable <- renderTable( bt(), digits=0, rownames=TRUE )
    i.qcc <- reactive( grepl( "^QCC", names(din()) ) )

    for( n in "Other" ) {
        output[[paste0(n,"table")]] <- renderTable({

            dd <- data.frame(
                Sample    = names(din()),
                DI        = din(),
                div       = div(),
                row.names = NULL
            )

            dd[!i.qcc(),]

        })
    }

    output$diPlot <- renderPlot({

        di <- req(di())

        t2 <- attr( di, "T2" )
        qres <- attr( di, "Qres" )

        names(t2) <- names(qres) <- names(di())

        di <- di()

        args <-
            list(
                t2, qres,

                main = basename( input_file()$name ),

                unit.scale = input$unit_scale,
                plot.names = input$plot_names,
                log.scale = input$log_scale,
                add.di.scale = input$add_di_scale,
                add.threshold = input$add_threshold,
                add.normal.backdrop = input$show_normals
            )

        do.call( plot_di, args )

    })
    output$downloadResults <- downloadHandler(

        filename = function() {
            sub( "\\.csv$", "-results.csv", input_file()$name, ignore.case=TRUE )
        },

        content = function(file) {
            d <- req(din())
            d <- d[ !grepl("QCC(30|29)",names(d)) ]
            b <- bt()
            bl <- bacteria.limits()
            colnames(b) <- probe.numbers( bl$Probe )
            dd <- cbind.data.frame( Sample=names(d), DI=d, b )
            write.csv( dd, file, row.names=FALSE )
        }

    )

    observeEvent( input_file(), {
        if( is.null(input_file()) ) {
            shinyjs::hideElement( id="download-button-column" )
        } else {
            shinyjs::showElement( id="download-button-column" )
        }
    })
    ## output$showDownload <- eventReactive( input_file(), {!is.null(input_file())}, ignoreInit=TRUE )

    ## Page 2: QC tables

    output$qcTableQCC30 <- renderTable(qctable("QCC30"), digits=1)
    output$qcTableQCC29 <- renderTable(qctable("QCC29"), digits=1)
    output$qcTableQCC23 <- renderTable(qctable("QCC23"), digits=1)
    output$qcTableQCC33 <- renderTable(qctable("QCC33"), digits=1)

    output$QCC23table <- renderTable(qcc("QCC23"))
    output$QCC33table <- renderTable(qcc("QCC33"))
    output$QCC29table <- renderTable(qcc("QCC29"))
    output$QCC30table <- renderTable(qcc("QCC30"))

    ## Page 3: DD
    plot_dd_qc_sample <- function( sname ) {
        pd <- req(plateData())
        rx <- sname
        if(grepl("^R",input$kitlot)){
           
	   
            pd$Platform=rep("lx200.RUOII",length(pd$Platform)) 	
            p<- try(
            
            plot_abundancy_qc(
                pd, start.from="file", kitlot=input$kitlot,
                sample_rx = rx, exact=TRUE,
                probenames=currentProbeAnnotation(),
                use.aa=TRUE,
                qc.check.qcc30=input$qc_all_filter,
                bacteria.table.revision="rev5"
            ) + ggtitle( sname )
        )
        if(inherits(p,"try-error")) {
            return(NULL)
        }
        p}
	else
	{ ##observeEvent(input$kitlot, {print(paste0("kitlot: ", input$kitlot))})
	    p <- try(
            plot_abundancy_qc(
                pd, start.from="file", kitlot=input$kitlot,
                sample_rx = rx, exact=TRUE,
                probenames=currentProbeAnnotation(),
                use.aa=TRUE,
                qc.check.qcc30=input$qc_all_filter,
                bacteria.table.revision="rev5"
            ) + ggtitle( sname )
        )
        if(inherits(p,"try-error")) {
            return(NULL)
        }
        p
    }}


    ## dd_qc_plots <- eventReactive( input$bc_file, {
    dd_qc_plots <- reactive( {

        pd <- req(plateData())
        set.qcc.indeces(pd, verbose=FALSE)
        l <- lapply(
            pd$Sample[ c(which(i.qcc23),which(i.qcc33)) ],
            plot_dd_qc_sample
        )
        names(l) <- pd$Sample[ c(which(i.qcc23),which(i.qcc33)) ]
        l <- l[ !sapply( l, is.null ) ]

        l

    })

    observeEvent( input$bc_file, {

        f <- input$bc_file$name

        if( is.character(f) && length(f) > 0 ) {
            try( {
                kl <- kitlot.from.filename(f)
                if(!is.na(kl)) {
                    updateSelectInput( session, "kitlot", selected=kl )
                }
            }, silent=TRUE )
        }

    })

    observeEvent( list(input$bc_file,input$kitlot,input$probe_labels), {
    ## observe({
        l <- req(dd_qc_plots())
        iwalk( l, ~{
            output_name <- paste0( "dd-qc-plot-", .y )
            output[[output_name]] <- renderPlot(.x)
        })
    })

    output$dd_qc_plots <- renderUI({

        l <- req( dd_qc_plots() )

        p_list <- imap( l, ~{

            tagList(
                plotOutput(
                    outputId = paste0("dd-qc-plot-",.y)
                ),
                br()
            )
        })

        tagList(p_list)

    })

    ## output$ddQcTables <- renderTable(ddQcTables())
    output$ddQcTables <- renderTable(ddQcTables())

    ## dd probe button text
    ## input$probe_labels %% 3 + 1
    observeEvent( input$probe_labels, {
        new.label <- as.character( probeAnnotations[ ga.utils::coalesce(input$probe_labels %% 3 + 1,1) ] )
        updateActionButton( session, "probe_labels", label=new.label )
    })
    output$probeButton <- renderUI({
        actionButton( "probe_labels", label=probeAnnotations[1] )
    })

    ## Page 4: Probe Annotations

    pr <- lx200.probes( include.technical=FALSE )

    probe.a <- data.frame(
        Probe = pr,
        PhylumNr = as.character( probe.numbers( pr ) ),
        Annotation = bacteria.names(pr )
    )

    output$ProbeAnnotations <- renderDataTable(
        probe.a,
        options = list(
            pageLength = 50
        )
    )

})
