library(shiny)

SERVER <- '/srv/shiny-server/DimaApp'
SERVERPROT <- '/srv/shiny-server/DimaApp/MeltonProteomics'
SERVERRNA <- '/srv/shiny-server/DimaApp/MeltonRNAseq'


Sweave(paste(SERVERPROT, 'ds_functions_21JAN17.Rnw', sep='/'))
Sweave(paste(SERVERPROT, 'ds_resources_21JAN17.Rnw', sep='/'))
Sweave(paste(SERVERPROT, 'ds_analysis_21JAN17.Rnw', sep='/'))
Sweave(paste(SERVERPROT, 'ds_analysis_A_2FEB17.Rnw', sep='/'))
Sweave(paste(SERVERRNA, 'ds_rnaseqFunctions_13MAR17.Rnw', sep='/'))
Sweave(paste(SERVERRNA, 'ds_rnaseqAnalysis_5MAR17.Rnw', sep='/'))



acclst <- list(s227=r227.df$Accession, s238=r238.df$Accession, s239=r239.df$Accession, s243=r243.df$Accession) # rows WITHOUT QUANTITATION REMOVED
xl <-  makeSubsetExprMatrix(1, adat.lst) # adat.lst in resources:label=msnsetfromsubsets, EACH data.frame HAS ADDITIONAL FEATURES COLUMN (don't confuse with acclst)
msl <- makeMSnSetfromSampleSubsets(xl, p.df) # p.df in resources:label=MSnSetcombined

isAcc <- function(idstring) {
    # test whether input ids are Uniprot accessions or symbols
    # returns TRUE if input is Uniprot accessions otherwise FALSE
    ids <- unlist(strsplit(unlist(strsplit(idstring, split=' ')), split=','))
    x <- grep('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', ids)
    return(any(x != 0))
}

## UI
ui <- fluidPage(
    titlePanel('Proteomic and Genomic Data Tool', windowTitle='ProteinAndGeneProfiles'),
    sidebarLayout(
        sidebarPanel(
            conditionalPanel(
                condition = "input.thetabs == 1 || input.thetabs == 2 || input.thetabs == 3 || input.thetabs == 4" ,
                textInput(inputId='geneprot.id', 'UniProt Accessions or Gene Symbols to Plot', value='D2IYK3 J7H3Y9 Q9H6I2 P52945 P78426 P01308 P01275'),
                checkboxInput('log', 'Abundance log transformed?', value=TRUE, width=NULL),
                verbatimTextOutput('info')),
            ##submitButton('Submit')),
            conditionalPanel(
                condition = "input.thetabs == 5",
                selectInput('contrast', 'Select a comparison and tabulate the expression of DE Proteins',
                            c('Make your choice'='', 'S1c vs S0c'=1, 'S2c vs S1c'=2, 'S3c vs S2c'=3,
                              'S4c vs S3c'=4, 'S5c vs S4c'=5, 'S3c vs Sc1'=7)),
                numericInput('threshold', 'Select a pValue  cutoff for DE proteins', value=0.05, min=0.05, max=0.2, width='250px'),
                checkboxInput('showheatmap', 'Show Heatmap?', value=FALSE, width=NULL)),
            conditionalPanel(
                condition = "input.thetabs == 6",
                selectInput('pathway', 'Select a pathway to display heatmap',
                            c('Choose pathway'='', 'FGF signaling'='afgf7', 'GSK3B signaling'='agsk3beta', 'TGFbeta signaling'='atgfbeta',
                              'Retinoid acid signaling'='arar', 'Baron, M et al. 2016. Cell'='abaronmetal', 'Pdx1'='pdx1_gm',
                              'Activin'='activin_gm', 'Hedgehog signaling'='ashh', 'BMP signaling'='abmp'))),
            conditionalPanel(
                condition = "input.thetabs == 7",
                textInput(inputId='clustof.id', 'Enter a UniProt Accession or Gene Symbol for Cluster Profile', value='B2R6T2')),
            ##submitButton('Submit'))
            ##verbatimTextOutput('info'))
            conditionalPanel(
                condition = "input.thetabs == 8",
                selectInput('selecdf', 'Select protein or gene set',
                            c('Choose a set'='', 'DE proteins'='mostlysig')))
        ),
        mainPanel(
            tabsetPanel(id='thetabs', selected=1, type='tabs',
                        tabPanel('Barplot', value=1, plotOutput('barplot', height='800px')),
                        ##tabPanel('Heatmap', value=2, plotOutput('heatmap', height='800px'), verbatimTextOutput('msg')),
                        tabPanel('Heatmap', value=2, fluidRow(column(3, offset=1, verbatimTextOutput('msg'))), plotOutput('heatmap')),
                        tabPanel('Cumulative Distribution', value=3, fluidRow(
                                                                         column(8, plotOutput('ecd1', width='800px', height='600px')),
                                                                         column(12, plotOutput('ecd2', width='800px', height='600px')),
                                                                         column(8, plotOutput('ecd3', width='800px', height='600px')),
                                                                         column(12, plotOutput('ecd4', width='800px', height='600px')))),
                        tabPanel('Parallel Plot', value=4, plotOutput('parallelplot', width='600px', height='800px')),
                        tabPanel('DE Proteins', value=5, fluidRow(
                                                             column(8, tableOutput('contrasttable')),
                                                             column(8, plotOutput('htmap', width='600px', height='600px')))),
                        tabPanel('Pathways', value=6, fluidRow(
                                                          column(8, plotOutput('pwymap1', width='400px', height='400px')),
                                                          column(12, plotOutput('pwymap2', width='400px', height='400px')),
                                                          column(8, plotOutput('pwymap3', width='400px', height='400px')),
                                                          column(12, plotOutput('pwymap4', width='400px', height='400px')))),
                        tabPanel('Cluster Profile', value=7, fluidRow(
                                                                 column(5, plotOutput('allhmapP', height='600px')),
                                                                 column(5, plotOutput('allhmapR', height='600px')),
                                                                 column(5, plotOutput('selclustP', height='600px')),
                                                                 column(5, plotOutput('selclustR', height='600px')))),
                        tabPanel('ECDF Selection', value=8, fluidRow(
                                                                column(8, plotOutput('ecdfspl1', width='600px', height='600px')),#, hover=hoverOpts(id='plothover'))),
                                                                #column(4, verbatimTextOutput("hoverinfo")),
                                                                column(8, plotOutput('ecdfspl2', width='600px', height='600px')),
                                                                column(8, plotOutput('ecdfspl3', width='600px', height='600px')),
                                                                column(8, plotOutput('ecdfspl4', width='600px', height='600px')),
                                                                column(8, plotOutput('ecdfspl5', width='600px', height='600px'))))
                        )
        )
    )
)



## Server
server <- function(input, output) {
    prot.id <- reactive({
        if (input$geneprot.id != '') {
            ## D031517: added toupper()
            gpval <- toupper(input$geneprot.id)
            ids <- unlist(strsplit(unlist(strsplit(gpval, split=' ')), split=','))
            x <- grep('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', ids)
            if (any(x != 0)) {
                prot.id <- ids
            } else {
                prot.id <- unlist(mget(ids, ealias2acc, ifnotfound=unlist(mget(ids, eBBsym2acc, ifnotfound=ids))))
            }
        }
    })

    symb <- reactive({
        symb <- unlist(mget(prot.id(), eacc2sym, ifnotfound=unlist(mget(prot.id(), eBBacc2sym, ifnotfound=prot.id()))))
    })

    # for Cluster Profile
    selclustprot.id <- reactive({
        if (input$clustof.id != '') {
            ## D031517: added toupper()
            cval = toupper(input$clustof.id)
            ids <- unlist(strsplit(unlist(strsplit(cval, split=' ')), split=','))
            x <- grep('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', ids)
            if (any(x != 0)) {
                selclustprot.id <- ids[1] # only one protein allowed
            } else {
                selclustprot.id <- unlist(mget(ids, ealias2acc, ifnotfound=unlist(mget(ids, eBBsym2acc, ifnotfound=ids))))[1]
            }
        }
    })

    clsymb <- reactive({
        clsymb <- unlist(mget(selclustprot.id(), eacc2sym, ifnotfound=unlist(mget(selclustprot.id(), eBBacc2sym, ifnotfound=selclustprot.id()))))
    })

    output$barplot <- renderPlot({
        output$info <- renderText({'Enter one or more UniProt Accessions'})

        if (length(prot.id()) != 0) {
            output$info <- renderText({symb()})
            
            z <- presentNSamples(prot.id(), acclst)
            df.lst <- genSubsetAbundMat(z, msl, input$log)
            ## here we have to scale abundances of 227 sample D21717
            #xmn <- mean(df.lst[['s238']][,1])
            #xdf <- df.lst[['s227']]
            #xdf <- apply(xdf, 2, function(x) {
            #    f <- mean(x)
            #    r <- xmn/f
            #    x * r
            #})
            #df.lst[['s227']] <- as.data.frame(xdf)
            p <- multiBarPlot(df.lst, msl)
            
            if (is.null(p)) {
                 htmlOutput({'Proteins not found!'}) # TO CORRECT: this probably doesn't do anything
            } else {
                 do.call(grid.arrange, c(p, ncol=1))
                
            }
        }
    })

    output$heatmap <- renderPlot({
        if (length(prot.id()) >= 2 ) { # cannot cluster with just one protein, no distance!!
            output$info <- renderText({symb()})

            m <- exprs(mss)
            m <- m[rownames(m) %in% prot.id(),]

            if (dim(m)[1] == 0) {
                output$msg <- renderText({'Protein not found!'})
            } else {
                sym <- unlist(mget(rownames(m), eacc2sym, ifnotfound=unlist(mget(rownames(m), eBBacc2sym, ifnotfound=rownames(m)))))
                sym <- ifelse(is.na(sym), '-', sym)
                #sym <- replReplacmnt(sym)
                df <- as.data.frame(m)
                rownames(df) <- sym
                useLevelplot_v1(df)
            }
        }
    })

    pcdf <- reactive({
        if (length(prot.id()) != 0) {
            output$info <- renderText({symb()})
        
            z <- presentNSamples(prot.id(), acclst)
            df.lst <- genSubsetAbundMat_v0(z, msl)
            pcdf <- multiEcdfPlot(df.lst, msl)
        }
    })
    
    output$ecd1 <- renderPlot({
         pcdf()[['s227']]
    })
    output$ecd2 <- renderPlot({
        pcdf()[['s238']]
    })
    output$ecd3 <- renderPlot({
        pcdf()[['s239']]
    })
    output$ecd4 <- renderPlot({
        pcdf()[['s243']]
    })


    output$contrasttable <- renderTable({
        if (!input$showheatmap) {
            if (input$contrast != '') {
                x <- topTable(fit2, coef=as.integer(input$contrast), number=Inf, p.value=input$threshold)
                addAnnotTt(x)
            }
        }
    })

    output$htmap <- renderPlot({
        if (input$showheatmap) {
            if (input$contrast != '') {
                x <- topTable(fit2, coef=as.integer(input$contrast), number=Inf, p.value=input$threshold)
                if (length(rownames(x)) >= 2 ) { # cannot cluster with just one protein, no distance!!
                    # modified: next line replaced, after next added, D22117
                    m <- fit$coef
                    colnames(m) <- c('S0c', 'S1c', 'S2c', 'S3c', 'S4c', 'S5c')
                    m <- m[rownames(m) %in% rownames(x),]
                
                    if (dim(m)[1] == 0) {
                        output$msg <- renderText({'Protein not found!'})
                    } else {
                        sym <- unlist(mget(rownames(m), eacc2sym, ifnotfound=unlist(mget(rownames(m), eBBacc2sym, ifnotfound=rownames(m)))))
                        sym <- ifelse(is.na(sym), '-', sym)
                        #sym <- replReplacmnt(sym)
                        df <- as.data.frame(m)
                        rownames(df) <- sym
                        useLevelplot_v1(df)
                    }
                } else {
                    output$msg <- renderText({'One protein is not enough'})
                }
            }
        }
    })



    
    ## modified: heatmap with each samples separately D22117

    pwydf <- reactive({
        if (input$pathway != '') {
            x <- get(input$pathway) # this is a character string

            dflst <- list()
            for (nm in names(msl)) {
                if (length(x) >= 2 ) { # cannot cluster with just one protein, no distance!!
                    m <- exprs(msl[[nm]])
                    m <- m[rownames(m) %in% x,]
                    
                    if (dim(m)[1] == 0) {
                        output$msg <- renderText({'Protein not found!'})
                    } else {
                        sym <- unlist(mget(rownames(m), eacc2sym, ifnotfound=unlist(mget(rownames(m), eBBacc2sym, ifnotfound=rownames(m)))))
                        sym <- ifelse(is.na(sym), '-', sym)
                        ##sym <- replReplacmnt(sym)
                        rownames(m) <- sym
                        df <- as.data.frame(m)
                        dflst[[nm]] <- df
                    }
                }
            }
            return(dflst)
        }
    })
    
    output$pwymap1 <- renderPlot({
        if (input$pathway != '') useLevelplot_v1(pwydf()[['s227']])
    })
    output$pwymap2 <- renderPlot({
        if (input$pathway != '') useLevelplot_v1(pwydf()[['s238']])
    })
    output$pwymap3 <- renderPlot({
        if (input$pathway != '') useLevelplot_v1(pwydf()[['s239']])
    })
    output$pwymap4 <- renderPlot({
        if (input$pathway != '') useLevelplot_v1(pwydf()[['s243']])
    })

    clobjlst.P <- reactive({
        f.m <- fit$coef
        f.df <- as.data.frame(f.m)
        colnames(f.df) <- c('S0c', 'S1c', 'S2c', 'S3c', 'S4c', 'S5c')

        sf.df <- scale(f.df)
        tsf.df <- t(sf.df)

        row.hc <- hclust(stats::dist(tsf.df))
        col.hc <- hclust(stats::dist(sf.df))

        return(list(fdf=f.df, sf=sf.df, tsf=tsf.df, rhc=row.hc, chc=col.hc))
    })

    clobjlst.R <- reactive({
        vfit <- vfit[vfit$genes$Acc %in% rownames(exprs(mss)),]
        rownames(vfit$coef) <- vfit$genes$Acc
        vf.m <- vfit$coef
        vf.df <- as.data.frame(vf.m)
        #colnames(vf.df) <- c('S0c', 'S1c', 'S2c', 'S3c', 'S4c', 'S5c')
        
        svf.df <- scale(vf.df)
        tsvf.df <- t(svf.df)
        
        row.hc <- hclust(stats::dist(tsvf.df))
        col.hc <- hclust(stats::dist(svf.df))
        
        return(list(fdf=vf.df, sf=svf.df, tsf=tsvf.df, rhc=row.hc, chc=col.hc))
    })

    
    output$allhmapP <- renderPlot({
        dd.col <- as.dendrogram(clobjlst.P()[['chc']])
        
        useLevelplot_v2(clobjlst.P()[['tsf']], dd.col, main='Proteins')
    })

    output$allhmapR <- renderPlot({
        dd.col <- as.dendrogram(clobjlst.R()[['chc']])
        
        useLevelplot_v2(clobjlst.R()[['tsf']], dd.col, main='Transcripts')
    })

    subclustdf <- reactive({
        df <- clustHeatSubcluster_v1(clobjlst.P()[['fdf']], clobjlst.P()[['chc']], nclust=50, prot=selclustprot.id())
        return(df)
    })
    output$selclustP <- renderPlot({
        ##clustHeatSubcluster(f.df, col.hc, nclust=50, prot='B2R6T2')
        ##clustHeatSubcluster(clobjlst.P()[['fdf']], clobjlst.P()[['chc']], nclust=50, prot=selclustprot.id())
        xdf <- subclustdf()
        sym <- unlist(mget(rownames(xdf), eacc2sym, ifnotfound=unlist(mget(rownames(xdf), eBBacc2sym, ifnotfound=rownames(xdf)))))
        sym <- sapply(sym, function(x) unlist(strsplit(x, split=';'))[1])
        rownames(xdf) <- sym
        useLevelplot_v1(xdf, scale=TRUE, pdf=FALSE, main='Proteins')
    })

    output$selclustR <- renderPlot({
        ##clustHeatSubcluster(f.df, col.hc, nclust=50, prot='B2R6T2')
        ##clustHeatSubcluster(clobjlst.R()[['fdf']], clobjlst.R()[['chc']], nclust=50, prot=selclustprot.id())
        xdf <- subclustdf()
        ydf <- clobjlst.R()[['fdf']]
        ydf <- ydf[rownames(ydf) %in% rownames(xdf),]

        sym <- unlist(mget(rownames(ydf), eacc2sym, ifnotfound=unlist(mget(rownames(ydf), eBBacc2sym, ifnotfound=rownames(ydf)))))
        sym <- sapply(sym, function(x) unlist(strsplit(x, split=';'))[1])
        rownames(ydf) <- sym

        useLevelplot_v1(ydf, scale=TRUE, pdf=FALSE, main='Transcripts')
    })

    newfitdf <- reactive({
        df <- as.data.frame(newfit$coef)
        colnames(df) <- c('S0c', 'S1c', 'S2c', 'S3c', 'S4c', 'S5c')
        return(df)
    })
    
    ## ECDF Selection
    output$ecdfspl1 <- renderPlot({
        plotECDF_v5(newfitdf(), 'S1c', rownames(c1), markers)
    })
    output$ecdfspl2 <- renderPlot({
        plotECDF_v5(newfitdf(), 'S2c', rownames(c2), markers)
    })
    output$ecdfspl3 <- renderPlot({
        plotECDF_v5(newfitdf(), 'S3c', rownames(c3), markers)
    })
    output$ecdfspl4 <- renderPlot({
        plotECDF_v5(newfitdf(), 'S4c', rownames(c4), markers)
    })
    output$ecdfspl5 <- renderPlot({
        plotECDF_v5(newfitdf(), 'S5c', rownames(c5), markers)
    })
    
    # for testing
    #output$contrasttable <- renderText({
    #    paste('Type is: ', class(get(input$contrast)))
    #})
    
    output$parallelplot <- renderPlot({
        if (length(prot.id()) != 0) {
            output$info <- renderText({symb()})

            z <- presentNSamples(prot.id(), acclst)
            df.lst <- genSubsetAbundMat(z, msl, input$log)
            p <- multiParallelPlot_v2(df.lst, msl, input$log)
            
            if (is.null(p)) {
                 output$info <- renderText({'Protein not  found!'})
                
            } else {
                 do.call(grid.arrange, c(p, ncol=1))
            }
        }
    })
}

shinyApp(ui=ui, server=server)
    
