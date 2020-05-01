library("shiny")
library("foreign")
library("ggplot2")

shinyServer(function(input, output) {
    
    ### Argument names:
    dataset<-reactive({
        if (is.null(input$datafile))
            return(NULL)                
        BioID<-read.csv(input$datafile$datapath)
        
        BioID=data.frame(BioID)
        BioID=BioID[,c('Bait','Prey','AvgSpec')]
        if (input$LC == "Whole cell"){
          Reference=read.csv('reference_data.csv')
        } else{
          Reference=read.csv('reference_data_mito.csv')
        }
        Reference=data.frame(Reference)
        Reference=Reference[,c('Bait','Prey','AvgSpec','Localization')]
        
        
        BaitU=unique(BioID$Bait)
        loc_reference=unique(Reference$Localization)
        loc_score_matrix=matrix(0,nrow=length(loc_reference),ncol=length(BaitU))
        rownames(loc_score_matrix)=loc_reference
        colnames(loc_score_matrix)=BaitU
        
        for (j in BaitU ){
            Bait_data=BioID[BioID$Bait==j,]
            
            loc_normalization=rep(1,length(loc_reference))
            names(loc_normalization)=loc_reference
            
            for (i in loc_reference){
                prey_set=Reference[Reference$Localization!=i,2]
                Reference_ff=Reference[Reference$Localization==i ,]
                Reference_l=Reference_ff[!is.element(Reference_ff$Prey,prey_set),]
                loc_normalization[i]=sum(Bait_data[is.element(Bait_data$Prey,Reference_l$Prey) ,3])
            }
            
            
            for (i in loc_reference){
                prey_set=Reference[Reference$Localization!=i,2]
                Reference_ff=Reference[Reference$Localization==i ,]
                Reference_l=Reference_ff[!is.element(Reference_ff$Prey,prey_set),]
                
                overlap_ref_spc=sum(Reference_l[is.element(Reference_l$Prey, Bait_data$Prey) ,3])+1
                spc_ref=sum(Reference_l[,3])
                #loc_score_matrix[i,j]=100*(overlap_ref_spc/spc_ref)*loc_normalization[i]/(sum(loc_normalization)-loc_normalization[i])
                loc_score_matrix[i,j]=100*(overlap_ref_spc/spc_ref)*loc_normalization[i]/sum(loc_normalization)
            }
            
        }
        loc_score_matrix_normalized=data.frame(loc_score_matrix)
        
        for (j in 1:dim(loc_score_matrix)[2])
            for (i in 1:dim(loc_score_matrix)[1])
                loc_score_matrix_normalized[i,j]=(loc_score_matrix[i,j]-min(loc_score_matrix[,j]))/(max(loc_score_matrix[,j])-min(loc_score_matrix[,j]))
        loc_score_matrix_normalized=cbind(Localization=rownames(loc_score_matrix_normalized),loc_score_matrix_normalized)
        # loc_score_matrix_normalized
    })
    # Argument selector:
    output$ArgSelect <- renderUI({
        #if (length(ArgNames())==0) return(NULL)
        
        selectInput("Bait","Bait:",colnames(dataset())[-1])
    })
    
    # Show table:
    Bait_r <- reactive({
        #if (is.null(input$file)) {
        # User has not uploaded a file yet
        return(input$Bait)
    })
    
    #BioID=dataset()
    # BioID=data.frame(BioID)
    # BioID=BioID[,c('Bait','Prey','AvgSpec')]
    output$plot <- renderPlot({
        dd=dataset()
        Bait_rr=Bait_r()
        ff=dd[,c('Localization',as.character(Bait_rr))]
        ggplot(ff, aes(x=ff[,1],
                       y=ff[,2], fill=ff[,2]))+
            geom_bar( stat='identity') +
            scale_fill_gradient2(low = "lightblue", mid = "white", high = "darkblue")+
            coord_polar()+ggtitle(Bait_rr)+
            theme(plot.title = element_text(size=20),
                  axis.title=element_blank(),axis.text.y =element_blank(), axis.ticks.y = element_blank(),
                  axis.text.x = element_text (vjust = 1, hjust=1,face=2, size=15),
                  panel.border = element_rect(colour="grey", fill=NA, size=1),
                  legend.title = element_blank(),
                  panel.background=element_blank(),
                  panel.grid.major.x   = element_line(colour = "grey",size=0.1))
        # ggsave("plot.pdf",plotInput())
    })
    plotInput = function() {
        dd=dataset()
        Bait_rr=Bait_r()
        ff=dd[,c('Localization',as.character(Bait_rr))]
        ggplot(ff, aes(x=ff[,1],
                       y=ff[,2], fill=ff[,2]))+
            geom_bar( stat='identity') +
            scale_fill_gradient2(low = "lightblue", mid = "white", high = "darkblue")+
            coord_polar()+ggtitle(Bait_rr)+
            theme(plot.title = element_text(size=20),
                  axis.title=element_blank(),axis.text.y =element_blank(), axis.ticks.y = element_blank(),
                  axis.text.x = element_text (vjust = 1, hjust=1,face=2, size=15),
                  panel.border = element_rect(colour="grey", fill=NA, size=1),
                  legend.title = element_blank(),
                  panel.background=element_blank(),
                  panel.grid.major.x   = element_line(colour = "grey",size=0.1))
    }
    
    output$table <- renderTable({
        
        
        
        return(dataset())
    })
    
    
    ### Download dump:
    # 
    # output$downloadFigure <- downloadHandler(
    #   
    #     #Bait_rr=Bait_r()
    #     filename = 'Figure.pdf',
    #     content = function(file) {
    #         device <- function(..., width, height) {
    #             grDevices::pdf(..., width = 8, height = 10)
    #         }
    #         ggsave(file, plot = plotInput(), device = device)
    #     }
    # )
    output$downloadFigure <- downloadHandler(
      #Bait_rr=Bait_r()
      filename = 'Figure.pdf',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::pdf(..., width = 8, height = 10)
        }
        ggsave(file, plot = plotInput(), device = "pdf")
      })
    
    output$downloadExample1 <- downloadHandler(
      filename = 'example1.csv',
      content = function(file) {
        file.copy("example1.csv", file)
      }
    )
    
    output$downloadOutput <- downloadHandler(
      filename = 'pubData_figures.pdf',
      content = function(file) {
        file.copy("pubData_figures.pdf", file)
      }
    )
    
    output$downloadExample2 <- downloadHandler(
      filename = 'example2.csv',
      content = function(file) {
        file.copy("example2.csv", file)
        }
    )
})