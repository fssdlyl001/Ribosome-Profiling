library(shiny)
library(ggplot2)
library(quantmod)
library(grid)
options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx8000m")

shinyServer(function(input, output) {
  Data = reactive({
    input$goButton
    if (input$goButton == 0) {return(NULL)}
    isolate({
      inFile = input$file1
      if (is.null(inFile))
        return(NULL)
      
      col=input$colour2
      lineweight=input$lineweight
      group=c()
      graph1=0
      top1=1
      if (length(inFile$datapath)==1) {
        polysome = read.csv(inFile$datapath, header=T, na.strings=c("NA","NaN", ""," "), skip=32, fileEncoding="latin1")
        polysome = na.omit(polysome[,1:2])
        colnames(polysome)=c("Position","Absorbance")
        if (length(which(polysome[,1]<0.03))>1) {
          firstcol=c(which(polysome[,1]<0.03))
          lastcol=c(firstcol[-1]-1,nrow(polysome))
          for (k in 1:length(firstcol)) {
            group=c(group,rep(k,(lastcol[k]-firstcol[k]+1)))
            if (max(polysome$Absorbance[-c(firstcol[k]:(firstcol[k]+217))])>1) {top=round(max(polysome$Absorbance[-c(firstcol[k]:(firstcol[k]+217))])+0.1,1)} else {top=1}
            top1=c(top1,top)
          }
        } else {
          group=1
          groupname=1
          if (is.null(input$n1)) {} else {
            groupname=eval(parse(text="input$n1"))
          }
          if (any(groupname=="")) {groupname[which(groupname=="")]=which(groupname=="")}
          if (max(polysome$Absorbance[-1:-218])>1) {top=round(max(polysome$Absorbance[-1:-218])+0.1,1);scale=0.1} else {top=1;scale=0.25}
          graph1=ggplot(polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=col),size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="Group")+ggtitle("Ribosome Profiling")
          graph1=graph1+scale_x_continuous(breaks=c(0,seq(0,56.5,5.65)))+scale_y_continuous(breaks=seq(0,top,scale))+coord_cartesian(xlim=c(0,56.5),ylim=c(0,top))+scale_colour_manual(values=col,labels=groupname)+theme(plot.title = element_text(face="bold",size=24),axis.title.x = element_text(face="bold", size=20),axis.title.y = element_text(face="bold", size=19),legend.title=element_text(face="bold",size=20),legend.text=element_text(face="bold",size=18),legend.key=element_rect(fill="white"),legend.key.height=unit(2,"cm"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
        }
      }
      if (length(inFile$datapath)>1) {
        polysome = data.frame()
        filenumber=length(inFile$datapath)
        number=1:filenumber
        a=0
        for (i in 1:filenumber) {
          datum = read.csv(inFile$datapath[i], header=T, na.strings=c("NA","NaN", ""," "), skip=32, fileEncoding="latin1")
          datum = na.omit(datum[,1:2])
          polysome=rbind(polysome,datum)
          num=number[i+a]
          if (length(which(datum[,1]<0.03))>1) {
            a=length(which(datum[,1]<0.03))-1
            number=1:(length(number)+a)
            firstcol=c(which(datum[,1]<0.03))
            lastcol=c(firstcol[-1]-1,nrow(datum))
            for (j in 1:length(firstcol)) {
              group=c(group,rep(num,lastcol[j]-firstcol[j]+1))
              num=num+1
              if (max(datum[-c(firstcol[j]:(firstcol[j]+217)),2])>1) {top=round(max(datum[-c(firstcol[j]:(firstcol[j]+217)),2])+0.1,1)} else {top=1}
              top1=c(top1,top)
            }
          }
          if (length(which(datum[,1]<0.03))==1) {
            group=c(group,rep(num,nrow(datum)))
            if (max(datum[-1:-218,2])>1) {top=round(max(datum[-1:-218,2])+0.1,1)} else {top=1}
            top1=c(top1,top)
          }
        }
      }
      normalize=c()
      if (length(group)>1) {
        polysome=data.frame(polysome,group)
        colnames(polysome)=c("Position","Absorbance","group")
        no=unique(group)
        normalize=polysome
        groupname=no
        if (is.null(input$n1)) {} else {
          groupname=c()
          for (i in unique(group)) {
            a=eval(parse(text=paste("input$n",i,sep="")))
            groupname=c(groupname,a)
          }
          if (any(groupname=="")) {groupname[which(groupname=="")]=which(groupname=="")}
        }
        Ab1=c()
        for (i in no) {
          Ab=subset(polysome,group==no[i],Absorbance)
          Ab=as.matrix(Ab)
          valley=findValleys(Ab)[1]
          if (valley>=90&valley<200) {
            range=Ab[90:1006]
            np=findPeaks(range)+88
          }
          if (valley<90|valley>=200) {
            range=Ab[218:1006]
            np=findPeaks(range)+216
          }
          if (Ab[np[2]]>Ab[np[1]]) {firstnp=np[2]} else {firstnp=np[1]}
          if (firstnp>340) Ab=append(Ab[-c(1:c(firstnp-340))],rep(min(Ab),firstnp-340))
          if (firstnp<340) Ab=append(rep(max(Ab),340-firstnp),Ab[-c(c(length(Ab)-339+firstnp):(length(Ab)+1))])
          Ab=Ab-rep(min(Ab),length(Ab))
          Ab1=c(Ab1,Ab)
        }
        normalize$Absorbance=Ab1
        top=max(top1)
        if (top > 1) {scale=0.1} else {scale=0.25}
      }
      polysome=as.data.frame(polysome)
      info=list(polysome=polysome,graph1=graph1,group=group,col=col,normalize=normalize,top=top,scale=scale,groupname=groupname,lineweight=lineweight)
      return(info)
    })
  })
  
  
  Value = reactive({
    input$goButton
    if (input$goButton == 0) {return(NULL)}
    polysome=Data()$polysome
    graph1=Data()$graph1
    group=Data()$group
    top=Data()$top
    scale=Data()$scale
    normalize=Data()$normalize
    groupname=Data()$groupname
    lineweight=input$lineweight
    isolate({
      graph2=0
      if (length(group)>1) {
        if (input$colour1=="classical") {
          pal=scale_colour_hue(labels=groupname,l=60,c=120)
        } else {
          if (input$similar1=="cold") {pal=scale_colour_hue(labels=groupname,h=c(150, 270),c=50)}
          if (input$similar1=="warm") {pal=scale_colour_hue(labels=groupname,h=c(0, 80))}
          if (input$similar1=="grey") {pal=scale_colour_grey(labels=groupname,start = 0, end = 0.8)}
        }
        
        graph1=ggplot(data=polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(group=factor(polysome$group),aes(colour=factor(polysome$group)),size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="Group")+ggtitle("Ribosome Profiling")
        graph1=graph1+scale_x_continuous(breaks=c(0,seq(0,56.5,5.65)))+scale_y_continuous(breaks=seq(0,top,scale))+coord_cartesian(xlim=c(0,56.5),ylim=c(0,top))+pal+theme(plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=12),axis.title.x = element_text(face="bold", size=20),axis.title.y = element_text(face="bold", size=19),legend.title=element_text(face="bold",size=14),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
        
        graph2=ggplot(data=normalize,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(group=factor(polysome$group),aes(colour=factor(polysome$group)),size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="Group")+ggtitle("Ribosome Profiling (Normalized)")
        graph2=graph2+scale_x_continuous(breaks=c(0,seq(0,56.5,5.65)))+coord_cartesian(xlim=c(0,56.5),ylim=c(0,1))+pal+theme(plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=12),axis.title.x = element_text(face="bold", size=20),axis.title.y = element_text(face="bold", size=19),legend.title=element_text(face="bold",size=14),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
      }
      info2=list(polysome=polysome,graph1=graph1,graph2=graph2,group=group,normalize=normalize)
      return(info2)
    })
  })
  
  output$insertname = renderUI({
    group=Value()$group
    if (input$check == TRUE) {
      text_list = lapply(unique(group), function(i) {
        textname = paste("n",i,sep="")
        textInput(textname, paste("Sample",i), "")
      })
      do.call(tagList, text_list)
    } else {return(NULL)}
  })
  
  
  output$plots = renderUI({
    group=Value()$group
    plot_output_list = lapply(unique(group), function(i) {
      plotname = paste("plot", i, sep="")
      plotOutput(plotname)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  observe({
    polysome=Value()$polysome
    group=Value()$group
    col=Data()$col
    groupname=Data()$groupname
    top=Data()$top
    scale=Data()$scale
    lineweight=input$lineweight
    if (length(group)>1) {
      for (i in unique(group)) {
        local({
          my_i = i
          plotname = paste("plot", my_i, sep="")
          singleribosome=subset(polysome,group==my_i,select=-group)
          singlename=paste("Ribosome Profiling","-",groupname[my_i],sep="")
          if (max(singleribosome$Absorbance[-1:-218])>1) {top2=round(max(singleribosome$Absorbance[-1:-218])+0.1,1);scale2=0.1} else {top2=1;scale2=0.25}
          
          output[[plotname]] = renderPlot({
            p1=ggplot(singleribosome,aes(x=Position,y=Absorbance))+geom_line(colour=col,size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)")+ggtitle(singlename)
            p1=p1+scale_x_continuous(breaks=c(0,seq(0,56.5,5.65)))+scale_y_continuous(breaks=seq(0,top2,scale2))+coord_cartesian(xlim=c(0,56.5),ylim=c(0,top2))+theme(plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", size=20),axis.title.y = element_text(face="bold", size=19),legend.title=element_text(face="bold",size=14),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),plot.margin = unit(c(1,1,1,1), "cm"))
            print(p1)
          },width=750,height=400)
          groupname=unique(group)
        })
      }
    } else {
      local({
        plotname="plot1"
        singlename=paste("Ribosome Profiling","-",groupname[1],sep="")
        output[[plotname]] = renderPlot({
          p2=ggplot(polysome,aes(x=Position,y=Absorbance))+geom_line(colour=col,size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)")+ggtitle(singlename)
          p2=p2+scale_x_continuous(breaks=c(0,seq(0,56.5,5.65)))+scale_y_continuous(breaks=seq(0,top,scale))+coord_cartesian(xlim=c(0,56.5),ylim=c(0,top))+theme(plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", size=20),axis.title.y = element_text(face="bold", size=19),legend.title=element_text(face="bold",size=14),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),plot.margin = unit(c(1,1,1,1), "cm"))
          print(p2)
        },width=750,height=400)
      })
    }
  })
  
  
  output$graph1 = renderPlot({
    if (input$goButton == 0) {return(NULL)} 
    if (is.null(input$file1))
      return(NULL)
    graph1 = Value()$graph1
    print(graph1)
  },width=800,height=400)
  
  output$graph2 = renderPlot({
    if (input$goButton == 0) {return(NULL)} 
    if (is.null(input$file1))
      return(NULL)
    group = Value()$group
    graph2 = Value()$graph2
    if (length(group)>1) {print(graph2)} else {return(NULL)}
  },width=800,height=400)
  
  
  output$table = renderTable({
    if (input$goButton == 0) {return(NULL)} 
    if (is.null(input$file1))
      return(NULL)
    table1 = Value()$polysome
    table1
  })
  
  output$information = renderTable({
    if (is.null(input$file1))
      return(NULL)
    information=read.csv(input$file1$datapath[1], header=F, na.strings=c("NA","NaN", ""," "), nrow=28, fileEncoding="latin1")
    infoname=as.character(information[1,1])
    information=data.frame(information[-1,])
    names(information)=infoname
    if (length(input$file1$datapath)>1) {
    for (i in 2:length(input$file1$datapath)) {
      inform=read.csv(input$file1$datapath[i], header=F, na.strings=c("NA","NaN", ""," "), nrow=28, fileEncoding="latin1")
      infoname=c(infoname,as.character(inform[1,1]))
      inform=data.frame(inform[-1,])
      information=data.frame(information, inform)
    }
    colnames(information)=infoname}
    return(information)
  })
})
