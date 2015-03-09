library(shiny)
library(ggplot2)
library(grid)
library(quantmod)
library(flux)
library(shinysky)
library(shinyBS)
library(snowfall)
options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx8000m")
#Be careful: group refers to sample, and groupname represents sample names.#

shinyServer(function(input, session, output) {
  sfInit(parallel=TRUE , cpus=2)
  Data = reactive({
    inFile = input$file1
    if (is.null(inFile))
      return(NULL)
    
    filenumber=length(inFile$datapath)
    group=c()
    polylist=as.list(filenumber)
    if (filenumber==1) {
      judge =  read.csv(inFile$datapath, header=F, na.strings=c("NA","NaN", ""," "), nrow=2, fileEncoding="latin1")
      if (grepl("Version",as.character(judge[2,1]))) {
        polysome = read.csv(inFile$datapath, header=T, na.strings=c("NA","NaN", ""," "), skip=32, fileEncoding="latin1")
      } else {
        if (grepl("Position",judge[1,1])) {polysome = read.csv(inFile$datapath, header=T, na.strings=c("NA","NaN", ""," "), fileEncoding="latin1")} else {return(NULL)}
      }
      xmax=round(max(polysome[,1])+0.05,1)
      polysome = na.omit(polysome[,1:2])
      colnames(polysome)=c("Position","Absorbance")
      if (length(which(polysome[,1]<0.03))>1) {
        firstcol=c(which(polysome[,1]<0.03))
        lastcol=c(firstcol[-1]-1,nrow(polysome))
        polylist=c(polylist,list(length(firstcol)-1))
        for (k in 1:length(firstcol)) {
          group=c(group,rep(k,(lastcol[k]-firstcol[k]+1)))
          polylist[[k]]=polysome[firstcol[k]:lastcol[k],]
        }
      } else {
        group=1
        polylist[[1]]=polysome
      }
    }
    
    if (filenumber>1) {
      polysome = data.frame()
      number=1:filenumber
      a=0
      xmax=0
      for (i in 1:filenumber) {
        judge =  read.csv(inFile$datapath[i], header=F, na.strings=c("NA","NaN", ""," "), nrow=2, fileEncoding="latin1")
        if (grepl("Version",as.character(judge[2,1]))) {skipno=32} else {skipno=0}
        
        datum = read.csv(inFile$datapath[i], header=T, na.strings=c("NA","NaN", ""," "), skip=skipno, fileEncoding="latin1")
        tempmax=max(datum[,1])
        if (xmax<tempmax) {xmax=tempmax}
        datum = na.omit(datum[,1:2])
        polysome=rbind(polysome,datum)
        num=number[i+a]
        if (length(which(datum[,1]<0.03))>1) {
          a=length(which(datum[,1]<0.03))-1
          number=1:(length(number)+a)
          firstcol=c(which(datum[,1]<0.03))
          lastcol=c(firstcol[-1]-1,nrow(datum))
          polylist=c(polylist,list(length(firstcol)-1))
          for (j in 1:length(firstcol)) {
            group=c(group,rep(num,lastcol[j]-firstcol[j]+1))
            polylist[[i+j-1]]=datum[firstcol[j]:lastcol[j],]
            num=num+1
          }
        }
        if (length(which(datum[,1]<0.03))==1) {
          polylist[[i]]=datum
          group=c(group,rep(num,nrow(datum)))
          if (max(datum[-1:-218,2])>1) {top=round(max(datum[-1:-218,2])+0.1,1)} else {top=1}
        }
      }
    }
    
    if (length(group)>1) {
      if (nrow(polysome) != length(group)) {return(NULL)}
      polysome=data.frame(polysome,group)
      colnames(polysome)=c("Position","Absorbance","group")
    }
    groupname=unique(group)
    polysome=as.data.frame(polysome)
    if (max(polysome$Absorbance[-1:-218])>1) {top=round(max(polysome$Absorbance[-1:-218])+0.1,1);scale=0.1} else {top=1;scale=0.25}
    xscale=round(xmax/10,2)
    
    info=list(polysome=polysome,polylist=polylist,group=group,top=top,scale=scale,groupname=groupname,xmax=xmax,xscale=xscale)
    return(info)
  })
  
  
  Normalize = reactive({
    if (is.null(Data()$groupname)) {return(NULL)} else{
      if (is.null(input$autonormalize)) {return(NULL)} else {
        if (input$autonormalize == 0) {return(NULL)}
      }
    }
    isolate({
    polysome=Data()$polysome
    polylist=Data()$polylist
    group=Data()$group
    groupname=Data()$groupname
    normalize=polysome
    
    if (length(group)==1) {
      mixlist=sfLapply(polylist,autonorm)
      
      normlist=sfLapply(mixlist,'[[',1)
      normalize$Absorbance=as.vector(sfSapply(mixlist,'[[',2))
      Calclist=sfLapply(mixlist,'[[',3)
    }
    
    if (length(group)>1) {
      mixlist=sfLapply(polylist,autonorm)

      normlist=sfLapply(mixlist,'[[',1)
      normalize$Absorbance=as.vector(unlist(sfSapply(mixlist,'[[',2)))
      Calclist=sfLapply(mixlist,'[[',3)
    }
    
    normalinfo=list(polysome=polysome,normlist=normlist,group=group,normalize=normalize,Calclist=Calclist,groupname=groupname)
    return(normalinfo)
    })
  })
  
  
  autonorm = function(polydata) {
    Ab=as.matrix(polydata[,2])
    polyrange=as.matrix(polydata[,1])
    
    basalmin=min(Ab[(length(Ab)*0.7):length(Ab)])
    fvalley=1
    Valleys=findValleys(Ab)-1
    if (length(Valleys)==0) {Valleys=1}
    
    if (Valleys[1]<=200) {
      low=which(Valleys<=200)
      if (any(Ab[Valleys[low]]<=basalmin)) {
        abnormal1=Valleys[which(Ab[1:200][Valleys]<=basalmin)]
        rmvalley1=which(Valleys<max(abnormal1)+30)
        fvalley=Valleys[last(rmvalley1)]
        Valleys=Valleys[-rmvalley1]
      }
    }
    
    if (max(Ab)-min(Ab)>0.4) {
      if (any(Ab[Valleys]>=max(Ab)-0.15)) {
        rmvalley2=which(Ab[Valleys]>=max(Ab)-0.15)
        Valleys=Valleys[-rmvalley2]
      }
    } else {
      if (any(Ab[Valleys]>=max(Ab)-0.05)) {
        rmvalley2=which(Ab[Valleys]>=max(Ab)-0.05)
        Valleys=Valleys[-rmvalley2]
      }
    }
    
    firstvalley=Valleys[1]
    range=Ab[firstvalley:length(Ab)]
    np=findPeaks(range)+firstvalley-1
    
    fprange=fvalley:firstvalley
    if (length(findPeaks(Ab[fprange]))==0) {fpeak=which.max(Ab[fprange]+fvalley-1)} else {fpeak=findPeaks(Ab[fprange])+fvalley-2}
    normpeak=last(fprange[which(Ab[fprange] >= min(Ab[fpeak])-0.05)])
    
    if (length(np)<3) {np=c(340,478,569,654,755,850,which(Ab==max(Ab[850:length(Ab)])))}
    
    if (length(which(Ab[np]>=Ab[firstvalley]+0.05))>4) {
      if (Ab[np[1]]<Ab[firstvalley]+0.05) {
        np=np[-1]
        firstvalley=Valleys[2]
      }
    }
    
    if (any(Ab[np]<basalmin)) {
      np=np[-which(Ab[np]<basalmin)]
    }
    
    rm=c()
    j=1
    while (j<=length(np)-2) {
      if (np[j+1]-np[j] <= 100) {
        rep=c(which(np>=np[j] & np<np[j]+100))
        retain=round(mean(rep))
        if (any(Ab[np[rep]]>=Ab[np[retain]]+0.01)){retain=round(mean(rep[which(Ab[np[rep]]>=Ab[np[retain]]+0.01)]))}
        rm=c(rm,rep[rep != retain])
        j=j+length(rep)
      } else {j=j+1}
    }
    if (length(rm)>0 & length(rm)<length(np)-4) {
      np=np[-rm]
    }
    
    
    if (any(Ab[np]<=basalmin)) {np[-which(Ab[np]<=basalmin)]}
    if (Ab[np[2]]>Ab[np[1]]+Ab[np[2]]/2) {
      np=np[-1]
      if (Ab[np[2]]>Ab[np[1]]+Ab[np[2]]/2){
        np=np[-1]
        if (Ab[np[2]]>Ab[np[1]]+Ab[np[2]]/2){
          np=np[-1]
        }
      }
    } else {
      if (Ab[np[3]]>Ab[np[1]]+Ab[np[3]]/2) {
        np=np[-1:-2]
      }
    }
    firstnp=np[1]
    
    stpeak=round(length(Ab)/5)
    if (firstnp>stpeak) {Ab=append(Ab[-c(1:c(firstnp-stpeak))],rep(min(Ab[firstnp:length(Ab)]),firstnp-stpeak))}
    if (firstnp<stpeak) {Ab=append(rep(max(Ab),stpeak-firstnp),Ab[-c(c(length(Ab)-stpeak+1+firstnp):(length(Ab)+1))])}
    np=np-firstnp+stpeak
    firstvalley=firstvalley-firstnp+stpeak
    normpeak=normpeak-firstnp+stpeak
    Ab=Ab-rep(basalmin,length(Ab))
    np3=np[3]
    
    normdata=data.frame(polyrange,Ab)
    colnames(normdata)=c("Position","Absorbance")
    
    return(list(normdata,Ab,list(Ab,normpeak,firstvalley,np3)))
  }
  
  
  
  CalcAUC = reactive({
    if (is.null(Normplot()$graph2)) {return(NULL)} else {
      if (is.null(input$showauc)) {return(NULL)} else{
        if (input$showauc==0) {return(NULL)}
      }
    }
    
    isolate({
      group=Normalize()$group
      Calclist=Normalize()$Calclist
      auclist=sfLapply(Calclist,calc.auc)
      
      if (length(group)==1) {
        MAUC=sfSapply(auclist,'[[',1)
        PAUC=sfSapply(auclist,'[[',2)
        Fm=sfSapply(auclist,'[[',3)
      }
      
      if (length(group)>1) {
        MAUC=sfSapply(auclist,'[[',1)
        PAUC=sfSapply(auclist,'[[',2)
        Fm=as.vector(unlist(sfSapply(auclist,'[[',3)))
      }
      
      aucinfo=list(MAUC=MAUC,PAUC=PAUC,Fm=Fm)
      return(aucinfo)
    })
  })
  
  calc.auc = function(auclist) {
    Ab=auclist[[1]]
    normpeak=auclist[[2]]
    firstvalley=auclist[[3]]
    np3=auclist[[4]]
    if (normpeak<=0) {normpeak=1}
    normmax=Ab[normpeak]
    fmname=paste("fmodel","i",sep="")
    if (normmax > Ab[firstvalley]) {
      x=(2*normpeak-firstvalley):firstvalley
      y=c(rev(Ab[(normpeak+1):firstvalley]),Ab[normpeak:firstvalley])
      sigma=1/(normmax*sqrt(2*pi))
      mu=normpeak
      normdata=data.frame(x,y)
      fm=nls(y~c * exp(-(x-a)^2/(2*b^2)) / (b*sqrt(2*pi)), data=normdata, control=nls.control(maxiter = 300), start = list(a=mu,b=sigma,c=2))
      deduct1=predict(fm,data.frame(x=firstvalley:np3))
      deduct2=predict(fm,data.frame(x=np3:length(Ab)))
      
      Fm=predict(fm,data.frame(x=1:length(Ab)))
      
      monorange=firstvalley:np3
      MAUC=auc(monorange,Ab[monorange])-auc(monorange,deduct1)
      polyrange=np3:length(Ab)
      PAUC=auc(polyrange,Ab[polyrange])-auc(polyrange,deduct2)
      
    } else {
      MAUC=auc(monorange,Ab[monorange])
      PAUC=auc(polyrange,Ab[polyrange])
      
      Fm=rep(0,length(Ab))
    }
    
    return(list(MAUC,PAUC,Fm))
  }
  
  
  
  Groupname = reactive({
    groupname=Data()$groupname
    if (is.null(groupname))
      return(NULL)
    group=Data()$group
    if (input$check == TRUE) {
      if (!is.null(input$n1)) {
        if (length(group)==1) {
          if (!is.null(input$n1)) {
            groupname=eval(parse(text="input$n1"))
          }
          if (groupname=="") {groupname=1}
        } else {
          groupname=c()
          groupname=sfLapply(unique(group), function(i) {
            eval(parse(text=paste("input$n",i,sep="")))
          })
          if (any(groupname=="")) {groupname[which(groupname=="")]=which(groupname=="")}
        }
      }
    }
    
    mygrouname=list(groupname=groupname)
    return(mygrouname)
  })
  
  
  Judge = reactive({
    if (is.null(input$file1)) {return(NULL)} else {
      if (is.null(Basicplot()$filepath)) {return(NULL)} else {
        if (!identical(Basicplot()$filepath,input$file1$datapath)) {return(NULL)} else {
          judgefile=T
        }
      }
    }
    return(list(judgefile=judgefile))
  })
  
  
  
  
  Basicplot = reactive({
    if (input$goButton == 0) 
      return(NULL)
    isolate({
    filepath=input$file1$datapath
    polysome=Data()$polysome
    polyjudge=polysome
    xmax=Data()$xmax
    xscale=Data()$xscale
    group=Data()$group
    top=Data()$top
    scale=Data()$scale    
    col=input$singlecol
    
    colorname=Temp()$colorname
    realgroup=Temp()$realgroup
    realgroup2=Temp()$realgroup2
    check2=input$check2
    groupingno=input$groupingno
    groupnumber=input$groupnumber
    groupnumber1=input$groupnumber1
    groupnumber2=input$groupnumber2
    linechoose=Temp()$linechoose
    colchoose=Temp()$colchoose
    colour1=input$colour1
    similar1=input$similar1
    if (is.null(similar1)) {colour1="classical"}
    singlegrname=NULL
    grname=Temp()$grname
    grname2=Temp()$grname2
    uniquegr=Temp()$uniquegr
    
    
    if (input$mode=='quick') {check2=FALSE} else {
    if (is.null(groupingno)) {
      check2=FALSE
    } else {
      if (groupingno=="single") {
        if (is.null(input$m1)) {check2=FALSE} else{
          if (input$m1=="") {check2=FALSE}
        }
      }
      if (groupingno=="double") {
        if (is.null(input$m11)|is.null(input$m22)) {check2=FALSE} else {
          if (input$m11==""|input$m22=="") {check2=FALSE}
        }
      }
      if (is.null(colchoose)) {check2=FALSE} else{
        if (colchoose[1]==0) {check2=FALSE}
      }
    }
    }
    
    if (input$mode=='quick'|input$check==FALSE) {groupname=Data()$groupname} else {
      if (!is.null(input$n1)) {
        groupname=Groupname()$groupname 
      } else {groupname=Data()$groupname}
    }
    
    graph1=NULL
    pal=NULL
    linestyle=NULL
    filltype=NULL
    legendguide=NULL
    beautify=NULL
    
    
    if (length(group) == 1) {
      group1=rep(1,nrow(polysome))
      polysome=data.frame(polysome,group1)
      beautify=list(labs(x="Position (mm)",y="Absorbance (254nm)",colour="Group"),scale_x_continuous(breaks=c(0,seq(0,xmax,xscale))),scale_y_continuous(breaks=seq(0,top,scale)),scale_colour_manual(values=col,labels=groupname),coord_cartesian(xlim=c(0,xmax),ylim=c(0,top)),theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black")))
      graph1=ggplot(data=polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(group1),size=group1))+ggtitle("Ribosome Profiling")
      graph1=graph1+beautify
      
    }
    
    if (check2 == FALSE) {
      if (length(group)>1) {
        if (colour1=="classical") {
          if (length(groupname)<=9){pal=scale_colour_brewer(labels=groupname,palette="Set1")} else{
            pal=scale_colour_hue(labels=groupname,l=55,c=120)
          }
        } else {
          getPalette1 = colorRampPalette(brewer.pal(9, "PuBu")[4:9])
          getPalette2 = colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
          getPalette3 = colorRampPalette(brewer.pal(9, "Greys")[4:9])
          if (similar1=="cold") {pal=scale_colour_manual(labels=groupname,values=getPalette1(length(groupname)))}
          if (similar1=="warm") {pal=scale_colour_manual(labels=groupname,values=getPalette2(length(groupname)))}
          if (similar1=="grey") {pal=scale_colour_manual(labels=groupname,values=getPalette3(length(groupname)))}
        }
        
        if (length(groupname)<=17) {legendno=1} else {legendno=2}
        legendguide=guides(colour = guide_legend(ncol=legendno))
        
        beautify=list(labs(x="Position (mm)",y="Absorbance (254nm)",colour="Group"),legendguide,scale_x_continuous(breaks=c(0,seq(0,xmax,xscale))),scale_y_continuous(breaks=seq(0,top,scale)),coord_cartesian(xlim=c(0,xmax),ylim=c(0,top)),pal,theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black")))
        graph1=ggplot(data=polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(group=factor(group),aes(colour=factor(group),size=group))+ggtitle("Ribosome Profiling")
        graph1=graph1+beautify
        
      }
      } else {
          polysome=Temp()$polysome
          uniquegr1=Temp()$uniquegr1
          uniquegr2=Temp()$uniquegr2
          
          if (!is.null(groupingno)) {
          if (groupingno=="single") {
            if (!is.null(groupnumber)) {
            if (!is.null(input$m1)) {
            if (!is.na(input$m1)& !input$m1=="") {
            if (length(group)>1) {
              if (length(groupname)<=12) {legendno=1} else {legendno=2}
              legendguide = guides(colour = guide_legend(ncol=legendno, order = 2), fill = guide_legend(order = 1))
              
              if (colour1=="classical") {
                pal=scale_colour_manual(values=colchoose, labels=groupname)
                filltype=scale_fill_manual(values=colorname, labels=unique(realgroup))
              }
              if (colour1=="similar") {
                getPalette1 = colorRampPalette(brewer.pal(9, "PuBu")[4:9])
                getPalette2 = colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
                getPalette3 = colorRampPalette(brewer.pal(9, "Greys")[4:9])
                if (similar1=="cold") {pal=scale_colour_manual(labels=groupname,values=getPalette1(length(groupname)))}
                if (similar1=="warm") {pal=scale_colour_manual(labels=groupname,values=getPalette2(length(groupname)))}
                if (similar1=="grey") {pal=scale_colour_manual(labels=groupname,values=getPalette3(length(groupname)))}
                filltype=scale_fill_manual(values=colorname, labels=unique(realgroup), guide='none')
                legendguide = guides(colour = guide_legend(ncol=legendno))
              }
              
              singlegrname=input$groupingname
              if (is.null(singlegrname)) {singlegrname="Group"} else {
                if (singlegrname=="") {singlegrname="Group"}
              }
              
              beautify=list(labs(x="Position (mm)",y="Absorbance (254nm)",colour="Sample",fill=singlegrname),legendguide,scale_x_continuous(breaks=c(0,seq(0,xmax,xscale))),scale_y_continuous(breaks=seq(0,top,scale)),coord_cartesian(xlim=c(0,xmax),ylim=c(0,top)),pal,filltype,theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=12),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black")))
              graph1=ggplot(data=polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(Sample),group=factor(Sample),size=Sample))+geom_rect(aes(xmin=-1,xmax=-1,ymin=-1,ymax=-1,fill=factor(Group)))+ggtitle("Ribosome Profiling")
              graph1=graph1+beautify
              
              if (!is.null(uniquegr)) {
              colnames(polysome)=c("Position","Absorbance","Sample",singlegrname)
              }
            }
           }
           }
           }
          } else {
            if (!is.null(groupnumber1) & !is.null(groupnumber2)) {
            if (!is.null(input$m11) & !is.null(input$m21)) { 
            if(!is.na(input$m11) & !input$m11=="" & !is.na(input$m21) & !input$m21==""){
            if (length(group)>1) {
              if (length(groupname)<=6) {legendno=1} else {legendno=2}
              legendguide = guides(colour = guide_legend(ncol=legendno, order = 3), linetype = guide_legend(order = 2), fill = guide_legend(order = 1))
              
              if (colour1=="classical") {
                pal=scale_colour_manual(values=colchoose, labels=groupname)
                filltype=scale_fill_manual(values=colorname, labels=unique(realgroup))
                linestyle=scale_linetype_manual(values=linechoose)
              } 
              if (colour1=="similar") {
                getPalette1 = colorRampPalette(brewer.pal(9, "PuBu")[4:9])
                getPalette2 = colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
                getPalette3 = colorRampPalette(brewer.pal(9, "Greys")[4:9])
                if (similar1=="cold") {pal=scale_colour_manual(labels=groupname,values=getPalette1(length(groupname)))}
                if (similar1=="warm") {pal=scale_colour_manual(labels=groupname,values=getPalette2(length(groupname)))}
                if (similar1=="grey") {pal=scale_colour_manual(labels=groupname,values=getPalette3(length(groupname)))}
                filltype=scale_fill_manual(values=colorname, guide = 'none')
                linestyle=scale_linetype_manual(values=linechoose)
                
                legendguide = guides(colour = guide_legend(ncol=legendno, order = 2), linetype = guide_legend(order = 1))
              }
              
              beautify=list(labs(x="Position (mm)",y="Absorbance (254nm)",colour="Sample",fill=grname,linetype=grname2),legendguide,scale_x_continuous(breaks=c(0,seq(0,xmax,xscale))),scale_y_continuous(breaks=seq(0,top,scale)),coord_cartesian(xlim=c(0,xmax),ylim=c(0,top)),pal,filltype,linestyle,theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=14),legend.text=element_text(face="bold",size=12),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black")))
              graph1=ggplot(data=polysome,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(Sample),group=factor(Sample),linetype=factor(Group2,levels=unique(uniquegr2)),size=Sample))+geom_rect(aes(xmin=-1,xmax=-1,ymin=-1,ymax=-1,fill=factor(Group1)))+ggtitle("Ribosome Profiling")
              graph1=graph1+beautify
              
              if (!is.null(uniquegr1) & !is.null(uniquegr2)) {
              colnames(polysome)=c("Position","Absorbance","Sample",grname,grname2)
              }
            }
          }
          }
          }
        }
       }
     }
    plotinfo1=list(polysome=polysome,filepath=filepath,polyjudge=polyjudge,graph1=graph1,group=group,groupname=groupname,check2=check2,singlegrname=singlegrname,realgroup=realgroup,realgroup2=realgroup2,groupingno=groupingno,beautify=beautify)
    return(plotinfo1)
    })
  })
  
  
  Normplot = reactive({
    if (is.null(Judge()$judgefile))
      return(NULL)
    if (is.null(Data()$groupname)) {return(NULL)} else{
      if (is.null(input$autonormalize)) {return(NULL)} else {
        if (input$autonormalize == 0) {return(NULL)}
      }
    }
    isolate({
      normalize=Normalize()$normalize
      if (is.null(normalize)) 
        return(NULL)
      graph2=NULL
      group=Data()$group
      beautify=Basicplot()$beautify
      groupingno=Basicplot()$groupingno
      check2=Basicplot()$check2
      singlegrname=Basicplot()$singlegrname
      grname=Temp()$grname
      grname2=Temp()$grname2
      
      if (length(group)==1) {
        group2=rep(1,nrow(normalize))
        normalize=data.frame(normalize,group2)
        graph2=ggplot(data=normalize,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(group2),size=group2))+ggtitle("Ribosome Profiling")
        graph2=graph2+beautify
      }
      
      if (check2 == FALSE) {
        if (length(group)>1) { 
          graph2=ggplot(data=normalize,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(group=factor(group),aes(colour=factor(group),size=group))+ggtitle("Ribosome Profiling (Normalized)")
          graph2=graph2+beautify
        }
      } else {
        uniquegr=Temp()$uniquegr
        uniquegr1=Temp()$uniquegr1
        uniquegr2=Temp()$uniquegr2
        realgroup=Basicplot()$realgroup
        
        if (!is.null(groupingno)) {
          if (groupingno=="single") {
            normalize=data.frame(normalize,realgroup)   
            colnames(normalize)=c("Position","Absorbance","Sample","Group")
            graph2=ggplot(data=normalize,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(Sample),group=factor(Sample),size=Sample))+geom_rect(aes(xmin=-1,xmax=-1,ymin=-1,ymax=-1,fill=factor(Group)))+ggtitle("Ribosome Profiling (Normalized)")
            graph2=graph2+beautify
          } else {
            realgroup2=Basicplot()$realgroup2
            normalize=data.frame(normalize,realgroup,realgroup2)
            colnames(normalize)=c("Position","Absorbance","Sample","Group1","Group2")
            graph2=ggplot(data=normalize,environment=environment(),aes(x=Position,y=Absorbance))+geom_line(aes(colour=factor(Sample),group=factor(Sample),linetype=factor(Group2,levels=unique(uniquegr2)),size=Sample))+geom_rect(aes(xmin=-1,xmax=-1,ymin=-1,ymax=-1,fill=factor(Group1)))+ggtitle("Ribosome Profiling (Normalized)")
            graph2=graph2+beautify
          }
       }
      }
      
      plotinfo2=list(graph2=graph2,normalize=normalize)
      return(plotinfo2)
    })
  })
  
  
  AUCplot = reactive({
    if (is.null(Normplot()$graph2)) {return(NULL)} else {
      if (is.null(input$showauc)) {return(NULL)} else{
        if (input$showauc==0) {return(NULL)}
      }
    }
    isolate({
      if (is.null(Judge()$judgefile))
        return(NULL)
      MAUC=CalcAUC()$MAUC
      PAUC=CalcAUC()$PAUC
      if (is.null(MAUC)|is.null(PAUC)) 
        return(NULL)
      group=Data()$group
      groupingno=input$groupingno
      groupnumber=input$groupnumber
      groupnumber1=input$groupnumber1
      groupnumber2=input$groupnumber2
      check2=Basicplot()$check2
      groupname=Basicplot()$groupname
      singlegrname=Basicplot()$singlegrname
      graph3=NULL
      graph4=NULL
      graph5=NULL
      graph6=NULL
      colorname=Temp()$colorname
      grname=Temp()$grname
      grname2=Temp()$grname2
      
      rcolorname=c("#0066CC","#CC3333","#669933","#993399","golden rod","#CCCCCC")
      rwidth=0.45
      
      if(check2==TRUE) {
        if (length(colorname)==2) {rcolor=c("#0033FF","#FF9900")} else {rcolor=rcolorname[1:length(colorname)]}
      }
      
      Ratio1=PAUC/MAUC
      Ratio2=MAUC/PAUC
        
      if (any(MAUC<=0)){
        na1=which(MAUC<=0)
        MAUC[na1]=0
        Ratio1[na1]=PAUC[na1]
      }
      if (any(PAUC<=0)){
        na2=which(PAUC<=0)
        PAUC[na2]=0
        Ratio2[na2]=MAUC[na2]
      }
      rgroupname=1:length(groupname)
      if (length(rgroupname)==length(Ratio1)) {
        RAUC1=data.frame(rgroupname,groupname,Ratio1)
        RAUC2=data.frame(rgroupname,groupname,Ratio2)
      }
      
      if (length(group) == 1) {
        if (!is.null(Ratio1) & !is.null(Ratio2)) {
          rnum=factor(c(0,1),labels=c("Monosome/Polysome","Polysome/Monosome"))
          Ratio=c(Ratio2,Ratio1)
          singleauc=data.frame(rnum,Ratio)
          
          rawnum=factor(c(0,1),labels=c("Monosome","Polysome"))
          AUC=c(MAUC,PAUC)
          rawauc=data.frame(rawnum,AUC)
          
          rmax=round(max(singleauc$Ratio))
          if (rmax < 4) {rscale=round(rmax+0.51,0); rdigit=0.5}
          if (rmax < 10 & rmax >= 4) {rscale=round(rmax+0.51,0); rdigit=1}
          if (rmax <= 50 & rmax >= 10) {rscale=round(rmax+5,-1); rdigit=5}
          if (rmax < 100 & rmax >= 50) {rscale=round(rmax+5,-1); rdigit=10}
          if (rmax < 1000 & rmax >= 100) {rscale=round(rmax+50,-2); rdigit=100}
          if (rmax < 10000 & rmax >= 1000) {rscale=round(rmax+500,-3); rdigit=1000}
          if (rmax < 100000 & rmax >= 10000) {rscale=round(rmax+5000,-3); rdigit=10000}
          if (rmax >= 100000) {rscale=round(rmax+rmax/5,-log10(rmax)+1); rdigit=round(rmax/10,-log10(rmax)+2)}
          
          rmax1=round(max(AUC))
          if (rmax1 < 4) {rscale1=round(rmax1+0.51,0); rdigit1=0.5}
          if (rmax1 < 10 & rmax >= 4) {rscale1=round(rmax1+0.51,0); rdigit1=1}
          if (rmax1 <= 50 & rmax1 >= 10) {rscale1=round(rmax1+5,-1); rdigit1=5}
          if (rmax1 < 100 & rmax1 >= 50) {rscale1=round(rmax1+5,-1); rdigit1=10}
          if (rmax1 < 1000 & rmax1 >= 100) {rscale1=round(rmax1+50,-2); rdigit1=100}
          if (rmax1 < 10000 & rmax1 >= 1000) {rscale1=round(rmax1+500,-3); rdigit1=1000}
          if (rmax1 < 100000 & rmax1 >= 10000) {rscale1=round(rmax1+5000,-3); rdigit1=10000}
          if (rmax1 >= 100000) {rscale1=round(rmax1+rmax1/5,-log10(rmax1)+1); rdigit1=round(rmax1/10,-log10(rmax1)+2)}
          
          
          graph3=ggplot(data=singleauc,aes(x=rnum,y=Ratio,fill=rnum))+geom_bar(stat="identity",width=0.3,show_guide=FALSE)+labs(x="Ratio Type",y="AUC Ratio")+ggtitle("Area Under Curve Ratio")
          graph3=graph3+scale_fill_manual(values=c("#0033FF","#FF9900"))+scale_y_continuous(breaks=seq(0,rscale,rdigit),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale))+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          
          graph4=ggplot(data=rawauc,aes(x=rawnum,y=AUC,fill=rawnum))+geom_bar(stat="identity",width=0.3,show_guide=FALSE)+labs(x="AUC Type",y="AUC")+ggtitle("Area Under Curve")
          graph4=graph4+scale_fill_manual(values=c("#0033FF","#FF9900"))+scale_y_continuous(breaks=seq(0,rscale1,rdigit1),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale1))+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          
        }
      }
      
      if (check2 == FALSE) {
        if (length(group)>1) {
          rawauc1=data.frame(rgroupname,groupname,MAUC)
          rawauc2=data.frame(rgroupname,groupname,PAUC)
          
          rmax1=round(max(RAUC1$Ratio1))
          if (rmax1 < 4) {rscale1=round(rmax1+0.51,0); rdigit1=0.5}
          if (rmax1 < 10 & rmax1 >= 4) {rscale1=round(rmax1+0.51,0); rdigit1=1}
          if (rmax1 <= 50 & rmax1 >= 10) {rscale1=round(rmax1+5,-1); rdigit1=5}
          if (rmax1 < 100 & rmax1 >= 50) {rscale1=round(rmax1+5,-1); rdigit1=10}
          if (rmax1 < 1000 & rmax1 >= 100) {rscale1=round(rmax1+50,-2); rdigit1=100}
          if (rmax1 < 10000 & rmax1 >= 1000) {rscale1=round(rmax1+500,-3); rdigit1=1000}
          if (rmax1 < 100000 & rmax1 >= 10000) {rscale1=round(rmax1+5000,-3); rdigit1=10000}
          if (rmax1 >= 100000) {rscale1=round(rmax1+rmax1/5,-log10(rmax1)+1); rdigit1=round(rmax1/10,-log10(rmax1)+2)}
          
          rmax2=round(max(RAUC2$Ratio2))
          if (rmax2 < 4) {rscale2=round(rmax2+0.51,0); rdigit2=0.5}
          if (rmax2 < 10 & rmax2 >= 4) {rscale2=round(rmax2+0.51,0); rdigit2=1}
          if (rmax2 <= 50 & rmax2 >= 10) {rscale2=round(rmax2+5,-1); rdigit2=5}
          if (rmax2 < 100 & rmax2 >= 50) {rscale2=round(rmax2+5,-1); rdigit2=10}
          if (rmax2 < 1000 & rmax2 >= 100) {rscale2=round(rmax2+50,-2); rdigit2=100}
          if (rmax2 < 10000 & rmax2 >= 1000) {rscale2=round(rmax2+500,-3); rdigit2=1000}
          if (rmax2 < 100000 & rmax2 >= 10000) {rscale2=round(rmax2+5000,-3); rdigit2=10000}
          if (rmax2 >= 100000) {rscale2=round(rmax2+rmax2/5,-log10(rmax2)+1); rdigit2=round(rmax2/10,-log10(rmax2)+2)}
          
          rmax3=round(max(MAUC))
          if (rmax3 < 4) {rscale3=round(rmax3+0.51,0); rdigit3=0.5}
          if (rmax3 < 10 & rmax3 >= 4) {rscale3=round(rmax3+0.51,0); rdigit3=1}
          if (rmax3 <= 50 & rmax3 >= 10) {rscale3=round(rmax3+5,-1); rdigit3=5}
          if (rmax3 < 100 & rmax3 >= 50) {rscale3=round(rmax3+5,-1); rdigit3=10}
          if (rmax3 < 1000 & rmax3 >= 100) {rscale3=round(rmax3+50,-2); rdigit3=100}
          if (rmax3 < 10000 & rmax3 >= 1000) {rscale3=round(rmax3+500,-3); rdigit3=1000}
          if (rmax3 < 100000 & rmax3 >= 10000) {rscale3=round(rmax3+5000,-3); rdigit3=10000}
          if (rmax3 >= 100000) {rscale3=round(rmax3+rmax3/5,-log10(rmax3)+1); rdigit3=round(rmax3/10,-log10(rmax3)+2)}
          
          rmax4=round(max(PAUC))
          if (rmax4 < 4) {rscale4=round(rmax4+0.51,0); rdigit4=0.5}
          if (rmax4 < 10 & rmax4 >= 4) {rscale4=round(rmax4+0.51,0); rdigit4=1}
          if (rmax4 <= 50 & rmax4 >= 10) {rscale4=round(rmax4+5,-1); rdigit4=5}
          if (rmax4 < 100 & rmax4 >= 50) {rscale4=round(rmax4+5,-1); rdigit4=10}
          if (rmax4 < 1000 & rmax4 >= 100) {rscale4=round(rmax4+50,-2); rdigit4=100}
          if (rmax4 < 10000 & rmax4 >= 1000) {rscale4=round(rmax4+500,-3); rdigit4=1000}
          if (rmax4 < 100000 & rmax4 >= 10000) {rscale4=round(rmax4+5000,-3); rdigit4=10000}
          if (rmax4 >= 100000) {rscale4=round(rmax4+rmax4/5,-log10(rmax4)+1); rdigit4=round(rmax4/10,-log10(rmax4)+2)}
          
          if (length(groupname)<=6) {
            rcolor=rcolorname[1:length(groupname)]
            rpal=scale_fill_manual(values=rcolor)
            if (length(groupname) == 2) {rwidth = 0.45}
            if (length(groupname) == 3) {rwidth = 0.55}
          } else {
            getPalette = colorRampPalette(brewer.pal(9, "Set1"))
            rpal=scale_fill_manual(values=getPalette(length(groupname)))
          }
          
          if (length(groupname) > 3) {rwidth = 0.75}
          
          graph3=ggplot(data=RAUC1,aes(x=factor(rgroupname),y=Ratio1,fill=factor(groupname)))+geom_bar(stat="identity",width=rwidth,weight=0.2)+labs(x="Group",y="AUC Ratio",fill="Group")+ggtitle("AUC Ratio (Polysome/Monosome)")
          graph3=graph3+scale_y_continuous(breaks=seq(0,rscale1,rdigit1),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale1))+rpal+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          
          graph4=ggplot(data=RAUC2,aes(x=factor(rgroupname),y=Ratio2,fill=factor(groupname)))+geom_bar(stat="identity",width=rwidth,weight=0.2)+labs(x="Group",y="AUC Ratio",fill="Group")+ggtitle("AUC Ratio (Monosome/Polysome)")
          graph4=graph4+scale_y_continuous(breaks=seq(0,rscale2,rdigit2),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale2))+rpal+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          
          graph5=ggplot(data=rawauc1,aes(x=factor(rgroupname),y=MAUC,fill=factor(groupname)))+geom_bar(stat="identity",width=rwidth,weight=0.2)+labs(x="Group",y="AUC",fill="Group")+ggtitle("AUC of Monosome")
          graph5=graph5+scale_y_continuous(breaks=seq(0,rscale3,rdigit3),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale3))+rpal+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          
          graph6=ggplot(data=rawauc2,aes(x=factor(rgroupname),y=PAUC,fill=factor(groupname)))+geom_bar(stat="identity",width=rwidth,weight=0.2)+labs(x="Group",y="AUC",fill="Group")+ggtitle("AUC of Polysome")
          graph6=graph6+scale_y_continuous(breaks=seq(0,rscale4,rdigit4),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale4))+rpal+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
        }
      } else {
        uniquegr=Temp()$uniquegr
        uniquegr1=Temp()$uniquegr1
        uniquegr2=Temp()$uniquegr2
        
        if (!is.null(groupingno)) {
          if (groupingno=="single") {
            if (!is.null(groupnumber)) {
              if (!is.null(uniquegr)) {
                RAUC1=data.frame(RAUC1,uniquegr)
                RAUC2=data.frame(RAUC2,uniquegr)
                colnames(RAUC1)=c("rgroupname","groupname","Ratio","uniquegr")
                colnames(RAUC2)=c("rgroupname","groupname","Ratio","uniquegr")
                AUC=rbind(RAUC1,RAUC2)
                ratiotype=rep(0:1,each=length(uniquegr))
                AUC=data.frame(AUC,ratiotype)
                rmean=with(AUC,tapply(Ratio,data.frame(ratiotype,uniquegr),mean))
                rsd=with(AUC,tapply(Ratio,data.frame(ratiotype,uniquegr),sd))
                if (length(colnames(rmean))==length(unique(uniquegr)) & all(colnames(rmean) %in% unique(uniquegr))) {} else {rmean=t(rmean)}
                if (length(colnames(rsd))==length(unique(uniquegr)) & all(colnames(rsd) %in% unique(uniquegr))) {} else {rsd=t(rsd)}
                order=colnames(rmean)
                rmean=c(rmean)
                rsd=c(rsd)
                rnum=rep(factor(c(0,1),labels=c("Polysome/Monosome","Monosome/Polysome")),length(order))
                rgroup=rep(order,each=2)
                auc=data.frame(rmean,rsd,rnum,rgroup)
                
                if (any(is.na(rmean))) {rmean[is.na(rmean)]=0}
                if (any(is.na(rsd))) {rsd[is.na(rsd)]=0}
                rmax=round(max(rmean+rsd))
                if (rmax < 4) {rscale=round(rmax+0.51,0); rdigit=0.5}
                if (rmax < 10 & rmax >= 4) {rscale=round(rmax+0.51,0); rdigit=1}
                if (rmax <= 50 & rmax >= 10) {rscale=round(rmax+5,-1); rdigit=5}
                if (rmax < 100 & rmax >= 50) {rscale=round(rmax+5,-1); rdigit=10}
                if (rmax < 1000 & rmax >= 100) {rscale=round(rmax+50,-2); rdigit=100}
                if (rmax < 10000 & rmax >= 1000) {rscale=round(rmax+500,-3); rdigit=1000}
                if (rmax < 100000 & rmax >= 10000) {rscale=round(rmax+5000,-3); rdigit=10000}
                if (rmax >= 100000) {rscale=round(rmax+rmax/5,-log10(rmax)+1); rdigit=round(rmax/10,-log10(rmax)+2)}
                
                if (length(rcolor) == 2) {rwidth = 0.45}
                if (length(rcolor) == 3) {rwidth = 0.55}
                if (length(rcolor) > 3) {rwidth = 0.75}
                
                graph3=ggplot(data=auc,aes(x=rnum,y=rmean,fill=factor(rgroup)))+geom_bar(stat="identity",position="dodge",width=rwidth,weight=0.2)+geom_errorbar(aes(ymin=rmean-rsd,ymax=rmean+rsd),width=0.1,position=position_dodge(rwidth))+labs(x="Ratio Type",y="AUC Ratio",fill=singlegrname)+ggtitle("Area Under Curve Ratio")
                graph3=graph3+scale_fill_manual(values=rcolor)+scale_y_continuous(breaks=seq(0,rscale,rdigit),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale))+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
              }
            }
          } else {
                if (!is.null(uniquegr1) & !is.null(uniquegr2)) {
                  RAUC1=data.frame(RAUC1,uniquegr1,uniquegr2)
                  RAUC2=data.frame(RAUC2,uniquegr1,uniquegr2)
                  rmean1=with(RAUC1,tapply(Ratio1,data.frame(uniquegr1,uniquegr2),mean))
                  rsd1=with(RAUC1,tapply(Ratio1,data.frame(uniquegr1,uniquegr2),sd))
                  rmean2=with(RAUC2,tapply(Ratio2,data.frame(uniquegr1,uniquegr2),mean))
                  rsd2=with(RAUC2,tapply(Ratio2,data.frame(uniquegr1,uniquegr2),sd))
                  
                  if (length(colnames(rmean1))==length(unique(uniquegr1)) & all(colnames(rmean1) %in% unique(uniquegr1))) {} else {rmean1=t(rmean1)}
                  if (length(colnames(rsd1))==length(unique(uniquegr1)) & all(colnames(rsd1) %in% unique(uniquegr1))) {} else {rsd1=t(rsd1)}
                  order1=colnames(rmean1)
                  order2=colnames(t(rmean1))
                  rnum1=rep(order1,each=length(order2))
                  rgroup1=rep(order2,length(order1))
                  rmean1=c(rmean1)
                  rsd1=c(rsd1)
                  
                  if (length(colnames(rmean2))==length(unique(uniquegr1)) & all(colnames(rmean2) %in% unique(uniquegr1))) {} else {rmean2=t(rmean2)}
                  if (length(colnames(rsd2))==length(unique(uniquegr1)) & all(colnames(rsd2) %in% unique(uniquegr1))) {} else {rsd2=t(rsd2)}
                  order3=colnames(rmean2)
                  order4=colnames(t(rmean2))
                  rnum2=rep(order3,each=length(order4))
                  rgroup2=rep(order4,length(order3))
                  rmean2=c(rmean2)
                  rsd2=c(rsd2)
                  
                  auc1=data.frame(rmean1,rsd1,rnum1,rgroup1)
                  auc2=data.frame(rmean2,rsd2,rnum2,rgroup2)
                  
                  if (any(is.na(rmean1))) {rmean1[is.na(rmean1)]=0}
                  if (any(is.na(rsd1))) {rsd1[is.na(rsd1)]=0}
                  if (any(is.na(rmean2))) {rmean2[is.na(rmean2)]=0}
                  if (any(is.na(rsd2))) {rsd2[is.na(rsd2)]=0}
                  
                  rmax1=round(max(rmean1+rsd1))
                  if (rmax1 < 4) {rscale1=round(rmax1+0.51,0); rdigit1=0.5}
                  if (rmax1 < 10 & rmax1 >= 4) {rscale1=round(rmax1+0.51,0); rdigit1=1}
                  if (rmax1 < 100 & rmax1 >= 50) {rscale1=round(rmax1+5,-1); rdigit1=10}
                  if (rmax1 < 1000 & rmax1 >= 100) {rscale1=round(rmax1+50,-2); rdigit1=100}
                  if (rmax1 < 10000 & rmax1 >= 1000) {rscale1=round(rmax1+500,-3); rdigit1=1000}
                  if (rmax1 < 100000 & rmax1 >= 10000) {rscale1=round(rmax1+5000,-3); rdigit1=10000}
                  if (rmax1 >= 100000) {rscale1=round(rmax1+rmax1/5,-log10(rmax1)+1); rdigit1=round(rmax1/10,-log10(rmax1)+2)}
                                    
                  rmax2=round(max(rmean2+rsd2))
                  if (rmax2 < 4) {rscale2=round(rmax2+0.51,0); rdigit2=0.5}
                  if (rmax2 < 10 & rmax2 >= 4) {rscale2=round(rmax2+0.51,0); rdigit2=1}
                  if (rmax2 <= 50 & rmax2 >= 10) {rscale2=round(rmax2+5,-1); rdigit2=5}
                  if (rmax2 < 100 & rmax2 >= 50) {rscale2=round(rmax2+5,-1); rdigit2=10}
                  if (rmax2 < 1000 & rmax2 >= 100) {rscale2=round(rmax2+50,-2); rdigit2=100}
                  if (rmax2 < 10000 & rmax2 >= 1000) {rscale2=round(rmax2+500,-3); rdigit2=1000}
                  if (rmax2 < 100000 & rmax2 >= 10000) {rscale2=round(rmax2+5000,-3); rdigit2=10000}
                  if (rmax2 >= 100000) {rscale2=round(rmax2+rmax2/5,-log10(rmax2)+1); rdigit2=round(rmax2/10,-log10(rmax2)+2)}
                  
                  if (length(unique(uniquegr2))==2) {rcolor=c("#0033FF","#FF9900")} else {rcolor=rcolorname[1:length(unique(linechoose))]}
                  if (length(order1) == 2) {rwidth = 0.45}
                  if (length(order1) == 3) {rwidth = 0.55}
                  if (length(order1) > 3) {rwidth = 0.75}
                  
                  lean=0
                  hmove=0.5
                  vmove=0
                  if (length(which(nchar(unique(rnum1))>10))>1 | max(nchar(unique(rnum1)))>12) {lean=30;hmove=0.95;vmove=1}
                  
                  graph3=ggplot(data=auc1,aes(x=factor(rnum1),y=rmean1,fill=factor(rgroup1),group=factor(rgroup1)))+geom_bar(stat="identity",position="dodge",width=rwidth,weight=0.2)+geom_errorbar(aes(ymin=rmean1-rsd1,ymax=rmean1+rsd1),width=0.1,position=position_dodge(rwidth))+labs(x=grname,y="AUC Ratio",fill=grname2)+ggtitle("AUC Ratio (Polysome/Monosome)")
                  graph3=graph3+scale_fill_manual(values=rcolor)+scale_y_continuous(breaks=seq(0,rscale1,rdigit1),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale1))+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",angle=lean,hjust=hmove,vjust=vmove,size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
                  
                  graph4=ggplot(data=auc2,aes(x=factor(rnum2),y=rmean2,fill=factor(rgroup2),group=factor(rgroup2)))+geom_bar(stat="identity",position="dodge",width=rwidth,weight=0.2)+geom_errorbar(aes(ymin=rmean2-rsd2,ymax=rmean2+rsd2),width=0.1,position=position_dodge(rwidth))+labs(x=grname,y="AUC Ratio",fill=grname2)+ggtitle("AUC Ratio (Monosome/Polysome)")
                  graph4=graph4+scale_fill_manual(values=rcolor)+scale_y_continuous(breaks=seq(0,rscale2,rdigit2),expand = c(0, 0))+coord_cartesian(ylim=c(0,rscale2))+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=24),legend.title=element_text(face="bold",size=17),legend.text=element_text(face="bold",size=16),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=20),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",angle=lean,hjust=hmove,vjust=vmove,size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
                }
          }
        }
      }
      plotinfo3=list(graph3=graph3,graph4=graph4,graph5=graph5,graph6=graph6)
      return(plotinfo3)
    })
  })
  
  
  Temp = reactive({
      reds=c("red","tomato","firebrick","light coral","indian red","pink","salmon","hot pink","dark red")
      yellows=c("orange","yellow","dark orange","gold","golden rod","#CC6600","dark golden rod","orange red")
      greens=c("green","forest green","chart reuse","dark green","lime green","light green","#808000","olive drab")
      blues=c("blue","cyan","dark blue","royal blue","deep sky blue","aqua marine","dodger blue","sky blue")
      blacks=c("black","grey","dim grey","dark grey","dark slate gray","light grey","#C0C0C0","#330000")
      purples=c("purple","deep pink","magenta","medium violet red","violet","medium orchid","dark magenta","dark orchid")
      coltype=list(reds,blues,greens,yellows,blacks,purples)
      colstyle=c("red","blue","green","yellow","black","purple")
      
      linetype=c("solid","dashed","dotdash","longdash","dotted","twodash")
      
      polysome=Data()$polysome
      group=Data()$group
      groupname=Data()$groupname      
      
      if (input$goButton == 0) {return(NULL)}
      isolate({
        groupingno=input$groupingno
        grname1=input$groupingname1
        grname2=input$groupingname2
        if (is.null(grname1)) {grname1="Group 1"} else {if (grname1=="") {grname1="Group 1"} } 
        if (is.null(grname2)) {grname2="Group 2"} else {if (grname2=="") {grname2="Group 2"} }    
        
        uniquegr=c()
        uniquegr1=c()
        uniquegr2=c()
        
        if (length(group)==1)
          return(NULL)
        if (input$mode == "quick")
          return(NULL)
        if (input$check2 == FALSE) {return(NULL)} else {
          if (is.null(groupingno)) {return(NULL)} else {
          if (groupingno=="single") {
            groupnumber=input$groupnumber
            if (is.null(groupnumber)) {return(NULL)} else {
              if (is.null(input$m1)) {return(NULL)} else {
            realgroup=group
            uniquegr=sfSapply(unique(group),function(i) {
              gr=eval(parse(text=paste("input$g",i,sep="")))
            })
            realgroup=as.character(uniquegr)[match(realgroup,unique(realgroup))]
            realgroup2=NULL
            
            truename=unique(realgroup)
            groupnumber=length(truename)
            
            if (groupnumber>6) {colchoose=0} else {
              colchoose=unique(group)
              for (i in 1:groupnumber) {
                whichcoltype=coltype[[i]]
                r=1
                for (j in unique(group)) {
                  if (!is.null(uniquegr[j]) & !is.null(truename[i])) {
                  if (!is.na(uniquegr[j]) & !is.na(truename[i])) {
                    if (uniquegr[j]==truename[i]) {colchoose[j]=whichcoltype[r]; r=r+1}
                    if (i!=1 & r>8) {r=1}
                    if (i==1 & r>9) {r=1}
                  } else {uniquegr=NULL}
                 } else {uniquegr=NULL}
                }
              }
            }
            
            colorname=colstyle[1:groupnumber]
            linechoose=0
            
            polysome=data.frame(polysome,realgroup)         
            colnames(polysome)=c("Position","Absorbance","Sample","Group")
            grname=0
            grname2=0
            }
          }
          } else {
            groupnumber1=input$groupnumber1
            groupnumber2=input$groupnumber2
            if (is.null(groupnumber1)|is.null(groupnumber2)) {return(NULL)} else {
            if (is.null(input$m11)|is.null(input$m21)) {return(NULL)} else {
            
            realgroup1=group
            realgroup2=group
            grlist=sfLapply(unique(group),function(i) {
              gr1=eval(parse(text=paste("input$g1",i,sep="")))
              gr2=eval(parse(text=paste("input$g2",i,sep="")))
              return(list(gr1,gr2))
            })
            uniquegr1=sfSapply(grlist,'[[',1)
            uniquegr2=sfSapply(grlist,'[[',2)
            realgroup1=as.character(uniquegr1)[match(realgroup1,unique(realgroup1))]
            realgroup2=as.character(uniquegr2)[match(realgroup2,unique(realgroup2))]
            realgroup=realgroup1
            tempgr=uniquegr1
            
            truename1=unique(realgroup1)
            groupnumber1=length(truename1)
            truename2=unique(realgroup2)
            groupnumber2=length(truename2)
            if (groupnumber1>=groupnumber2) {
              if (groupnumber1>6) {colchoose=0} else {
                colchoose=unique(group)
                for (i in 1:groupnumber1) {
                  whichcoltype=coltype[[i]]
                  r=1
                  for (j in unique(group)) {
                    if (!is.null(uniquegr1[j]) & !is.null(truename1[i])) {
                    if (!is.na(uniquegr1[j]) & !is.na(truename1[i])) {
                      if (uniquegr1[j]==truename1[i]) {colchoose[j]=whichcoltype[r]; r=r+1}
                      if (i!=1 & r>8) {r=1}
                      if (i==1 & r>9) {r=1}
                    } else {uniquegr1=NULL}
                   } else {uniquegr1=NULL}
                  }
                }
              }
              colorname=colstyle[1:groupnumber1]
              
              if (groupnumber2>6) {linechoose=0} else {
                  linechoose=unique(linetype[1:groupnumber2][match(uniquegr2,unique(uniquegr2))])
                  }
              
              grname=grname1
            } else {
              realgroup=realgroup2
              realgroup2=realgroup1
              if (groupnumber2>6) {colchoose=0} else {
                colchoose=unique(group)
                for (i in 1:groupnumber2) {
                  whichcoltype=coltype[[i]]
                  r=1
                  for (j in unique(group)) {
                    if (!is.null(uniquegr2[j]) & !is.null(truename2[i])) {
                    if (!is.na(uniquegr2[j]) & !is.na(truename2[i])) {
                      if (uniquegr2[j]==truename2[i]) {colchoose[j]=whichcoltype[r]; r=r+1}
                      if (i!=1 & r>8) {r=1}
                      if (i==1 & r>9) {r=1}
                    } else {uniquegr2=NULL}
                   } else {uniquegr2=NULL}
                  }
                }
              }
              colorname=colstyle[1:groupnumber2]
              
              if (groupnumber1>6) {linechoose=0} else {
                linechoose=unique(linetype[1:groupnumber1][match(uniquegr1,unique(uniquegr1))])
              }
              
              uniquegr1=uniquegr2
              uniquegr2=tempgr
              grname=grname2
              grname2=grname1
            }
            
            polysome=data.frame(polysome,realgroup,realgroup2)
            colnames(polysome)=c("Position","Absorbance","Sample","Group1","Group2")
          }
        }
      }
      }
      }
      
      info3=list(polysome=polysome,realgroup=realgroup,realgroup2=realgroup2,colchoose=colchoose,colorname=colorname,linechoose=linechoose,grname=grname,grname2=grname2,uniquegr=uniquegr,uniquegr1=uniquegr1,uniquegr2=uniquegr2)
      return(info3)  
    })  
  })

  
  
  textInput2=function (inputId, label, value = "",...) 
  {
    tagList(tags$label(label, `for` = inputId), tags$input(id = inputId, type = "text", value = value,...))
  }
  
  
  output$insertname = renderUI({
    group=Data()$group
    if (input$check == TRUE) {
      text_list = sfLapply(unique(group), function(i) {
        textname = paste("n",i,sep="")
        div(style="display:inline-block",textInput2(inputId = textname, label = paste("Sample",i), class="input-mini"))
      })
      do.call(tagList, text_list)
    } else {return(NULL)}
  })
  
  
  output$singlecolor = renderUI({
    group=Data()$group
    if (length(group) != 1)
      return(NULL)
    div(class="row",div(class="span1",""),div(class="span6",div(class="span5",helpText(h4("Single curve color:"))),div(class="span3",selectInput("singlecol", "", choices=c("Blue" = "blue", "Red" = "red", "Orange" = "orange", "Green"="green", "Black" = "black")))))
  })
  
  
  output$groupingnumber = renderUI({
    if (input$check2 == TRUE) {
      radioButtons("groupingno","",list("Single-grouping"="single","Double-groupings"="double"))
    } else {return(NULL)}
  })
  
  
  output$groupno = renderUI({
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      isolate({
        if (groupingno=="double") {          
          return(list(h5("Group Information:"),uiOutput("groupno2"),uiOutput("typegroup2"),uiOutput("grouping2"),uiOutput("groupno3"),uiOutput("typegroup3"),uiOutput("grouping3")))
        } else {return(list(uiOutput("groupno1"),uiOutput("typegroup1"),uiOutput("grouping1")))}
      })
    } else {return(NULL)}
  })
  
  
  output$groupno1 = renderUI({
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      isolate({
        no=list(div(style="display:inline-block",numericInput("groupnumber",h5("Grouping Number:"),"",min=1,max=100,step=1)), rep(list(div(style="display:inline-block",helpText(" "))),3), div(style="display:inline-block",textInput2("groupingname",h5("Grouping Name:"), class="input-small")))
      })
    } else {return(NULL)}
  })
  
  output$groupno2 = renderUI({
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      isolate({
        no1=list(div(style="display:inline-block",numericInput("groupnumber1",h5("Grouping 1 Number:"),"",min=1,max=100,step=1)), rep(list(div(style="display:inline-block",helpText(" "))),3), div(style="display:inline-block",textInput2("groupingname1",h5("Grouping 1 Name:"), class="input-small")))
      })
    } else {return(NULL)}
  })
  
  output$groupno3 = renderUI({
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      isolate({
        no2=list(div(style="display:inline-block",numericInput("groupnumber2",h5("Grouping 2 Number:"),"",min=1,max=100,step=1)), rep(list(div(style="display:inline-block",helpText(" "))),3), div(style="display:inline-block",textInput2("groupingname2",h5("Grouping 2 Name:"), class="input-small")))
      })
    } else {return(NULL)}
  })
  
  
  
  output$typegroup1 = renderUI({
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
        if (groupingno=="single") {
          groupnumber=input$groupnumber
          if (!is.numeric(groupnumber)) {return(NULL)}
            text_list2 = sfLapply(1:groupnumber, function(i) {
              textname2 = paste("m",i,sep="")
              div(style="display:inline-block",textInput2(inputId = textname2, label = paste("Group",i), class="input-small"))
            })
            do.call(tagList, text_list2)
        } else {return(NULL)}
    } else {return(NULL)}
  })
  
  output$typegroup2 = renderUI({
    if (input$check2 == TRUE) {
    groupingno=input$groupingno
    if (is.null(groupingno)) {return(NULL)}
    groupnumber1=input$groupnumber1
    if (!is.numeric(groupnumber1)) {return(NULL)}
      text_list4 = sfLapply(1:groupnumber1, function(i) {
        textname4 = paste("m1",i,sep="")
        div(style="display:inline-block",textInput2(inputId = textname4, label = paste("Grouping-1-",i,sep=""), class="input-small"))
      })
      do.call(tagList, text_list4)
    } else {return(NULL)}
  })
  
  output$typegroup3 = renderUI({
    if (input$check2 == TRUE) {
    groupingno=input$groupingno
    if (is.null(groupingno)) {return(NULL)}
    groupnumber2=input$groupnumber2
    if (!is.numeric(groupnumber2)) {return(NULL)}
      text_list5 = sfLapply(1:groupnumber2, function(i) {
        textname5 = paste("m2",i,sep="")
        div(style="display:inline-block",textInput2(inputId = textname5, label = paste("Grouping-2-",i,sep=""), class="input-small"))
      })
      do.call(tagList, text_list5)
  } else {return(NULL)}
  })
  
  
  
  output$grouping1 = renderUI({
    group=Data()$group
    groupname=Data()$groupname
    if (input$check==TRUE) {
      if (!is.null(input$n1)) {
        groupname=Groupname()$groupname 
      }
    }
    
    if (!is.null(groupname)) {
      if (!all(groupname %in% Data()$groupname)) {groupname=paste(1:length(groupname),": ",groupname,sep="")}
    }
    
    if (input$check2 == TRUE) {
    groupingno=input$groupingno
    if (is.null(groupingno)) {return(NULL)}
    if (groupingno=="double") {return(NULL)}
      groupnumber=input$groupnumber
      if (is.null(groupnumber)) {return(NULL)}
      if (is.null(groupname)) {return(NULL)}
      if (!is.numeric(groupnumber)|groupnumber>length(groupname)) {return(NULL)}
        groupchoice=c()
        for (i in 1:groupnumber) {
          groupinput=eval(parse(text=paste("input$m",i,sep="")))
          if (is.null(groupinput)) {groupinput="Ungrouped"}
          groupchoice=c(groupchoice,groupinput)
        }
        select_list = sfLapply(unique(group), function(i) {
          groupings = paste("g",i,sep="")
          div(style="display:inline-block",selectInput(groupings,groupname[i],choices=groupchoice,width="110px"))
        })
        do.call(tagList, select_list)
    } else {return(NULL)}
  })
  
  output$grouping2 = renderUI({
    group=Data()$group
    groupname=Data()$groupname
    if (input$check==TRUE) {
      if (!is.null(input$n1)) {
        groupname=Groupname()$groupname 
      }
    }
    
    if (!is.null(groupname)) {
      if (!all(groupname %in% Data()$groupname)) {groupname=paste(1:length(groupname),": ",groupname,sep="")}
    }
    
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      if (groupingno=="double") {
        groupnumber1=input$groupnumber1
        if (is.null(groupnumber1)) {return(NULL)}
        if (is.null(groupname)) {return(NULL)}
        if (is.numeric(groupnumber1) & groupnumber1<=length(groupname)) {
          groupchoice2=c()
          for (i in 1:groupnumber1) {
            groupinput2=eval(parse(text=paste("input$m1",i,sep="")))
            if (is.null(groupinput2)) {groupinput2="Ungrouped"}
            groupchoice2=c(groupchoice2,groupinput2)
          }
          select_list2 = sfLapply(unique(group), function(i) {
            groupings = paste("g1",i,sep="")
            div(style="display:inline-block",selectInput(groupings,groupname[i],choices=groupchoice2,width="110px"))
          })
          do.call(tagList, select_list2)
        } else {return(NULL)}
      } else {return(NULL)}
      } else {return(NULL)}
  })
  
  output$grouping3 = renderUI({
    group=Data()$group
    groupname=Data()$groupname
    if (input$check==TRUE) {
      if (!is.null(input$n1)) {
        groupname=Groupname()$groupname 
      }
    }
    
    if (!is.null(groupname)) {
      if (!all(groupname %in% Data()$groupname)) {groupname=paste(1:length(groupname),": ",groupname,sep="")}
    }
    
    if (input$check2 == TRUE) {
      groupingno=input$groupingno
      if (is.null(groupingno)) {return(NULL)}
      if (groupingno=="double") {
        groupnumber2=input$groupnumber2
        if (is.null(groupnumber2)) {return(NULL)}
        if (is.null(groupname)) {return(NULL)}
        if (is.numeric(groupnumber2) & groupnumber2<=length(groupname)) {
          groupchoice3=c()
          for (i in 1:groupnumber2) {
            groupinput3=eval(parse(text=paste("input$m2",i,sep="")))
            if (is.null(groupinput3)) {groupinput3="Ungrouped"}
            groupchoice3=c(groupchoice3,groupinput3)
          }
          select_list3 = sfLapply(unique(group), function(i) {
            groupings = paste("g2",i,sep="")
            div(style="display:inline-block",selectInput(groupings,groupname[i],choices=groupchoice3,width="110px"))
          })
          do.call(tagList, select_list3)
        } else {return(NULL)}
      } else {return(NULL)}
    } else {return(NULL)}
  })

  
  output$singlesetting = renderUI({
    div(class="row",div(class="span1", ''),div(class="span4", checkboxInput("checkfit",h5("Show normal fitting curve"),FALSE)),div(class="span1", ''),div(class="span5",div(class="span4",helpText(h4("Curve color:"))),div(class="span4",selectInput("colour2", "", choices=c("Blue" = "blue", "Red" = "red", "Orange" = "orange", "Green"="green", "Black" = "black")))))
  })
  
  observe({
    polysome=Data()$polysome
    if (is.null(polysome)) 
      return(NULL)
    group=Basicplot()$group
    top=Data()$top
    scale=Data()$scale
    lineweight=input$lineweight
    col=input$colour2
    groupname=Basicplot()$groupname
    xmax=Data()$xmax
    xscale=Data()$xscale
    plotchoose=input$singleplotchoose
    fm=NULL
    if (is.null(groupname)) {groupname=Data()$groupname}
    if (is.null(input$checkfit))
      return(NULL)
    if (input$checkfit==T) {
      polysome=Normalize()$normalize
      if (is.null(polysome)) 
        return(NULL)
      Fm=polysome
      fm=CalcAUC()$Fm
        if (is.null(fm)) {
          addfit=0
        } else {
        Fm$Absorbance=fm
        addfit=rep(c(0,1),each=length(CalcAUC()$Fm))
        polysome=rbind(polysome,Fm)
      }
    } else {
      addfit=0
    }
    
    polysome=data.frame(polysome,addfit)
    
    if (length(group)>1) {
      panelrow=ceiling(length(groupname)/4)
      local({
        output$latticeplot = renderPlot({
          if (is.null(Data()$groupname))
            return(NULL)
          theme.novpadding =
            list(layout.heights =
                   list(top.padding = 0,
                        main.key.padding = 0,
                        key.axis.padding = 0,
                        axis.xlab.padding = 0,
                        xlab.key.padding = 0,
                        key.sub.padding = 0,
                        bottom.padding = 0),
                 layout.widths =
                   list(left.padding = 0,
                        key.ylab.padding = 0,
                        ylab.axis.padding = 0,
                        axis.key.padding = 0,
                        right.padding = 0),
                 strip.background=list(col="lightblue"))
          xyplot(Absorbance~Position|factor(group,labels=groupname),groups=addfit,data=polysome,type="l",lwd=lineweight,strip=T,lty=c(1,2),col.line=c(col,"black"),aspect=0.5,layout=c(4,panelrow),as.table=T,par.settings=theme.novpadding)
        })
        


        output$singleplot = renderPlot({
          if (is.null(plotchoose))
            return(NULL)
          if (is.na(plotchoose))
            return(NULL)
          singleribosome=subset(polysome,group==plotchoose,select=-group)
          singlename=paste("Ribosome Profiling","-",groupname[plotchoose],sep="")
          if (max(singleribosome$Absorbance,na.rm=T)>1) {top2=round(max(singleribosome$Absorbance,na.rm=T)+0.1,1);scale2=0.1} else {top2=1;scale2=0.25}
          
          if (is.null(fm)) {
            p1=ggplot(data=singleribosome,aes(x=Position,y=Absorbance))+geom_line(colour=col,size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="",linetype="")+ggtitle(singlename)
            p1=p1+scale_x_continuous(breaks=c(0,seq(0,xmax,xscale)))+scale_y_continuous(breaks=seq(0,top2,scale2))+coord_cartesian(xlim=c(0,xmax),ylim=c(0,top2))+theme(plot.margin = unit(c(0.5,0.7,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.title=element_text(face="bold",size=14),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
          } else {
            p1=ggplot(data=singleribosome,aes(x=Position,y=Absorbance,group=addfit,colour=factor(addfit),linetype=factor(addfit)))+geom_line(size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="",linetype="")+ggtitle(singlename)
            p1=p1+scale_colour_manual(breaks=c(0,1),values=c(col,"black"),labels=c("ribosome profiling curve","normal fitting curve"))+scale_linetype_manual(breaks=c(0,1),values=c("solid","longdash"),labels=c("ribosome profiling curve","normal fitting curve"))+scale_x_continuous(breaks=c(0,seq(0,xmax,xscale)))+scale_y_continuous(breaks=seq(0,top2,scale2))+coord_cartesian(xlim=c(0,xmax),ylim=c(0,top2))+theme(plot.margin = unit(c(0.5,0.7,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
          }
          print(p1)
        },width=750,height=400)
      })
        
    } else {
      local({
        output$singleplot = renderPlot({
          singlename=paste("Ribosome Profiling","-",groupname[1],sep="")
          if (is.null(fm)) {
            p2=ggplot(data=polysome,aes(x=Position,y=Absorbance))+geom_line(colour=col,size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="",linetype="")+ggtitle(singlename)
            p2=p2+scale_x_continuous(breaks=c(0,seq(0,xmax,xscale)))+scale_y_continuous(breaks=seq(0,top,scale))+coord_cartesian(xlim=c(0,xmax),ylim=c(0,top))+theme(plot.margin = unit(c(0.5,0.7,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.title=element_text(face="bold",size=14),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y = element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          } else {
            p2=ggplot(data=polysome,aes(x=Position,y=Absorbance,group=addfit,colour=factor(addfit),linetype=factor(addfit)))+geom_line(size=lineweight)+labs(x="Position (mm)",y="Absorbance (254nm)",colour="",linetype="")+ggtitle(singlename)
            p2=p2+scale_colour_manual(breaks=c(0,1),values=c(col,"black"),labels=c("ribosome profiling curve","normal fitting curve"))+scale_linetype_manual(breaks=c(0,1),values=c("solid","longdash"),labels=c("ribosome profiling curve","normal fitting curve"))+scale_x_continuous(breaks=c(0,seq(0,xmax,xscale)))+scale_y_continuous(breaks=seq(0,top,scale))+coord_cartesian(xlim=c(0,xmax),ylim=c(0,top))+theme(plot.margin = unit(c(0.5,0.7,0.5,0.5), "cm"), plot.title = element_text(face="bold",size=23),axis.title.x = element_text(face="bold", vjust=-0.4, size=20),axis.title.y = element_text(face="bold", vjust=1.2, size=19),legend.key=element_rect(fill="white"),axis.text.x=element_text(colour="black",face="bold",size=16),axis.text.y=element_text(colour="black",face="bold",size=16),panel.background=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
          }
          print(p2)
        },width=750,height=400)
      })
    }
  })
  
  output$lattice=renderUI({
    groupname=Basicplot()$groupname
    if (length(groupname)==1)
      return(NULL)
    panelrow=ceiling(length(groupname)/4)
    plotheight=paste(panelrow*(200-panelrow*20),"px",sep="")
    return(plotOutput("latticeplot",width="100%",height=plotheight))
  })
  
  
  observe({
    if (is.null(input$checkfit))
      return(NULL)
    if (input$checkfit==F)
      return(NULL)
    if (!is.null(CalcAUC()$Fm))
      return(NULL)
    showshinyalert(session, "shinyalert1", paste("Please normalize and caculate AUC first.", "<button type='button' class='btn btn-success'>Okay</button>"), 
                   styleclass = "success")
  })
  
  
  output$plotchoose = renderUI({
    groupname=Data()$groupname
    if (is.null(groupname)) 
      return(NULL)
    if (length(groupname)==1) 
      return(NULL)
    div(class="row",div(class="span1", ''),div(class="span3",numericInput("singleplotchoose",h5("Sample Number:"),"",min=1,max=length(groupname),step=1)))
  })
  
  
  output$totalcolour = renderUI({
    if (input$colour1 == "classical") {return(NULL)}
    ccselect=selectInput('similar1','Select hue of similar colors:',choices=c("Cold" = "cold", "Warm" = "warm", "Grey" = "grey"))
    if (!is.null(input$check2)) {
      if (input$check2==TRUE) {
        ccselect=NULL
      }
    }
    return(ccselect)
  })
  
  
  output$NormalizeDIY = renderUI({
    if (is.null(Judge()$judgefile))
      return(NULL)
    graph1 = Basicplot()$graph1
    if (is.null(graph1))
      return(NULL)
    div(class="row",div(class="span1", ''),div(class="span7",uiOutput("Autonormalize"),uiOutput("DIY")))
  })
  
  output$Autonormalize = renderUI({
    if (is.null(Judge()$judgefile))
      return(NULL) 
    graph1 = Basicplot()$graph1
    if (is.null(graph1))
      return(NULL)
    div(class="span4",actionButton("autonormalize", "Auto Normalization",styleclass="success"))
  }) 
  
  output$DIY = renderUI({  
    if (is.null(input$autonormalize)) {return(NULL)} else {
      if (input$autonormalize == 0) {return(NULL)}
    }
    isolate({
      div(class="span6",bsButton("DIYbutton", "Open Customization Window", style="info"))
    })
  })
  
  
  
  observe({
    if (is.null(Normplot()$graph2))
      return(NULL)
    if (input$autonormalize == 0) 
      return(NULL)
    groupname=Basicplot()$groupname
    local({
      output$DIYModal = renderUI({
        mybsModal("DIYwindow", "Graph Customization", trigger = "DIYbutton",
                  tags$div(class = "row-fluid",
                           tags$div(class = "span3 well control-panel",
                                    textInput("motitle", h5("Graph Title:"), "Ribosome Profiling"),
                                    sliderInput("lineweight2",h5("Line weight:"),min=0,max=1.8,value=0.9,step=0.05),
                                    checkboxGroupInput("removesample", h5("Remove Sample:"), as.character(groupname), inline=TRUE),
                                    helpText(h5("Move curve:")),
                                    uiOutput("movesliders"),
                                    busyIndicator(text = "Plotting...",wait = 0)
                           ),
                           tags$div(class = "span9",
                                    plotOutput("flexplot")
                           )
                  )
        )
      })
    })
  })
  
  toggleModal(session,"DIYwindow",data=NULL)
  
  output$movesliders = renderUI({
    groupname=Basicplot()$groupname
    if (is.null(groupname))
      return(NULL)
    xmax=Data()$xmax
    move_output_list = sfLapply(1:length(groupname), function(i) {
      slidername = paste("move", i, sep="")
      sliderInput(slidername,groupname[i],min=-xmax/10,max=xmax/10,0,step=0.01)
    })
    do.call(tagList, move_output_list)
  })
  
  
  output$Showauc = renderUI({
    if (is.null(Data()$groupname))
      return(NULL)
    if (is.null(Normplot()$graph2))
      return(NULL)
    div(class="row",div(class="span1",""),div(class="span3",actionButton("showauc","Show AUC",styleclass = "warning")))
  })
  
  
  output$graph1 = renderPlot({
    if (is.null(Judge()$judgefile))
      return(NULL)
    graph1 = Basicplot()$graph1
    if (is.null(graph1))
      return(NULL)
    lineweight=input$lineweight
    graph1 = graph1+scale_size(range=c(lineweight, lineweight), guide=FALSE)
    
    if (is.null(graph1))
      return(NULL)
    print(graph1)
  },width=800,height=400)
  
  
  
  output$graph2 = renderPlot({
    if (is.null(Judge()$judgefile))
      return(NULL)
    if (is.null(input$autonormalize)) 
      return(NULL)
    graph2 = Normplot()$graph2
    if (is.null(graph2))
      return(NULL)
    
    lineweight=input$lineweight
    graph2 = graph2+scale_size(range=c(lineweight, lineweight), guide=FALSE)
    
    print(graph2)
  },width=800,height=400)
  
  
  
  observe({
    if (is.null(Judge()$judgefile))
      return(NULL)
    if (is.null(input$autonormalize)) 
      return(NULL)
    if (input$autonormalize == 0) 
      return(NULL)
    
    group=Data()$group
    if (is.null(group))
      return(NULL)
    groupname=Basicplot()$groupname
    check2=Basicplot()$check2
    beautify=Basicplot()$beautify
    lineweight2=input$lineweight2
    
    num=1:length(groupname)
    remove=input$removesample
    if (!is.null(remove)) {
      rmmatch=match(remove,groupname)
      if (any(is.na(rmmatch))) {return(NULL)} else {
        num=num[-rmmatch]
      }
    }
    
    normlist=Normalize()$normlist
    if (is.null(normlist)) 
      return(NULL)
    if (length(normlist) != length(groupname))
      return(NULL)
    
    if (length(group)==1|check2 == FALSE) {
        multilines=sfLapply(num,function(i){        
        move=eval(parse(text=paste("input$move",i,sep="")))
        if (is.null(move)) {move=0}
        polysome1=normlist[[i]]
        polysome1$Position=polysome1$Position+move
        polysome1$group=i
        geom_line(data=polysome1,aes(x=Position,y=Absorbance,colour=factor(group)),size=lineweight2)
      })
    } else {
      groupingno=input$groupingno
      normalize=Normplot()$normalize
      if (is.null(normalize)) 
        return(NULL)      
      uniquegr2=Temp()$uniquegr2
      realgroup=Temp()$realgroup
      
      Group=unique(realgroup)
      rect=data.frame(Group)
      beautify=list(beautify,geom_rect(data=rect,aes(xmin=-1,xmax=-1,ymin=-1,ymax=-1,fill=factor(Group))))
      if (!is.null(groupingno)) {
        if (groupingno=="single") {
          multilines=sfLapply(num,function(i){
            move=eval(parse(text=paste("input$move",i,sep="")))
            if (is.null(move)) {move=0}
            polysome1=subset(normalize,Sample==i)
            polysome1$Position=polysome1$Position+move
            geom_line(data=polysome1,aes(x=Position,y=Absorbance,colour=factor(Sample)),size=lineweight2)
          })
        } else {
          normalize$uniquegr2=uniquegr2
          multilines=sfLapply(num,function(i){
            move=eval(parse(text=paste("input$move",i,sep="")))
            if (is.null(move)) {move=0}
            polysome1=subset(normalize,Sample==i)
            polysome1$Position=polysome1$Position+move
            geom_line(data=polysome1,aes(x=Position,y=Absorbance,colour=factor(Sample),linetype=factor(Group2,levels=unique(uniquegr2))),size=lineweight2)
          })
        }
      }
    }
    
    flexplot = ggplot(environment = environment())+multilines+ggtitle(input$motitle)+beautify
    
    output$flexplot = renderPlot({
      if (is.null(Judge()$judgefile))
        return(NULL)
      if (is.null(input$autonormalize)) 
        return(NULL)
      if (input$autonormalize == 0) 
        return(NULL)
      if (length(num) == 0)
        return(NULL)
      print(flexplot)
    },width=900,height=450)
  })
  
  
  
  
  output$graph3 = renderPlot({
    if (is.null(input$showauc)) {return(NULL)} else {
      if (input$showauc == 0) {return(NULL)} else {
        if (is.null(CalcAUC()$MAUC)|is.null(CalcAUC()$PAUC)) {return(NULL)}
      }
    }
    isolate({
      graph3 = AUCplot()$graph3
      if (is.null(graph3))
        return(NULL)
      print(graph3)
    })
  },width=600,height=400)
  
  output$graph4 = renderPlot({
    if (is.null(input$showauc)) {return(NULL)} else {
      if (input$showauc == 0) {return(NULL)} else {
        if (is.null(CalcAUC()$MAUC)|is.null(CalcAUC()$PAUC)) {return(NULL)}
      }
    }
    isolate({
      graph4 = AUCplot()$graph4
      if (is.null(graph4))
        return(NULL)
      print(graph4)
    })
  },width=600,height=400)
  
  output$graph5 = renderPlot({
    if (is.null(input$showauc)) {return(NULL)} else {
      if (input$showauc == 0) {return(NULL)} else {
        if (is.null(CalcAUC()$MAUC)|is.null(CalcAUC()$PAUC)) {return(NULL)}
      }
    }
    isolate({
      graph5 = AUCplot()$graph5
      if (is.null(graph5))
        return(NULL)
      print(graph5)
    })
  },width=600,height=400)
  
  output$graph6 = renderPlot({
    if (is.null(input$showauc)) {return(NULL)} else {
      if (input$showauc == 0) {return(NULL)} else {
        if (is.null(CalcAUC()$MAUC)|is.null(CalcAUC()$PAUC)) {return(NULL)}
      }
    }
    isolate({
      graph6 = AUCplot()$graph6
      if (is.null(graph6))
        return(NULL)
      print(graph6)
    })
  },width=600,height=400)
  
  
  output$table = renderDataTable({
    if (is.null(input$file1))
      return(NULL)
    table1 = Basicplot()$polysome
    if (is.null(table1)) {
      table1 = Data()$polysome
    }
    table1
  })
  
  
  output$information = renderTable({
    if (is.null(input$file1))
      return(NULL)
    isolate({
      judge =  read.csv(input$file1$datapath[1], header=F, na.strings=c("NA","NaN", ""," "), nrow=2, fileEncoding="latin1")
      if (grepl("Version",as.character(judge[2,1]))) {
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
      } else {information=NULL}
      return(information)
    })
  })
  
  
  
  mybsModal <- function(id, title, trigger, ..., href) {
    mo <- tags$div(class = "modal sbs-modal hide fade", id = id, style = "width: 1200px; height :500px; left: 350px; top: 30px",
                   "data-trigger" = trigger,
                   tags$div(class = "modal-header",
                            tags$button(Type = "button", class = "close",
                                        "data-dismiss" = "modal", HTML("&times;")),
                            tags$h3(title)),
                   body <- tags$div(class = "modal-body"),
                   tags$div(class = "modal-footer",
                            tags$div(class="row-fluid",
                                     tags$div(class="span6","   "),
                                     tags$div(class="span6",
                                              tags$div(class="span4",""),
                                              tags$div(class="span2",actionButton("apply", "Apply",style="info")),
                                              tags$div(class="span4",actionButton("getimage", "Get High Quality Graph", style="success")),
                                              tags$div(class="span2",tags$a(href = "#", class = "btn", "data-dismiss" = "modal", "Close"))
                                     )
                            )
                   )
    )
    
    if(!missing(href)) {
      mo <- addAttribs(mo, "data-remote" = href)
    } else {
      mo$children[[2]] <- tagAppendChildren(mo$children[[2]], list = list(...))
    }
    
    return(mo)
  }

  sfStop()
})
