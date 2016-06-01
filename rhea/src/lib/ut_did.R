## file name: ut_did.R
## last modified: October 6, 2015
## author: Ju-Sung Lee, University of Twente

f.run.did <- function(mode="1a",load,robj=F,obj.path="robj/",csv.path="csv/",
                      png=F,png.path="png/",print=F,draw.mode="gui") {

  old.mgp = par("mgp");
  old.mar = par("mar");
  par("mgp"= c(2,.6,0));
  if ("1a" %in% mode) {
    ## buyer spike
    if (missing(load)) load = !exists("x.bs.200");
    if (load) if (robj) {
      load.obj("x.bs.200",path=obj.path);
    } else {
      if (print) catnf("reading",ps(csv.path,"x_bs_200.csv"),"...");
      x.bs.200 <<- read.csv(ps(csv.path,"x_bs_200.csv"));
    }
    if (png) png(ps(png.path,"did0.png"),pointsize=18);
    fig.nB.did(data=x.bs.200,print=print);
    if (png) dev.off();
  }
  if ("1b" %in% mode) {
    if (missing(load)) load = !exists("x.bs.200");
    if (load) if (robj) {
      load.obj("x.bs.200",path=obj.path);
    } else {
      if (print) catnf("reading",ps(csv.path,"x_bs_200.csv"),"...");
      x.bs.200 <<- read.csv(ps(csv.path,"x_bs_200.csv"));
    }
    if (png) png(ps(png.path,"did1aa.png"),pointsize=18);
    fig.nBS.did(data=x.bs.200);
    if (png) dev.off();
  }
  if ("1c.0" %in% mode) {
    ## phase space: trades/buyers
    if (missing(load)) load = !exists("x.yM.10.200");
    if (load) if (robj) {
      load.obj("x.yM.10.200");
    } else {
      catnf("reading",ps(csv.path,"x_yM_10_200.csv"));
      x.yM.10.200 <<- read.csv(ps(csv.path,"x_yM_10_200.csv"));
    }
    if (png) png(ps(png.path,"did5.png"),pointsize=18);
    fig.yM.traj.did(data=x.yM.10.200,x="nT",y="nB",s=1,fx=I,fy=I,
                    xlab="Number of Trades",ylab="Number of Buyers");
    if (png) dev.off();
  }
  if ("1c.1" %in% mode) {
    ## phase space: coastal price/buyers
    if (missing(load)) load = !exists("x.yM.10.200");
    if (load) if (robj) {
      load.obj("x.yM.10.200");
    } else {
      if (print) catnf("reading",ps(csv.path,"x_yM_10_200.csv"));
      x.yM.10.200 <<- read.csv(ps(csv.path,"x_yM_10_200.csv"));
    }
    if (png) png(ps(png.path,"did6.png"),pointsize=18);
    fig.yM.traj.did(data=x.yM.10.200,x="avPriceCoastFront",y="nB",s=1,fx=I,fy=I,
                    xlab="Average Price Coastal Front",
                    ylab="Numbers of Buyers");
    if (png) dev.off();
  }
  if ("1c.2" %in% mode) {
    ## phase space: ticks/price coastal
    if (missing(load)) load = !exists("x.yM.10.200");
    if (load) if (robj) {
      load.obj("x.yM.10.200");
    } else {
      if (print) catnf("reading",ps(csv.path,"x_yM_10_200.csv"));
      x.yM.10.200 <<- read.csv(ps(csv.path,"x_yM_10_200.csv"));
    }
    if (png) png(ps(png.path,"did7.png"),pointsize=18);
    fig.yM.traj.did(data=x.yM.10.200,x="t",y="avPriceCoastFront",s=1,fx=I,fy=I,
                    xlab="Ticks",ylab="Average Price of Coastal Front")
    if (png) dev.off();
  }
  if ("1d" %in% mode) {
    ## spikes/mortgage
    if (missing(load)) load = !exists("x.yM2.200");
    if (load) if (robj) {
      load.obj("x.yM2.200",path=obj.path);
    } else {
      if (print) catnf("reading",ps(csv.path,"x_yM2_200.csv"),"...");
      x.yM2.200 <<- read.csv(ps(csv.path,"x_yM2_200.csv"));
    }
    if (png) png(ps(png.path,"did1.png"),pointsize=18);
    fig.yM2.did(data=x.yM2.200,lwd=.25,FOS=NULL,print=print);
    if (png) dev.off();
  }

  if ("2" %in% mode) {
    if (missing(load)) load = !exists("x.bs.coef.200");
    if (load) if (robj) {
      load.obj("x.bs.coef.200",path=obj.path);
    } else {
      catnf("reading",ps(csv.path,"x_bs_coef_200.csv"),"...");
      x.bs.coef.200 <<- read.csv(ps(csv.path,"x_bs_coef_200.csv"));
    }
    if (png) png(ps(png.path,"did3.png"),pointsize=18);
    plot.coef(data=x.bs.coef.200,Y="dBS",off.div=.75,
              ylab="Buyers - Sellers",scale=F,ylim=c(-500,620));
    if (png) dev.off();
  }

  if ("3" %in% mode) {
    ## depreciation/mds trajectories
    if (missing(load)) load = !exists("xc0.s1"); 
    if (load) f.load.deprec(path=csv.path,print=T); 
    fig.deprec.rdist.did(png=png,png.path=png.path,print=print,
                         mode=draw.mode);
  }

  if ("4" %in% mode) {
    if (missing(load)) load = !exists("x.sr.14.dm1");
    if (load) if (robj) {
      load.obj(c("x.sr.14.dm1","x.sr.20.dm1"),path=obj.path);
    } else {
      if (print) catnf("reading",ps(csv.path,"x_sr_14_dm1.csv"),"...");
      x.sr.14.dm1 <<- read.csv(ps(csv.path,"x_sr_14_dm1.csv"));
      if (print) catnf("reading",ps(csv.path,"x_sr_20_dm1.csv"),"...");
      x.sr.20.dm1 <<- read.csv(ps(csv.path,"x_sr_20_dm1.csv"));
    }
    l.lm.yM2.fos.did(do.res=F,lr=3,r=0,mode=0)
  }
  par("mgp"=old.mgp);
  par("mar"=old.mar);
}

l.lm.yM2.fos.did <- function(data=x.fos.dm1.200,
                             y="avPriceCoastFront",latex=T,mode=0,YM=c(10,30),
                             T=60,scale=1000,fn=I,FOS=c(.14,.20),tvar="t",
                             form="RH+Ins+G+Ins:RH+Ins:G+Ins:RH:G",
                             ## xvar=c("RH","Ins","G"),
                             yvar=c("yM","fos"),
                             lr=3,r=0,do.res=F,xlab=NULL,
                             yvar.fmt=c("%d","%.2f"),
                             ret=F,print=T,combine.pval=F,
                             outfile="",
                             ...) {
  ## Example 4
  ## data=rbind(cbind(fos=.14,x.sr.14.dm1),cbind(fos=.20,x.sr.20.dm1));
  if (outfile != "") sink(outfile);
  data[[y]] = fn(data[[y]]/scale);
  ## FOS = c(.14,.20);
  ## data = subset(data,t==T);
  data = data[data[[tvar]] %in% T,];
  k = 1;
  resl = list();
  lab = c();
  ## forml = list();
  ## forml[[1]] = ps(y,"~",xvar[1],"+",xvar[2],"+",xvar[3],"+",xvar[2],":",xvar[1],"+", xvar[2],":",xvar[3],"+",xvar[2],":",xvar[1],":",xvar[3]);
  ## forml[[2]] = ps(y,"~",xvar[1],"+",xvar[2],"+",xvar[3],"+",xvar[2],":",xvar[1],"+",xvar[1],":",xvar[3],"+",xvar[2],":",xvar[1],":",xvar[3]);
  form = ps(y,"~",form);
  for (i in 1:length(YM)) {
    for (j in 1:length(FOS)) {
      ## data1 = subset(data,yM == YM[i] & fos == FOS[j]);
      data1 = data[data[[yvar[1]]] == YM[i] & data[[yvar[2]]] == FOS[j],];
      if (F) {
        if (mode == 0) 
          resl[[k]] = lm(data1[[y]]~RH+Ins+G+Ins:RH+Ins:G+Ins:RH:G,data=data1);
        if (mode == 1)
          resl[[k]] = lm(data1[[y]]~RH+Ins+G+Ins:RH+Ins:G+RH:G+Ins:RH:G,data=data1);
      } else {
        ## catnf(forml[[mode+1]]);
        ## resl[[k]] = lm(as.formula(forml[[mode+1]]),data=data1);
        ## catnf(form);
        resl[[k]] = lm(as.formula(form),data=data1);
      }
      k = k + 1;
      ## lab = c(lab,sprintf("\\T{\\scriptsize %s}",ps("\\ttf{yM}=",YM[i],", \\ttf{fos}=",sprintf("%.2f",FOS[j]))));
      lab = c(lab,ps(yvar[1],"=",sprintf(yvar.fmt[1],YM[i]),",",yvar[2],"=",
        sprintf(yvar.fmt[2],FOS[j])));
    }
  }
  if (!latex) return(resl);
  names(resl) = 1:4;
  if (is.null(xlab)) {
    if (F) {
      xn = sprintf("\\T{%s}",c("Intercept"));
      xn = c(xn,sapply(c("Insurance","Adaptive","Growth"),function(a)
        sprintf("%s\\T{%s}",substr(a,1,1),substr(a,2,nchar(a, keepNA = FALSE)))));
      if (mode == 0) 
        xn = c(xn,"I \\stimes A","I \\stimes G","I \\stimes A \\stimes G");
      if (mode == 1) 
        xn = c(xn,"I \\stimes A","I \\stimes G","A \\stimes G","I \\stimes A \\stimes G");
    }
    xn = c("Intercept","I(nsurance)","A(daptive)","G(rowth)",
      "IxA","IxG","IxAxG");
  } else {
    xn = xlab;
  }
  z = f.fmt.lm(resl=resl,use.first=1,xn=xn,lr=lr,r=r,do.res=do.res,...);
  z[1,] = c("\\T{Predictor}",c(rbind(lab,"")));
  z[1,] = gsub("\\\\ttf{([a-zA-Z]*)}","\\1",z[1,],perl=T)
  z[1,] = gsub("\\\\T{(.*)}","\\1",z[1,],perl=T)
  z[1,] = gsub("\\\\scriptsize ","",z[1,],perl=T)
  z[1,] = gsub(", ",",",z[1,],perl=T)
  zs = strsplit(z[1,],",");
  z1 = sapply(zs,function(a) a[1]);
  z2 = sapply(zs,function(a) if (length(a) == 1) "" else a[2]);
  z1 = ifelse(is.na(z1),"",z1);
  z2 = ifelse(is.na(z2),"",z2);
  z = rbind(z1,z2,z[-1,]);

  z[,1] = gsub(" [\\]stimes ","x",z[,1]);
  z[,1] = gsub("\\\\T{(.*)}","(\\1)",z[,1],perl=T)
  z[,1] = gsub("(Intercept)","Intercept",z[,1],fixed=T);
  z[,1] = gsub("\\mathcal{L}","Adj-R^2",z[,1],fixed=T);
  z[,-1] = gsub("\\\\s{(.*)}","(\\1)",z[,-1],perl=T)
  z[,-1] = gsub("\\\\p{(.*)}","\\1",z[,-1],perl=T)
  z[,1] = realign.col(z[,1]);
  for (i in 2:ncol(z)) {
    pre = i %% 2 == 0;
    z[,i] = realign.col(z[,i],pre=pre);
  }
  rownames(z) = NULL;
  if (print) {
    catlm(z,do.array=T,align.header=T,align=T,header=F,simple=T,
          use.colnames=F,
          sep="",end="");
  }
  if (outfile != "") sink(NULL);
  if (ret) {
    if (combine.pval) { ## doesn't work well cos kable eliminates whitespace
      z1 = z[,1];
      for (i in seq(2,ncol(z),2)) {
        z1 = cbind(z1,c(ps(z[,i],align.col(z[,i+1]))));
      }
      rownames(z1) = NULL;
      colnames(z1) = z1[1,];
      z = z1;
    }
    ## colnames(z) = z[1,];
    return(z);
  }
}

fig.yM2.did <- function(data=x.yM2.200,ext="",pause=F,lwd=1,quality=75,
                        YM=seq(10,40,by=10),yM.fn=f.yM2.did,labx=60,
                        fn="fig_BS_yM2_grid",
                        ## FOS=c(.14,.20),
                        FOS=NULL,
                        print=F,
                        SR=NULL,mfrow=c(8,4),height=5.5*2,width=4*2,ps=14,
                        vlwd=1,vlty=2,sd.mult=-1,mode=0,S=0:24,
                        path="../tex/rhea/fig/",
                        yvar=c("nB","nS"),
                        zvar=c("RH","Ins","G","yM"),
                        R=0:1,Is=0:1,G=0:1,svar="s",RH,Ins,
                        outfile="",
                        ylab="Buyers/Sellers",
                        xlab="Ticks"
                        ) {
  ## Example 1d
  data = data[data[[svar]] %in% S,];
  if (!missing(RH)) R = RH;
  if (!missing(Ins)) Is = Ins;
  do.matrix = TRUE;
  if (length(R) == 1) do.matrix = FALSE;
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  if (length(data) != 2 && !is.null(FOS)) {
    data0 = data;
    data = list(subset(data0,fos==FOS[1]),subset(data0,fos==FOS[2]));
  }
  yM.fn.str = deparse(substitute(yM.fn));

  old.mfrow = par("mfrow");
  old.mar = par("mar");
  old.oma = par("oma");
  old.mgp = par("mgp");
  if (do.matrix) {
    par(mfrow=mfrow,mar=c(.05,.05,.05,.05),oma=c(.1,.1,.1,.1));
  } else {
    par(mar=c(2.5,2.5,.1,.1));
    par(mgp=c(2,.6,0));
  }
  axlab = !do.matrix;
  if (yM.fn.str == "f.yM2.did") { 
    yM.fn(data=data,pause=pause,axlab=axlab,lwd=lwd,YM=YM,labx=labx,
          print=print,yvar=yvar,zvar=zvar,svar=svar,R=R,Is=Is,G=G,
          xlab=xlab,ylab=ylab);
  } else {
    yM.fn(data=data,pause=pause,axlab=axlab,lwd=lwd,YM=YM,labx=labx,FOS=FOS,
          print=print,SR=SR,mfrow=mfrow,vlwd=vlwd,vlty=vlty,sd.mult=sd.mult,
          mode=mode,S=S,print=print,ylab=ylab,xlab=xlab);
  }
  par(mfrow=old.mfrow,mar=old.mar,oma=old.oma);
  if (outfile != "") dev.off();
  par("mfrow"=old.mfrow);
  par("mar"=old.mar);
  par("oma"=old.oma);
  par("mgp"=old.mgp);
}

f.yM2.did <- function(data=x.yM2.200,pause=T,axlab=T,lwd=1,R=0:1,
                      YM=seq(10,40,by=10),labx=60,print=T,yvar=c("nB","nS"),
                      Is=0:1,G=0:1,
                      zvar=c("RH","Ins","G","yM"),
                      svar="s",
                      ylab=ylab,xlab=xlab,do.matrix=do.matrix) {
  k = 1;
  for (r in R) {
    for (i in Is) {
      for (g in G) {
        for (ym in YM) { 
          ## if (ym == 10 && r == 0 && i == 0 && g == 0) {
          if (ym == YM[1] && r == R[1] && i == Is[1] && g == G[1]) {
            legend = list(pos="topleft",text=paste(r,i,g,", yM =",ym))
          ## } else if (ym == 10) {
          } else if (ym == YM[1]) {
            legend = list(pos="topleft",text=paste(r,i,g))
          } else if (r == R[1] && i == Is[1] && g == G[1]) {
            legend = list(pos="topleft",text=paste("yM =",ym));
          } else legend = NULL;
          data1 = data[data[[zvar[1]]] == r & data[[zvar[2]]] == i & data[[zvar[3]]] == g & data[[zvar[4]]] == ym,];
          ## plot.nB(subset(data,RH==r & Ins==i & G==g & yM == ym),yM=ym,
          plot.nB(data1,yM=ym,
                  sd.mult=-1,rev=T,include.s=T,use.min=T,axlab=axlab,
                  lwd=lwd,print=F,legend=legend,labx=labx,yvar=yvar,svar=svar,
                  ylab=ylab,xlab=xlab);
          if (pause) readline(paste(k,":",ym,r,i,g))
          else {
            if (print) catnf(paste(k,":",ym,r,i,g));
          }
          k = k + 1;
        }
      }
    }
  }
}

fig.yM.traj.did <- function(data=x.yM.10.200,s=0,x="nT",y="nB",fx=I,fy=I,
                            xt="t",Ts=NULL,ret=F,xlab="",ylab="",ext="",
                            height=5,width=6,ps=18,tag="_yM10_s1",
                            mar=c(3,3,.5,.5),
                            outfile="",...) {
  ## Example 1c
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  old.mar = par("mar");
  par("mar"=mar);
  old.mgp = par("mgp");
  par("mgp"=c(2,.6,0));
  seed = s;
  data = subset(data,s == seed & RH == 1 & Ins == 0 & t > 0);
  if (!is.null(Ts)) {
    data = subset(data,t %in% Ts);
  }
  x = fx(data[[x]]);
  y = fy(data[[y]]);
  plot(x,y,col=8,type='l',xlab=xlab,ylab=ylab,lwd=2);
  points(x,y,col=jet.cols(data$t),pch=16);
  par("mar"=old.mar);
  par("mgp"=old.mgp);
  if (outfile != "") dev.off();
  if (ret) return(y);
}

f.aggr.deprec <- function() {
  z = c();
  for (d in 0:2) {
    for (s in 1:4) {
      catnf(d,s);
      var = ps("xc",d,".s",s);
      z = rbind(z,cbind(d=d,s=s,get(var)$r));
    }
  }
  return(z);
}

f.load.deprec <- function(skip=F,path="../did/src/",print=F,old=F) {
  if (!old) {
    assign("x.deprec.60",read.csv("csv/x_deprec_60.csv"),envir=.GlobalEnv);
  } else {
    for (d in 0:2) {
      for (s in 1:4) {
        var = ps("xc",d,".s",s);
        if (skip && exists(var)) next;
        if (print) catnf("reading from",ps(path,"test_trPr_s",s,"_d",d,"/"));
        data = f.deprec.trPr(seed18=T,nlogo=F,ret=T,tick=60,
          path0=ps(path,"test_trPr_s",s,"_d",d,"/"),
          print=print);
        ## print=!print);
        assign(var,data,envir=.GlobalEnv);
      }
    }
  }
}

fig.deprec.rdist.did <- function(data=x.deprec.60,
                                 cex=1,png=F,png.path="png/",
                                 mode="gui",S=c(1,2,4),
                                 D=0:2,print=F,...) {
  ## Example 3
  old.mar = par("mar");
  mar=c(2.0,1.8,.1,.5);
  old.mgp = par("mgp");
  par("mar"=mar);
  par("mgp"=c(2,.6,0))
  for (d in D) {
    for (s in S) {
      if (png) {
        png(ps(png.path,"did2_d",d,"_s",s,".png"),pointsize=18);
        par("mar"=mar);
        par("mgp"=c(2,.6,0))
      }
      plot.deprec.rdist(data=data,dm=d,s=s,lwd=2,mode="rel",ret=FALSE,cex=cex,
                        ...);
      if (mode == "png") {
        dev.off();
        catnf(ps("dm=",d,",seed=",s));
      }
      if (mode == "rmd") {
        if (print) catnf(ps("dm=",d,",seed=",s));
      }
      if (mode == "gui") {
        readline(ps("dm=",d,",seed=",s,", Press <Enter> to continue"));
      }
    }
  }
  par("mar"=old.mar);
  par("mgp"=old.mgp);
}

plot.deprec.rdist <- function(data=x.deprec.60,fn=I,lwd=2,cex=1.5,dm=0,s=1,
                              mode="rel",ret=F,cex.axis=1,lab.line=-2,
                              ignore=c(1,2,3,26-0:4),scalingfn="cmdscale",
                              sseed=1,
                              outfile="",width=7,height=7) {
  ## Example 3
  ## 100215: modified so it's just one dataset
  ## 083015: previously mode="scale"
  ## if (missing(data)) {
  ##  data = get(ps("xc",dm,".s",s),envir=.GlobalEnv)$r;
  ## }
  ## x = data[,c(-1,-(ncol(data)-0:4))];
  ## x = subset(data,d==dm&data$s==s);
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  old.mar = par("mar");
  old.mgp = par("mgp");
  par("mar"=c(2.0,1.8,.1,.5))
  par("mgp"=c(2,.6,0))
  x = data[data$d==dm & data$s==s,];
  ## catnf(nrow(x));
  x = x[,-ignore];
  if (mode == "scale") xs = scale(x);
  if (mode == "rel") {
    xs = apply(x,2,function(a) a/a[1]);
  }
  xd = dist(xs);
  ## return(xd);
  if (is.character(scalingfn)) {
    if (scalingfn == "cmdscale") {
      coord = cmdscale(d=xd,k=2);
    }
    if (scalingfn == "isoMDS") {
      set.seed(sseed);
      coord = isoMDS(d=xd,k=2,tol=1e-10,maxit=10000,trace=F)$points;
    }
    if (scalingfn == "monoMDS") {
      set.seed(sseed);
      coord = monoMDS(d=xd,y=cmdscale(d=xd,k=2),k=2,scaling=F,pc=F,smin=1e-10,sfgrmin=1e-10,sratmax=1-1e-10,niter=10000)$points;
    }
  } else {
    coord = scalingfn(dist=xd,k=2)$points;
  }
  ## dc = as.matrix(dist(coord));
  ## xm = as.matrix(xd);
  ## catnf(sum((dc[upper.tri(dc)] - xm[upper.tri(xm)])^2));
  ## coord = fn(coord);
  n = nrow(coord);
  plot(coord[-1,],type="p",cex.axis=cex.axis);
  cols = jet.cols(0:n);
  for (i in 3:nrow(coord)) {
    lines(coord[c(i-1,i),1],coord[c(i-1,i),2],col=cols[-1][i],lwd=lwd);
  }
  bgs = ifelse(abs((0:60)[-1]-30) <= 17,"black","white");
  shadowtext(coord[-1,1],coord[-1,2],labels=(0:60)[-1],
             col=cols[-1],
             bg=bgs,
             r=.07,cex=cex);
  dm.col = c(2,3,1);
  mtext(ps("dm = ",dm,", seed = ",s),line=lab.line,col=dm.col[dm+1],
        cex=cex*1.5);
  par("mar"=old.mar);
  par("mgp"=old.mgp);
  if (outfile != "") dev.off();
  if (ret) coord else invisible();
}

f.deprec.trPr <- function(tick=1,nlogo=F,form=T,negot=T,diff=F,
                          seed18=T,ret=F,path0,fillCF=F,mode=0,s=1,
                          do.mkPr=T,print=F) {
  if (nlogo) {
    if (!seed18) {
      x = read.csv("../did/test_trPr/stats.csv");
      y = read.csv("../did/test_trPr/trade.csv");
      z = read.csv("d:/usr/proj/ut/did/src/in/dell/runs_deprec0.d/Adaptive.07.0/parcels.csv");
    } else {
      if (missing(path0)) path = "../did/test_trPr_s18_d0/";
      x = read.csv(ps(path0,"stats.csv"));
      T.fn = ps(path0,ps("RS ",s,"-R Adaptive-GA 0.5-T.csv"));
      R.fn = ps(path0,ps("RS ",s,"-R Adaptive-GA 0.5-R.csv"));
      if (!file.exists(ps("ls \"",T.fn,"\""))) {
        T.fn = ps(path0,ps("RS_",s,"-R_Adaptive-Ins_false-T.csv"));
        R.fn = ps(path0,ps("RS_",s,"-R_Adaptive-Ins_false-R.csv"));
      }
      y = read.csv(T.fn);
      z = read.csv("../did/src/runs_deprec0/Adaptive.17.0/parcels.csv");
      r = read.csv(R.fn);
      a = read.csv(ps(path0,"formAskPrice.csv"),header=T);
      ro = read.csv(ps(path0,"R_Output.csv"));
      dp = NULL;
    }
  } else {
    if (!seed18) {
      x = read.csv("../did/src/test_trPr/stats.csv");
      y = read.csv("../did/src/test_trPr/trade.csv");
      z = read.csv("../did/src/test_trPr/parcels.csv");
      a = read.csv("../did/src/test_trPr/formAskPrice.csv",header=F);
      d = read.csv("../did/src/test_trPr/deprec.csv",header=F);
    } else {
      if (missing(path0)) path0 = "../did/src/runs_deprec0/Adaptive.17.0/";
      if (print) catnf(path0);
      x = read.csv(ps(path0,"stats.csv"));
      y = read.csv(ps(path0,"trade.csv"));
      z = read.csv(ps(path0,"parcels.csv"));
      a = read.csv(ps(path0,"formAskPrice.csv"),header=F);
      dp = read.csv(ps(path0,"deprec.csv"),header=F);
      r = read.csv(ps(path0,"realtors.csv"));
      ro = read.csv(ps(path0,"R_Output.csv"));
    }
  }
  if (!nlogo) {
    old = FALSE;
    if (ncol(a) == 12) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","pr","pr0","coastalFront",
                "deprec","pTime","newHome","age","depPr");
    } else if (ncol(a) == 11) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","pr","pr0","coastalFront",
                "deprec","pTime","newHome","age");
    } else if (ncol(a) == 10) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","pr","pr0","coastalFront",
                "deprec","pTime","newHome");
    } else if (ncol(a) == 9) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","pr","pr0","coastalFront",
                "deprec","pTime");
      old = TRUE; 
    } else if (ncol(a) == 8) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","pr","coastalFront",
                "deprec","pTime");
      old = TRUE; 
    } else  if (ncol(a) == 6) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","coastalFront","deprec");
    } else if (ncol(a) == 5) {
      colnames(a) = c("ticks","trPID","trPr","mkPr","coastalFront");
    } else {
      colnames(a) = c("ticks","trPID","trPr","coastalFront");
    }
  }
  if (fillCF) {
  }
  if (mode == 0) {
    z1 = z[z$residential. == "true" & z$coastalFront == 1,];
    za = z[z$residential. == "true",];
    base = sum(z1$askPr);
    price = z1$askPr;      names(price) = z1$trPID;
    pricea = za$askPr;     names(pricea) = za$trPID;
    newHomea = za$newHome; names(newHomea) = za$trPID;
    agea = za$age;         names(agea) = za$trPID;
    
    ma = mc = d = nha = mag = c();
    price0 = price;
    for (i in 1:tick) {
      a1 = a[a$ticks==i & a$coastalFront==1,];
      aa = a[a$ticks==i,];
      if (form) {
        price[ps(a1$trPID)] = a1$trPr;
        pricea[ps(aa$trPID)] = aa$trPr;
        if (!old) {
          newHomea[ps(aa$trPID)] = aa$newHome;
          agea[ps(aa$trPID)] = aa$age;
        }
      }
      mc = cbind(mc,price);
      ma = cbind(ma,pricea);
      nha = cbind(nha,newHomea);
      mag = cbind(mag,agea);
      y1 = y[y$ticks==i & y$coastalFront==1,];
      if (diff) {
        for (j in 1:nrow(y1)) {
          w = which(a1$trPID == y1$trPID[j]);
          if (length(w) != 1) catnf("oops =",i,j,w,y1$trPID[j]);
          pa = a1$trPr[w];
          py = y1$trPr[j];
          if (pa != py) {
            d = rbind(d,c(i,j,pa,py,py-pa));
          }
        }
        w = z1$trPID[!(z1$trPID %in% a1$trPID)];
        u = !(a1$trPID %in% y1$trPID);
      }
      if (negot) price[ps(y1$trPID)] = y1$trPr;
      if (F) {
        catnf(x$avPriceCoastFront[i+1],mean(price),
              x$avPriceCoastFront[i+1]-mean(price),
              (sum(a1$trPr)+sum(d[d[,1]==i,5])+sum(price0[ps(w)]))/nrow(z1),
              (sum(a1$trPr[u])+sum(y1$trPr)+sum(price0[ps(w)]))/nrow(z1)
              );
        
      }
    }
    colnames(ma) = colnames(mc) = 1:tick;
    if (diff) return(d);
    av = NULL;
  }
  if (do.mkPr) {
    for (i in 1:max(a$ticks)) {
      aa = a[a$ticks == i,];
      mkPr = aa$mkPr;
      pTime = aa$pTime;
      names(mkPr) = aa$trPID;
      names(pTime) = aa$trPID;
      yy = y[y$ticks==i,];
      y$mkPr[y$ticks==i] = mkPr[ps(yy$trPID)];
      y$pTime[y$ticks==i] = pTime[ps(yy$trPID)];
    }
  }
  if (mode == 1) {
    price = z$askPr;
    names(price) = z$trPID;
    fc = z$coastalFront == 1 & z$residential. == "true";
    if (tick == -1) tick = length(table(y$ticks));
    av = numeric(tick);
    for (i in 1:tick) {
      y1 = y[y$ticks==i,];
      if (negot) price[ps(y1$trPID)] = y1$trPr;
      av[i] = mean(price[fc]);
    }
    if (ret) return(list(x=x,y=y,z=z,a=a,av=av));
  }
  a$delta = a$pr - a$mkPr;
  y$delta = y$trPr - y$mkPr;
  return(list(x=x,y=y,z=z,a=a,av=av,d=dp,r=r,ma=ma,mc=mc,nha=nha,age=mag,ro=ro));
}

fig.coef.did <- function(...) {
  plot.coef(...);
}

plot.coef <- function(data=x.bs.coef.200,Y="avPriceCoastFront",off.div=2.5,
                      len=15,## ylim=c(100,300),
                      ylim=NULL,
                      ## scale=T,
                      rel=T,
                      ylab=NULL,seed=1,
                      nbcs = c(.60,.65,.70,.75,.80),
                      RH=0:1,Ins=0:1,yvar=c("RH","Ins","nBC"),
                      legend=c("-A -I","-A +I","+A -I","+A +I"),
                      tvar="t",xlab="Year",
                      outfile="",width=7,height=7
                      ) {
  ## Example 2
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  old.mar = par("mar");
  par("mar"=c(2.5,2.5,.1,.1));
  old.mgp = par("mgp");
  par("mgp"=c(2,.6,0))

  if (is.null(data)) data = f.comp.bh(plot=F,comp.rh=F,ret=T);
  if (is.null(ylab)) ylab = "% @ Time 1";
  if (is.null(ylim)) {
    ## ylim = c(100,250);
    X = sort(unique(data[[tvar]]));
    Xf = round(seq(0,max(X),length.out=len))+1;
    X = X/2;
    xlim = range(X);
    xlim[2] = xlim[2]+2;
    
    ylim = c()
    for (k in 1:length(nbcs)) {
      nbc = nbcs[k];
      l = 1;
      for (r in RH) {
        for (i in Ins) {
          ## f = data$I == i & data$RH == r & data$nBC == as.numeric(nbc);
          f = data[[yvar[2]]] == i & data[[yvar[1]]] == r & data[[yvar[3]]] == as.numeric(nbc);
          x1 = data[f,];
          y = tapply(x1[,Y],x1$t,mean);
          y = y[Xf];
          ylim = range(y,ylim);
        }
      }
    }
  }
  X = sort(unique(data[[tvar]]));
  Xf = round(seq(0,max(X),length.out=len))+1;
  X = X/2;
  xlim = range(X);
  xlim[2] = xlim[2]+2;
  Y2 = c("");
  ## cols = jet.cols(0:4);
  cols = jet.cols(0:(length(nbcs)-1));
  plot(0,0,ylim=ylim,xlim=xlim,type="n",xlab="",ylab="");
  ## xlab="Year";
  mtext(xlab,side=1,line=1.5);
  mtext(ylab,side=2,line=1.5);
  grid(lwd=1.5,col=8);
  ## nbcs = c(".60",".65",".70",".75",".80");
  z = z.pch = z.col = c();
  for (k in 1:length(nbcs)) {
    nbc = nbcs[k];
    l = 1;
    for (r in RH) {
      for (i in Ins) {
        ## f = data$I == i & data$RH == r & data$nBC == as.numeric(nbc);
        f = data[[yvar[2]]] == i & data[[yvar[1]]] == r & data[[yvar[3]]] == as.numeric(nbc);
        x1 = data[f,];
        for (j in 1:length(Y2)) {
          yvar1 = Y;
          y = tapply(x1[,yvar1],x1$t,mean);
          y = y[Xf];
          if (rel) {
            lines(X[Xf]+(l-1)/off.div,y/y[1]*100,lty=1,type='b',lwd=2,
                  pch=NA_integer_,
                  col=cols[k]);
            z = rbind(z,cbind(X[Xf]+(l-1)/off.div,y/y[1]*100));
            z.pch = c(z.pch,rep(l,length(y)));
            z.col = c(z.col,rep(cols[k],length(y)));
          } else {
            lines(X[Xf]+(l-1)/off.div,y,lty=1,type='l',lwd=2,
                  col=cols[k]);
            z = rbind(z,cbind(X[Xf]+(l-1)/off.div,y));
            z.pch = c(z.pch,rep(l,length(y)));
            z.col = c(z.col,rep(cols[k],length(y)));
          }
        }
        l = l + 1;
      }
    }
  }
  set.seed(seed);
  o = sample(nrow(z)); 
  points(z[o,1],z[o,2],pch=z.pch[o],col=z.col[o],lwd=2);
  pos = if (Y == "nD" && ylim[1] == -350) "bottomleft" else "topleft"
  legend(pos,
         ## legend=c("-A -I","-A +I","+A -I","+A +I"),
         legend=legend,
         lty=0,
         pch=1:4,
         lwd=2,
         bg=0);
  par("mgp"=old.mgp);
  par("mar"=old.mar);
  if (outfile != "") dev.off();
}

fig.nB.did <- function(data=x.bs.200,print=F,
                       outfile="",width=7,height=7,...
                       ) {
  ## Example 1a
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  old.mar = par("mar");
  old.mgp = par("mgp"); 
  par("mgp"=c(2,.6,0));
  par("mar"=c(2.5,2.5,.1,.1));
  plot.nB(data=data,lwd=.7,axlab=T,print=print,...);
  par("mar"=old.mar);
  par("mgp"=old.mgp);
  if (outfile != "") dev.off();
}

plot.nB <- function(data=x.bs.200,mode=0,ret=F,lwd=1,crit=1300,yM=30,
                    seeds, ## 082015
                    by1=F,rev=F,y0=500,S=0:24,subIns=F,sd.mult=0,include.s=F,
                    use.min=F,axlab=T,print=F,legend=NULL,labx=60,vlwd=1,
                    vlty=2,annotate=F,
                    yvar=c("nB","nS"),ylab="Buyers/Sellers",
                    xlab="Ticks",
                    svar='s',
                    tvar='t',do.matrix=F) {
  if (!missing(seeds)) S = seeds;
  if (is.null(S)) S = sort(unique(data[[svar]]));
  ## data = subset(data,t > 0); 
  data = data[data[[tvar]] > 0,];
  if (subIns) data = data[data$RH==1&data$Ins==0,];
  ylim = c(y0,range(data[[yvar[1]]])[2]);
  if (include.s) { 
    ylim = c(y0,range(data[[yvar[2]]],data[[yvar[1]]])[2]);
  }
  if (use.min) {
    ylim = c(min(data[[yvar[2]]],data[[yvar[1]]]),ylim[2]);
    if (print) catnf(ylim);
  }
  old.xaxt = par("xaxt"); old.yaxt = par("xaxt"); 
  if (!axlab) {
    par(xaxt="n",yaxt="n");
    xlab = ""; ylab = "";
  }
  ## with(data,plot(t,nB,type='n',ylab="",xlab="",ylim=ylim));
  ## plot(data[[tvar]],data[[yvar[1]]],type='n',ylab=ylab,xlab=xlab,ylim=ylim);
  ## axes labels come down below
  plot(data[[tvar]],data[[yvar[1]]],type='n',ylab="",xlab="",ylim=ylim);
  w = c();
  if (mode == 0) {
    if (sd.mult > 0) { 
      u = mean(data[[yvar[1]]]); s = sd(data[[yvar[1]]]);
      crit = u + s*sd.mult;
    } else if (sd.mult == -1) {
      crit = max(data[[yvar[2]]]); 
    }
    for (i in 1:length(S)) {
      x1 = data[data[[svar]]==S[i],];
      if (any(x1[[yvar[1]]] > crit)) w = c(w,i);
    }
  }
  if (mode == 1) {
    if (F) {
      d = c();
      for (i in 1:length(S)) {
        x1 = subset(data,s==S[i]);
        d = c(d,max(x1[[yvar[1]]])-min(x1[[yvar[1]]]));
      }
      u = mean(d);
      s = sd(d);
    }
    for (i in 1:length(S)) {
      ## x1 = subset(data,s==S[i]);
      x1 = data[data[[svar]]==S[i],];
      if (F && which.max(x1$nB) < 100) {
        catnf(S[i],x1$RH[1],x1$Ins[1],x1$G[1],max(x1$nB)-min(x1$nB),
              which.max(x1$nB));
      }
      if (max(x1[[yvar[1]]])>(sd.mult*sd(x1[[yvar[1]]][x1$t>60])+mean(x1[[yvar[1]]][x1$t>60]))) {
        w = c(w,i);
      }
    }
  }
  nw = length(w);
  if (length(w) > 0) {
    if (print) catnf("w =",S[w]);
    Is = c( w, (1:length(S))[-w] ); 
  } else {
    Is = 1:length(S);
  }
  for (i in 1:length(S)) {
    x1 = data[data[[svar]]==S[i],];
    ## with(x1,lines(t,nS,col="#FFCCCC",type='l',lwd=lwd));
    lines(x1[[tvar]],x1[[yvar[2]]],col="#FFCCCC",type='l',lwd=lwd);
  }
  range = 1:length(Is)
  if (rev) range = rev(range);
  for (ii in range) {
    i = Is[ii];
    x1 = data[data[[svar]]==S[i],];
    if (length(w) > 0 && i %in% w) col=ii else col=8;
    ## with(x1,{
    ##  lines(t,nB,type='l',col=col,lwd=lwd);
    ## });
    ## lines(x1$t,x1[[yvar[1]]],type='l',col=col,lwd=lwd);
    lines(x1[[tvar]],x1[[yvar[1]]],type='l',col=col,lwd=lwd);
    if (by1) readline(S[Is[ii]]);
  }
  abline(h=crit,col=8,lwd=1,lty=3);
  abline(v=yM*2,col=3,lwd=vlwd,lty=vlty);
  abline(v=60,col=2,lwd=vlwd,lty=vlty); 
  if (axlab) {
    ## ylab = "Buyers/Sellers";
    mtext(xlab,side=1,line=1.5);
    mtext(ylab,side=2,line=1.5);
  }
  if (!is.null(legend)) {
    shadowtext(labx,ylim[1]+(ylim[2]-ylim[1])*.925,legend$text,
               col=1,bg="white",r=.08,pos=4);
  }
  if (annotate) { 
    if (nw > 0) {
      shadowtext(180,ylim[1]+(ylim[2]-ylim[1])*.055,
                 ps(nw,"/",length(S)),
                 col=1,bg="white",r=.08)
    }
  }
  par("xaxt"=old.xaxt); par("yaxt"=old.yaxt); 
  if (ret) return(data);
}

fig.nBS.did <- function(data=NULL,ext="",lwd=1,abs=T,only=NULL,cex=.3,
                        scale=1000,S=NULL,yvar=c("nB","nS"),tvar="t",
                        svar="s",
                        xlab="Ticks",ylab="Buyers/Sellers",
                        outfile="",width=7,height=7,...) {
  ## Example 1b
  ## 083015: updated for DiD/rmd
  ## 111214: plot the Buyer/Seller for each seed
  ## scale = scale.val;
  ## abs = is.abs;
  if (outfile != "") pdf(file=outfile,width=width,height=height);
  only = S;
  old.mar = par("mar");
  old.mgp = par("mgp");
  par("mgp"=c(2,.6,0))
  par("mar"=c(1,1,.05,.01))
  fn <- function(data=NULL,only=NULL) {
    do.matrix = FALSE;
    if (is.null(only) || all(only == 0:24)) {
      par(mfrow=c(7,4),mar=c(.05,.05,.05,.05),oma=c(1.7,1.5,1,.1));
      ## par(mfrow=c(9,3),
      ## layout(matrix(1:27,byrow=T,nrow=9),widths=rep(1/3,3),heights=rep(1/9,9));
      ## layout(matrix(1:27,byrow=T,nrow=9),
      ##        widths=rep(lcm(2),3),heights=rep(lcm(1.2),9));
      par(mar=c(.05,.05,.05,.05),oma=c(1.7,1.5,1,.1));
      do.matrix = TRUE;
    }
    if (!do.matrix) {
      par("mar"=c(3,3,.05,.01))
    }
    if (is.null(data)) {
      data = x.bs.200[x.bs.200$RH==1&x.bs.200$Ins==0,];
    }
    if (!abs) scale = 1;
    ## if (abs) ylim0 = with(data[data$t > 0,],range(c(nB,nS)/scale));
    if (abs) {
      data1 = data[data[[tvar]] > 0,];
      ylim0 = range(c(data1[[yvar[1]]],data1[[yvar[2]]])/scale);
    }
    Is = if (!is.null(only)) only else 0:24;
    for (i in Is) {
    ## for (i in 0:23) {
      ## par(mar=c(.1,.1,.1,.1));
      yaxt = if (i %in% (c(1,5,9,13,17,21,25)-1)) "s" else "n";
      xaxt = if (i %in% (c(22:25)-1)) "s" else "n";
      if (do.matrix) par(xaxt=xaxt,yaxt=yaxt);
      ## x = data[data$s==i,][-1,];
      x = data[data[[svar]]==i,][-1,];
      ylim = range(c(x[[yvar[1]]],x[[yvar[2]]]));
      if (!abs) ylim0 = ylim;
      ## plot(x[[tvar]],x[[yvar[1]]]/scale,type=if (!is.null(only)) 'o' else 'l',
      plot(x[[tvar]],x[[yvar[1]]]/scale,type='l',
           ylim=ylim0,
           xlab=if (do.matrix) "" else xlab,
           ylab=if (do.matrix) "" else ylab,
           xpd=T,
           lwd=lwd);
      lines(x[[tvar]],x[[yvar[2]]]/scale,col=2,lwd=lwd);
      ## xlab = "Ticks";
      ## ylab = "Buyers/Sellers";
      ## catnf(nrow(x));
      ## w = with(x,which(nB[-1] - nB[-200] > 0));
      w = rep(0,nrow(x));
      ## d = nB[-1] - nB[-200];
      for (j in 1:(nrow(x)-1)) {
        for (k in (j+1):nrow(x)) {
          ## if (x$nB[k] > x$nB[k-1]) {
          if (x[[yvar[1]]][k] > x[[yvar[1]]][k-1]) {
            w[j] = w[j] + 1;
          } else {
            break;
          }
        }
      }
      w = ifelse(w > 3,w,0);
      ## catnf(w);
      if (F) {
      with(x,points(t,nB/scale,
                    ## cex=log(.2*w),
                    ## cex=log(.3*w),
                    cex=cex,
                    col=ifelse(w == 0,0,3),pch=1 ));
      }
      if (F) {
      points(x$t,x$nB/scale,
             ## cex=log(.2*w),
             ## cex=log(.3*w),
             cex=cex,
             col=ifelse(w == 0,1,3),pch=1);
      }
      f = which(w > 0);
      ## points(x$t[f],x$nB[f]/scale,
      points(x[[tvar]][f],x[[yvar[1]]][f]/scale,
             ## cex=log(.2*w),
             ## cex=log(.3*w),
             cex=cex,
             col=3,pch=1);
      
      text(0-10,(ylim0[2]-ylim0[1])*.9+ylim0[1],ps("s = ",i+1),pos=4);
      ## readline(i);
    }
  }
  fn(data=data,only=only);
  par("mgp"=old.mgp);
  par("mar"=old.mar);
  ## fig.call(fn,use.call=T,args=list(data=data,only=only),ext="",fonts="serif",
  ##         ps=ps,
  ##         mar=c(0,0,0,0),fn=ps("../tex/fig/fig_BS_seeds_",
  ##                          if (abs) "raw" else "rel"),
  ##         height=5.5,width=4);
  ## par(mar=c(2,2,2,2));
  ## par(mar=c(.1,2,.1,.1));
  ## par(mfrow = c(3,3));
  ## par("oma"=c(.9,.9,.9,.9));
  ## par("oma"=c(0,0,0,0));
  ## layout(matrix(c(2,0,1,3),2,2,byrow=T));
  ## layout.show(2);
  ## fn();
  if (outfile != "") dev.off();
}

l.details.did <- function(dlab=c("x.yM2.200","x.yM.10.200","x.bs.200",
                            "x.bs.coef.200")) {
  ## 100615: provide table of details of the data sets
  ## data also = list(x.deprec.60)
  add.to.dict = function(var,tag,dict) {
    if (is.null(dict[[var]])) {
      dict[[var]] = tag;
    } else {
      if (dict[[var]] != tag) catnf(var,"has a conflict!");
    }
    return(dict);
  }
  if (length(dlab) > 1) {
    dict.verbose = matrix(c(
      "Ins","Insurance",
      "RH","Adaptive hedonics",
      "G","Income growth",
      "yM","Mortgage years",
      "nBC","New buyers coefficient",
      "s","Random number generator seed",
      "t","Time step (aka tick)",
      "avPriceFP0","Average price of non-flood risk properties",
      "avPriceFP100","Average price of properties with 1/100 flood risk",
      "avPriceFP500","Average price of properties with 1/500 flood risk",
      "avPriceCoastFront","Average price of coastal properties",
      "avPriceNowFP0","Average price of non-flood risk properties (just traded)",
      "avPriceNowFP100","Average price of properties with 1/100 risk (just traded)",
      "avPriceNowFP500","Average price of properties with 1/500 risk (just traded)",
      "avPriceNowCF","Average price of coastal properties (just traded)",
      "nTradedNowFP0","Number of traded properties that are non-flood risk",
      "nTradedNowFP100","Number of just traded properties that are 1/100 risk",
      "nTradedNowFP500","Number of just traded properties that are 1/500 risk",
      "nTradedNowCF","Number of just traded properties that are coastal",
      "nT","Number of parcels just traded",
      "nB","Current number of buyers",
      "nS","Current number of sellers",
      "nNC","Number of newcomers",
      "nOS","Number of old sellers",
      "nOB","Number of old buyers",
      "nBL","Number of buyers who left",
      "dBS","Difference between number of buyers and sellers",
      "nD","Difference between number of buyers and sellers"
      ),ncol=2,byrow=T);
    mode = 0;
  } else if (length(dlab) == 1 && dlab == "x.deprec.60") {
    dict.verbose = matrix(c(
      "d","Depreciation Mode",
      "s","Random number generator seed",
      "time","Time step or tick",
      "intercept","Intercept",
      "bathrooms","Number of bathrooms",
      "bathrooms.2","(Number of bathrooms)^2",
      "age","Age of the house on the property",
      "age.2","Age^2",
      "SQFT","Square footage of the house",
      "SQFT.2","SQFT^2",
      "lotsize","Size of the property",
      "lotsize.2","lotsize^2",
      "newhome","Is the house new?",
      "postFirm","Property is post-FIRM (Flood Insurance Rate Map)",
      "FP100","Is the parcel at 1/100 flood risk?",
      "FP500","Is the property at 1/500 flood risk?",
      "coastalFront","Is the property coastal front?",
      "log.distAmen.","log(Distance to beach)",
      "log.distCBD.","log(Distance to center of business)",
      "log.distHwy.","log(Distance to nearest highway)",
      "log.distPark.","log(Distance to nearest park)",
      ## "town1M","Is the parcel in Beaufort",
      ## "town15B","Is the parcel in Morehead",
      "R2","Adjusted R^2 of the regression"
      ## "regCode","Regression status code"
      ## "i","??"
      ),ncol=2,byrow=T);
    mode = 1;
  }
  rownames(dict.verbose) = dict.verbose[,1];
  data = list();
  for (i in 1:length(dlab)) {
    if (!exists(dlab[i])) {
      catnf("reading",dlab[i]);
      assign(dlab[i],read.csv(ps("csv/",gsub("[.]","_",dlab[i]),".csv")),
             envir=.GlobalEnv);
    }
    data[[i]] = get(dlab[i],envir=.GlobalEnv);
  }
  ## catnf(length(data));
  names(data) = dlab;
  labs = c();
  dict = list();
  for (i in 1:length(data)) {
    ## catnf(dlab[i],colnames(data[[i]]));
    labs = c(labs,colnames(data[[i]]));
    for (j in 1:ncol(data[[i]])) {
      lab = colnames(data[[i]])[j];
      x = data[[i]][,j];
      if (all(x %in% c(0,1))) {
        tag = "integer (binary)";
      } else if (all(x == floor(x))) {
        tag = "integer";
      } else {
        tag = "real";
      }
      dict = add.to.dict(lab,tag,dict);
    }
  }
  labs = unique(labs);
  ## u = which(labs == 's');
  ## v = which(labs == 'nBC');
  ## m = length(labs);
  ## labs = c(labs[1:(u-1)],labs[v],labs[(u:length(labs))[-which(u:m == v)]]);
  put.before <- function(src,pos,labs) {
    p = which(labs == pos);
    s = which(labs == src);
    la = labs[1:(p-1)];
    lb = labs[p:(s-1)];
    lc = labs[(s+1):length(labs)];
    c(la,labs[s],lb,lc);
  }
  if (mode == 0) {
    labs = put.before('nBC','s',labs);
    labs = put.before('nT','nB',labs);
  }
  if (mode == 1) {
    labs = labs[-(length(labs):(length(labs)-4))];
  }
  ## return(labs);
  z = c();
  for (i in 1:length(labs)) {
    ## catnf(i);
    z = rbind(z,c(labs[i],dict[[labs[i]]],dict.verbose[labs[i],2]));
    ## z = rbind(z,c(labs[i],dict[[labs[i]]]));
  }
  ## return(dict.verbose);
  colnames(z) = c("Variable","Data Type","Description");
  z;
}


## utility functions

jet.cols <- function(x,n=100,mode=0,alpha=1) {
  # 020614: added pmax
  if (mode == 0) jet.col2(n,alpha)[scale1(x,n)]
  else if (mode == 1) jet.col2(n,alpha)[scale2(x,n)];
}

jet.col2 <- function (n = 100, alpha = 1) 
{
    red <- c(0, 0, 0, 255, 255, 128)
    green <- c(0, 0, 255, 255, 0, 0)
    blue <- c(143, 255, 255, 0, 0, 0)
    x.from <- c(0, seq(0.125, 1, by = 0.25), 1)
    x.to <- seq(0, 1, length.out = n)
    expand <- function(col) approx(x = x.from, y = col, xout = x.to)$y
    return(rgb(expand(red), expand(green), expand(blue), maxColorValue = 255, 
        alpha = alpha * 255))
}

scale1 <- function(x,n=100) {
  # scale values from either 0.0 - 1.0 or 1-n
  range = maxn(x)-minn(x);
  normed = (x - minn(x))/range;
  if (is.null(n))
    normed
  else
    pmin(floor(normed*(n-1))+1,n,na.rm=T);
}

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
	theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
	
	xy <- xy.coords(x,y)
	xo <- r*strwidth('A')
	yo <- r*strheight('A')
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
		  labels, col=bg, ... )
	}
	text(xy$x, xy$y, labels, col=col, ... )
}

catlm <- function(x,header=T,fill=T,simple=F,align=F,empty=F,
                  align.header=T,
                  align.row=T,
		  do.array=F,breaks=NULL,
		  sep = ' & ',
		  end = ' \\\\',
		  spaces=T,
		  ownline=F,
		  text=F,
                  use.colnames=F,offset=0,use.rownames=F,
                  na.string="--",text1=F,
                  lm=T,aspec,align.codes=c('r','r'),
                  hline=T,pre0) {
  x2 = c();
  x0 = x;
  if (use.colnames) {
    x = rbind(colnames(x),x);
  }
  if (use.rownames) {
    if (use.colnames) col = c("",rownames(x0)) else col = rownames(x0);
    x = cbind(col,x);
  }
  if (text1) {
    x[,1] = sprintf("\\T{%s}",x[,1]);
  }
  if (align) {
    for (i in 1:ncol(x)) {
      x[,i] = ifelse(is.na(x[,i]),na.string,x[,i]);
    }
    if (align.header) {
      pre = if (!missing(pre0)) pre0 else TRUE;
      pre0 = pre; 
      for (i in 1:ncol(x)) {
        if (i == 1 && !align.row) {
          x2 = cbind(x2,align.col(x[,i],pre=F));
        } else {
          x2 = cbind(x2,align.col(x[,i],pre=pre0));
        }
      }
    } else {
      xh = x[1,];
      pre = if (!missing(pre0)) pre0 else TRUE;
      pre0 = pre; 
      for (i in 1:ncol(x)) {
        if (i == 1 && !align.row) {
          x2 = cbind(x2,align.col(x[-1,i],pre=F));
        } else {
          x2 = cbind(x2,align.col(x[-1,i],pre=pre0));
        }
      }
      x2 = rbind(xh,x2);
    }
    x = x2;
  }
  if (!simple) {
    if (!do.array) cat('\\begin{center} ');
    if (fill) {
      if (!do.array) {
        cats('\\begin{tabular*}{1\\textwidth}{@{\\extracolsep{\\fill}}');
        catsn('r|',rep('r',ncol(x)-1),'}');
      } else {
	cats('\\[ \\begin{array}');
        if (missing(aspec)) {
          catsn('{l',rep(' r@{}l',(ncol(x)-1)/2),'}');
        } else {
          catsn("{",aspec,"}");
        }
      }
    } else {
      if (missing(aspec)) {
        catsn("\\begin{tabular}{",align.codes[1],"|",rep(align.codes[2],ncol(x)-1),"}");
      } else {
        catsn("\\begin{tabular}{",aspec,"}");
      }
    }
  }
  if (hline && !simple) catnf("\\hline"); 
  if (is.null(breaks)) {
    for (i in 1:nrow(x)) {
      if (empty) cat(sep);
      if ((ownline == F || i == 1) && lm) {
	cat(x[i,],sep=sep);
      } else {
	if (i %% 2 == 0 || i == nrow(x) || !lm) {
          word = x[i,1];
          word = gsub("^[ ]*","",word);
          cats(word);catsn();
        }
        else if (length(grep("IC",x[i,1])) > 0 ||
                 length(grep("Adj",x[i,1])) > 0) {
          cats(x[i,1]);catsn();
        }
	cat('',x[i,-1],sep=sep);
      }
      catsn(end);
      if (header && i == 1) catn('\\hline');
    }
  } else {
    sep0 = sep;
    len = nchar(x[2,1], keepNA = FALSE)+3+offset;
    at.end = F; 
    for (i in 1:nrow(x)) {
      for (j in 1:ncol(x)) {
	sep = if (j == ncol(x)) end else if (spaces) sep0 else '&';
	
	
	if (i == 1 && header) {
	  if (j == 1 || j %% 2 == 0) {
	      cat(sprintf(if (text) "\\text{%s}" else "%s",
			  gsub(" ","",x[i,j])),sep);
	    next;
	  }
	}
	if (ownline == T && j == 1) {
	  if (i %% 2 == 0 || at.end) {
	    catnf(sprintf(if (text) "\\text{%s}" else "%s",
			  
			  gsub("[ ]+$","",x[i,j])),sep);
	  } else {
	    cat(sep);
	  }
	  next;
	}
	cats(x[i,j],sep);
	if (j %in% breaks) {
	  catn();
	  if (j > 1) cat(paste(rep(if (spaces) " " else "",len),collapse=""));
	  if (ownline==T && i %% 2 ==1) {
	    cat(" ");
	  }
        }
	if (j == ncol(x)) catn();
      }
      if (header && i == 1) catn('\\hline');
      if (length(grep("Adj",x[i,1])) > 0) at.end = T;
    }
  }
  if (hline && !simple) catnf("\\hline");
  if (!simple) {
    if (!do.array) {
      if (fill) { cat('\\end{tabular*} '); }
      else { cat('\\end{tabular} '); }
      catn('\\end{center}');
    } else {
      catn('\\end{array} \\]');
    }
  }
}

realign.col <- function(x,pre=F) {
  M = max(nchar(x, keepNA = FALSE));
  x = gsub("[ ]*$","",x); 
  align.col(x,pre=pre);
}

align.col <- function(x,pre=T) {
  M = 0;
  for (i in 1:length(x)) {
    M = max(nchar(x[i], keepNA = FALSE),M);
  }
  y = character(length(x));
  for (i in 1:length(x)) {
    m = nchar(x[i], keepNA = FALSE);
    if (pre) {
      y[i] = pastes(paste(rep(" ",M-m),collapse=""),x[i]);
    } else {
      y[i] = pastes(x[i],paste(rep(" ",M-m),collapse=""));
    }
  }
  y;
}

f.fmt.lm <- function(x=NULL,resl,b,s,y='\\T{Predictor}',
		     yn,xn,
		     r=3,n,
		     unalign.header=F,l=NULL,p=NULL,do.signif=F,
		     do.glm=F,lr=0,use.first=NULL,
		     do.ll=T,do.n=T,nlab='n',
                     do.glm.p=F,lname,math=F,
                     r2,do.res=F,do.nls=F,mult=1,mult.mode='y',
                     do.summary=T,do.newton=F,do.aic=F,do.bic=F,aic,bic,
                     resl2,lm.replace,sci.big=F,
                     r3 
                     ) {

  if (missing(resl) && !missing(b)) {
    if (missing(yn)) yn = colnames(b);
    if (missing(xn)) xn = rownames(b);
  }
  if (do.res) {
    b = resl[[1]]$coefficients[,1];
  }
  if (missing(r2) && !do.glm & !missing(b)) r2 = rep(r,ncol(b));
  if (missing(lname)) lname = '\\mathcal{L}';
  if (is.list(mult)) {
    mult = unlist(mult);
  }
  if (!missing(resl)) {
    b = s = p = l = n = aic = bic = c();
    yn = names(resl);
    coef1 = coef(resl[[1]]); 
    if (is.matrix(coef1)) {
      if(missing(xn)) xn = rownames(coef1);
      m = nrow(coef1);
    } else {
      if (missing(xn)) xn = names(coef1);
      m = length(coef1);
    }
    if (mult.mode == 'y') {
      if (length(mult) < length(resl)) mult = rep(mult,length(resl));
    }
    if (mult.mode == 'x') {
      if (length(mult) < length(xn)) mult = rep(mult,m);
    }
    not.same = FALSE;
    if (length(resl) > 1) {
      coefl = list();
      for (i in 1:length(resl)) coefl[[i]] = coef(resl[[i]]);
      for (i in c(use.first,(1:length(resl))[-use.first])) {
        if (length(coefl[[i]]) != length(coefl[[use.first]])) not.same = TRUE;
      }
      if (not.same) {
        lab = c();
        for (i in c(use.first,(1:length(coefl))[-use.first])) {
          lab = c(lab,names(coefl[[i]]));
        }
        lab = unique(lab);
      }
    }
    for (i in 1:length(resl)) {
      if (is.null(resl[[i]]$adj.r.squared)) {
        sres = summary(resl[[i]]);
      } else sres = resl[[i]];
      coef = sres$coefficients;
      if (mult.mode == 'y')
        coef[-1,1] = coef[-1,1]*mult[i];
      if (mult.mode == 'x') {
        coef[,1] = coef[,1]*mult;
        coef[,2] = coef[,2]*mult;
      }
      if (not.same) {
        b = cbind(b,coef[,1][lab]);
        s = cbind(s,coef[,2][lab]);
        p = cbind(p,coef[,4][lab]);
        rownames(b) = rownames(s) = rownames(p) = lab;
        xn = lab;
      } else {
        b = cbind(b,coef[,1]);
        s = cbind(s,coef[,2]);
        p = cbind(p,coef[,4]);
      }
      if (!do.nls) {
        l = c(l,sres$adj.r.squared);
        n = if (!is.null(sres[['n']])) c(n,sres[['n']]) else c(n,nrow(resl[[i]]$model));
        if (do.aic) aic = c(aic,AIC(resl[[i]]));
        if (do.bic) bic = c(bic,BIC(resl[[i]]));
      } else {
        if (!missing(resl2)) {
          sres = summary(resl2[[i]]);
          l = c(l,sres$adj.r.squared);
          if (!missing(lm.replace) && lm.replace[[1]] == i) {
            f = lm.replace[[2]];
            s[f,i] = sres$coefficients[,2];
            p[f,i] = sres$coefficients[,4];
          }
        }
        n = c(n,length(resl[[i]]$model[[1]]));
        if (do.aic) aic = c(aic,AIC(resl[[i]]));
        if (do.bic) bic = c(bic,BIC(resl[[i]]));
      }
    }
  } else if (do.glm && !is.null(x)) {
    b = c(); s = c(); l = c(); n = c();
    aic = bic = c();
    if (!is.null(x[[1]]$aic) && !do.aic)
      lname = '\\text{AIC}'
    else lname = '\\text{Adj-}R^2';
    if (is.null(use.first)) use.first = 1;
    names.fn = if (is.matrix(x[[use.first]]$coefficients)) rownames else names;
    if (missing(xn) || is.null(xn)) {
      if (is.null(use.first)) xn = names.fn(x[[1]]$coefficients)
      else xn = names.fn(x[[use.first]]$coefficients);
    }
    M = length(xn);
    if (is.list(x)) yn = names(x);
    xlab0 = NULL;
    for (i in 1:length(x)) {
      if (!is.null(x[[i]][['n']])) n = c(n,x[[i]][['n']])
      else n = c(n,nrow(x[[i]]$model)); 
      res = if (do.summary) summary(x[[i]]) else x[[i]];
      coef = res$coefficients;
      m = if (is.matrix(coef)) nrow(coef) else length(coef);
      if (!is.null(x[[i]]$aic) && !do.aic) {
        l = c(l,round(x[[i]]$aic));
      } else {
        l = c(l,round(res$adj.r.squared,lr));
      }
      xlab = names.fn(x[[i]]$coefficients);
      if (mult.mode == 'y') {
        if (length(mult) < length(x)) mult = rep(mult,length(x));
        coef[-1,1] = coef[-1,1]*mult[i];
      }
      if (mult.mode == 'x') {
        coef[,1] = coef[,1]*mult[xlab];
        coef[,2] = coef[,2]*mult[xlab];
      }
      if (!is.null(use.first)) {
	if (is.null(xlab0)) xlab0 = names.fn(x[[use.first]]$coefficients);
	
	b0 = numeric(length(xlab0));
	
	
	
	b0[which(xlab0 %in% xlab)] = coef[,1];
	b0[which(!(xlab0 %in% xlab))] = NA;
	s0 = numeric(length(xlab0));
	
	s0[which(xlab0 %in% xlab)] = coef[,2];
	s0[which(!(xlab0 %in% xlab))] = NA;
      } else {
	b0 = c(res$coefficients[,1],rep(NA,M-m));
	s0 = c(res$coefficients[,2],rep(NA,M-m));
      }
      b = cbind(b,b0);
      s = cbind(s,s0);

      if (do.aic) aic = c(aic,x[[i]]$aic);
      if (do.bic) bic = c(bic,x[[i]]$bic);
    }
  } else if (do.newton) { 
    if (!is.null(x)) {
      b = c(); s = c(); l = c();
      for (i in 1:length(x)) {
	b = cbind(b,x[[i]][[1]]);

	s = cbind(s,sqrt(diag(x[[i]][[2]])));
	l = c(l,x[[i]][[3]]);
      }
    } else {
      if (is.null(n)) n = rep(800+144+79,ncol(b));
      if (length(n) == 1) n = rep(n,ncol(b));
    }
  } else {
    if (mult.mode == 'x') {
      b0 = b; s0 = s;
      b = s = c();
      catnf(mult);
      for (i in 1:ncol(b0)) {
        b = cbind(b,b0[,i]*mult);
        s = cbind(s,s0[,i]*mult)
      }
    }
  }
  z = matrix('',nrow=nrow(b)*2+1+as.numeric(do.n)+as.numeric(do.ll) +
    as.numeric(do.aic) + as.numeric(do.bic),
    ncol=ncol(b)*2+1); 
  if (math) {
    z[1:nrow(b)*2,1] = sprintf("\\T{%s}",xn);
  } else {
    z[1:nrow(b)*2,1] = xn;
  }
  if (do.aic) z[nrow(z)-3,1] = "\\T{AIC}";
  if (do.bic) z[nrow(z)-2,1] = "\\T{BIC}";
  if (do.ll) z[nrow(z)-1,1] = lname;
  if (do.n) z[nrow(z),1] = nlab;
  z[1,1] = y;
  z[1,1:ncol(b)*2] = yn;
  if (!do.signif) { # as in significant digits
    fmt.b = sprintf("%%.%df",r);
    fmt.s = sprintf("\\s{%%.%df}",r);
    fmt.na = sprintf("%%s");
  } else {
    fmt.b = sprintf("%%s");
    fmt.s = sprintf("\\s{%%s}");
  }
    
  fmt.l = sprintf("%%.%df",lr);

  for (j in 1:ncol(b)) {
    k = 2;
    for (i in 1:nrow(b)) {
      if (is.na(b[i,j])) {
	z[i*2,j*2] = z[i*2+1,j*2] = "";
	next;
      } 
      if (!do.signif) {
        if (sci.big) {
          z[i*2,j*2] = sprintf(fmt.b,if (is.na(b[i,j])) "" else b[i,j]);
          z[i*2+1,j*2] = sprintf(fmt.s,if (is.na(s[i,j])) "" else s[i,j]);
          if (abs(b[i,j]) >= 1000) {
            str = format(signif(b[i,j],2),scientific=T);
            str = gsub("[+]","",str);
            str = gsub("e0","e",str);
            z[i*2,j*2] = str;
          }
          if (abs(s[i,j]) >= 1000) {
            str = format(signif(s[i,j],2),scientific=T);
            str = gsub("[+]","",str);
            str = gsub("e0","e",str);
            z[i*2+1,j*2] = str;
          }
        } else {
          z[i*2,j*2] = sprintf(fmt.b,if (is.na(b[i,j])) "" else b[i,j]);
          z[i*2+1,j*2] = sprintf(fmt.s,if (is.na(s[i,j])) "" else s[i,j]);
        }
      } else {
        this.r = if (missing(r3)) r2[j] else r3[i]; 
        z[i*2,j*2] = sprintf(fmt.b,
           if (is.na(b[i,j])) "" else round(b[i,j],this.r));
        z[i*2+1,j*2] = sprintf(fmt.s,
           if (is.na(s[i,j])) "" else round(s[i,j],this.r));
      }
      pval = if (is.null(p))
	if (do.glm || do.glm.p)
	  2*pnorm(-abs(b[i,j]/s[i,j]))
	    else 2*pt(abs(b[i,j])/s[i,j],n[j],lower.tail=F)
		   else p[i,j];
      z[i*2,j*2+1] = sprintf("\\p{%s}",
         if (is.na(pval)) "" else f.p2stars(pval,latex=T));
    }
    if (do.aic) {
      if (aic[j] > 1e5) {
        str = format(signif(a[j],2),scientific=T);
        str = gsub("[+]","",str);
        str = gsub("e0","e",str);
        z[nrow(z)-3,j*2] = str;
      } else {
        z[nrow(z)-3,j*2] = sprintf("%.0f",round(aic[j],0));
      }
    }
    if (do.bic) z[nrow(z)-2,j*2] = sprintf("%.0f",round(bic[j],0));
    if (do.ll) {
      if (F && l[j] > 1e5) {
        str = format(signif(l[j],2),scientific=T);
        str = gsub("[+]","",str);
        str = gsub("e0","e",str);
        z[nrow(z)-1,j*2] = str;
      } else {
        z[nrow(z)-1,j*2] = sprintf(fmt.l,l[j]);
      }
    }
    if (do.n) z[nrow(z),j*2] = n[j];
  }
  for (i in 1:ncol(z)) {
    if (unalign.header) {
      z[-1,i] = align.col(z[-1,i],pre=i > 1);
    } else {
      z[,i] = align.col(z[,i],pre=i > 1);
    }
  }
  z;
}

f.p2stars <- function(p,align=T,latex=F) {

  caret = if (latex) "\\ca" else "^"

  p.star = if (align) c("   ","***","** ","*  ",
			if (latex) caret else pastes(caret,"  "),"   ") else 
    c("","***","**","*",caret,"");

  ifelse(is.na(p),p.star[1],
	 ifelse(p<.001,p.star[2],
		ifelse(p<.01,p.star[3],
		       ifelse(p<.05,p.star[4],
			      ifelse(p<.1,p.star[5],p.star[6])))));
}

maxn <- function(...) {
  return(if (all(is.na(...))) NA else max(...,na.rm=T));
}

minn <- function(...) {
  return(if (all(is.na(...))) NA else min(...,na.rm=T));
}

pastes <- function(...) {
  paste(...,sep='');
}

catnf <- function(...) {
  cat(...,'\n'); flush.console();invisible();
}

catsn <- function(...) {
  cat(...,'\n',sep='');
}

ps <- function(...) {
  pastes(...);
}

printf <- function(...) {
  print(...); flush.console();
}

load.obj <- function(list,x,
                     path,
                     raw=F,
                     global=T,
                     envir,
                     print=F,
                     ret=F) {
  if (!missing(list)) x = list;
  if (inherits(try(is.character(x),silent=T),"try-error"))
    x = deparse(substitute(x));
  if (exists('globals') && !is.null(globals$obj.path) && missing(path)) {
    path = globals$obj.path;
  }
  for (i in x) {
    fn = if (!raw) pastes(path,'.Robj.',i) else i;
    if (print) catnf('loading',fn);
    if (global) {
      load(file=fn,envir=.GlobalEnv);
    } else if (!missing(envir)) {
      load(file=fn,envir=envir);
    } else {
      load(file=fn);
      return(get(i));
    }
  }
}

exists2 <- function(x) {
  ## do not search R packages
  (x %in% ls(envir=.GlobalEnv));
}

# utrr <- function() source("did/ut_did.R");
utrr <- function() source("ut_did.R");

render2 <- 
function (input, output_format = NULL, output_file = NULL, output_dir = NULL, 
    output_options = NULL, intermediates_dir = NULL, runtime = c("auto", 
        "static", "shiny"), clean = TRUE, params = NULL, envir = parent.frame(), 
    quiet = FALSE, encoding = getOption("encoding")) 
{
    rmarkdown:::perf_timer_start("render")
    if (identical(output_format, "all")) {
        output_format <- enumerate_output_formats(input, envir, 
            encoding)
        if (is.null(output_format)) 
            output_format <- "html_document"
    }
    if (is.character(output_format) && length(output_format) > 
        1) {
        outputs <- character()
        for (format in output_format) {
            output <- render(input = input, output_format = format, 
                output_file = NULL, output_dir = output_dir, 
                output_options = output_options, intermediates_dir = intermediates_dir, 
                runtime = runtime, clean = clean, params = params, 
                envir = envir, quiet = quiet, encoding = encoding)
            outputs <- c(outputs, output)
        }
        return(invisible(outputs))
    }
    required_pandoc <- "1.12.3"
    if (!pandoc_available(required_pandoc)) {
        stop("pandoc version ", required_pandoc, " or higher ", 
            "is required and was not found.", call. = FALSE)
    }
    intermediates <- c()
    on.exit(if (clean) unlink(intermediates, recursive = TRUE), 
        add = TRUE)
    if (!is.null(intermediates_dir)) {
        if (!dir_exists(intermediates_dir)) 
            dir.create(intermediates_dir, recursive = TRUE)
        intermediates_dir <- normalizePath(intermediates_dir, 
            winslash = "/")
    }
    intermediates_loc <- function(file) {
        if (is.null(intermediates_dir)) 
            file
        else file.path(intermediates_dir, file)
    }
    if (!is.null(output_dir)) {
        if (!dir_exists(output_dir)) 
            dir.create(output_dir, recursive = TRUE)
        output_dir <- normalizePath(output_dir, winslash = "/")
    }
    original_input <- normalizePath(input, winslash = "/")
    if (grepl(rmarkdown:::.shell_chars_regex, basename(input))) {
        input_no_shell_chars <- intermediates_loc(file_name_without_shell_chars(basename(input)))
        if (file.exists(input_no_shell_chars)) {
            stop("The name of the input file cannot contain the special shell ", 
                "characters: ", .shell_chars_regex, " (attempted to copy to a ", 
                "version without those characters '", input_no_shell_chars, 
                "' ", "however that file already exists)", call. = FALSE)
        }
        file.copy(input, input_no_shell_chars, overwrite = TRUE)
        intermediates <- c(intermediates, input_no_shell_chars)
        input <- input_no_shell_chars
    }
    oldwd <- setwd(dirname(tools::file_path_as_absolute(input)))
    on.exit(setwd(oldwd), add = TRUE)
    input <- basename(input)
    knit_input <- input
    knit_output <- intermediates_loc(rmarkdown:::file_with_meta_ext(input, 
        "knit", "md"))
    intermediates <- c(intermediates, knit_output)
    utf8_input <- intermediates_loc(rmarkdown:::file_with_meta_ext(input, 
        "utf8", "md"))
    intermediates <- c(intermediates, utf8_input)
    md_input <- identical(tolower(tools::file_ext(input)), "md")
    if (identical(tolower(tools::file_ext(input)), "r")) {
        spin_input <- intermediates_loc(file_with_meta_ext(input, 
            "spin", "R"))
        file.copy(input, spin_input, overwrite = TRUE)
        intermediates <- c(intermediates, spin_input)
        spin_rmd <- knitr::spin(spin_input, knit = FALSE, envir = envir, 
            format = "Rmd")
        intermediates <- c(intermediates, spin_rmd)
        knit_input <- spin_rmd
        metadata <- paste("\n", "---\n", "title: \"", input, 
            "\"\n", "author: \"", Sys.info()[["user"]], "\"\n", 
            "date: \"", date(), "\"\n", "---\n", sep = "")
        if (!identical(encoding, "native.enc")) 
            metadata <- iconv(metadata, to = encoding)
        cat(metadata, file = knit_input, append = TRUE)
    }
    input_lines <- rmarkdown:::read_lines_utf8(knit_input, encoding)
    yaml_front_matter <- rmarkdown:::parse_yaml_front_matter(input_lines)
    if (!rmarkdown:::is_output_format(output_format)) {
        output_format <- rmarkdown:::output_format_from_yaml_front_matter(input_lines, 
            output_options, output_format)
        output_format <- rmarkdown:::create_output_format(output_format$name, 
            output_format$options)
    }
    pandoc_to <- output_format$pandoc$to
    run_citeproc <- rmarkdown:::citeproc_required(yaml_front_matter, input_lines)
    if (is.null(output_file)) 
        output_file <- rmarkdown:::pandoc_output_file(input, output_format$pandoc)
    if (!is.null(output_dir)) {
        output_file <- file.path(output_dir, basename(output_file))
    }
    output_dir <- dirname(output_file)
    files_dir <- file.path(output_dir, rmarkdown:::knitr_files_dir(basename(output_file)))
    files_dir <- pandoc_path_arg(files_dir)
    cache_dir <- NULL
    knit_meta <- NULL
    if (tolower(tools::file_ext(input)) %in% c("r", "rmd", "rmarkdown")) {
        optk <- knitr::opts_knit$get()
        on.exit(knitr::opts_knit$restore(optk), add = TRUE)
        optc <- knitr::opts_chunk$get()
        on.exit(knitr::opts_chunk$restore(optc), add = TRUE)
        hooks <- knitr::knit_hooks$get()
        on.exit(knitr::knit_hooks$restore(hooks), add = TRUE)
        rmarkdown:::knit_meta_reset()
        on.exit(rmarkdown:::knit_meta_reset(), add = TRUE)
        knitr::render_markdown()
        knitr::opts_chunk$set(tidy = FALSE, error = FALSE)
        knitr::opts_knit$set(rmarkdown.pandoc.to = pandoc_to)
        knitr::opts_knit$set(rmarkdown.keep_md = output_format$keep_md)
        knitr::opts_knit$set(rmarkdown.version = 2)
        if (packageVersion("knitr") < "1.5.23") {
            local({
                hook_source = knitr::knit_hooks$get("source")
                knitr::knit_hooks$set(source = function(x, options) {
                  hook_source(strip_white(x), options)
                })
            })
        }
        figures_dir <- paste(files_dir, "/figure-", pandoc_to, 
            "/", sep = "")
        knitr::opts_chunk$set(fig.path = figures_dir)
        cache_dir <- rmarkdown:::knitr_cache_dir(input, pandoc_to)
        knitr::opts_chunk$set(cache.path = cache_dir)
        cache_dir <- gsub("/$", "", cache_dir)
        if (!is.null(output_format$knitr)) {
            knitr::opts_knit$set(as.list(output_format$knitr$opts_knit))
            knitr::opts_chunk$set(as.list(output_format$knitr$opts_chunk))
            knitr::knit_hooks$set(as.list(output_format$knitr$knit_hooks))
        }
        runtime <- match.arg(runtime)
        if (identical(runtime, "auto")) {
            if (!is.null(yaml_front_matter$runtime)) 
                runtime <- yaml_front_matter$runtime
            else runtime <- "static"
        }
        knitr::opts_knit$set(rmarkdown.runtime = runtime)
        if (!is.null(yaml_front_matter$params)) {
            if (packageVersion("knitr") < "1.10") 
                stop("knitr >= 1.10 required to use rmarkdown params")
            knit_params <- mark_utf8(knitr::knit_params(input_lines))
            default_params <- list()
            for (param in knit_params) default_params[[param$name]] <- param$value
            if (!is.null(params)) {
                if (!is.list(params) || (length(names(params)) != 
                  length(params))) 
                  stop("render params argument must be a named list")
                invalid_params <- setdiff(names(params), names(default_params))
                if (length(invalid_params) > 0) 
                  stop("render params not declared in yaml: ", 
                    paste(invalid_params, sep = ", "))
            }
            params <- merge_lists(default_params, params)
            if (!exists("params", envir = envir, inherits = FALSE)) {
                assign("params", params, envir = envir)
                lockBinding("params", envir)
                on.exit({
                  do.call("unlockBinding", list("params", envir))
                  remove("params", envir = envir)
                }, add = TRUE)
            }
            else {
                stop("params object already exists in knit environment ", 
                  "so can't be overwritten by render params", 
                  call. = FALSE)
            }
        }
        env <- environment(render)
        do.call("unlockBinding", list("metadata", env))
        on.exit({
            env$metadata <- list()
            lockBinding("metadata", env)
        }, add = TRUE)
        env$metadata <- yaml_front_matter
        rmarkdown:::perf_timer_start("knitr")
        input <- knitr::knit(knit_input, knit_output, envir = envir, 
            quiet = quiet, encoding = encoding)
        rmarkdown:::perf_timer_stop("knitr")
        rmd_warnings <- rmarkdown:::knit_meta_reset(class = "rmd_warning")
        for (rmd_warning in rmd_warnings) {
            message("Warning: ", rmd_warning)
        }
        knit_meta <- rmarkdown:::knit_meta_reset()
        if (!(rmarkdown:::is_pandoc_to_html(output_format$pandoc) || identical(tolower(tools::file_ext(output_file)), 
            "html"))) {
            if (rmarkdown:::has_html_dependencies(knit_meta)) {
                stop("Functions that produce HTML output found in document targeting ", 
                  pandoc_to, " output.\nPlease change the output type ", 
                  "of this document to HTML.", call. = FALSE)
            }
            if (!identical(runtime, "static")) {
                stop("Runtime '", runtime, "' is not supported for ", 
                  pandoc_to, " output.\nPlease change the output type ", 
                  "of this document to HTML.", call. = FALSE)
            }
        }
    }
    if (output_format$clean_supporting && (is.null(cache_dir) || 
        !rmarkdown:::dir_exists(cache_dir))) 
        intermediates <- c(intermediates, files_dir)
    input_text <- rmarkdown:::read_lines_utf8(input, encoding)
    writeLines(input_text, utf8_input, useBytes = TRUE)
    rmarkdown:::perf_timer_start("pre-processor")
    if (!is.null(output_format$pre_processor)) {
        extra_args <- output_format$pre_processor(yaml_front_matter, 
            utf8_input, runtime, knit_meta, files_dir, output_dir)
        output_format$pandoc$args <- c(output_format$pandoc$args, 
            extra_args)
    }
    if (!is.null(intermediates_dir) && !is.null(output_format$intermediates_generator)) {
        intermediates <- c(intermediates, output_format$intermediates_generator(original_input, 
            encoding, intermediates_dir))
    }
    rmarkdown:::perf_timer_stop("pre-processor")
    if (!is.null(yaml_front_matter$bibliography)) {
        output_format$pandoc$args <- c(output_format$pandoc$args, 
            rbind("--bibliography", pandoc_path_arg(yaml_front_matter$bibliography)))
    }
    rmarkdown:::perf_timer_start("pandoc")
    if (output_format$pandoc$keep_tex) {
        pandoc_convert(utf8_input, pandoc_to, output_format$pandoc$from, 
            file_with_ext(output_file, "tex"), run_citeproc, 
            output_format$pandoc$args, !quiet)
      }
    ## printf(output_format$pandoc$args)
    opa = output_format$pandoc$args
    if (opa[5] == "--latex-engine")
      output_format$pandoc$args[6] = "c:\\bin\\x64\\MIKTEX~1.9\\miktex\\bin\\x64\\pdflatex"
    ## printf(output_format$pandoc$args)
    pandoc_convert(utf8_input, pandoc_to, output_format$pandoc$from, 
        output_file, run_citeproc, output_format$pandoc$args, 
        !quiet)
    if (!is.null(intermediates_dir)) {
        intermediate_output <- file.path(intermediates_dir, basename(output_file))
        if (file.exists(intermediate_output)) {
            file.rename(intermediate_output, output_file)
        }
    }
    rmarkdown:::perf_timer_stop("pandoc")
    rmarkdown:::perf_timer_start("post-processor")
    if (!is.null(output_format$post_processor)) 
        output_file <- output_format$post_processor(yaml_front_matter, 
            utf8_input, output_file, clean, !quiet)
    if (!quiet) {
        message("\nOutput created: ", relative_to(oldwd, output_file))
    }
    rmarkdown:::perf_timer_stop("post-processor")
    rmarkdown:::perf_timer_stop("render")
    if (output_format$keep_md && !md_input) {
        md <- c(md_header_from_front_matter(yaml_front_matter), 
            partition_yaml_front_matter(input_text)$body)
        writeLines(md, file_with_ext(output_file, "md"), useBytes = TRUE)
    }
    invisible(tools::file_path_as_absolute(output_file))
}
