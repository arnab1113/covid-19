input1<-read.csv(file.choose(),header=T)
input2<-read.csv(file.choose(),header=T)


incidence<-function(input,district,dty){
  
  n<-dim(input)[1]
  
  dates<-unique(input$Date.Announced)
  d<-length(dates)
  
  newcases1<-function(i){
    length(seq(1:n)[input$Date.Announced==dates[i] & input$Detected.District==district])}
  
  newcases2<-function(i){
    rowid<-seq(1:n)[input$Date.Announced==dates[i] & input$Detected.District==district]
    sum(input$Num.Cases[rowid])
  }
  
  if(dty==1){vec<-cbind(dates,sapply(seq(1:d),newcases1))
  }else{
    vec<-cbind(dates,sapply(seq(1:d),newcases2))
  }
  return(vec)
}

incid1=incidence(input1,"Indore",1)
incid2=incidence(input2,"Indore",2)

incid=as.numeric(rbind(incid1,incid2)[,2])

day1<-min(seq(1:length(incid))[incid>0])

estR=estimate_R(incid[day1:length(incid)],method="parametric_si",
                config=list(t_start=2:42,t_end=5:45,n1=500,mean_si=3.96,std_si=4.75,
                            n2=100,seed=1,mcmc_control=list(init.pars=NULL,burnin=10000,thin=1000,seed=1)))


