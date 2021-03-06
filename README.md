# Probability-Models
---
---
title: "Flu trend forcast"
author: "Ruyi Ding"
date: "April 12, 2016"
output: pdf_document
---
# data preparation
```{r}
setwd("/Users/dingruyi/Documents/Wustl SP2016/500T")

data7=read.csv("newdata.csv")
data7[,4]=as.numeric(data7[,4])
n=dim(data7)[1]
```

# Beta Bionomial Model
```{r}
bbll = function(par) {
    a=exp(par[1])
    b=exp(par[2])
    sum=0
    for(i in 1:n){
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        result=sum(log(N/D)+lgamma(a+data7[i,3])+lgamma(b+data7[i,2]
               -data7[i,3])-lgamma(a+b+data7[i,2])-log(beta(a,b)))
        sum=sum+result
    }
    return(sum)}
 

lprob = function(par) {
    a=exp(par[1])
    b=exp(par[2])
    
    lprobs = vector("numeric",n) 
    for(i in 1:n){
            N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
            D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
            ND=N/D
          c=1 
       for( j in 1:length(ND)) {
           c=c*ND[j]
       }
        lprobs[i]= c*exp(lgamma(a+data7[i,3])+lgamma(b+data7[i,2]-data7[i,3])-
                             lgamma(a+b+data7[i,2]))/(beta(a,b))
    }
    return(lprobs)}

flurate=function(par){
    a=exp(par[1])
    b=exp(par[2])
    lprobs=vector("numeric",n)
    for(i in 1:n){
        lprobs[i]=(a+data7[i,3])/(a+b+data7[i,2])
    }
    return (lprobs)
}


par.init=c(0,0)
result = optim(par=par.init,fn=bbll,control=list(fnscale=-1))
lprob(par = result$par)
data7$averageflurate = flurate(par = result$par)

library("dplyr")
temp = data7 %>% mutate(nextweek = averageflurate*data7[,2]) %>% arrange(desc(averageflurate)) 

write.csv(temp,"temp.csv")

```

# mixed logit model
```{r}
data7=read.csv("newdata.csv")

n=dim(data7)[1]

# +0.000001 or rescale?

set.seed(914141)

Z = rnorm(10000)

# param [mean,var]
LL.mixedlogit = function(par) {
    mu.draws = par[1] + Z*sqrt(exp(par[2]))
    
    p.draws = 1/(1+exp(-mu.draws));
    notp.draws = 1-p.draws;
    
    llsum = 0;
    
    n=dim(data7)[1]
    for (i in 1:n) {
        x = data7[i,3]
        m = data7[i,2]
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        lprobp = (p.draws)^x*(notp.draws)^(m-x)
        llsum = llsum + sum(log(N/D)) + log(mean(lprobp)+0.000001)  
    }
    return(llsum)
}

probs.mixedlogit = function(par) {
    mu.draws = par[1] + Z*sqrt(exp(par[2]))
    
    p.draws = 1/(1+exp(-mu.draws));
    notp.draws = 1-p.draws;
    n=dim(data7)[1]
    lprobs = vector("numeric",n) 
    for(i in 1:n){
        
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        ND=N/D
        c=1 
        for( j in 1:length(ND)) {
            c=c*ND[j]
        }
        x = data7[i,3]
        m = data7[i,2]
        lprobs[i]= c*mean((p.draws^x)*(notp.draws^(m-x)))
    }
    
    return(lprobs)
}

par.init=c(0,0)
result = optim(par=par.init,fn=LL.mixedlogit,control=list(fnscale=-1))

probs.mixedlogit(result$par)
```

# latent class model + mixed logit model
```{r}
data7=read.csv("newdata.csv")

# data7 = data7[data7[,3]<100 & data7[,2]<500,]

# par [mu1,mu2,pi]
LL.mixedlogit = function(par) {
    mu1=par[1]
    mu2=par[2]
    pi=1/(1+exp(par[3]))
    p.draws1 = 1/(1+exp(-mu1));
    notp.draws1 = 1-p.draws1;
    p.draws2 = 1/(1+exp(-mu2));
    notp.draws2 = 1-p.draws2;
    
    llsum = 0;
    
    n=dim(data7)[1]
    for (i in 1:n) {
        x = data7[i,3]
        m = data7[i,2]
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        lprobp = (p.draws1^x)*(notp.draws1^(m-x))*pi+(p.draws2^x)*(notp.draws2^(m-x))*(1-pi)
        llsum = llsum + sum(log(N/D)) +log(lprobp+0.000001)
    }
    return(llsum)
}

probs.mixedlogit = function(par) {
    mu1=par[1]
    mu2=par[2]
    pi=1/(1+exp(par[3]))
    p.draws1 = 1/(1+exp(-mu1));
    notp.draws1 = 1-p.draws1;
    p.draws2 = 1/(1+exp(-mu2));
    notp.draws2 = 1-p.draws2;
    n=dim(data7)[1]
    lprobs = vector("numeric",n);
    for(i in 1:n){
        
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        ND=N/D
        c=1 
        for( j in 1:length(ND)) {
            c=c*ND[j]
        }
        x = data7[i,3]
        m = data7[i,2]
        lprobs[i]= c*(p.draws1^x)*(notp.draws1^(m-x))*pi+
            c*(p.draws2^x)*(notp.draws2^(m-x))*(1-pi)
    }

    return(lprobs)
}

result = optim(c(0,0,0),LL.mixedlogit,control=list(fnscale=-1))

probs.mixedlogit(result$par)
```

# mixed logit regression
## consider density
```{r}
data7=read.csv("newdata.csv")
# data8[,2]=data8[,2]/100000
data7[,4]=as.numeric(data7[,4])
# par [mu,beta1]
LL.mixedlogit = function(par) {
    mu=par[1]
    beta1=par[2]
    
    llsum = 0;
    
    n=dim(data7)[1]
    for (i in 1:n) {
        x = data7[i,3]
        m = data7[i,2]
        density = data7[i,4]
        newmu = mu+beta1*density
        p.draws = 1/(1+exp(-newmu))
        notp.draws = 1-p.draws;
        N=seq(from=data7[i,2],to=data7[i,3]+1,by=-1)
        D=seq(from=data7[i,2]-data7[i,3],to=1,by=-1)
        lprob = (p.draws^x)*(notp.draws^(m-x))
        llsum = llsum + sum(log(N/D)) + log(lprob+0.000001)   
    }
    return(llsum)
}

probs.mixedlogit = function(par) {
    mu=par[1]
    beta1=par[2]
   
     n=dim(data8)[1]
     lprobs = vector("numeric",n);

    for (i in 1:n) {
        x = data8[i,3]
        m = data8[i,2]
        e = data8[i,4]
        p = data8[i,5]
        newmu = mu+beta1*e+beta1*p
        p.draws = 1/(1+exp(-newmu))
        notp.draws = 1-p.draws;
        lprob = choose(m,x)*(p.draws^x)*(notp.draws^(m-x))
        lprobs[i] = lprob
        }
    
    return(lprobs)
}

result = optim(c(0,0),LL.mixedlogit,control=list(fnscale=-1))

probs.mixedlogit(result$par)
```


```{r}
data9=read.csv("newpo.csv")
n=dim(data9)[1]
m=dim(data9)[2]
pll = function(par,index){
    data=data9[,index]
    r=exp(par[1])
    a=exp(par[2])
    sum=0
    for(i in 1:n){
        x=data[i]
        N=rep(1,x)
        D=seq(from=x,to=1,by=-1)
        result=sum(log(N/D)+lgamma(r+x)-lgamma(r)+r*log(a/(a+1))
                   +x*log(1/(1+a)))
        sum=sum+result
    }
    return(sum)
}

weeklydata=c(0)
for(i in 1:(m-1)){
    result = optim(par=c(0,0),index=i+1,pll,control=list(fnscale=-1))
    r=exp(result$par[1])
    a=exp(result$par[2])
    weeklydata[i]=r/a
}

write.csv(weeklydata,"temp2.csv")
```
