library(mvtnorm)
library(ggplot2)
library(dplyr)


#####Globle parameters
pi <- 0.4
delta1 <- 3
delta2 <- 3
sigma <- 5
type1.error <- 0.05
beta <- 0.8
alpha.min <- 0.0001

###FUNCTION: calculate sample size and return to optimal alpha

IfPowerSatisfy <- function(pi.1,delta1,delta2,sigma,N,type1.error,beta,alpha.min){
        ###functions
        #mean and covariate matrix under alternative
        ConvertToMuSigma <- function(pi.1,delta1,delta2,sigma,N){
                pi.1 <- pi.1
                delta.1 <- delta1
                delta.2 <- delta2
                sigma <- sigma
                #compute Sigma
                pi.2 <- 1-pi.1
                Sigma <- matrix(c(1,sqrt(pi.1),sqrt(pi.2),sqrt(pi.1),1,0,sqrt(pi.2),0,1),nrow = 3,ncol = 3)
                row.names(Sigma) <- c(0,1,2)
                colnames(Sigma) <- c(0,1,2)
                #compute mu
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                N.0 <- N #total sample size
                N.1 <- N.0*pi.1
                N.2 <- N.0*pi.2
                #alternative under 3 different power condition
                #(,1,1)
                mu.0 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                #(,1,0)
                delta.1 <- delta1
                delta.2 <- 0
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                mu.1 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                #(,0,1)
                delta.1 <- 0
                delta.2 <- delta2
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                mu.2 <- c(delta.0/(2*sigma/sqrt(N.0)),delta.1/(2*sigma/sqrt(N.1)),delta.2/(2*sigma/sqrt(N.2)))
                return(list(mu = cbind(mu.0,mu.1,mu.2), Sigma = Sigma))
        }
        #convert to C
        alpha.matrix <- function(alpha.star){
                C <- matrix(0,nrow = 7,ncol = 3)
                type1.error <- sum(alpha.star)
                C[1,1] = C[2,2] = C[3,3] = type1.error
                #test for (c,1)
                C[4,1] = alpha.star[1]
                C[4,2] = alpha.star[2]+alpha.star[3]
                #test for (c,2)
                C[5,1] = alpha.star[1]
                C[5,3] = alpha.star[2]+alpha.star[3]
                #test for (1,2)
                C[6,2] = alpha.star[2]+1/2*alpha.star[1]
                C[6,3] = alpha.star[3]+1/2*alpha.star[1]
                #test for (c,1,2)
                C[7,1] = alpha.star[1]
                C[7,2] = alpha.star[2]
                C[7,3] = alpha.star[3]
                rownames(C) <- c("0","1","2","0,1","0,2","1,2","0,1,2")
                colnames(C) <- c("0","1","2")
                return(C)
        }
        #calculate power
        PowerMB <- function(C,mu,Sigma,j){
                C <- C
                mu <- mu #a matrix under different alternative
                Sigma <- Sigma
                ####Function: local rejection region for H_F
                RejRegionHFBon <- function(F){ 
                        for (i in 1:nrow(C)){
                                if(row.names(C)[i] == F){
                                        c_j <- C[i,]
                                        break;
                                }
                        }
                        alpha.current <- c_j
                        u <- numeric(3)
                        index <- as.numeric(unlist(strsplit(F,",")))
                        index <- index+1
                        for (i in 1:length(index)){
                                u[index[i]] <- qnorm(alpha.current[index[i]],mean = 0, sd = 1,lower.tail = FALSE)
                        }
                        return(u)
                }
                
                ####Function: rejection region for rejecting a single test by closed testing procedure, i.e reject all local tests          ####containing this single test
                RejRegionSingle <- function(j){
                        cut.point <- NULL #store local rejection region, 4 rows
                        for (i in 1:nrow(C)){#local rejection regions
                                if(j %in% unlist(strsplit(rownames(C)[i],","))){
                                        F <- rownames(C)[i]
                                        cut.point[[i]] <- RejRegionHFBon(F)
                                }
                        }
                        #combine all rejection region together
                        cut.point <- data.frame(t(matrix(unlist(cut.point),nrow = 3))) #equivalent with A vector!
                        colnames(cut.point) <- c("0","1","2")
                        
                        ####Function:Calculate all rejection cubes of intersection of unions
                        RejCube <- function(cut.point){ 
                                index <- as.numeric(j)+1
                                index.others <- setdiff(c(1,2,3),index)
                                represent.point <- NULL
                                #construct for test 1
                                a <- cut.point[,index]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- (a[2]+a[3])/2
                                temp[4] <- (a[3]+a[4])/2
                                temp[5] <- a[4]+1
                                represent.point[[index]] <- temp 
                                #construct for test 0
                                index <- index.others[1]
                                a <- cut.point[,index][cut.point[,index]!=0]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- a[2]+1
                                represent.point[[index]] <- temp
                                #constrct for test 2
                                index <- index.others[2]
                                a <- cut.point[,index][cut.point[,index]!=0]
                                a <- a[order(a)]
                                temp <- NULL
                                temp[1] <- a[1]-1
                                temp[2] <- (a[1]+a[2])/2
                                temp[3] <- a[2]+1
                                represent.point[[index]] <- temp
                                
                                IfPointInLocalRegion <- function(p,i){ #i is the row index in cut.point
                                        index <- which(cut.point[i,]!=0)
                                        if(any(p[index]>cut.point[i,][index]))#OR
                                                return(TRUE)
                                        else return(FALSE)
                                }
                                
                                IfPointInRegion <- function(p){
                                        temp <- 1
                                        for (i in 1:nrow(cut.point)){
                                                temp <- temp * IfPointInLocalRegion(p,i)
                                        }
                                        if (temp==1) return(TRUE)
                                        else return(FALSE)
                                }
                                
                                cube <- NULL
                                x <- cut.point[,1][cut.point[,1]!=0]
                                x <- x[order(x)]
                                x <- append(append(x,Inf),-Inf,after=0)
                                y <- cut.point[,2][cut.point[,2]!=0]
                                y <- append(append(y,Inf),-Inf,after=0)
                                y <- y[order(y)]
                                z <- cut.point[,3][cut.point[,3]!=0]
                                z <- z[order(z)]
                                z <- append(append(z,Inf),-Inf,after=0)
                                
                                id <- 1
                                for(i in 1:length(represent.point[[1]])){#for loop
                                        for(j in 1:length(represent.point[[2]])){
                                                for(k in 1:length(represent.point[[3]])){
                                                        p <- c(represent.point[[1]][i],represent.point[[2]][j],represent.point[[3]][k])
                                                        if(IfPointInRegion(p)){
                                                                cube[[id]] <- rbind(c(x[i],x[i+1]),c(y[j],y[j+1]),c(z[k],z[k+1]))
                                                                id <- id+1
                                                        }
                                                }
                                        }
                                }#end of for
                                return(cube)
                        }
                        return(RejCube(cut.point))
                }
                
                ####Function: compute the exact power for rejecting test j under alternative mu_j
                PowerRejSingle <- function(j,mu,Sigma){
                        shape <- RejRegionSingle(j) #rejection region
                        number.of.cube <- length(shape)
                        area <- 0
                        for (i in 1:number.of.cube){
                                lower <- shape[[i]][,1]
                                upper <- shape[[i]][,2]
                                area <- area + pmvnorm(lower = lower, upper = upper, mean = mu[,(as.numeric(j)+1)], sigma = Sigma,algorithm=GenzBretz(abseps = 10^-10 ,maxpts=10^5))
                        }
                        return(as.numeric(area))
                }
                
                return(PowerRejSingle(j,mu,Sigma))
        } 
        
        mu <- ConvertToMuSigma(pi.1,delta1,delta2,sigma,N)$mu
        Sigma <- ConvertToMuSigma(pi.1,delta1,delta2,sigma,N)$Sigma
        type1.error <- type1.error/2
        
        #search over all possible alpha.star values
        for(alpha.star.1 in seq(alpha.min,type1.error,by=0.005)){
                for(alpha.star.2 in seq(alpha.min,type1.error-alpha.star.1,by=0.005)){
                        alpha.star.0=type1.error-alpha.star.1-alpha.star.2
                        alpha.star <- c(alpha.star.0,alpha.star.1,alpha.star.2)
                        C <- alpha.matrix(alpha.star)
                        if((PowerMB(C,mu,Sigma,"0") > beta) & (PowerMB(C,mu,Sigma,"1") > beta) & (PowerMB(C,mu,Sigma,"2") > beta)) return(list(satisfy=TRUE,alpha=alpha.star))
                }
        }
        return(list(satisfy=FALSE,alpha=NULL))
        
}
##example
##if N=250
#IfPowerSatisfy(0.4,3,3,5,250,0.05,0.8,0.0001)
##if N=251
#IfPowerSatisfy(0.4,3,3,5,251,0.05,0.8,0.0001)

#SampleSizeBon <- function(pi.1,delta.1,delta.2,sigma,alpha,beta){
        alpha <- alpha/2
        alpha <- alpha/3 #bon adjustment, one-sided test
        pi.1 <- pi.1
        pi.2 <- 1-pi.1
        delta.0 <- pi.1*delta.1+pi.2*delta.2
        delta.1 <- delta.1
        delta.2 <- delta.2
        sigma <- sigma
        n0.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.0^2)
        n1.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.1^2)
        n2.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.2^2)
        #satisfy all 3 power conditions
        return(ceiling(max(n0.req,n1.req/pi.1,n2.req/pi.2)))
}
#example
#SampleSizeBon(0.4,3,3,5,0.05,0.8)

#####FUNCTION: binary search for sample size and opt alpha
SampleSizeMB <- function(pi.1,delta.1,delta.2,sigma,alpha,beta,alpha.min){
        SampleSizeBon <- function(pi.1,delta.1,delta.2,sigma,alpha,beta){
                alpha <- alpha/2
                alpha <- alpha/3 #bon adjustment, one-sided test
                pi.1 <- pi.1
                pi.2 <- 1-pi.1
                delta.0 <- pi.1*delta.1+pi.2*delta.2
                delta.1 <- delta.1
                delta.2 <- delta.2
                sigma <- sigma
                n0.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.0^2)
                n1.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.1^2)
                n2.req <- (qnorm(1-alpha) + qnorm(beta))^2*sigma^2*4/(delta.2^2)
                #satisfy all 3 power conditions
                return(ceiling(max(n0.req,n1.req/pi.1,n2.req/pi.2)))
        }
        #search from size.bon
        n.upper <- floor(SampleSizeBon(pi.1,delta.1,delta.2,sigma,alpha,beta)) 
        n.lower <- 0
        n.current <- ceiling((n.upper +n.lower)/2)
        #current mu and Sigma
        repeat{
                temp <- IfPowerSatisfy(pi.1,delta.1,delta.2,sigma,n.current,alpha,beta,alpha.min)
                if(temp$satisfy){
                        n.upper <- n.current
                        n.current <- ceiling((n.upper +n.lower)/2)
                        alpha.current <- temp$alpha
                }
                if(!temp$satisfy){
                        n.lower <- n.current
                        n.current <- ceiling((n.upper +n.lower)/2)
                }
                if(n.upper-n.lower <= 1) break
        }
        
        return(list(sample.size=n.upper,alpha.opt=alpha.current))
}
#example
#SampleSizeMB(0.4,3,3,5,0.05,0.8,0.0001)

set.seed(1)
design <- SampleSizeMB(pi,delta1,delta2,sigma,type1.error,beta,alpha.min)


###Simulation setting
##test statistics (z statistics) is the same, to calculate power, see if the 3 dimension z statistic fall in the rejection region
##for each simulation, 3 scenarios are generated for 3 alternatives, must satisfy all 3 power condition

###FUNCTION: calculate rejection region
RejRegion <- function(alpha.star,j){
        ###functions
        #convert to C
        alpha.matrix <- function(alpha.star){
                C <- matrix(0,nrow = 7,ncol = 3)
                type1.error <- sum(alpha.star)
                C[1,1] = C[2,2] = C[3,3] = type1.error
                #test for (c,1)
                C[4,1] = alpha.star[1]
                C[4,2] = alpha.star[2]+alpha.star[3]
                #test for (c,2)
                C[5,1] = alpha.star[1]
                C[5,3] = alpha.star[2]+alpha.star[3]
                #test for (1,2)
                C[6,2] = alpha.star[2]+1/2*alpha.star[1]
                C[6,3] = alpha.star[3]+1/2*alpha.star[1]
                #test for (c,1,2)
                C[7,1] = alpha.star[1]
                C[7,2] = alpha.star[2]
                C[7,3] = alpha.star[3]
                rownames(C) <- c("0","1","2","0,1","0,2","1,2","0,1,2")
                colnames(C) <- c("0","1","2")
                return(C)
        }
        #calculate rejection region
        RejRegionSingle <- function(C,j){
                C <- C
                RejRegionHFBon <- function(F){ 
                        for (i in 1:nrow(C)){
                                if(row.names(C)[i] == F){
                                        c_j <- C[i,]
                                        break;
                                }
                        }
                        alpha.current <- c_j
                        u <- numeric(3)
                        index <- as.numeric(unlist(strsplit(F,",")))
                        index <- index+1
                        for (i in 1:length(index)){
                                u[index[i]] <- qnorm(alpha.current[index[i]],mean = 0, sd = 1,lower.tail = FALSE)
                        }
                        return(u)
                }
                cut.point <- NULL #store local rejection region, 4 rows
                for (i in 1:nrow(C)){#local rejection regions
                        if(j %in% unlist(strsplit(rownames(C)[i],","))){
                                F <- rownames(C)[i]
                                cut.point[[i]] <- RejRegionHFBon(F)
                        }
                }
                #combine all rejection region together
                cut.point <- data.frame(t(matrix(unlist(cut.point),nrow = 3))) #equivalent with A vector!
                colnames(cut.point) <- c("0","1","2")
                
                ####Function:Calculate all rejection cubes of intersection of unions
                RejCube <- function(cut.point){ 
                        index <- as.numeric(j)+1
                        index.others <- setdiff(c(1,2,3),index)
                        represent.point <- NULL
                        #construct for test 1
                        a <- cut.point[,index]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- (a[2]+a[3])/2
                        temp[4] <- (a[3]+a[4])/2
                        temp[5] <- a[4]+1
                        represent.point[[index]] <- temp 
                        #construct for test 0
                        index <- index.others[1]
                        a <- cut.point[,index][cut.point[,index]!=0]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- a[2]+1
                        represent.point[[index]] <- temp
                        #constrct for test 2
                        index <- index.others[2]
                        a <- cut.point[,index][cut.point[,index]!=0]
                        a <- a[order(a)]
                        temp <- NULL
                        temp[1] <- a[1]-1
                        temp[2] <- (a[1]+a[2])/2
                        temp[3] <- a[2]+1
                        represent.point[[index]] <- temp
                        
                        IfPointInLocalRegion <- function(p,i){ #i is the row index in cut.point
                                index <- which(cut.point[i,]!=0)
                                if(any(p[index]>cut.point[i,][index]))#OR
                                        return(TRUE)
                                else return(FALSE)
                        }
                        
                        IfPointInRegion <- function(p){
                                temp <- 1
                                for (i in 1:nrow(cut.point)){
                                        temp <- temp * IfPointInLocalRegion(p,i)
                                }
                                if (temp==1) return(TRUE)
                                else return(FALSE)
                        }
                        
                        cube <- NULL
                        x <- cut.point[,1][cut.point[,1]!=0]
                        x <- x[order(x)]
                        x <- append(append(x,Inf),-Inf,after=0)
                        y <- cut.point[,2][cut.point[,2]!=0]
                        y <- append(append(y,Inf),-Inf,after=0)
                        y <- y[order(y)]
                        z <- cut.point[,3][cut.point[,3]!=0]
                        z <- z[order(z)]
                        z <- append(append(z,Inf),-Inf,after=0)
                        
                        id <- 1
                        for(i in 1:length(represent.point[[1]])){#for loop
                                for(j in 1:length(represent.point[[2]])){
                                        for(k in 1:length(represent.point[[3]])){
                                                p <- c(represent.point[[1]][i],represent.point[[2]][j],represent.point[[3]][k])
                                                if(IfPointInRegion(p)){
                                                        cube[[id]] <- rbind(c(x[i],x[i+1]),c(y[j],y[j+1]),c(z[k],z[k+1]))
                                                        id <- id+1
                                                }
                                        }
                                }
                        }#end of for
                        return(cube)
                }
                return(RejCube(cut.point))
        }
        
        C <- alpha.matrix(alpha.star)
        rej.region <- RejRegionSingle(C,j)
        rej.region
}

###FUNCTION: if z statistics fall in the rejection region
IfInRejRegion <- function(rej.region,z){
        IfinCube <- function(i){
                cube <- rej.region[[i]]
                cube1 <-cube[1,]
                cube2 <-cube[2,]
                cube3 <-cube[3,]
                if((z[1]>=cube1[1])&(z[1]<=cube1[2])&(z[2]>=cube2[1])&(z[2]<=cube2[2])&(z[3]>=cube3[1])&(z[3]<=cube3[2])) return(TRUE)
                else return(FALSE)
        } #judge if in the single cube
        for (i in 1:length(rej.region)){
                if (IfinCube(i)) return(TRUE) #if in one of the single cube, return
                else next() #if not, judge for next cube
        }
        #if not in any of the cubes, return false
        return(FALSE)
}

##example
# alpha.opt <- IfPowerSatisfy(0.4,3,3,5,251,0.05,0.8,0.0001)
# z <- c(1,2,3)
# rej.region <- RejRegion(alpha.opt,"0")
# IfInRejRegion(rej.region,z)


###simulation
###major time consuming is calculate N and alpha.opt, not simulation 

N <- design$sample.size #total
alpha.opt <- design$alpha.opt

N1 <- N*pi #total sub1
N2 <- N*(1-pi)
n1 <- ceiling(N1/2) #sample in each arm in sub 1
n2 <- ceiling(N2/2)#sample in each arm in sub 2

#fix reject region
rej.region0 <- RejRegion(alpha.opt,"0")
rej.region1 <- RejRegion(alpha.opt,"1")
rej.region2 <- RejRegion(alpha.opt,"2")

S <- 1000
power <- c(0,0,0)

for(i in 1:S){
        #alter for "0": (,1,1)
        ysub1_t <- rnorm(n1,delta1,sigma)
        ysub1_c <- rnorm(n1,0,sigma)
        ysub2_t <- rnorm(n2,delta2,sigma)
        ysub2_c <- rnorm(n2,0,sigma)
        z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
        if(IfInRejRegion(rej.region0,z)) power[1] <- power[1]+1
        
        #alter for "1":(,1,0)
        ysub1_t <- rnorm(n1,delta1,sigma)
        ysub1_c <- rnorm(n1,0,sigma)
        ysub2_t <- rnorm(n2,0,sigma)
        ysub2_c <- rnorm(n2,0,sigma)
        z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
        if(IfInRejRegion(rej.region1,z)) power[2] <- power[2]+1
        
        #alter for "2":(,0,1)
        ysub1_t <- rnorm(n1,0,sigma)
        ysub1_c <- rnorm(n1,0,sigma)
        ysub2_t <- rnorm(n2,delta2,sigma)
        ysub2_c <- rnorm(n2,0,sigma)
        z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
        if(IfInRejRegion(rej.region2,z)) power[3] <- power[3]+1
        
}

power/S

###sensitivity analysis

sens <- matrix(0,nrow = 21,ncol = 3)
count <- 0
for (p in seq(pi-0.1,pi+0.1,by=0.01)){
        count <- count+1
        N1 <- N*p #total sub1
        N2 <- N*(1-p)
        n1 <- ceiling(N1/2) #sample in each arm in sub 1
        n2 <- ceiling(N2/2)#sample in each arm in sub 2
        
        #simulation
        S <- 10000
        power <- c(0,0,0)
        
        for(i in 1:S){
                #alter for "0": (,1,1)
                ysub1_t <- rnorm(n1,delta1,sigma)
                ysub1_c <- rnorm(n1,0,sigma)
                ysub2_t <- rnorm(n2,delta2,sigma)
                ysub2_c <- rnorm(n2,0,sigma)
                z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
                if(IfInRejRegion(rej.region0,z)) power[1] <- power[1]+1
                
                #alter for "1":(,1,0)
                ysub1_t <- rnorm(n1,delta1,sigma)
                ysub1_c <- rnorm(n1,0,sigma)
                ysub2_t <- rnorm(n2,0,sigma)
                ysub2_c <- rnorm(n2,0,sigma)
                z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
                if(IfInRejRegion(rej.region1,z)) power[2] <- power[2]+1
                
                #alter for "2":(,0,1)
                ysub1_t <- rnorm(n1,0,sigma)
                ysub1_c <- rnorm(n1,0,sigma)
                ysub2_t <- rnorm(n2,delta2,sigma)
                ysub2_c <- rnorm(n2,0,sigma)
                z <- c((pi*(mean(ysub1_t)-mean(ysub1_c))+(1-pi)*(mean(ysub2_t)-mean(ysub2_c)))/sqrt(4*sigma^2/N),(mean(ysub1_t)-mean(ysub1_c))/sqrt(4*sigma^2/N1),(mean(ysub2_t)-mean(ysub2_c))/sqrt(4*sigma^2/N2))
                if(IfInRejRegion(rej.region2,z)) power[3] <- power[3]+1
                
        }
        sens[count,] <- power/S
}

sens <- data.frame(cbind(seq(pi-0.1,pi+0.1,by=0.01),sens)) %>% rename(p=X1,power_reject_0=X2,power_reject_1=X3,power_reject_2=X4)


ggplot(sens) + geom_line(aes(p,power_reject_0),color="black")
ggplot(sens) + geom_line(aes(p,power_reject_1),color="black")
ggplot(sens) + geom_line(aes(p,power_reject_2),color="black")

#plot(seq(pi-0.1,pi+0.1,by=0.01),sens[,1],xlab = "p",ylab = "power for reject 0",type="l")
#plot(seq(pi-0.1,pi+0.1,by=0.01),sens[,2],xlab = "p",ylab = "power for reject 1",type="l")
#plot(seq(pi-0.1,pi+0.1,by=0.01),sens[,3],xlab = "p",ylab = "power for reject 2",type="l")










