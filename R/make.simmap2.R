#### MAKESIMMAP2 FUNCTIONS ####

make.simmap2 <- function (tree, x, model="SYM", nsim=1, monitor = FALSE, rejmax = NULL, rejint = 1000000, ...){
  if (inherits(tree, "multiPhylo")) {
    ff <- function(yy, x, model, nsim, ...) {
      zz <- make.simmap(yy, x, model, nsim, ...)
      if (nsim > 1) 
        class(zz) <- NULL
      return(zz)
    }
    if (nsim > 1) 
      mtrees <- unlist(sapply(tree, ff, x, model, nsim, 
                              ..., simplify = FALSE), recursive = FALSE)
    else mtrees <- sapply(tree, ff, x, model, nsim, ..., 
                          simplify = FALSE)
    class(mtrees) <- c("multiSimmap", "multiPhylo")
  } else {
    if (hasArg(pi)) 
      pi <- list(...)$pi
    else pi <- "equal"
    if (hasArg(message)) 
      pm <- list(...)$message
    else pm <- TRUE
    if (hasArg(tol)) 
      tol <- list(...)$tol
    else tol <- 0
    if (hasArg(Q)) 
      Q <- list(...)$Q
    else Q <- "empirical"
    if (hasArg(burnin)) 
      burnin <- list(...)$burnin
    else burnin <- 1000
    if (hasArg(samplefreq)) 
      samplefreq <- list(...)$samplefreq
    else samplefreq <- 100
    if (hasArg(vQ)) 
      vQ <- list(...)$vQ
    else vQ <- 0.1
    prior <- list(alpha = 1, beta = 1, use.empirical = FALSE)
    if (hasArg(prior)) {
      pr <- list(...)$prior
      prior[names(pr)] <- pr
    }
    if (!inherits(tree, "phylo")) 
      stop("tree should be object of class \"phylo\".")
    if (!is.matrix(x)) 
      xx <- to.matrix(x, sort(unique(x)))
    else xx <- x
    xx <- xx[tree$tip.label, ]
    xx <- xx/rowSums(xx)
    tree <- bt <- reorder.phylo(tree, "cladewise")
    if (!is.binary(bt)) 
      bt <- multi2di(bt, random = FALSE)
    N <- Ntip(tree)
    m <- ncol(xx)
    root <- N + 1
    sim <- 0

    if (is.character(Q) && Q == "empirical") {
      XX <- getPars(bt, xx, model, Q = NULL, tree, tol, 
                    m, pi = pi, args = list(...))
      L <- XX$L
      Q <- XX$Q
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "empirical")
      mtrees <- replicate(nsim,
                          smap2(tree, x, N, m, root, L, Q, pi, logL, rejmax, rejint, monitor, sim),
                          simplify = FALSE)
    }
    else if (is.character(Q) && Q == "mcmc") {
      if (prior$use.empirical) {
        qq <- fitMk(bt, xx, model)$rates
        prior$alpha <- qq * prior$beta
      }
      get.stationary <- if (pi[1] == "estimated") 
        TRUE
      else FALSE
      if (pi[1] %in% c("equal", "estimated")) 
        pi <- setNames(rep(1/m, m), colnames(xx))
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      XX <- mcmcQ(bt, xx, model, tree, tol, m, burnin, 
                  samplefreq, nsim, vQ, prior, pi = pi)
      L <- lapply(XX, function(x) x$L)
      Q <- lapply(XX, function(x) x$Q)
      logL <- lapply(XX, function(x) x$loglik)
      pi <- if (get.stationary) 
        lapply(Q, statdist)
      else if (pi[1] == "fitzjohn") 
        lapply(XX, function(x) x$pi)
      else lapply(1:nsim, function(x, y) y, y = pi)
      if (pm) 
        printmessage(Reduce("+", Q)/length(Q), Reduce("+", 
                                                      pi)/length(pi), method = "mcmc")
      mtrees <- if (nsim > 1) 
        mapply(smap2, L = L, Q = Q, pi = pi, logL = logL, 
               MoreArgs = list(tree = tree, x = x, N = N, 
                               m = m, root = root, rejmax = rejmax,
                               rejint = rejint, monitor = monitor, sim = sim), SIMPLIFY = FALSE)
      else list(smap2(tree = tree, x = x, N = N, m = m, 
                      root = root, rejmax = rejmax, rejint = rejint, monitor = monitor, sim = sim,
                      L = L[[1]], Q = Q[[1]], pi = pi[[1]], logL = logL[[1]]))
    }
    else if (is.matrix(Q)) {
      XX <- getPars(bt, xx, model, Q = Q, tree, tol, m, 
                    pi = pi, args = list(...))
      L <- XX$L
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "fixed")
      
      mtrees <- replicate(nsim,
                          c(sim <<- sim + 1,
                            if(monitor == TRUE){print(paste("simulation", sim, "of", nsim, sep = " "))},
                            smap2(tree, x, N, m, root, L, Q, pi, logL, rejmax, rejint, monitor, sim)),
                          simplify = FALSE)
    }
    if (nsim==1) {
      mtrees <- mtrees[[1]]
      class(mtrees) <- c("simmap","phylo")
    } else {
      #Change class of list
      class(mtrees) <- c("multiSimmap", "multiPhylo")
      
      #Change class of individual objects
      for(i in 1:nsim){
        class(mtrees[[i]]) <- c("simmap","phylo")
      }
    }
  }
  (if (hasArg(message)) 
    list(...)$message
    else TRUE)
  if ((if (hasArg(message)) 
    list(...)$message
    else TRUE) && inherits(tree, "phylo")) 
    message("Done.")
  return(mtrees)
}

smap2 <- function(tree,x,N,m,root,L,Q,pi,logL,rejmax,rejint,monitor,sim){
  # create the map tree object
  mtree<-tree
  mtree$maps<-list()
  mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",
                                                                    tree$edge[,2],sep=""),colnames(L)))
  # now we want to simulate the node states & histories by pre-order traversal
  NN<-matrix(NA,nrow(tree$edge),2) # our node values
  NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]/
                                             sum(L[as.character(root),])) # assign root
  
  #Create list for failed tips
  fail <- list()
  fail.count <- 0
  
  for(j in 1:nrow(tree$edge)){
    # conditioned on the start value, assign end value of node (if internal)
    p<-EXPM(Q*tree$edge.length[j])[as.numeric(NN[j,1]),]*L[as.character(tree$edge[j,2]),]
    NN[j,2]<-rstate(p/sum(p))
    NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
    # now simulate on the branches
    accept <- FALSE
    counter <- 0
    while(!accept){
      map<-sch(as.numeric(NN[j,1]),tree$edge.length[j],Q)
      if (counter == rejmax){
        if (monitor == TRUE){
          print(paste("sim",sim,": branch", j, "of", nrow(tree$edge),
                      "has exceeded the rejection limit of", format(rejmax, scientific = FALSE),
                      "and will be skipped", sep = " "))
        }
        map <- NA
        fail.count <- fail.count + 1
        fail[[fail.count]] <- j
        accept = TRUE
      } else
        if(names(map)[length(map)]==NN[j,2]){
          if (monitor == TRUE){
            print(paste("sim",sim,": branch", j, "of", nrow(tree$edge),
                        "ACCEPTED after", format(counter, scientific = FALSE),
                        "total rejections", sep = " "))
          }
          accept = TRUE
        } else {
          counter <- counter + 1
          if ((counter/rejint) %% 1 == 0){
            if (monitor == TRUE){  
              print(paste("sim",sim,": branch", j, "of", nrow(tree$edge),
                          "rejected with", format(counter, scientific = FALSE),
                          "total rejections", sep = " "))
            }
          }
        }
    }
    mtree$maps[[j]]<-map
    for(k in 1:length(mtree$maps[[j]])){
      mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-
        mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
    }
  }
  for(L in 1:length(mtree$maps)){
    if(is.na(mtree$maps[[L]][1]) == TRUE){
      named.edge <- mtree$edge.length[L]
      names(named.edge) <- "fail"
      mtree$maps[[L]] <- named.edge
    }
  }
  mtree$Q<-Q
  mtree$logL<-logL
  mtree$fail <- unlist(fail)
  if(!inherits(mtree,"simmap")) class(mtree)<-c("simmap",setdiff(class(mtree),"simmap"))
  attr(mtree,"map.order")<-"right-to-left"
  return(mtree)
}

mcmcQ<-function(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior,pi,args=list()){
  update<-function(x){
    x<-abs(x+rnorm(n=np,mean=0,sd=sqrt(vQ)))
    return(x)
  }
  # get model matrix
  if(is.character(model)){
    rate<-matrix(NA,m,m)
    if(model=="ER"){
      np<-rate[]<-1
      diag(rate)<-NA
    }
    if(model=="ARD"){
      np<-m*(m-1)
      rate[col(rate)!=row(rate)]<-1:np
    }
    if(model=="SYM") {
      np<-m*(m-1)/2
      sel<-col(rate)<row(rate)
      rate[sel]<-1:np
      rate<-t(rate)
      rate[sel]<-1:np
    }
  } else {
    if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
    if(ncol(model)!=m) 
      stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
    rate<-model
    np<-max(rate)
  }
  # burn-in
  p<-rgamma(np,prior$alpha,prior$beta)
  Q<-matrix(c(0,p)[rate+1],m,m)
  diag(Q)<--rowSums(Q,na.rm=TRUE)
  yy<-getPars(bt,xx,model,Q,tree,tol,m,pi=pi,args=args)
  cat("Running MCMC burn-in. Please wait....\n")
  flush.console()
  for(i in 1:burnin){
    pp<-update(p)
    Qp<-matrix(c(0,pp)[rate+1],m,m)
    diag(Qp)<--rowSums(Qp,na.rm=TRUE)
    zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi,args)
    p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
                  yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
    if(p.odds>=runif(n=1)){
      yy<-zz
      p<-pp
    }
  }
  # now run MCMC generation, sampling at samplefreq
  cat(paste("Running",samplefreq*nsim,"generations of MCMC, sampling every",samplefreq,"generations.\nPlease wait....\n\n"))
  flush.console()
  XX<-vector("list",nsim)
  for(i in 1:(samplefreq*nsim)){
    pp<-update(p)
    Qp<-matrix(c(0,pp)[rate+1],m,m)
    diag(Qp)<--rowSums(Qp,na.rm=TRUE)
    zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi,args)
    p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
                  yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
    if(p.odds>=runif(n=1)){
      yy<-zz
      p<-pp
    }
    if(i%%samplefreq==0){
      Qi<-matrix(c(0,p)[rate+1],m,m)
      diag(Qi)<--rowSums(Qi,na.rm=TRUE)
      XX[[i/samplefreq]]<-getPars(bt,xx,model,Qi,tree,tol,m,TRUE,pi=pi,args)
    }
  }
  return(XX)
}

statdist<-function(Q){
  foo<-function(theta,Q){
    Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
    sum((Pi%*%Q)^2)
  }
  k<-nrow(Q)
  if(nrow(Q)>2){ 
    fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
    return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
  } else {
    fit<-optimize(foo,interval=c(0,1),Q=Q)
    return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
  }
}

getPars <- function(bt,xx,model,Q,tree,tol,m,liks=TRUE,pi,args=list()){
  if(!is.null(args$pi)) args$pi<-NULL
  args<-c(list(tree=bt,x=xx,model=model,fixedQ=Q,output.liks=liks,pi=pi),args)
  obj<-do.call(fitMk,args)
  N<-length(bt$tip.label)
  pi<-obj$pi
  II<-obj$index.matrix+1
  lvls<-obj$states
  if(liks){
    L<-obj$lik.anc
    rownames(L)<-N+1:nrow(L)
    if(!is.binary(tree)){
      ancNames<-matchNodes(tree,bt)
      L<-L[as.character(ancNames[,2]),]
      rownames(L)<-ancNames[,1]
    }
    L<-rbind(xx,L)
    rownames(L)[1:N]<-1:N
  } else L<-NULL	
  Q<-matrix(c(0,obj$rates)[II],m,m,dimnames=list(lvls,lvls))
  if(any(rowSums(Q,na.rm=TRUE)<tol)){
    message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",tol,"\n"))
    ii<-which(rowSums(Q,na.rm=TRUE)<tol)
    for(i in 1:length(ii)) Q[ii[i],setdiff(1:ncol(Q),ii[i])]<-tol/(ncol(Q)-1)
  }
  diag(Q)<--rowSums(Q,na.rm=TRUE)
  return(list(Q=Q,L=L,loglik=logLik(obj),pi=pi))
}

printmessage <- function(Q,pi,method){
  if(method=="empirical"||method=="fixed")
    cat("make.simmap is sampling character histories conditioned on\nthe transition matrix\n\nQ =\n")
  else if(method=="mcmc"){
    cat("make.simmap is simulating with a sample of Q from\nthe posterior distribution\n")
    cat("\nMean Q from the posterior is\nQ =\n")
  }
  print(Q)
  if(method=="empirical") cat("(estimated using likelihood);\n")
  else if(method=="fixed") cat("(specified by the user);\n")
  cat("and (mean) root node prior probabilities\npi =\n")
  if(is.list(pi)) pi<-Reduce("+",pi)/length(pi)
  print(pi)
  flush.console()
}

sch <- function(start,t,Q){
  tol<-t*1e-12
  dt<-setNames(0,start)
  while(sum(dt)<(t-tol)){
    s<-as.numeric(names(dt)[length(dt)])
    dt[length(dt)]<-if(-Q[s,s]>0) rexp(n=1,rate=-Q[s,s]) else t-sum(dt)
    if(sum(dt)<(t-tol)){
      dt<-c(dt,0)
      if(sum(Q[s,][-match(s,colnames(Q))])>0)
        names(dt)[length(dt)]<-rstate(Q[s,][-match(s,colnames(Q))]/sum(Q[s,][-match(s,colnames(Q))]))
      else names(dt)[length(dt)]<-s
    } else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
  }
  return(dt)
}

EXPM <- function(x,...){
  e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
  dimnames(e_x)<-dimnames(x)
  e_x
}

