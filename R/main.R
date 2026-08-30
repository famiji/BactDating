#' Main function of BactDating package to date the nodes of a phylogenetic tree
#' @param tree Tree wih branches measured in unit of substitutions
#' @param date Sampling dates for the leaves of the tree
#' @param initMu Initial rate of substitutions per genome (not per site), or zero to use root-to-tip estimate
#' @param initAlpha Initial coalescent time unit
#' @param initSigma Initial std on per-branch substitution rate (only used in relaxed or mixed-relaxed models)
#' @param updateMu Whether or not to update the substitution rate
#' @param updateAlpha Whether or not to update the coalescent time unit
#' @param updateSigma Whether or not to per-branch substitution rate (only used in relaxed or mixed-relaxed models)
#' @param updateRoot Whether or not to use the root finding algorithm
#' @param nbIts Number of MCMC iterations to perform
#' @param thin Thining interval between recorded MCMC samples
#' @param useCoalPrior Whether or not to use a coalescent prior on the tree
#' @param model Which model to use (poisson or negbin or strictgamma or relaxedgamma or mixedgamma or arc or carc or mixedcarc)
#' @param useRec Whether or not to use results from previous recombination analysis
#' @param minbralen Minimum branch length of the phylogenetic tree (in number of substitutions)
#' @param showProgress Whether or not to show a progress bar
#' @param tuning Whether or not to use tuning
#' @param fast Whether or not to use the optimized implementation (C++ node-date sweep, incremental coalescent prior, localized root moves; statistically identical to the default loop)
#' @param exact Only used when fast=TRUE: force the optimized implementation to stay bit-identical to the default loop, at the cost of the incremental prior and localized root moves (slower than fast=TRUE,exact=FALSE but still faster than fast=FALSE)
#' @return Dating results
#' @export
bactdate = function(tree, date, initMu = NA, initAlpha = NA, initSigma = NA, updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T, nbIts = 10000, thin=ceiling(nbIts/1000), useCoalPrior = T,  model = 'arc', useRec = F, minbralen = 0.1, showProgress = F, tuning = T, fast = T, exact = F)
{
  if (sum(tree$edge.length)<5) warning('Warning: input tree has small branch lengths. Make sure branch lengths are in number of substitutions (NOT per site).\n')
  #Rerranging of dates, if needed
  if (is.matrix(date)) {
    rangedate=date
    if (!is.null(rownames(rangedate))) rangedate=cbind(findDates(tree,rangedate[,1]),findDates(tree,rangedate[,2]))
    date=(rangedate[,1]+rangedate[,2])/2
  } else rangedate=NULL
  if (is.vector(date) && !is.null(names(date))) date=findDates(tree,date)
  if (is.null(rangedate)) rangedate=cbind(date,date)

  #Initial rooting of tree, if needed
  if (useRec==T && is.null(tree$unrec)) stop("To use recombination, the proportion of unrecombined needs to be input.")
  if (length(which(tree$edge[,1]==Ntip(tree)+1))==3) tree=unroot(tree)
  if (is.rooted(tree)==F) tree=initRoot(tree,date,useRec=useRec)

  #If the initial rate was not specified, start with the rate implied by root-to-tip analysis
  if (is.na(initMu)) {
    r=suppressWarnings(unname(roottotip(tree,date,showFig=F)$rate))
    if (is.na(r) || r<0) r=1
    initMu=r
  }
  mu=initMu
  if (is.na(initSigma)) initSigma=mu
  sigma=initSigma

  #Select prior function
  prior=function(...) return(0)
  if (useCoalPrior) prior=coalpriorC else updateAlpha=F

  #Selection likelihood function
  #if (!is.element(model,c('poisson','poissonR','negbin','negbinR','arc','arcR'))) tree$edge.length=pmax(tree$edge.length,minbralen)
  if (is.element(model,c('strictgamma','strictgammaR','poisson','poissonR'))) {updateSigma=F;sigma=0}
  if (model == 'mixedgamma' || model == 'mixedcarc') sigma=0
  if (model == 'poisson') likelihood=function(tab,mu,sigma) return(likelihoodPoissonC(tab,mu))
  if (model == 'poissonR') likelihood=function(tab,mu,sigma) return(likelihoodPoisson(tab,mu))
  if (model == 'negbin') likelihood=function(tab,mu,sigma) return(likelihoodNegbinC(tab,mu,sigma))
  if (model == 'negbinR') likelihood=function(tab,mu,sigma) return(likelihoodNegbin(tab,mu,sigma))
  if (model == 'strictgamma') likelihood=function(tab,mu,sigma) return(likelihoodGammaC(tab,mu))
  if (model == 'strictgammaR') likelihood=function(tab,mu,sigma) return(likelihoodGamma(tab,mu))
  if (model == 'relaxedgamma'||model == 'mixedgamma') likelihood=function(tab,mu,sigma) return(likelihoodRelaxedgammaC(tab,mu,sigma))
  if (model == 'relaxedgammaR') likelihood=function(tab,mu,sigma) return(likelihoodRelaxedgamma(tab,mu,sigma))
  if (model == 'arc') likelihood=function(tab,mu,sigma) return(likelihoodArcC(tab,mu,sigma))
  if (model == 'arcR') likelihood=function(tab,mu,sigma) return(likelihoodArc(tab,mu,sigma))
  if (model == 'carc'||model == 'mixedcarc') likelihood=function(tab,mu,sigma) return(likelihoodCarcC(tab,mu,sigma))
  if (model == 'carcR') likelihood=function(tab,mu,sigma) return(likelihoodCarc(tab,mu,sigma))
  if (model == 'null') {updateMu=0;likelihood=function(tab,mu,sigma) return(0)}
  if (!exists('likelihood')) stop('Unknown model.')

  #Map model to an integer code for the optimized C++ internal-date update
  modelcode=0
  if (model=='poisson') modelcode=1
  if (model=='negbin') modelcode=2
  if (model=='strictgamma') modelcode=3
  if (model=='relaxedgamma'||model=='mixedgamma') modelcode=4
  if (model=='arc') modelcode=5
  if (model=='carc'||model=='mixedcarc') modelcode=6

  #Deal with missing dates
  rangedate[which(is.na(rangedate[,1])),1]=min(rangedate[,1],na.rm = T)
  rangedate[which(is.na(rangedate[,2])),2]=max(rangedate[,2],na.rm = T)
  misDates=which(rangedate[,1]<rangedate[,2])
  date[misDates]=(rangedate[misDates,1]+rangedate[misDates,2])/2
  orderedleafdates=sort(date,decreasing = T)

  #Create table of nodes (col1=name,col2=observed substitutions on branch above,col3=date,col4=father,col5=unrecombined proportion, only if useRec=T)
  n = Ntip(tree)
  tab = matrix(NA, n + tree$Nnode, ifelse(useRec,5,4))
  nrowtab = nrow(tab)
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    tab[i, 1] = i
    tab[i, 4] = tree$edge[r, 1]
    tab[i, 2] = tree$edge.length[r]
    if (useRec) tab[i,5]=tree$unrec[r]
    if (i <= n)
      tab[i, 3] = date[i]
  }
  tab[n+1,2]=minbralen#store minbralength above root

  #Create initial tree
  tree = reorder.phylo(tree, 'postorder')
  for (r in 1:(nrow(tree$edge))) {
    i = tree$edge[r, 2]
    if (i <= n)
      next
    children = which(tab[, 4] == i)
    if (length(children) != 0)
      tab[i, 3] = min(tab[children,3]-0.01,mean(tab[children, 3]-tab[children,2]/mu))
  }
  i = n + 1
  children = which(tab[, 4] == i)
  tab[i, 1] = i
  tab[i, 3] = min(tab[children, 3]-tab[children,2]/mu)

  #Sample Alpha if no initial point provided
  if (is.na(initAlpha)) {
    if (fast) {
      su=coalAlphaSumC(tab[1:n,3],tab[(n+1):nrowtab,3])
    } else {
      s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
      k=cumsum(2*(s$ix<=n)-1)
      difs=s$x[1:(nrowtab-1)]-s$x[2:nrowtab]
      su=sum(k[1:(nrowtab-1)]*(k[1:(nrowtab-1)]-1)*difs)
    }
    initAlpha=1/rgamma(1,shape=n+0.001-1,scale=2000/(su*1000+2))
  }
  alpha=initAlpha
  initHeight=max(tab[,3])-min(tab[,3])

  #Proposal distributions and tuning
  sdMu=0.1*ifelse(initMu==0,1e-3,initMu)
  sdSigma=0.1*ifelse(initSigma==0,1e-3,initSigma)
  sdDates=0.05*ifelse(initHeight==0,1,initHeight)
  acceptance_target=0.234
  delta=1/acceptance_target/(1-acceptance_target)

  #Start of MCMC algorithm
  l = likelihood(tab, mu, sigma)
  orderednodedates=sort(tab[(n+1):nrow(tab),3],method='quick',decreasing = T)
  p = prior(orderedleafdates,orderednodedates,alpha)
  record = matrix(NA, floor(nbIts / thin), nrow(tab)*3 + 6)
  colnames(record)<-c(rep(NA,nrow(tab)*3),'likelihood','mu','sigma','alpha','prior','root')
  if (showProgress) pb <- utils::txtProgressBar(min=0,max=nbIts,style = 3)
  children=vector("list", max(tab[,4],na.rm = T))
  for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
  curroot=findCurrootC(tree$edge,children[[n+1]],n+1)

  if (fast && modelcode>0 && !showProgress) {
    #The entire MCMC loop runs in C++ (mcmcLoopC): mode 1 keeps the trajectory
    #bit-identical to the default loop, mode 2 accumulates mathematically
    #identical increments; RNG consumption matches the R loop draw by draw
    res=mcmcLoopC(tab,tree$edge,n,nbIts,thin,orderedleafdates,orderednodedates,misDates,rangedate,mu,sigma,alpha,l,p,sdMu,sdSigma,sdDates,updateMu,updateSigma,updateAlpha,updateRoot==T,updateRoot=='branch',model=='mixedgamma'||model=='mixedcarc',tuning,useCoalPrior,modelcode,useRec,minbralen,ifelse(exact,1L,2L),record)
    tab=res$tab;l=res$l;p=res$p;mu=res$mu;sigma=res$sigma;alpha=res$alpha;record=res$record
  } else {
  for (i in 1:nbIts) {
    #Record
    if (i %% thin == 0) {
      if (showProgress) utils::setTxtProgressBar(pb, i)
      record[i / thin, 1:nrow(tab)] = tab[, 3]
      record[i / thin, (1:nrow(tab))+nrow(tab)]=tab[, 4]
      record[i / thin, (1:nrow(tab))+2*nrow(tab)]=tab[, 2]
      record[i / thin, 'likelihood'] = l
      record[i / thin, 'mu'] = mu
      record[i / thin, 'sigma'] = sigma
      record[i / thin, 'alpha'] = alpha
      record[i / thin, 'prior'] = p
      record[i / thin, 'root'] = curroot
    }

    if (updateMu == T) {
      #MH move using Gamma(1e-3,1e3) prior
      mu2=abs(rnorm(1,mu,sdMu))
      l2=likelihood(tab,mu2,sigma)
      mh=l2-l+dgamma(mu2,shape=1e-3,scale=1e3,log=T)-dgamma(mu,shape=1e-3,scale=1e3,log=T)
      if (log(runif(1))<mh) {l=l2;mu=mu2}
      if (tuning) sdMu=sdMu*exp(delta/(i+1)*(min(1,exp(mh))-acceptance_target))
    }

    if (updateSigma == T && sigma>0) {
      #MH move using Gamma(1e-3,1e3) prior
      sigma2=abs(rnorm(1,sigma,sdSigma))
      l2=likelihood(tab,mu,sigma2)
      mh=l2-l+dgamma(sigma2,shape=1e-3,scale=1e3,log=T)-dgamma(sigma,shape=1e-3,scale=1e3,log=T)
      if (log(runif(1))<mh) {l=l2;sigma=sigma2}
      if (tuning) sdSigma=sdSigma*exp(delta/(i+1)*(min(1,exp(mh))-acceptance_target))
    }

    if (model == 'mixedgamma' || model == 'mixedcarc') {
      #Reversible-jump move
      if (sigma==0) {
        sigma2=rexp(1)
        qratio=sigma2#-dexp(sigma2,1,log=T)
        pratio=dgamma(sigma2,shape=1e-3,scale=1e3,log=T)
      } else {
        sigma2=0
        qratio=-sigma#dexp(sigma,1,log=T)
        pratio=-dgamma(sigma,shape=1e-3,scale=1e3,log=T)
      }
      l2=likelihood(tab,mu,sigma2)
      if (log(runif(1))<l2-l+pratio+qratio) {l=l2;sigma=sigma2}
    }

    if (updateAlpha == T) {
      #Gibbs move using inverse-gamma prior
      if (fast) {
        su=coalAlphaSumC(orderedleafdates,orderednodedates)
      } else {
        s <- sort.int(tab[, 3], method='quick',decreasing = T, index.return = TRUE)
        k=cumsum(2*(s$ix<=n)-1)
        difs=s$x[1:(nrowtab-1)]-s$x[2:nrowtab]
        su=sum(k[1:(nrowtab-1)]*(k[1:(nrowtab-1)]-1)*difs)
      }
      alpha=1/rgamma(1,shape=n+0.001-1,scale=2000/(su*1000+2))
      p=prior(orderedleafdates,orderednodedates,alpha)
    }

    #MH to update internal dates
    rn=rnorm(nrow(tab)-n,0,sdDates)
    if (fast && modelcode>0) {
      upd=nodeDatesUpdateC(tab,n,orderedleafdates,orderednodedates,rn,mu,sigma,l,p,alpha,minbralen,sdDates,i,tuning,useCoalPrior,modelcode,useRec,!exact)
      l=upd$l;p=upd$p;sdDates=upd$sdDates
    } else for (j in c((n + 1):nrow(tab))) {
      r=rn[j-n]
      old=tab[j,3]
      new=old+r
      if (r<0&&(!is.na(tab[j,4])&&new<tab[tab[j,4],3])) {mh=-Inf} else#can't be older than father
      if (r>0&&new>min(tab[children[[j]],3])) {mh=-Inf} else {#can't be younger than sons
      if (j>(n+1)) {mintab=rbind(tab[children[[j]],],tab[tab[j,4],],tab[j,]);mintab[,4]=c(4,4,NA,3)}
      else {mintab=rbind(tab[children[[j]],],tab[j,]);mintab[,4]=c(3,3,NA)}
      mintab[3,2]=minbralen
      l2=l-likelihood(mintab,mu,sigma)
      mintab[nrow(mintab),3]=new
      l2=l2+likelihood(mintab,mu,sigma)
      changeinorderedvec(orderednodedates,old,new)
      p2 = prior(orderedleafdates, orderednodedates, alpha)
      mh = l2 - l + p2 - p
      if (log(runif(1)) < mh)
      {l = l2; p = p2;tab[j,3]=new}
      else
        changeinorderedvec(orderednodedates,new,old)
      }
      if (tuning) {
        ii=(i-1)*(nrow(tab)-n)+j
        sdDates=sdDates*exp(delta/(ii+1)*(min(1,exp(mh))-acceptance_target))
      }
    }

    #MH to update missing leaf dates
    for (j in misDates) {
      old = tab[j, 3]
      new = rnorm(1,old,(rangedate[j,2]-rangedate[j,1])*0.05)
      if (new-old<0&&(!is.na(tab[j,4])&&new<tab[tab[j,4],3])) next#can't be older than father
      if (new>rangedate[j,2]||new<rangedate[j,1]) next#stay within prior range
      mintab=rbind(tab[j,],tab[tab[j,4],])
      mintab[,4]=c(2,NA)
      mintab[2,2]=minbralen
      l2=l-likelihood(mintab,mu,sigma)
      mintab[1,3]=new
      l2=l2+likelihood(mintab,mu,sigma)
      changeinorderedvec(orderedleafdates,old,new)
      p2 = prior(orderedleafdates, orderednodedates, alpha)
      if (log(runif(1)) < l2 - l + p2 - p)
      {l = l2; p = p2;tab[j,3]=new}
      else
        changeinorderedvec(orderedleafdates,new,old)
    }

    if (updateRoot == T || updateRoot == 'branch') {
      #Move root on current branch
      if (fast && modelcode>0 && !exact) {
        #The root is always node n+1; only the two side branches' substitution
        #counts change, so only their two likelihood terms need updating
        sides=children[[n+1]]
        old=tab[sides,2]
        oldlocal=localTermsC(tab,sides,mu,sigma,minbralen,modelcode,useRec)
        r=runif(1)
        tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
        l2=l-oldlocal+localTermsC(tab,sides,mu,sigma,minbralen,modelcode,useRec)
        if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
      } else if (fast && modelcode>0) {
        #exact mode: full likelihood recompute (bit-identical), but the root
        #is always node n+1 so which() is avoided
        sides=children[[n+1]]
        old=tab[sides,2]
        r=runif(1)
        tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
        l2=likelihood(tab,mu, sigma)
        if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
      } else {
        root=which(is.na(tab[,4]))
        sides=which(tab[,4]==root)
        old=tab[sides,2]
        r=runif(1)
        tab[sides,2]=c(sum(old)*r,sum(old)*(1-r))
        l2=likelihood(tab,mu, sigma)
        if (log(runif(1))<l2-l) l=l2 else tab[sides,2]=old
      }
  }

    if (updateRoot == T) {
      #Move root branch
      if (fast && modelcode>0 && !exact) {
        #Only the terms of nodes a, left and right are affected by this move
        sides=children[[n+1]]
        if (tab[sides[1],3]<tab[sides[2],3]) {left=sides[1];right=sides[2]} else {left=sides[2];right=sides[1]}
        if (left>n) {
          ab=children[[left]]
          if (runif(1)<0.5) {a=ab[1];b=ab[2]} else {a=ab[2];b=ab[1]}
          aff=sort(c(a,left,right))
          oldlocal=localTermsC(tab,aff,mu,sigma,minbralen,modelcode,useRec)
          olda2=tab[a,2];oldleft2=tab[left,2];oldright2=tab[right,2]
          olda4=tab[a,4];oldright4=tab[right,4]
          if (useRec) oldleft5=tab[left,5]
          tab[a,4]=n+1
          tab[right,4]=left
          r=runif(1)
          tab[a,2]=olda2*r
          tab[left,2]=olda2*(1-r)
          tab[right,2]=oldright2+oldleft2
          if (useRec) tab[left,5]=tab[a,5]
          l2=l-oldlocal+localTermsC(tab,aff,mu,sigma,minbralen,modelcode,useRec)
          if (log(runif(1))<l2-l+log(ifelse(olda2==0,1,olda2)/ifelse(tab[right,2]==0,1,tab[right,2])))
            {l=l2
            children=vector("list", max(tab[,4],na.rm = T))
            for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
            curroot=findCurrootC(tree$edge,children[[n+1]],n+1)
          } else {
            tab[a,4]=olda4;tab[right,4]=oldright4
            tab[a,2]=olda2;tab[left,2]=oldleft2;tab[right,2]=oldright2
            if (useRec) tab[left,5]=oldleft5
          }
        }
      } else if (fast && modelcode>0) {
      #Move root branch (exact mode): full likelihood recompute (bit-identical),
      #but which() scans are replaced by the children structure
      sides=children[[n+1]]
      if (tab[sides[1],3]<tab[sides[2],3]) {left=sides[1];right=sides[2]} else {left=sides[2];right=sides[1]}
      if (left>n) {
        oldtab=tab
        ab=children[[left]]
        if (runif(1)<0.5) {a=ab[1];b=ab[2]} else {a=ab[2];b=ab[1]}
        tab[a,4]=n+1
        tab[right,4]=left
        r=runif(1)
        tab[a,2]=oldtab[a,2]*r
        tab[left,2]=oldtab[a,2]*(1-r)
        tab[right,2]=oldtab[right,2]+oldtab[left,2]
        if (useRec) tab[left,5]=tab[a,5]
        l2=likelihood(tab,mu,sigma)
        if (log(runif(1))<l2-l+log(ifelse(oldtab[a,2]==0,1,oldtab[a,2])/ifelse(tab[right,2]==0,1,tab[right,2])))
          {l=l2
          children=vector("list", max(tab[,4],na.rm = T))
          for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
          curroot=findCurrootC(tree$edge,children[[n+1]],n+1)
        } else tab=oldtab
      }
      } else {
      #Move root branch
      root=which(is.na(tab[,4]))
      sides=which(tab[,4]==root)
      if (tab[sides[1],3]<tab[sides[2],3]) {left=sides[1];right=sides[2]} else {left=sides[2];right=sides[1]}
      if (left>n) {
        oldtab=tab
        ab=which(tab[,4]==left)
        if (runif(1)<0.5) {a=ab[1];b=ab[2]} else {a=ab[2];b=ab[1]}
        tab[a,4]=root
        tab[right,4]=left
        r=runif(1)
        tab[a,2]=oldtab[a,2]*r
        tab[left,2]=oldtab[a,2]*(1-r)
        tab[right,2]=oldtab[right,2]+oldtab[left,2]
        if (useRec) tab[left,5]=tab[a,5]
        l2=likelihood(tab,mu,sigma)
        if (log(runif(1))<l2-l+log(ifelse(oldtab[a,2]==0,1,oldtab[a,2])/ifelse(tab[right,2]==0,1,tab[right,2])))
          {l=l2
          children=vector("list", max(tab[,4],na.rm = T))
          for (i in 1:nrow(tab)) if (!is.na(tab[i,4])) children[[tab[i,4]]]=c(children[[tab[i,4]]],i)
          curroot=findCurrootC(tree$edge,children[[n+1]],n+1)
        } else tab=oldtab
      }
      }
    }

  }#End of MCMC loop
  }
  if (showProgress) close(pb)

  #Output
  inputtree = tree
  bestroot = as.numeric(names(sort(table(record[floor(nrow(record) / 2):nrow(record),'root']),decreasing=T)[1]))
  bestrows = intersect(floor(nrow(record) / 2):nrow(record),which(record[,'root']==bestroot))
  meanRec = colMeans(record[bestrows, ,drop=F])
  for (i in 1:nrow(tree$edge)) {
    tree$edge[i,1]=record[bestrows[1],nrow(tab)+tree$edge[i,2]]
    tree$edge.length[i] = meanRec[tree$edge[i, 2]] - meanRec[tree$edge[i, 1]]
    tree$subs[i] = meanRec[nrow(tab)*2+tree$edge[i,2]]
  }

  #Calculate DIC
  meantab=tab
  meantab[,3]=meanRec[1:nrow(tab)]
  meantab[,4]=meanRec[(1:nrow(tab))+nrow(tab)]
  dic=-2*likelihood(meantab,meanRec['mu'],meanRec['sigma'])+var(-2*record[bestrows,'likelihood'])

  tree$root.time = unname(meanRec[n + 1])#max(date)-max(leafDates(tree))
  #Remove pre-existing order attribute as the tree is no longer ordered in any particular way
  attributes(tree)$order <- NULL
  attributes(inputtree)$order <- NULL
  CI = matrix(NA, nrow(tab), 2)
  for (i in 1:nrow(tab)) {
    s=sort(record[bestrows,i])
    CI[i,1]=s[ceiling(length(s)*0.025)]
    CI[i,2]=s[max(1,floor(length(s)*0.975))]
  }
  out = list(
    inputtree = inputtree,
    tree = tree,
    model = model,
    record = record,
    rootdate = unname(meanRec[n + 1]),
    rootprob = length(bestrows)/length(floor(nrow(record) / 2):nrow(record)),
    CI = CI,
    dic = dic
  )
  class(out) <- 'resBactDating'
  if (model=='mixedgamma'||model=='mixedcarc') out$pstrict=length(which(record[,'sigma']==0))/nrow(record)
  return(out)
}

#' Perform Bayesian comparison between two models based on DIC values
#' @param res1 An output from the bactdate function
#' @param res2 Another output from the bactdate function
#' @return Prints out results of model comparison
#' @export
modelcompare = function(res1,res2) {
  dic1=res1$dic
  dic2=res2$dic
  dif=dic2-dic1
  cat(sprintf('The first model has DIC=%.2f and the second model has DIC=%.2f.\n',dic1,dic2))
  if (dif>10) cat('Model 1 is definitely better.\n')
  if (dif>5 && dif<10) cat('Model 1 is slightly better.\n')
  if (abs(dif)<5) cat('The difference is not significant.\n')
  if (dif< -5 && dif>-10) cat('Model 2 is slightly better.\n')
  if (dif< -10) cat('Model 2 is definitely better.\n')
}
