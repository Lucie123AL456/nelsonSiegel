### Fonctions pour le modele de Nelson-Siegel

##########################################

phi=function(x){
  y=(1-exp(-x))/x
  return(y)
}

##########################################

#Cas ou param = [mu_1,mu_2,mu_3,tau] # vecteur de dim 4
NS1=function(param,maturite){ 
  #Representation parametrique usuelle de NS (param[4]=tau)
  x=maturite/param[4]
  R=param[1]+param[2]*phi(x)+param[3]*(phi(x)-exp(-x))
  return(R)
}

##########################################

#Cas ou par = [mu_1,mu_2,mu_3] # vecteur de dim 3
NS1bis=function(par,tau,maturite){
  #Representation parametrique usuelle de NS avec taux0 fixe "a la main"
  param=c(par,tau)
  R=NS1(param,maturite)
  return(R)
}

##########################################

ajusteNelsonSiegel=function(dataTaux,fixeTau0){
  # A function written by Diethelm Wuertz / modifiee a la marge par FPL
  # Description:
  #    Fit the Yield Curve by the Nelson-Siegel Method
  # Details:
  #    This function finds a global solution. The start values for the
  #    betas are solved exactly as a function of tau using OLS.
  # Copyright:
  #    Diethelm Wuertz, (c) 2004 fBonds
  # Source:
  #    Partial copy from 'fBonds' from 'Rmetrics' (unpublished).
  
  Maturity=maturites  ##maturités des taux considerees
  Yield=dataTaux[1,] ##taux a une date donnee suivant les maturites
  
  # Find Optimal Start Solution by OLS of beta's vs. Yields:
  n = length(Maturity)
  gmin = 1.0e99
  if (is.null(fixeTau0)){tailleBoucle=n}else{tailleBoucle=1}   # fixeTau0 est fixe dans les parametres donc on n effectura que 1 boucle
  for (i in 1:tailleBoucle) {
    if (is.null(fixeTau0)){tau=Maturity[i]}else{tau=fixeTau0}
    x = Maturity/tau # Pour chaque maturites considerees on va calculer le rapport T/alpha (cf formule) que l'on note x 
    a = matrix(rep(NA, times = 9), nrow = 3)
    a[1,1] = 1
    a[1,2] = a[2,1] = mean((1-exp(-x))/x)
    a[1,3] = a[3,1] = mean((1-exp(-x))/x - exp(-x))
    a[2,2] = mean( ((1-exp(-x))/x)^2 )
    a[2,3] = a[3,2] = mean(((1-exp(-x))/x)*((1-exp(-x))/x-exp(-x)))
    a[3,3] = mean(((1-exp(-x))/x - exp(-x))^2) # a matrice pour realiser la regression
    b = c(
      mean ( Yield ),
      mean ( Yield *  ((1-exp(-x))/x)),
      mean ( Yield * (((1-exp(-x))/x - exp(-x))))) # b est le vecteur dans le systeme lineaire a*beta=b
    beta = solve(a, b) # beta sont les coefficient mu de la formule de nelson siegel
    #print("Voici les coefficients mu initiaux : ")
    #print(beta)
    yfit = beta[1] + beta[2]*phi(x) + beta[3]*(phi(x)-exp(-x))  #beta[1] + beta[2]*exp(-x) + beta[3]*x*exp(-x)
    #print(yfit)
    fmin = sum( (Yield-yfit)^2 ) # erreur entre les taux observes du fichier et les taux calcules via la formule de NS
    #print(fmin)
    if (fmin < gmin) { # si l erreur trouvee est inf a 1.e99 alors on garde c est OK
      gmin = fmin
      gvec = c(beta, tau) #concatenation entre beta qui est le vecteur des parametres estimes avec la regression et taux qui est fixeTau0
    }
  }
  #print(gvec)
  
  if (is.null(fixeTau0)){
    
    # Function to be optimized:       
    func=function(x){sum((Yield-NS1(x,Maturity))^2)}
    # Optimize:
    fit = nlminb(objective = func, start = gvec)
    #fit$start = gvec
    #names(fit$par) = c("mu1", "mu2", "mu3", "tau")
    #Recalcul des coefficient r,l et c e partir des beta(i) : l=beta(1), l-r=-beta(2), c=beta(3)
    #fit$parAdjust=c(fit$par["mu1"]+fit$par["mu2"],fit$par["mu1"],fit$par["mu3"],fit$par["tau"])
    #names(fit$parAdjust) = c("tauxCourt", "tauxLong", "convexite","tau")
    #fit$res=NS1(fit$par,Maturity)
    
  }else{
    # Function to be optimized:       
    func=function(x){sum((Yield-NS1bis(x,fixeTau0,Maturity))^2)}
    # Optimize:
    fit = nlminb(objective = func, start = beta)
    #fit$start = beta
    fit$par=c(fit$par,fixeTau0)
  }
  
  names(fit$par) = c("mu1", "mu2", "mu3","tau")
  #Recalcul des coefficient r,l et c e partir des beta(i) : l=beta(1), l-r=-beta(2), c=beta(3)
  fit$parAdjust=c(fit$par["mu1"]+fit$par["mu2"],fit$par["mu1"],fit$par["mu3"],fit$par["tau"])
  names(fit$parAdjust) = c("tauxCourt", "tauxLong", "convexite","tau")
  fit$res=NS1(fit$par,Maturity)
  
  return(fit)
}

##########################################

estimationNelsonSiegel=function(date, dataTaux){
  dataTauxuni = dataTaux[date,]
  #print(dataTauxuni)
  return(ajusteNelsonSiegel(dataTauxuni,fixeTau0)$parAdjust)
}

##########################################

