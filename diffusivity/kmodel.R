## Gomi et al. (2013) model for conductivity in the core.

el.si<- list( list(z=1,  x=22.5,name='Si') )      ## Si as light element
el.c <- list( list(z=1,  x=30,  name='C') )       ## C as light element
el.o <- list( list(z=0.5,x=23.2,name='O') )       ## O as light element
el.s <- list( list(z=0.5,x=19.4,name='S') )       ## S as light element

el <- el.si          ## Change this assignment for different light element

res <- function(f, T, zx=list(), theta0=417, gamma=1.52){
   ## f - V/V0
   ## T - absolute temperature
   ## zx - vector of impurities, each a list:
   ##   z= valence charge, x= concentration (atomic fraction)

   theta <- theta0*exp(-gamma*log(f))
   rho <- function(f, f1, f2, f3){f1*(f2 - f)^f3}

   tfac <- integrate(
      function(z){z^5/((exp(z)-1)*(1-exp(-z)))},
      0, theta/T
   )$value * (T/theta)^5

   id.fe <- rho(f, 5.26e-9, 1.24, -3.21)*tfac
   id.si <- rho(f, 3.77e-8, 1.48, -3.10)
   imp <- 0; for(el in zx){
      imp <- imp + el$z^2*el$x
   }
   id.fe + id.si*imp
}

k <- function(id,T,f, L=2.44e-8){
   ## id - ideal resistivity
   ## T - absolute temperature
   ## f - V/V0

   sat <- 1.68e-6 * f^(1/3)

   (L*T)*(1/id + 1/sat)
}

## EOS functions; generic, could be replaced by more specific ones.

eos <- function(p,t,meos){
   ## Murnaghan equation of state
   ## p - pressure in Pa
   ## t - absolute temperature
   ## meos - material properties list for EOS
   rho.p <- meos$rho0*(meos$k.p * p/meos$k + 1)^(1/meos$k.p)
   al.pt <- exp(-(meos$al*exp(-meos$al.p/meos$k*p))*(t - meos$t0))
   # cat('p,t,rho.p,al.pt',p,t,rho.p,al.pt,'\n')
   rho.p*al.pt
}

eos.al <- function(p,meos){
   meos$al*exp(-meos$al.p/meos$k*p)
}

eos.tad <- function(p,t0,p0,meos){
   ## Return adiabatic lapse from t0, p0 to p
   ## t0 - absolute T at P=p0
   ## p0 - P0 pressure (Pa)
   ## p - pressure (Pa)

   t0*exp(integrate(function(p){
      eos.al(p,meos)/eos(p,t0,meos)
   }, p0, p)$value / meos$Cp)
}

# Calculate alpha such that value at 135 GPa is 1e-5, 360 GPa is 0.5e-5
rho0.met <- 7010; k.met <- 130         # Literature values
al.p <- k.met*log(2)/(360-135); al <- 1e-5/exp(-al.p*135/k.met)
meos.met <- list(
   t0=1812, rho0=rho0.met, k=k.met*1e9, k.p=4, al=al, al.p=al.p,
   Cp=800, cond=140
)

fval <- function(P, T, meos=meos.met){
   ## Return f = V/V0 for material
   ## P - pressure in GPa (note)
   ## T - absolute temperature
   eos(0,300,meos)/eos(P*1e9,T,meos)
} 

## Test to calculate variation of k along an adiabat rooted at the CMB with
##   temperature T=3750 K and P=135 GPa to ICB P=330 GPa

p <- seq(135,330,5)

t <- sapply(p, function(p)eos.tad(p*1e9, 3750, 135e9, meos.met)) ## adiabat

print(meos.met)

plot(p, t, type='l', xlab='P (GPa)', ylab='T (K)')  ## Plot T along adiabat

kvals <- sapply(1:length(p),function(i){            ## Calculate k along adiabat
   f <- fval(p[i], t[i])
   rho <- res(f, t[i], el)
   k(rho, t[i], f)
})

dev.new(); spar <- par('mar'); par(mar=0.1+c(5,5,4,2)) ## Extra margin room

plot(p, kvals, type='l',                            ## Plot k along adiabat
   xlab='P (GPa)',
   ylab=expression(k~(W~m^{-1}~K^{-1}))
)
text(p[1]+0.1*diff(range(p)), mean(kvals),
   sprintf("   %s (%.1f at.%%)",el[[1]]$name,el[[1]]$x)
)

par(spar)
