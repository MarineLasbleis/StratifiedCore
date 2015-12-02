## Thermal model of evolution of core temperature.  Consists of a parameterized
## scheme in the mantle side of the CMB that sets the heat flux, coupled with
## an analytic solution to the heat equation in the core with the specified
## flux boundary condition at the surface.

flux <- function(rc, rb, Tc, Tb, Ta, ka, Cpa, rhoa, ga, alpha,
   E=530e3/8.3144, beta=1/3, Rac=1100, eta0=8.8e10){
   ## Returns inverse of boundary layer thickness (km)
   ## rc, rb - core, base of magma ocean radii (km)
   ## Tc, Tb, Ta - temperature (K) at CMB, MOB, bottom of adiabat at CMB (km)
   ## ka, Cpa, rhoa, alpha - thermal conductivity, heat capacity, density,
   ##   expansivity of solid silicate
   ## ga - grav. acc. at CMB
   kappa <- ka/(rhoa*Cpa)
   Ra <- alpha*rhoa*ga*(rb-rc)^3*1e9*(Tc-Tb)/(kappa*eta0*exp(2*E/(Tb+Ta)))
   ka*(Tc-Tb)/(rb-rc)*(Ra/Rac)^beta*1e-3 ## flux
}

CJcoeff.lin0l <- function(m,l,h,k,f){
   ## For solution on linear interval 0 <= x <= l, radiation boundary cond.
   ##   Carslaw & Jaeger (1959) ss 3.9
   stopifnot(m>0)
   h1 <- h[1]; h2 <- h[2]; k1 <- k[1]; k2 <- k[2]
   beta <- Zi <- Zn <- rep(0,m)
   for(n in 1:m){
      beta[n] <- uniroot(
	 function(b){(k1*k2*b^2 - h1*h2)*sin(b*l) - b*(k1*h1+k2*h1)*cos(b*l)},
	 pi/2*c(2*(n-1),2*n-1)
      )$root
      b <- beta[n]; b2 <- beta[n]^2
      c <- sqrt(
	 2*(k2^2*b2) /
	   ( (k1^2*b2 + h1^2)*(l*(k2^2*b2 + h2^2) + k2*h2) +
	      k1*h1*(k2^2*b2 + h2^2)
	   )
      )
      Zn[n] <- c
      Zi[n] <- integrate(
	 function(x){f(x) * c*(k1*b*cos(b*x) + h1*sin(b*x))},
	 0, l
      )$value
   }
   list(h=h1, k=k1, beta=beta, Zn=Zn, Zi=Zi)
}

CJcoeff.lin <- function(m,l,h,f,lim=c(0,1)){
   ## For solution on linear interval -l <= x <= l, radiation boundary cond.
   ## f(x) is an even function of x (initial temperaure profile)
   ##   Carslaw & Jaeger (1959) ss 3.10
   stopifnot(m>0)
   allok <- 'OK'
   a <- Z <- rep(0,m)
   for(n in 1:m){
      a[n] <- uniroot(
	 function(b){b*sin(b*l) - h*cos(b*l)},
	 pi*c(n-1,n - 10*.Machine$double.eps)
      )$root
      cc <- h^2 + a[n]^2
      res <- integrate(
	 function(x){f(x) * cos(a[n]*x)},
	 l*lim[1], l*lim[2],
	 subdivisions=200,
	 stop.on.error=FALSE
      )
      if(res$message != 'OK'){
         allok <- res; break
      }
      Z[n] <- res$value * 2*cc/(cc*l + h)
   }
   list(h=h, a=a, Z=Z, OK=allok)
}

CJcoeff.sph.flux <- function(m,a,h,f,lim=c(0,1)){
   ## Solution in spherical geometry, 0 <= r <= a
   ## From Caretto notes eqs 132, 129.  h is h/k here; k is conductivity.
   stopifnot(m>0)
   allok <- 'OK'
   ah <- a*h; al <- Z <- rep(0,m)
   if (1-ah > 0) {
      bnd <- c(-2,-1)+ah
   } else if (1-ah < 0) {
      bnd <- c(-1, 0) ## - 10*.Machine$double.eps
   } else {
      bnd <- c(-1.1,-0.9)
   }
   for(n in 1:m){
      al[n] <- uniroot(
	 function(al){a*al*cos(a*al) - (1 - ah)*sin(a*al)},
	 pi/2 * (2*n+bnd) / a
      )$root
      res <- integrate( function(r){r*f(r) * sin(al[n]*r)},
	 a*lim[1], a*lim[2],
         stop.on.error=FALSE, subdivisions=1000)
      if (res$message != 'OK'){
         #cat(paste('n',n,'value',res$value,'abs.error',res$abs.error,
	 #          'subdivisions',res$subdivisions),'\n')
         # print(res); stopifnot(res$message == 'OK')
         allok <- res; break
      }
      Zi <- res$value
      Z[n] <- Zi / (a/2 - (sin(al[n]*a)*cos(al[n]*a))/(2*al[n]))
   }
   if (abs(al[1]) < 10*.Machine$double.eps) cn[1] <- 0
   list(a=a, al=al, Z=Z, OK=allok)
}

CJcoeff.sph <- function(m,a,h,c0){
   ## Solution in spherical geometry, 0 <= r <= a
   ## Crank 6.3.4; c0 is constant initial concentration, h is h/D
   stopifnot(m>0)
   L <- a*h; bt <- rep(0,m)
   for(n in 1:m){
      if (1-L > 0) {
	 bnd <- c(-2,-1)+L
      } else if (1-L < 0) {
	 bnd <- c(-1, 0) ## - 10*.Machine$double.eps
      } else {
	 bnd <- c(-1.1,-0.9)
      }
      bt[n] <- uniroot(
	 function(al){al*cos(al) - (1 - L)*sin(al)},
	 pi/2 * (2*n+bnd)
      )$root
   }
   list(a=a, bt=bt, L=L, c0=c0)
}

CJsoln.lin <- function(x.all,t,kap,coef){
   with(coef, {
      sapply(x.all,function(x){
         sum(exp(-kap*t*a^2) * cos(a*x) * Z)
      })
   })
}

CJsoln.lin0l <- function(x,t,kap,coef){
   ## Solution on interval 0 <= x <= 1
   Zn <- coef$Zn; Zi <- coef$Zi; beta <- coef$beta
   sum(
      Zi*Zn*(coef$k*beta*cos(beta*x) + coef$h*sin(beta*x)) * exp(-kap*t*beta^2)
   )
}

CJsoln.sph.flux <- function(r,t,kap,coef){
   with(coef, sapply(r, function(r){
      if (r <= 0) s <- al else s <- sin(al*r)/r
      sum(s * Z * exp(-kap*t*al^2))
   }))
}

CJsoln.sph <- function(r,t,kap,coef){
   with(coef, {
      c0*2*L*a*sapply(r, function(r){
         if (r>0) s <- sin(bt*r/a)/(r*sin(bt)) else s <- bt/(a*sin(bt))
	 sum(s * exp(-kap*t*(bt/a)^2) / (bt^2 + L*(L-1)))
      })
   })
}
