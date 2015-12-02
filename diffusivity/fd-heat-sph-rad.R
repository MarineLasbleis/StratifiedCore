##Crank-Nicholson finite difference solution to the 1D heat equation T(x,t)
##  Testing out understanding of FD solution methods
##  This implements a solution on the spherical domain 0 <= r <= 1 with an init.
##  condition and with a flux boundary condition dT/dr = h*T(1,t>0)

if (!exists('CJsoln.sph')) source('tmodel.R')

init<-function(x){rep(1,length(x))}           ## Initial condition T(x,0) = 1

eps <- 0.0                                    ## Factor controlling adiabat
                                              ## alpha^2*g/(3 Cp)

h <- 0.5

true <- function(x,t){
   coef <- CJcoeff.sph(200,1,h,init(0))       ## True solution: coefficients
   CJsoln.sph(x,t,1,coef)
}

tridagsoln<-function(a,b,c,r)
{	n<-length(b)
	gam<-u<-vector(length=n)
	bet<-b[1]
	u[1]<-r[1]/bet
	for (j in 2:n) {
		gam[j]<-c[j-1]/bet
		bet<-b[j]-a[j]*gam[j]
		u[j]<-(r[j]-a[j]*u[j-1])/bet
	}
	for (j in (n-1):1) u[j]<-u[j]- gam[j+1]*u[j+1]
	u
}

xmin<-0
xmax<-1
N.xsteps<-150
dx<-(xmax-xmin)/N.xsteps
xgrid<-dx*0:N.xsteps + xmin
Ngrid <- length(xgrid)
Neq <- Ngrid                ## Number of equations is number of grid pts

taumax<-0.5
N.timesteps<-20
dtau<-taumax/N.timesteps
taugrid<-(taumax/N.timesteps)*0:N.timesteps

r<-dtau/dx^2

w<-matrix(0,nrow=Ngrid,ncol=(N.timesteps+1))
u<-vector(length=Ngrid)

ai <- -r*(1:Neq - 1)/1:Neq
bi <- rep(2 + 2*r,Neq)
ci <- -r*(1:Neq + 1)/1:Neq

di <- vector(length=Neq)

ai[1] <- 0; bi[1] <- 2+6*r; ci[1] <- -6*r ## central point
ai[Neq] <- -2*r; bi[Neq] <- 2+2*r*(1 + dx*h*(Neq+1)/Neq); ci[Neq] <- 0

u <- init(xgrid)

for (j in 1:N.timesteps){
   di[1] <- (2-6*r)*u[1] + 6*r*u[2]       ## central point
   di[Neq] <- (2*r*u[Neq-1] + (2 - 2*r*(1 + dx*h*(Neq+1)/Neq))*u[Neq]) * (
             1 / (1 + eps*Neq*dx*u[Neq])
   )


   ix <- 1 + 1:(Neq-2)

   di[ix] <- (r*(ix-1)/ix*u[ix-1] + (2-2*r)*u[ix] + r*(ix+1)/ix*u[ix+1]) * (
             1 / (1 + eps*ix*dx*u[ix])
   )
     

   w[,j+1] <- u <- tridagsoln(ai,bi,ci,di)
}

plot(xgrid,true(xgrid,taumax),type='l',
   xlab="r",ylab=bquote(f(r, t=.(taumax))),
   ylim=c(0,1),lwd=2
)
lines(xgrid,w[,(N.timesteps+1)], lty=2)
lines(xgrid,init(xgrid), lty=3)
for(j in 2:N.timesteps)lines(xgrid, w[,j], lty=2)

## Plot spherical solution too
scof <- CJcoeff.sph(200,1,h,init(0))
for(tau in taugrid)lines(
   xgrid, CJsoln.sph(xgrid,tau,1,scof), col='blue'
)
lines(xgrid, CJsoln.sph(xgrid,taumax,1,scof), col='red')

lines(c(0,0.2),rep(0.1,2),lwd=2)
   text(0.2,0.1,'Analytic spherical',adj=c(-0.1,0.5))
lines(c(0,0.2),rep(0.05,2),lty=2); text(0.2,0.05,'Finite diff.',adj=c(-0.1,0.5))

print(sqrt(sum((w[,(N.timesteps+1)]-true(xgrid,taumax))^2)*dx))
