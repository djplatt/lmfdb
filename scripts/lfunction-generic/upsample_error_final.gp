/* assume nus=0 */

w(t,q,r,mu,nu,t0,H)=sqrt(q)*zeta(1.5)^r*2^((3-r)/2)*Pi^((1-r)/2)*(1+mu+nu+abs(t))^(r*(mu+1)/2)*exp(r*Pi*nu/4-(t-t0)^2/(2*H^2));

trunc_err(q,r,mu,nu,t0,H,N,B)={local(w0,w1);w0=w(t0+N/2/B,q,r,mu,nu,t0,H);w1=w(t0+(N+1)/2/B,q,r,mu,nu,t0,H);return(2*w0/(Pi*N*(1-w1/w0)));};

P(d,M,t0,H,mu,nu)={local(n);n=d*(M+mu)/2;return(2^(n-M+2)*Pi^(-M+1/2)+exp(d*Pi*nu/4)*H*(sqrt(Pi)*(M+1+mu+nu)^n+2^n*((H*sqrt(2))^n*gamma((n+1)/2)+sqrt(Pi)*t0^n)))};

Ierr(d,q,M,t0,H,mu,nu,B)=return(4*q^(M/2)*zeta(M)^d*exp(M^2/(2*H^2))*exp(-2*Pi*M*B)*P(d,M,t0,H,mu,nu));

max_r=12;
max_q=2^64;
mu=8;nu=8;
max_err=2^-200;
max_div_err=2^-35;
M=2;

num_zeros=10;

div_Ierr(d,q,M,t0,H,mu,nu,B,nz)=Ierr(d,q,M,t0,H,mu,nu,B)/(M-1/2)^nz;
div_trunc_err(q,r,mu,nu,t0,H,N,B,nz)=trunc_err(q,r,mu,nu,t0,H,N,B)/(N/(2*B))^nz;

NN=2^16; /* sampling rate in final iFFT */

for(r=1,max_r,t0=64/r;A=NN/64/8*r;B=A/2;this_H=solve(H=0.005,1,Ierr(r,max_q,M,t0,H,mu,nu,B)-max_err);this_H=round(this_H*1024*1024)/(1024*1024);this_N=solve(N=10,1000,trunc_err(max_q,r,mu,nu,t0,this_H,N,B)-max_err);this_N=ceil(this_N);printf("r: %d N: %d H: %d / %d %6e %6e\n",r,this_N,this_H*1024*1024,1024*1024,Ierr(r,max_q,M,t0,this_H,mu,nu,B)+trunc_err(max_q,r,mu,nu,t0,this_H,this_N,B),div_Ierr(r,max_q,M,t0,this_H,mu,nu,B,num_zeros)+div_trunc_err(max_q,r,mu,nu,t0,this_H,this_N,B,num_zeros)));



/*quit;*/
