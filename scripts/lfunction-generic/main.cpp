// generic l-function calculator


// modform function files look like
// <N> <wt> <char> <ver> <file name>
// elliptic curve/ genus 2 curves
// <N> <smalljac spec>

#include "slint.h" // thanks Bober/NTL


// set up the L_comp_t instance
// called read_comp coz I expected to read these parameters in
bool read_comp(L_comp_t *L, FILE *cfile, uint64_t prec)
{
  // only called once, so use init/clear not statics
  arb_t tmp;
  arb_init(tmp);

  L->N=1<<11; // size of convolution FFT
#ifdef HI_PREC
  L->NN=1<<20; // size of output FFT
#else
  L->NN=1<<16;
#endif
  printf("Set FFT length to %lu. Not checking whether this is sufficient. Maybe need FFT Length > 2*(# G values used).\n",L->N);

  L->G=(acb_t *)malloc(sizeof(acb_t)*L->N);
  for(uint64_t n=0;n<L->N;n++)
    acb_init(L->G[n]);

  L->eta=0.0;
  arb_init(L->delta);
  arb_const_pi(L->delta,prec);
  arb_mul_2exp_si(L->delta,L->delta,-1); // pi/2
  arb_set_d(tmp,1.0-L->eta);
  arb_mul(L->delta,L->delta,tmp,prec); // (1-eta)pi/2

  L->w=(acb_t *)malloc(sizeof(acb_t)*L->N/2); // twiddles for little FFT
  if(!L->w)
    {
      printf("Error allocating memory for w. Exiting.\n");
      return(0);
    }

  for(uint64_t i=0;i<L->N/2;i++)
    acb_init(L->w[i]);
  acb_initfft(L->w,L->N,prec); // set twiddles for little FFT

  printf("Setting iFFT length to %lu\n",L->NN);
  fflush(stdout);

  L->ww=(acb_t *)malloc(sizeof(acb_t)*L->NN/2); // twiddles for big FFT
  if(!L->ww)
    {
      printf("Error allocating memory for ww. Exiting.\n");
      return(0);
    }

  for(uint64_t i=0;i<L->NN/2;i++)
    acb_init(L->ww[i]);
  acb_initfft(L->ww,L->NN,prec); // set twiddles for big FFT

  // space for the zeros once we isolate them
  L->zeros[0]=(arb_t *)malloc(sizeof(arb_t)*MAX_ZEROS);
  L->zeros[1]=(arb_t *)malloc(sizeof(arb_t)*MAX_ZEROS);
  for(uint64_t z=0;z<MAX_ZEROS;z++)
    {
      arb_init(L->zeros[0][z]);
      arb_init(L->zeros[1][z]);
    }

  arb_init(L->zero_prec);
  arb_zero(L->zero_prec);
  arb_set_ui(tmp,1);
  arb_mul_2exp_si(tmp,tmp,-L->OP_ACC-1);
  arb_add_error(L->zero_prec,tmp); // +/- 2^-{OP_ACC+1}
  printf("Zero error set to ");arb_printd(L->zero_prec,20);printf("\n");
  arb_clear(tmp);
  return(true);
} /* read_comp */



// each line consists of <i> <k> <m1> <e1> <m2> <e2>
// so Gs[i][k]=m1*2^e1 +/- m2*2^e2
bool read_Gs(L_comp_t *Lc, L_family_t *Lf, FILE *infile, uint64_t prec)
{

  int64_t fi,i,j;uint64_t fk;

  for(i=Lf->low_i,j=0;i<=Lf->hi_i;i++,j++)
    for(uint64_t k=0;k<Lf->max_K;k++)
      {
	if(fscanf(infile,"%ld %lu",&fi,&fk)!=2)
	  return(false);
	if(fi!=i)
	  return(false);
	if(fk!=k)
	  return(false);
	if(!read_arb(Lf->Gs[k][j],infile))
	  return(false);
      }
  return(true);
}

uint64_t set_X(uint64_t r,double *mus,double one_over_B)
{
  double max_mu=mus[0];
  for(uint64_t j=1;j<r;j++)
    if(mus[j]>max_mu)
      max_mu=mus[j];
  printf("max mu_j = %10.8e\n",max_mu);
  max_mu+=2.5;
  double max_t=1.0/(OUTPUT_RATIO*one_over_B);
  printf("t0 = %10.8e\n",max_t);
  double X2=max_t*max_t-max_mu*max_mu;
  if(X2<=25.0)
    {
      printf("Error setting X (eq 4-10). A mu_j was too large. Exiting.\n");
      return(0);
    }
  return(sqrt(X2));
}

bool read_family(L_family_t *L_family, L_comp_t *Lc, L_error_t *Le, FILE *ffile, uint64_t prec)
{
  arb_t tmp;
  arb_init(tmp);

  L_family->r=read_r(ffile);
  if(!L_family->r)
    return(false);
  printf("r=%lu\n",L_family->r);
  L_family->mus=(double *) malloc(sizeof(double)*L_family->r);
  if(!L_family->mus)
    {
      printf("Error allocating memory for mus in read_family. Exiting.\n");
      exit(0);
    }
  L_family->nus=(arb_t *) malloc(sizeof(arb_t)*L_family->r);
  if(!L_family->nus)
    {
      printf("Error allocating memory for nus in read_family. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=0;i<L_family->r;i++)
    {
      arb_init(L_family->nus[i]);
    }

  for(uint64_t i=0;i<L_family->r;i++)
    if(fscanf(ffile,"%lf",L_family->mus+i)!=1)
      {
	printf("Error reading mu from family file. Exiting.\n");
	//printf("mu[%lu] was ",i);acb_printd(L_family->mus[i],10);printf("\n");
	exit(0);
      }
    else
      printf("mu[%lu]=%10.8e\n",i,L_family->mus[i]);



  // Lemma 5.2 defines mu
  acb_init(L_family->mu);
  acb_zero(L_family->mu);
  for(uint64_t i=0;i<L_family->r;i++)
    {
      arb_set_d(tmp,L_family->mus[i]);
      acb_add_arb(L_family->mu,L_family->mu,tmp,prec);
    }
  acb_div_ui(L_family->mu,L_family->mu,L_family->r,prec);
  arb_set_d(tmp,0.5);
  arb_sub(acb_realref(L_family->mu),acb_realref(L_family->mu),tmp,prec);

  // lemma 3 defines nu_j and nu
  for(uint64_t i=0;i<2;i++)
    {
      arb_set_d(tmp,L_family->mus[i]-0.5);
      arb_mul_2exp_si(L_family->nus[i],tmp,-1);
    }
  for(uint64_t i=2;i<L_family->r;i++)
    {
      arb_set_d(tmp,L_family->mus[i]-1.0);
      arb_mul_2exp_si(L_family->nus[i],tmp,-1);
    }
  arb_init(L_family->nu);
  arb_set(L_family->nu,L_family->nus[0]);
  for(uint64_t i=1;i<L_family->r;i++)
    arb_add(L_family->nu,L_family->nu,L_family->nus[i],prec);
  arb_div_ui(L_family->nu,L_family->nu,L_family->r,prec);
  arb_mul_2exp_si(L_family->nu,L_family->nu,1);
  arb_set_d(tmp,0.5);
  arb_add(L_family->nu,L_family->nu,tmp,prec);

  if(!read_C_alpha_line(ffile))
    {
      printf("Error reading C,alpha line. Exiting.\n");
      exit(0);
    }
  printf("Assuming Ramanujan bound.\n");
  //printf("C=%10.8e alpha=%10.8e\n",L_family->C,L_family->alpha);

  // alpha must be >=1/r so use 1/r
  arb_init(L_family->alpha);
  arb_set_d(L_family->alpha,1.0);
  arb_div_ui(L_family->alpha,L_family->alpha,L_family->r,prec);

  // c defined in Lemma 4  
  arb_init(L_family->c);
  arb_add(L_family->c,L_family->nu,L_family->alpha,prec);
  arb_set_d(tmp,0.5);
  arb_add(L_family->c,L_family->c,tmp,prec);

  L_family->one_over_B=read_gap(ffile);
  if(L_family->one_over_B==0.0)
    {
      printf("Error reading 1/B. Exiting.\n");
      exit(0);
    }
  arb_init(L_family->B);
  arb_set_d(L_family->B,L_family->one_over_B);
  arb_inv(L_family->B,L_family->B,prec);
  printf("B=");arb_printd(L_family->B,10);printf("\n");
  arb_init(L_family->two_pi_by_B);
  arb_set_d(L_family->two_pi_by_B,L_family->one_over_B*2.0);
  arb_const_pi(tmp,prec);
  arb_mul(L_family->two_pi_by_B,L_family->two_pi_by_B,tmp,prec);

  printf("File has 2pi/B=");arb_printd(L_family->two_pi_by_B,10);printf("\n");

  L_family->X=set_X(L_family->r,L_family->mus,L_family->one_over_B);
  if(L_family->X==0)
    {
      printf("Error setting X in read_family, Exiting.\n");
      exit(0);
    }
  printf("X for eq 4-10 set to %lu\n",L_family->X);

  arb_init(L_family->imint);
  arb_zero(L_family->imint);
  arb_t sig,a,b,im_tmp;
  arb_init(sig);
  arb_init(a);
  arb_init(b);
  arb_init(im_tmp);
  arb_set(sig,L_family->B);
  arb_div_ui(a,sig,OUTPUT_RATIO,prec);
  arb_div_ui(im_tmp,sig,TURING_RATIO,prec);
  arb_add(b,a,im_tmp,prec);
  arb_mul_2exp_si(a,a,-1); // /2 for change of variable in integral
  arb_mul_2exp_si(b,b,-1);

  
  if(L_family->r<=MAX_TURING_R)
    {
      for(uint64_t j=0;j<L_family->r;j++)
	{
	  arb_set_d(sig,0.25+L_family->mus[j]/2.0);
	  
	  if(!im_int(im_tmp,sig,a,b,prec))
	    {
	      printf("Error computing Int of Im loggamma. Exiting.\n");
	      exit(0);
	    }
	  arb_add(L_family->imint,L_family->imint,im_tmp,prec);
	}
      arb_mul_2exp_si(L_family->imint,L_family->imint,2); // x2 for conj x2 for change of var
      printf("Im int = ");arb_printd(L_family->imint,10);printf("\n");
    }
  else
    printf("Degree too high (%lu). Not using Turing's method.\n",L_family->r);
  
  arb_clear(sig);
  arb_clear(a);
  arb_clear(b);
  arb_clear(im_tmp);


  L_family->low_i=read_low_i(ffile);
  if(L_family->low_i==BAD_64)
    {
      printf("Error reading low_i. Exiting.\n");
      exit(0);
    }
  L_family->hi_i=read_hi_i(ffile);
  if(L_family->hi_i==BAD_64)
    {
      printf("Error reading hi_i. Exiting.\n");
      exit(0);
    }
  L_family->G_len=L_family->hi_i-L_family->low_i+1;

  printf("We have Gs from %ld to %lu inclusive.\n",L_family->low_i,L_family->hi_i);

  L_family->max_K=read_max_K(ffile);
  if(!L_family->max_K)
    {
      printf("Error reading max_K. Exiting.\n");
      exit(0);
    }

  Lc->K=L_family->max_K;
  printf("Number of differentials in file K = %lu\n",L_family->max_K);
  
  //Lc->sks=(arb_t *)malloc(sizeof(arb_t)*MAX_M);
  //for(uint64_t m=0;m<MAX_M;m++)
  //arb_init(Lc->sks[m]);
  
  arb_init(Le->eq59);
  if(!read_eq59(Le->eq59,ffile))
    {
      printf("Error reading equation 5-9 error. Exiting.\n");
      exit(0);
    }

  printf("Error for equation 5-9 (Taylor truncation) = ");
  arb_printd(Le->eq59,10);
  printf("\n");


  L_family_t *Lf=L_family;



  Lc->A=Lc->NN*Lf->one_over_B;
  printf("A=%10.8e\n",Lc->A);
  arb_init(Lc->arb_A);
  arb_set_d(Lc->arb_A,Lc->A);
  arf_init(Lc->arf_A);
  arf_set_d(Lc->arf_A,Lc->A);

  arb_init(Lc->one_over_A);
  arb_inv(Lc->one_over_A,Lc->arb_A,prec);
  printf("1/A = ");arb_printd(Lc->one_over_A,20);printf("\n");
  arf_init(Lc->arf_one_over_A);
  arf_ui_div(Lc->arf_one_over_A,1,Lc->arf_A,prec,ARF_RND_NEAR);
  // allocate vector of vectors to hold G^(k)(i/B)
  Lf->Gs=(arb_t **)malloc(sizeof(arb_t *)*Lc->K);
  for(uint64_t k=0;k<Lc->K;k++)
    {
      Lf->Gs[k]=(arb_t *)malloc(sizeof(arb_t)*Lf->G_len);
      for(uint64_t n=0;n<Lf->G_len;n++)
	arb_init(Lf->Gs[k][n]);
    }

  // now we read the G^(k)(2 pi i/B)
  // read them all for now, we'll select later
  if(!read_Gs(Lc,Lf,ffile,prec))
    {
      printf("Error reading Gs. Exiting.\n");
      exit(0);
    }

  Lc->kres=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->skm=(acb_t **)malloc(sizeof(acb_t *)*Lc->K);
  for(uint64_t k=0;k<Lc->K;k++)
    {
      Lc->skm[k]=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
      for(uint64_t n=0;n<Lc->N;n++)
	acb_init(Lc->skm[k][n]);
    }

  Lc->res=(acb_t *)malloc(sizeof(acb_t)*Lc->NN);
  for(uint64_t n=0;n<Lc->N;n++)
    acb_init(Lc->kres[n]);
  for(uint64_t n=0;n<Lc->NN;n++)
    acb_init(Lc->res[n]);

  arb_init(L_family->ftwiddle_error);
  if(!init_ftwiddle_error(L_family,Lc,prec))
    {
      printf("Error initialising ftwiddle error. Exiting.\n");
      exit(0);
    }


  arb_clear(tmp);
  return(true);
} /* read_family */

bool init_func(L_func_t *L)
{
  arb_init(L->one_over_root_N);
  arb_init(L->sum_ans);
  acb_init(L->epsilon);
  acb_init(L->epsilon_sqr);
  arb_init(L->ftwiddle_error);
  L->ans=(acb_t *)malloc(sizeof(acb_t)*8192);
  if(!L->ans)
    {
      printf("Error allocating memory for a_n's. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=0;i<8192;i++)
    acb_init(L->ans[i]);
  L->allocated_M=8192;

  return(true);
}

int64_t calc_m(uint64_t i, double two_pi_by_B, double dc)
{
  double x=log((double) i/dc)/two_pi_by_B;
  return(round(x));
}

// sks=log(m/sqrt{N})-u_m
void comp_sks(arb_t sks, uint64_t m, int64_t ms, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
  arb_mul_si(tmp1,Lf->two_pi_by_B,ms,prec);
  arb_mul_ui(tmp2,Lfu->one_over_root_N,m+1,prec);
  arb_log(tmp3,tmp2,prec);
  arb_sub(sks,tmp3,tmp1,prec);
}

// read Dirichlet coefficients, a_n
// normalise them if required
// then construct skm[k][m] vectors for convolutions
bool read_func(L_func_t *L, L_comp_t *Lc, L_family_t *Lf, double normalisation, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,n,sks,normal;
  static acb_t coeff,ctmp;
  if(!init)
    {
      init=true;
      acb_init(coeff);acb_init(ctmp);
      arb_init(tmp);arb_init(tmp1);arb_init(n);arb_init(sks);
      arb_init(normal);
    }

  arb_set_d(normal,-normalisation-0.5);
  arb_sqrt_ui(tmp,L->N,prec);
  arb_inv(L->one_over_root_N,tmp,prec);
  complete_ftwiddle_error(Lf,L,prec);
  printf("Final Ftwiddle Error set to ");arb_printd(L->ftwiddle_error,10);printf("\n");

  L->dc=sqrt((double) L->N);
  Lc->M0=ceil(L->dc/100);
  printf("M0 set to %lu.\n",Lc->M0);
  L->M=L->dc*exp(2*M_PI*(Lf->hi_i+0.5)*Lf->one_over_B);
  printf("M computed from hi_i = %lu\n",L->M);

  /*
  // try to reduce M a bit.
  arb_zero(tmp);
  arb_set_d(tmp1,MAX_M_ERROR);
  uint64_t old_M;
  while(true)
    {
      old_M=L->M;
      L->M/=1.05;
      M_error(n,tmp,L->M,Lf,L,prec); // > 0
      arb_sub(tmp1,tmp1,n,prec);
      if(!arb_is_positive(tmp1))
	{
	  L->M=old_M;break;
	}
    }
  printf("M now set to %lu\n",L->M);
    */
  if(L->M>L->allocated_M)
    {
      printf("Need more space for Dirichlet coefficients.\n");
      for(uint64_t i=0;i<L->allocated_M;i++)
	acb_clear(L->ans[i]);
      while(L->allocated_M<L->M)
	{
	  L->allocated_M<<=1;
	  if(L->allocated_M==0) // L->M was huge!
	    L->allocated_M=L->M;
	}
      L->ans=(acb_t *)realloc(L->ans,sizeof(acb_t)*L->allocated_M);
      if(!L->ans)
	{
	  printf("Attempt to (re-)allocate memory for Dirichlet coefficients failed. Exiting.\n");
	  exit(0);
	}
      for(uint64_t i=0;i<L->allocated_M;i++)
	acb_init(L->ans[i]);
      printf("re-allocated enough memory for Dirichlet coefficients.\n");
    }
  fflush(stdout);

  // call the routine to produce the Euler Product polynomials
  // and put the Dirichlet co-efficients into L->ans
  printf("Getting coefficients ");fflush(stdout);system("date");
  if(!do_lpolys(L,Lf,Lc,prec))
    {
      printf("Error getting Dirichlet coefficients. Skipping.\n");
      return(false);
    }
  printf("Got our coefficients ");fflush(stdout);system("date");
  // when we get here, the unnormalised Lc->M dirichlet coefficients are in L->ans[0]..[M-1]
  /*
  for(uint64_t m=0;m<32;m++)
    {
      printf("ans[%lu]=",m+1);acb_printd(L->ans[m],20);printf("\n");
    }
  */
  // use the first M0 of them 
  arb_zero(L->sum_ans);
  for(uint64_t m=0;m<Lc->M0-1;m++)
    {
      arb_set_ui(n,m+1);
      arb_pow(tmp1,n,normal,prec); // n^(-normalisation-1/2)
      acb_mul_arb(L->ans[m],L->ans[m],tmp1,prec);
      acb_abs(tmp1,L->ans[m],prec);
      arb_add(L->sum_ans,L->sum_ans,tmp1,prec);
    }

  for(uint64_t k=0;k<Lc->K;k++)
    for(uint64_t n=0;n<Lc->N;n++)
      acb_zero(Lc->skm[k][n]);

  double two_pi_by_B=2.0*M_PI*Lf->one_over_B;
  uint64_t one_pc=L->M/100,pc=(100*(Lc->M0-1))/L->M;
  for(uint64_t m=Lc->M0-1;m<L->M;m++)
    {
      int64_t ms=calc_m(m+1,two_pi_by_B,L->dc);
      int64_t b=(-ms)%Lc->N;
      arb_set_ui(n,m+1);
      arb_pow(tmp1,n,normal,prec); // n^(-normalisation-1/2)
      acb_mul_arb(coeff,L->ans[m],tmp1,prec);
      acb_abs(tmp1,coeff,prec);
      arb_add(L->sum_ans,L->sum_ans,tmp1,prec);
      acb_add(Lc->skm[0][b],Lc->skm[0][b],coeff,prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^0
      comp_sks(sks,m,ms,Lc,Lf,L,prec);
      //printf("m=%lu sks=",m);arb_printd(sks,20);printf("\n");
      acb_mul_arb(ctmp,coeff,sks,prec);
      acb_add(Lc->skm[1][b],Lc->skm[1][b],ctmp,prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^1
      arb_set(tmp1,sks);
      for(uint64_t k=2;k<Lc->K;k++)
	{
	  arb_mul(tmp1,tmp1,sks,prec);		
	  acb_mul_arb(ctmp,coeff,tmp1,prec);
	  acb_add(Lc->skm[k][b],Lc->skm[k][b],ctmp,prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^k
	}
    }

  printf("sum |an/sqrt(n)|=");arb_printd(L->sum_ans,10);printf(" ");fflush(stdout);system("date");
  return(true);
} /* read_func */


// just check our G values go down far enough
bool finalise_comp(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, FILE *ffile, uint64_t prec)
{
  double two_pi_by_B=Lf->one_over_B*2*M_PI;
  Lc->offset=calc_m(1,two_pi_by_B,Lfu->dc);
  if(Lc->offset<Lf->low_i)
    {
      printf("G values need to go down to %ld. We only have down to %ld. Skipping.\n",
	     Lc->offset,Lf->low_i);
      return(false);
    }

  return(true);
} /* finalise_comp */


// do the k convolutions, summing results in res.
bool do_convolves(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static acb_t tmp;
  static arb_t tmp1;
  if(!init)
    {
      init=true;
      acb_init(tmp);
      arb_init(tmp1);
    }
  //int64_t n0=round(log(Lc->M0/Lfu->dc)/(Lf->one_over_B*2*M_PI));
  int64_t n0=Lf->low_i;
  printf("Taking G values from %ld to %ld\n",n0,Lf->hi_i);

  // let's do k=0
  //printf("Doing k=%lu\n",0);fflush(stdout);
  for(uint64_t n=0;n<Lc->N;n++)
    acb_zero(Lc->G[n]);
  /*
  for(uint64_t m=Lc->M0-1;m<Lfu->M;m++)
    {
      uint64_t b=bucket(m,Lc,Lf);			
      acb_add(Lc->skm[b],Lc->skm[b],Lfu->ans[m],prec);
      //printf("Put coeff a_%lu (",m+1);acb_printd(Lfu->ans[m],20);printf(") into bucket %lu\n",b);

    }
  */
  // just copy those G we actually need
  // i.e. from hi_i down to u_m=round(log(1/sqrt{conductor})*B/2/Pi)
  // forget that, copy them all
  for(int64_t n=n0,n2=n0-Lf->low_i;n<=Lf->hi_i;n++,n2++)
    {
      int64_t n1=n%Lc->N;
      //printf("copying G value ");arb_printd(Lf->Gs[0][n2],20);printf("\n");
      arb_set(acb_realref(Lc->G[n1]),Lf->Gs[0][n2]);
      //arb_zero(acb_imagref(Lc->G[n1]));
    }
  
  
  acb_convolve(Lc->res,Lc->skm[0],Lc->G,Lc->N,Lc->w,prec);
  //printf("res[0] after convolve = ");acb_printd(Lc->res[0],10);printf("\n");
  //printf("res[-1]=");acb_printd(Lc->res[Lc->N-1],10);printf("\n");
  
  for(uint64_t k=1;k<Lc->K;k++)
    {
      //printf("Doing k=%lu\n",k);fflush(stdout);
      for(uint64_t n=0;n<Lc->N;n++)
	{
	  //acb_zero(Lc->skm[n]);
	  acb_zero(Lc->G[n]); 
	}
      for(int64_t n=n0,n2=n0-Lf->low_i;n<=Lf->hi_i;n++,n2++)
	{
	  int64_t n1=n%Lc->N;
	  arb_set(acb_realref(Lc->G[n1]),Lf->Gs[k][n2]);
	}
      acb_convolve(Lc->kres,Lc->skm[k],Lc->G,Lc->N,Lc->w,prec);
      
      for(int64_t n=0;n<=Lc->N/2;n++)
	acb_add(Lc->res[n],Lc->res[n],Lc->kres[n],prec);    
      // we need [N-1] to compute epsilon when F_hat(0)=0  
      acb_add(Lc->res[Lc->N-1],Lc->res[Lc->N-1],Lc->kres[Lc->N-1],prec);
      //printf("res[0] after convolve = ");acb_printd(Lc->res[0],10);printf("\n");
      //printf("res[-1]=");acb_printd(Lc->res[Lc->N-1],10);printf("\n");
    }

  return(true);
} /* do_convolves */

// handle the coefficients from m=1 to M0-1
// that weren't convolved
bool finish_convolves(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static acb_t tmp,tmp2;
  static arb_t tmp1,sks;
  if(!init)
    {
      init=true;
      acb_init(tmp);acb_init(tmp2);
      arb_init(tmp1);arb_init(sks);
    }
  if(Lc->M0==1) return(true);

  printf("Finishing off convolutions...\n");
  printf("F_hat[-1]=");acb_printd(Lc->res[Lc->N-1],30);printf("\n");

  double two_pi_by_B=2.0*M_PI*Lf->one_over_B;  
  for(uint64_t m=0;m<Lc->M0-1;m++)
    {
      //printf("Doing m=%lu\n",m);
      int64_t ms=calc_m(m+1,two_pi_by_B,Lfu->dc);
      for(int64_t n=-1;;n++)
	{
	  int64_t nn=ms+n;
	  if(nn>Lf->hi_i) // run out of G values
	    break;
	  acb_mul_arb(tmp2,Lfu->ans[m],Lf->Gs[0][nn-Lf->low_i],prec);
	  //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
	  acb_add(Lc->res[n%Lc->N],Lc->res[n%Lc->N],tmp2,prec);
	}
      comp_sks(sks,m,ms,Lc,Lf,Lfu,prec);
      for(int64_t n=-1;;n++)
	{
	  int64_t nn=ms+n;
	  if(nn>Lf->hi_i) // run out of G values
	    break;
	  acb_mul_arb(tmp,Lfu->ans[m],sks,prec);
	  acb_mul_arb(tmp2,tmp,Lf->Gs[1][nn-Lf->low_i],prec);
	  //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
	  acb_add(Lc->res[n%Lc->N],Lc->res[n%Lc->N],tmp2,prec);
	}
      arb_set(tmp1,sks);
      for(uint64_t k=2;k<Lc->K;k++)
	{
	  arb_mul(tmp1,tmp1,sks,prec);
	  acb_mul_arb(tmp,Lfu->ans[m],tmp1,prec); // an/sqrt(n)(log(m/sqrt(N))-um)^k
	  for(int64_t n=-1;;n++)
	    {
	      int64_t nn=ms+n;
	      if(nn>Lf->hi_i) // run out of G values
		break;
	      acb_mul_arb(tmp2,tmp,Lf->Gs[k][nn-Lf->low_i],prec);
	      //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
	      acb_add(Lc->res[n%Lc->N],Lc->res[n%Lc->N],tmp2,prec);
	    }
	}
    }

  printf("Convolutions finished.\n");

  return(true);
} /* finish_convolves */

// do the final iFFT (with implicit upsample
// multiply output by 2Pi/B
// add the ftwiddle error
bool final_ifft(L_comp_t *Lc, L_func_t *Lfu, L_family_t *Lf, L_error_t *Le, uint64_t prec)
{
  printf("Going to iFFT.\n");fflush(stdout);
  acb_ifft(Lc->res,Lc->NN,Lc->ww,prec);
  printf("iFFT done.\n");fflush(stdout);
  for(uint64_t n=0;n<Lc->NN;n++) // only need to go to end of Turing + upsample region
    {
      arb_mul(acb_realref(Lc->res[n]),acb_realref(Lc->res[n]),Lf->two_pi_by_B,prec);
      arb_add_error(acb_realref(Lc->res[n]),Lfu->ftwiddle_error);
    }
  return(true);
} /* final_ifft */

// normalise both sides by exp(pit/4)
// copy relevant portions to values
// need output, turing +/- the upsampling width
// ensure that Lambda(0)>0 or Lambda(delta)>0
void pit4(L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Lu, L_func_t *Lfu, uint64_t prec)
{

  static bool init=false;
  static arb_t g,t;
  if(!init)
    {
      init=true;
      arb_init(g);
      arb_init(t);  
    }

  bool negate_me=false;
  if(arb_is_negative(acb_realref(Lc->res[0])))
    negate_me=true;
  else
    {
      if(arb_contains_zero(acb_realref(Lc->res[0]))
	 &&(arb_is_negative(acb_realref(Lc->res[1]))))
	negate_me=true;
    }

  int64_t nn=-Lu->N*Lu->stride*2; // this is where we start from
  arb_mul_si(t,Lc->one_over_A,nn,prec);
  for(uint64_t n=0;n<Lu->no_values;n++,nn++,arb_add(t,t,Lc->one_over_A,prec))
    {
      exp_pi_t_4(g,t,prec);
      arb_mul(Lu->values[0][n],acb_realref(Lc->res[nn%Lc->NN]),g,prec);
      arb_mul(Lu->values[1][n],acb_realref(Lc->res[(-nn)%Lc->NN]),g,prec);
    }
  if(negate_me)
    for(uint64_t n=0;n<Lu->no_values;n++)
      {
	arb_neg(Lu->values[0][n],Lu->values[0][n]);
	arb_neg(Lu->values[1][n],Lu->values[1][n]);
      }
}

// a version of exp that takes more care
// exploits fact that exp is monotonic
void safe_exp(arb_t res, arb_t x, int64_t prec)
{
  static bool init=false;
  static fmpz_t l,r,e;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init=true;
      fmpz_init(l);
      fmpz_init(r);
      fmpz_init(e);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  arb_get_interval_fmpz_2exp(l,r,e,x);
  arb_set_fmpz_2exp(tmp,l,e);
  arb_exp(tmp1,tmp,prec);
  arb_set_fmpz_2exp(tmp,r,e);
  arb_exp(tmp2,tmp,prec);
  arb_union(res,tmp1,tmp2,prec);
}


void simple_sinc(arb_t res, arb_t pi_x, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_sin(tmp,pi_x,prec);
  arb_div(res,tmp,pi_x,prec);
}

// sin(pi x)/(pi x)
void safe_sinc(arb_t res, arb_t x, uint64_t prec)
{
  static arb_t pi,tmp1,tmp,tmp3,tmp2,sinc_err;
  static bool init=false;
  static int64_t facts[5]={-6,120,-5040,362880,-3991680};
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp);
      arb_init(sinc_err);
      arb_set_d(sinc_err,1.6e-10);
    }
  arb_mul(tmp,x,pi,prec);
  //printf("in my_sinc_interval with x = ");arb_printd(x,20);printf(" sin(pi*x) = ");arb_printd(sin_x,20);printf("\n");
  arb_sub_ui(tmp1,tmp,1,prec);
  if(!arb_is_negative(tmp1))
    {
      simple_sinc(res,tmp,prec);
      return;
    }
  arb_add_ui(tmp1,tmp,1,prec);
  if(!arb_is_positive(tmp1))
    {
      simple_sinc(res,tmp,prec);
      return;
    }
  arb_set_ui(res,1);
  arb_sqr(tmp1,tmp,prec);
  arb_div_si(tmp2,tmp1,facts[0],prec);
  arb_add(res,res,tmp2,prec); // 1-x^2/6
  arb_sqr(tmp,tmp1,prec); // x^4
  arb_div_si(tmp2,tmp,facts[1],prec);
  arb_add(res,res,tmp2,prec); // 1-x^2/6+x^4/120
  arb_mul(tmp2,tmp1,tmp,prec); // x^6
  arb_div_si(tmp3,tmp2,facts[2],prec);
  arb_add(res,res,tmp3,prec); // 1-x^2/6+x^4/120-x^6/7!
  arb_sqr(tmp2,tmp,prec); // x^8
  arb_div_si(tmp3,tmp2,facts[3],prec);
  arb_add(res,res,tmp3,prec); // 1-x^2/6+x^4/120-x^6/7!+x^8/9!
  arb_mul(tmp,tmp2,tmp1,prec); // x^10
  arb_div_si(tmp3,tmp,facts[4],prec);
  arb_add(res,res,tmp,prec); // 1-x^2/6+x^4/120-x^6/7!+x^8/9!-x^10/11!
  arb_mul(tmp2,tmp,tmp1,prec); // x^12
  arb_mul(tmp,tmp2,sinc_err,prec);
  arb_add_error(res,tmp);
  //printf("new sinc of ");arb_printd(x,20);printf(" returning ");arb_printd(res,20);printf("\n");
}

// sinc(x)=sin(pi*x)/(pi*x)
// expsinc=exp(-dt^2/(2H^2)*sinc(A*dt)
void safe_expsinc(arb_t res, arb_t A, arb_t inv_2H2, arb_t dt, int64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,tmp2,tmp3,tmp4,pi_A;

  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(pi_A);
      arb_const_pi(tmp1,prec);
      arb_mul(pi_A,tmp1,A,prec);
    }
  //printf("delta = ");arb_printd(dt,20);printf("\n");
  arb_sqr(tmp,dt,prec);
  //printf("delta^2 = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp1,tmp,inv_2H2,prec);
  //printf("-delta^2/2H^2 = ");arb_printd(tmp1,10);printf("\n");
  safe_exp(tmp,tmp1,prec);
  //printf("exp term = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp2,dt,A,prec);
  safe_sinc(tmp1,tmp2,prec);
  //printf("sinc term = ");arb_printd(tmp1,10);printf("\n");
  arb_mul(res,tmp1,tmp,prec);
  //printf("exp sinc=");arb_printd(res,20);printf("\n");
}

// sinc(x)=sin(pi*x)/(pi*x)
// expsinc=exp(-dt^2/(2H^2)*sinc(A*dt)
// handles dt contains 0
void expsinc(arb_t res, arb_t A, arb_t inv_2H2, arb_t dt, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,tmp4,pi_A;

  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(pi_A);
      arb_const_pi(tmp1,prec);
      arb_mul(pi_A,tmp1,A,prec);
    }
  arb_mul(tmp1,dt,dt,prec); // ^2
  arb_mul(tmp2,tmp1,inv_2H2,prec); // -dt^2/(2H^2)
  arb_exp(tmp1,tmp2,prec); // exp bit
  arb_mul(tmp2,dt,pi_A,prec);
  arb_sin(tmp3,tmp2,prec);
  arb_div(tmp4,tmp3,tmp2,prec);
  arb_mul(res,tmp4,tmp1,prec);
  //printf("expsinc(");arb_printd(dt,20);printf(")=");arb_printd(res,20);printf("\n");
}

#define OLD_DIV
void setup_expsincs(L_comp_t *Lc, L_upsample_t *Lu, int64_t prec)
{
#ifndef NEW_DIV // only need these data for NEW_DIV
  return;
#endif
  arb_t err,delta,tmp;
  arb_init(err);arb_init(delta);arb_init(tmp);
  arb_set(err,Lc->one_over_A);
  arb_mul_2exp_si(err,err,-LOG_UPS-1);
  printf("************Check that computing expsincs to NN/4 is sufficient.\n");
  for(uint64_t n=1;n<DIV_UPSAMPLE_RATE;n++)
    {
      Lu->upsample_expsincs[n-1]=(acb_t *)malloc(sizeof(acb_t)*Lc->NN);
      if(!Lu->upsample_expsincs[n-1])
	{
	  printf("Failed to allocate memory for upsample_expsincs[%lu]. Exiting.\n",n);
	  exit(0);
	}
      for(uint64_t i=0;i<Lc->NN;i++)
	acb_init(Lu->upsample_expsincs[n-1][i]);

      arb_mul_ui(delta,err,2*n-1,prec);
      arb_add_error(delta,err);
      arb_set(tmp,delta);
      for(uint64_t nn=0;nn<=Lc->NN/4;nn++)
	{
	  safe_expsinc(acb_realref(Lu->upsample_expsincs[n-1][nn]),Lc->arb_A,Lu->inv_2H2,tmp,prec);
	  arb_add(tmp,tmp,Lc->one_over_A,prec);
	}
      arb_sub(tmp,delta,Lc->one_over_A,prec);
      for(uint64_t nn=0;nn<=Lc->NN/4;nn++)
	{
	  safe_expsinc(acb_realref(Lu->upsample_expsincs[n-1][Lc->NN-1-nn]),Lc->arb_A,Lu->inv_2H2,tmp,prec);
	  arb_sub(tmp,tmp,Lc->one_over_A,prec);
	}
      acb_fft(Lu->upsample_expsincs[n-1],Lc->NN,Lc->ww,prec);
    }
  arb_clear(err);arb_clear(delta);arb_clear(tmp);
}



// shift all the values in values[]
// to be 1/(2*DIV_UPSAMPLE_RATE) to the right into values_div[]
// this avoids problems with central zeros when dividing out zeros
void shift_values(L_comp_t *Lc, L_upsample_t *Lu, uint64_t side, int64_t prec)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      setup_expsincs(Lc,Lu,prec);
      Lu->expsincs=(acb_t *)malloc(sizeof(acb_t)*Lc->NN);
      if(!Lu->expsincs)
	{
	  printf("Ran out of memory allocating Lu.expsincs. Exiting.\n");
	  exit(0);
	}
      for(uint64_t n=0;n<Lc->NN;n++)
	acb_init(Lu->expsincs[n]);
      arb_t delta,tmp;
      arb_init(delta);arb_init(tmp);
      arb_div_ui(delta,Lc->one_over_A,2*DIV_UPSAMPLE_RATE,prec);
      arb_neg(tmp,delta);
      for(uint64_t nn=0;nn<Lc->NN/4;nn++)
	{
	  expsinc(acb_realref(Lu->expsincs[nn]),Lc->arb_A,Lu->inv_2H2,tmp,prec);
	  arb_add(tmp,tmp,Lc->one_over_A,prec);
	}
      arb_add(tmp,delta,Lc->one_over_A,prec);
      arb_neg(tmp,tmp);
      for(uint64_t nn=0;nn<Lc->NN/4;nn++)
	{
	  expsinc(acb_realref(Lu->expsincs[Lc->NN-1-nn]),Lc->arb_A,Lu->inv_2H2,tmp,prec);
	  arb_sub(tmp,tmp,Lc->one_over_A,prec);
	}
      acb_fft(Lu->expsincs,Lc->NN,Lc->ww,prec);
      arb_clear(tmp);arb_clear(delta);
    }
  //printf("Shift_values called.\n");
  for(uint64_t n=0;n<Lu->no_values;n++)
    {
      arb_set(acb_realref(Lc->res[n]),Lu->values[side][n]);
      arb_zero(acb_imagref(Lc->res[n]));
    }
  for(uint64_t n=Lu->no_values;n<Lc->NN;n++)
    acb_zero(Lc->res[n]);

  acb_convolve1(Lc->res,Lc->res,Lu->expsincs,Lc->NN,Lc->ww,prec);
  
  // add the upsampling error and put into values_div, but only offset once
  for(uint64_t n=0;n<Lu->no_values_div;n++)
    arb_add(Lu->values_div[side][n],acb_realref(Lc->res[n+Lu->N*Lu->stride]),Lu->upsampling_error,prec);
  /*
  for(uint64_t n=0;n<Lu->no_values_div-2*Lu->N*Lu->stride;n++)
    {
      printf("values[%lu] ",n);arb_printd(Lu->values_start[side][n],20);
      printf(" values_div[] ");arb_printd(Lu->values_div_start[side][n],20);
      printf("\n");
    }
  */
}

int main(int argc, char **argv)
{
  printf("Command line:- %s",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  fflush(stdout);
  if(argc!=6)
    {
      printf("Usage:- %s <type> <prec> <family file> <func file names>  <target acc>.\n",argv[0]);
      return(0);
    }
  L_family_t L_family;
  L_family.ltype=atoi(argv[1]);
  print_ltype(&L_family);
  printf("Setting up L Family...\n");
  //system("date");
  uint64_t arb_prec=atol(argv[2]);
  FILE *ffile=fopen(argv[3],"r");
  if(!ffile)
    {
      printf("Failed to open family file %s. Exiting.\n",argv[3]);
      return(0);
    }
  FILE *fnamefile=fopen(argv[4],"r");
  if(!fnamefile)
    {
      printf("Failed to open L function list file %s. Exiting.\n",argv[4]);
      return(0);
    }
  L_comp_t L_comp;
  L_comp.OP_ACC=atol(argv[5]);
  printf("Operating accuracy set to %ld\n",L_comp.OP_ACC);
  FILE *cfile;
  L_func_t L_func;
  self_dual_p(&L_func,&L_family);
  L_error_t L_errors;
  L_upsample_t L_upsample;

  if(!read_comp(&L_comp,cfile,arb_prec))
    {
      printf("Error reading computation file. Exiting.\n");
      return(0);
    }
  fflush(stdout);
  
  
  if(!read_family(&L_family,&L_comp,&L_errors,ffile,arb_prec))
    {
      printf("Error reading L function family file. Exiting.\n");
      return(0);
    }
  fflush(stdout);
  
  //print_family(&L_family);


  if(!init_func(&L_func))
    {
      printf("Error initialising L-func structure. Exiting.\n");
      return(0);
    }
  fflush(stdout);


  if(!setup_upsampling(&L_upsample,&L_comp,&L_family,arb_prec))
    {
      printf("Error setting up for upsampling. Exiting.\n");
      return(0);
    }
  fflush(stdout);

  // iterate over a set of L functions in the same family
  arb_t tc,tc1; // Turing estimate
  arb_init(tc);arb_init(tc1);
  arb_t tmp,tmp1,t,delta_t;
  acb_t ctmp;
  acb_init(ctmp);
  arb_init(t);arb_init(delta_t);arb_init(tmp);arb_init(tmp1);
  arb_set_d(tmp,L_family.one_over_B*OUTPUT_RATIO*NUM_OUTPUT_POINTS);
  arb_inv(delta_t,tmp,arb_prec);
  fflush(stdout);
  system("date");
  printf("Starting...............................................\n");
  fflush(stdout);
  while(parse_line(fnamefile,&L_func,&L_family))
    {
      L_func.wt=get_weight(&L_func,&L_family); // a noop in the case of Modforms

      double normalisation=((double)L_func.wt-1.0)/2.0;
      printf("Normalisation for coefficients = n^(-%10.8e)\n",normalisation);

      if(!read_func(&L_func,&L_comp,&L_family,normalisation,arb_prec))
	  continue;

      if(!finalise_comp(&L_comp, &L_family, &L_func, ffile, arb_prec))
	{
	  printf("Error finalising computation parameters. Skipping.\n");
	  continue;
	}

      printf("Doing convolutions ");fflush(stdout);system("date");

      if(!do_convolves(&L_comp,&L_family,&L_func,arb_prec))
	{
	  printf("Error doing computation. Skipping.\n");
	  continue;
	}

      if(!finish_convolves(&L_comp,&L_family,&L_func,arb_prec))
	{
	  printf("Error finishing convolutions. Skipping.\n");
	  continue;
	}
      printf("Finished convolutions ");fflush(stdout);system("date");

      if(!do_pre_iFFT_errors(&L_comp,&L_family,&L_func,&L_errors,arb_prec))
	{
	  printf("Error doing pre-iFFT error terms. Skipping.\n");
	  continue;
	}
      printf("pre iFFT errors done ");fflush(stdout);system("date");

      if(!final_ifft(&L_comp,&L_func,&L_family,&L_errors,arb_prec))

	{
	  printf("Error doing final IFFT. Skipping.\n");
	  continue;
	}

      printf("Final IFFT done ");fflush(stdout);system("date");

      pit4(&L_comp,&L_family,&L_upsample,&L_func,arb_prec);

      if(L_family.r>MAX_TURING_R)
	{
	  shift_values(&L_comp, &L_upsample,0,arb_prec);
	  if(!L_func.self_dual_p)
	    shift_values(&L_comp, &L_upsample,1,arb_prec);
	}

      int64_t zeros_found=find_zeros(&L_comp,&L_family,&L_upsample,&L_func,0,arb_prec);
      printf("zeros_found=%ld\n",zeros_found);
      fflush(stdout);
      
      if(zeros_found==0)
	{
	  printf("Error finding zeros. Skipping.\n");
	  continue;
	}

      int64_t conj_zeros=0;

      if(!L_func.self_dual_p)
	{
	  conj_zeros=find_zeros(&L_comp,&L_family,&L_upsample,&L_func,1,arb_prec);
	  if(conj_zeros==0)
	    {
	      printf("Error finding zeros. Skipping.\n");
	      continue;
	    }
	  printf("conj zeros_found=%ld\n",conj_zeros);
	  fflush(stdout);
	}

      int64_t total_zeros;
      if(L_func.self_dual_p)
	{
	  if(zeros_found<0) // there was a central zero
	    {
	      zeros_found=-zeros_found;
	      total_zeros=zeros_found*2;
	      for(uint64_t i=0;arb_is_zero(L_comp.zeros[0][i]);i++,total_zeros--);
	    }
	  else	
	    total_zeros=zeros_found*2;
	      
	  if(zeros_found<OUTPUT_ZEROS)
	    {
	      printf("Not enough zeros for our purposes. Skipping.\n");
	      continue;
	    }
	}
      else
	{
	  if(zeros_found<0)
	    {
	      zeros_found=-zeros_found;
	      total_zeros=zeros_found+conj_zeros;
	      for(uint64_t i=0;arb_is_zero(L_comp.zeros[0][i]);i++,total_zeros--);
	    }
	  else
	    total_zeros=zeros_found+conj_zeros;

	  if((zeros_found<OUTPUT_ZEROS)||(conj_zeros<OUTPUT_ZEROS))
	    {
	      printf("Not enough zeros for our purposes. Skipping.\n");
	      continue;
	    }
	}
      
      if(L_family.r<=MAX_TURING_R) // use Turing's method
	{
	  if(!turing_count(tc,&L_comp,&L_family,&L_upsample,&L_func,arb_prec))
	    {
	      printf("Error computing Turing count. Skipping.\n");
	      continue;
	    }
	  fflush(stdout);
	  arb_sub_ui(tc1,tc,total_zeros+2,arb_prec); // can only miss pairs of zeros?
	  if(!arb_is_negative(tc1)) // need more zeros
	    {
	      printf("We only found %ld zeros.\nWe expected <=",total_zeros);
	      arb_printd(tc,10);printf(" Skipping.\n");
	      /*
	      for(uint64_t n=0;n<=L_comp.NN/OUTPUT_RATIO;n++)
		{
		  printf("%10.8e ",n/L_comp.A);
		  arb_printd(L_upsample.values_off[0][n],20);
		  printf("\n");
		}
	      */

	      continue;
	    }
	  printf("We found the %ld zeros we expected. RH is OK!\n",total_zeros);
	}
      else // can't use Turing's method
	{
	  if(!old_div_zeros(&L_comp,&L_family,&L_upsample,&L_func,0,OUTPUT_ZEROS,arb_prec))
	    continue;
	  if(!L_func.self_dual_p)
	    {
	      if(!old_div_zeros(&L_comp,&L_family,&L_upsample,&L_func,1,OUTPUT_ZEROS,arb_prec))
		continue;
	    }
	}
      
      // everything appears to have worked. Let's output it all
      char out_fname[BUFF_LEN];
      FILE *outfile;
      
      if(!make_out_fname(out_fname,&L_family,&L_func,""))
	{
	  printf("Error constructing file name for output. Skipping.\n");
	  continue;
	}
      if(!write_data(out_fname,&L_family,&L_func,&L_comp,&L_upsample,0,arb_prec))
	continue;

      if(!L_func.self_dual_p)
	{
	  if(!make_out_fname(out_fname,&L_family,&L_func,"conj."))
	    {
	      printf("Error constructing file name for output. Skipping.\n");
	      continue;
	    }
	  if(!write_data(out_fname,&L_family,&L_func,&L_comp,&L_upsample,1,arb_prec))
	    continue;
	}
      printf("L function processed OK.\n");
      fflush(stdout);
      system("date");

    }
  printf("Computation completed.\n");
  return(0);
}
