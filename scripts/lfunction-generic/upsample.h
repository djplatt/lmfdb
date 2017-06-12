// upsampling parameters from upsample_error_final.gp on lmfdb4 code/generic
#define MAX_R (12)
#define MAX_MU (8.0)

uint64_t Ns[]={143,198,249,296,342,386,430,473,515,557,599,641};
uint64_t Hs[]={60994,40040,31947,27352,24297,22077,20373,19016,17899,16958,16152,15450};
double errs[]={8.062039e-61,9.905012e-61,7.547808e-61,1.043349e-60,7.601862e-61,9.261187e-61,8.876395e-61,9.698565e-61,1.126559e-60,1.324721e-60,9.377305e-61,1.090669e-60};
double div_errs[]={7.321928e-62,4.944110e-60,1.130228e-59,8.479909e-59,1.035862e-58,3.390292e-58,2.818711e-58,4.496419e-58,1.379638e-57,2.340042e-57,2.663942e-57,2.074413e-57};

bool setup_upsampling(L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, uint64_t prec)
{
  if(Lf->r>MAX_R)
    {
      printf("Can only cope with degree <=%lu at the moment. Exiting.\n",MAX_R);
      exit(0);
    }
  for(uint64_t j=0;j<Lf->r;j++)
    if(Lf->mus[j]>MAX_MU)
      {
	printf("Can only handle mu_j upto %f. We have %f. Exiting\n",MAX_MU,Lf->mus[j]);
	exit(0);
      }
  arb_init(Lu->H);
  arb_set_ui(Lu->H,Hs[Lf->r-1]);
  arb_mul_2exp_si(Lu->H,Lu->H,-20);
  printf("Upsampling H set to ");arb_printd(Lu->H,10);printf("\n");
  arb_init(Lu->inv_2H2); // -1/2H^2
  arb_inv(Lu->inv_2H2,Lu->H,prec); // 1/H
  arb_mul(Lu->inv_2H2,Lu->inv_2H2,Lu->inv_2H2,prec); // 1/H^2
  arb_mul_2exp_si(Lu->inv_2H2,Lu->inv_2H2,-1); // 1/2H^2
  arb_neg(Lu->inv_2H2,Lu->inv_2H2); // -1/2H^2
  printf("-1/2H^2 = ");arb_printd(Lu->inv_2H2,20);printf("\n");
  Lu->N=Ns[Lf->r-1];
  printf("Upsampling N set to %lu\n",Lu->N);
  arb_init(Lu->upsampling_error);
  arb_set_d(Lu->upsampling_error,errs[Lf->r-1]);
  arb_init(Lu->div_upsampling_error);
  arb_set_d(Lu->div_upsampling_error,div_errs[Lf->r-1]);
  printf("Upsampling error set to ");arb_printd(Lu->upsampling_error,10);printf("\n");
  printf("Division upsampling error set to ");arb_printd(Lu->div_upsampling_error,10);printf("\n");
#ifdef HI_PREC
  Lu->stride=1<<4;
#else
  Lu->stride=1;
#endif
  Lu->no_values=Lc->NN/OUTPUT_RATIO+Lc->NN/TURING_RATIO+Lu->N*4*Lu->stride+1;
  Lu->no_values_div=Lc->NN/OUTPUT_RATIO+2*Lu->N*Lu->stride+1;
  Lu->values[0]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values);
  Lu->values[1]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values);
  Lu->values_off[0]=Lu->values[0]+Lu->N*Lu->stride*2;
  Lu->values_off[1]=Lu->values[1]+Lu->N*Lu->stride*2;
  Lu->values_div[0]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values_div);
  Lu->values_div[1]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values_div);
  for(uint64_t n=0;n<Lu->no_values;n++)
    {
      arb_init(Lu->values[0][n]);
      arb_init(Lu->values[1][n]);
    }
  for(uint64_t n=0;n<Lu->no_values_div;n++)
    {
      arb_init(Lu->values_div[0][n]);
      arb_init(Lu->values_div[1][n]);
    }
  Lu->values_start[0]=Lu->values[0]+Lu->N*2;
  Lu->values_start[1]=Lu->values[1]+Lu->N*2;
  Lu->values_div_start[0]=Lu->values_div[0]+Lu->N;
  Lu->values_div_start[1]=Lu->values_div[1]+Lu->N;

  return(true);
}

// sin(pi x)/(pi x)
void my_sinc(arb_t res, arb_t x, uint64_t prec)
{
  static arb_t pi,tmp,tmp1;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_mul(tmp1,x,pi,prec);
  arb_sin(tmp,tmp1,prec);
  arb_div(res,tmp,tmp1,prec);
}

// sinc x = sin(pi x)/(pi x)
// on entry sin_x = sin(pi x)
void new_my_sinc(arb_t res, arb_t x, arb_t sin_x, uint64_t prec)
{
  static arb_t pi,tmp1;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(tmp1);
    }
  arb_mul(tmp1,x,pi,prec);
  arb_div(res,sin_x,tmp1,prec);
}

// sinc x = sin(pi x)/(pi x)
// on entry pi x = pi*x, sin_x = sin(pi*x)
void new_my_sinc1(arb_t res, arb_t pix, arb_t sin_x, uint64_t prec)
{
  arb_div(res,sin_x,pix,prec);
}

// sin(pi x)/(pi x)
void my_sinc_interval(arb_t res, arb_t x, arb_t sin_x, uint64_t prec)
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
      new_my_sinc1(res,tmp,sin_x,prec);
      return;
    }
  arb_add_ui(tmp1,tmp,1,prec);
  if(!arb_is_positive(tmp1))
    {
      new_my_sinc1(res,tmp,sin_x,prec);
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



int64_t nearest_n(arb_ptr diff, arb_ptr t0, arb_t A, uint64_t prec)
{
  static arb_t tmp,tmp1;
  static fmpz_t fmpz_tmp;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      fmpz_init(fmpz_tmp);
    }
  arb_mul(tmp,t0,A,prec);
  arb_get_mid_arb(tmp1,tmp);
  arb_floor(tmp,tmp1,prec);
  if(arb_get_unique_fmpz(fmpz_tmp,tmp)==0)
    {
      printf("t0 = ");arb_printd(t0,20);printf(" tmp = ");arb_printd(tmp,20);printf("\n");
      printf("Error rounding to int in upsample routines.\n");
      return(BAD_64);
    }
  int64_t res=fmpz_get_si(fmpz_tmp); // this is going to be 1st point in left tail
  arb_set_si(tmp,res);
  arb_div(tmp,tmp,A,prec);
  arb_sub(diff,tmp,t0,prec); // will be -ve n pts to left of t0
  return(res);
}

// n points to f(n/A)
// f(t) is what we are trying to estimate
// delta = t-n/A
// note that A is Lu->stride times smaller
void do_point(arb_ptr res, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, int64_t n, arb_t t, arb_t delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2;
  static acb_t s;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(s);
      arb_set_d(acb_realref(s),0.5);
    }
  //printf("In do_point with t= ");arb_printd(t,10);printf("\n");
  //printf("delta = ");arb_printd(delta,10);printf("\n");
  arb_mul(tmp,delta,delta,prec);
  arb_mul(tmp1,tmp,Lu->inv_2H2,prec);
  arb_exp(tmp,tmp1,prec);
  //printf("exp term = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp2,delta,A,prec);
  my_sinc(tmp1,tmp2,prec);
  //printf("sinc term = ");arb_printd(tmp1,10);printf("\n");
  arb_mul(tmp2,tmp1,tmp,prec);

  /*
  arb_set(acb_imagref(s),t);
  absQ(tmp,s,Lf,Lfu->N,prec);
  arb_sqrt(tmp1,tmp,prec);
  arb_div(tmp,tmp2,tmp1,prec);
  */
  arb_set(tmp,tmp2);

  //printf("W term = ");arb_printd(Lu->values[side][n],10);printf("\n");
  arb_mul(res,tmp,Lu->values[side][n],prec);
  //printf("Upsample contribution from %ld (delta = ",n);arb_printd(delta,10);
  //printf(" was ");arb_printd(res,10);printf("\n");
}

void new_do_point(arb_ptr res, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2;
  static acb_t s;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(s);
      arb_set_d(acb_realref(s),0.5);
    }
  //printf("\nIn do_point with t= ");arb_printd(t,10);printf("\n");
  //printf("delta = ");arb_printd(delta,10);printf("\n");
  arb_mul(tmp,delta,delta,prec);
  arb_mul(tmp1,tmp,Lu->inv_2H2,prec);
  arb_exp(tmp,tmp1,prec);
  //printf("exp term = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp2,delta,A,prec);
  new_my_sinc(tmp1,tmp2,sin_delta,prec);
  //printf("sinc term = ");arb_printd(tmp1,10);printf("\n");
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,Lu->values[side][n],prec);
  //printf("Upsample contribution from %ld (delta = ",n);arb_printd(delta,10);
  //printf(" was ");arb_printd(res,10);printf("\n");
}

// a version of exp that takes more care
// exploits fact that exp is monotonic
void my_exp_interval(arb_t res, arb_t x, int64_t prec)
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

void do_point_div(arb_ptr res, L_upsample_t *Lu, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2;
  static acb_t s;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(s);
      arb_set_d(acb_realref(s),0.5);
    }
  //printf("\nIn do_point_div with t= ");arb_printd(t,10);printf("\n");
  //printf("delta = ");arb_printd(delta,10);printf("\n");
  arb_sqr(tmp,delta,prec);
  //printf("delta^2 = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp1,tmp,Lu->inv_2H2,prec);
  //printf("-delta^2/2H^2 = ");arb_printd(tmp1,10);printf("\n");
  my_exp_interval(tmp,tmp1,prec);
  //printf("exp term = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp2,delta,A,prec);
  my_sinc_interval(tmp1,tmp2,sin_delta,prec);
  //printf("sinc term = ");arb_printd(tmp1,10);printf("\n");
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,Lu->values_div[side][n],prec);
  //printf("(Int) Upsample contribution from %ld (delta = ",n);arb_printd(delta,10);
  //printf(" was ");arb_printd(res,10);printf("\n");
}

// estimate f(t0) by upsampling off Lu->values
bool upsample_stride(arb_ptr res, arb_ptr t0, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,sin_diff,neg_sin_diff;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(step);
      arb_init(t);
      arb_init(t1);
      arb_init(t_delta);
      arb_mul_ui(step,Lc->one_over_A,Lu->stride,prec);
      arb_init(A);
      arb_div_ui(A,Lc->arb_A,Lu->stride,prec);
      arb_init(diff);
      arb_init(this_diff);
      arb_init(term);
      arb_init(sin_diff);
      arb_init(neg_sin_diff);
    }

  //printf("upsampling with t0=");arb_printd(t0,100);printf("\n");
  int64_t n=nearest_n(diff,t0,Lc->arb_A,prec);
  if(n==BAD_64)
    return(false);
  arb_mul(neg_sin_diff,diff,Lc->arb_A,prec);
  arb_sin_pi(sin_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  //printf("Nearest n = %ld\n",n);
  arb_mul_ui(t_delta,Lc->one_over_A,Lu->stride,prec);
  arb_mul_si(t1,Lc->one_over_A,n,prec);
  arb_set(t,t1);
  int64_t nn=n+Lu->N*Lu->stride*2,nn1=nn; // offset into values

  //printf("nearest point is n=%lu\n",n);
  //printf("delta to t0 is ");arb_printd(diff,20);printf("\n");
  // do nearest point
  new_do_point(res,Lu,Lc,Lf,Lfu,nn,t,diff,sin_diff,A,side,prec);
  //printf("nearest point contributed ");arb_printd(res,20);printf("\n");
  arb_set(this_diff,diff);
  // do Lu->N points to left of nearest
  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_sub(this_diff,this_diff,step,prec);
      arb_sub(t,t,t_delta,prec);
      nn-=Lu->stride;
      if(count&1)
	new_do_point(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,sin_diff,A,side,prec);
      else
	new_do_point(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,neg_sin_diff,A,side,prec);

      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
    }

  arb_set(t,t1); // point to middle again
  arb_set(this_diff,diff);
  nn=nn1;

  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_add(this_diff,this_diff,step,prec);
      arb_add(t,t,t_delta,prec);
      nn+=Lu->stride;
      if(count&1)
	new_do_point(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,sin_diff,A,side,prec);
      else
	new_do_point(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,neg_sin_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
    }
  arb_add_error(res,Lu->upsampling_error);
  //exit(0);
  //printf("upsampled value is ");arb_printd(res,100);printf("\n");
  return(true);
}

// upsample off Lu->values_div, but with a radius on t (and therefore delta)
// here we use every point in values_div to keep the error term simple
bool upsample_stride_div(arb_ptr res, arb_ptr t0, uint64_t len, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,t_error,diff_error;
  static arb_t sin_diff,neg_sin_diff;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(step);
      arb_init(t);
      arb_init(t1);
      arb_init(t_delta);
      arb_mul_ui(step,Lc->one_over_A,Lu->stride,prec);
      arb_init(A);
      arb_div_ui(A,Lc->arb_A,Lu->stride,prec);
      arb_init(diff);
      arb_init(this_diff);
      arb_init(sin_diff);
      arb_init(neg_sin_diff);
      arb_init(term);
      arb_init(t_error);
      arb_init(diff_error);
      arb_mul_2exp_si(t_error,Lc->one_over_A,-1);
      arb_set_d(diff_error,Lu->stride*2);
      arb_inv(diff_error,diff_error,prec);
    }

  //printf("upsampling with t0=");arb_printd(t0,100);printf("\n");
  int64_t n=nearest_n(diff,t0,Lc->arb_A,prec);
  if(n==BAD_64)
    return(false);
  //printf("Nearest n = %ld diff = ",n);arb_printd(diff,10);printf("\n");
  arb_mul(neg_sin_diff,diff,Lc->arb_A,prec);
  arb_sin_pi(sin_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  arb_mul_ui(t_delta,Lc->one_over_A,Lu->stride,prec);
  arb_mul_si(t1,Lc->one_over_A,2*n+1,prec);
  arb_mul_2exp_si(t1,t1,-1);
  arb_add_error(t1,t_error);
  arb_set(t,t1);
  int64_t nn=n+Lu->N*Lu->stride,nn1=nn; // offset into values

  //printf("nearest point is n=%lu\n",n);
  //printf("delta to t0 is ");arb_printd(diff,20);printf("\n");
  // do nearest point
  do_point_div(res,Lu,nn,t,diff,sin_diff,A,side,prec);
  //printf("nearest point contributed ");arb_printd(res,20);printf("\n");
  arb_set(this_diff,diff);
  for(uint64_t count=0;;count++)
    {
      arb_sub(this_diff,this_diff,step,prec);
      arb_sub(t,t,t_delta,prec);
      nn-=Lu->stride;
      if(nn<0) break;
      do_point_div(term,Lu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
    }

  arb_set(t,t1); // point to middle again
  arb_set(this_diff,diff);
  nn=nn1;

  for(uint64_t count=0;;count++)
    {
      arb_add(this_diff,this_diff,step,prec);
      arb_add(t,t,t_delta,prec);
      nn+=Lu->stride;
      if(nn>len) break;
      do_point_div(term,Lu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
    }
  arb_add_error(res,Lu->div_upsampling_error);
  //printf("upsampled value is ");arb_printd(res,100);printf("\n");
  return(true);
}

// delta=n/A-t
// sin_delta=sin(Pi*A*(n/A-t))
void do_point_f_dash(arb_ptr res, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t cos_delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2,tmp3;
  static arb_t pi,A_pi_delta,A_pi;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(A_pi);
      arb_mul(A_pi,Lc->arb_A,pi,prec);
      arb_init(A_pi_delta);
    }
  //printf("In do_point_f_dash with n/A-t= ");arb_printd(delta,20);printf("\n");
  arb_mul(A_pi_delta,delta,A_pi,prec);
  arb_div(tmp,sin_delta,A_pi_delta,prec); // sin(A pi delta)/A pi delta
  arb_sub(tmp1,tmp,cos_delta,prec); // sin(A pi delta)/A pi delta-cos(A pi delta)
  arb_div(tmp,tmp1,delta,prec);
  arb_mul(tmp1,delta,delta,prec);
  arb_mul(tmp2,tmp1,Lu->inv_2H2,prec);
  arb_exp(tmp1,tmp2,prec);
  //printf("exp() =");arb_printd(tmp1,20);printf("\n)");
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,Lu->values[side][n],prec);
  //printf("do_point_f_dash returning ");arb_printd(res,20);printf("\n");
}

void do_point_f_dash_div(arb_ptr res, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t cos_delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2,tmp3;
  static arb_t pi,A_pi_delta,A_pi;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(A_pi);
      arb_mul(A_pi,Lc->arb_A,pi,prec);
      arb_init(A_pi_delta);
    }
  //printf("In do_point_f_dash with n/A-t= ");arb_printd(delta,20);printf("\n");
  arb_mul(A_pi_delta,delta,A_pi,prec);
  arb_div(tmp,sin_delta,A_pi_delta,prec); // sin(A pi delta)/A pi delta
  arb_sub(tmp1,tmp,cos_delta,prec); // sin(A pi delta)/A pi delta-cos(A pi delta)
  arb_div(tmp,tmp1,delta,prec);
  arb_mul(tmp1,delta,delta,prec);
  arb_mul(tmp2,tmp1,Lu->inv_2H2,prec);
  arb_exp(tmp1,tmp2,prec);
  //printf("exp() =");arb_printd(tmp1,20);printf("\n)");
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,Lu->values_div[side][n],prec);
  //printf("do_point_f_dash returning ");arb_printd(res,20);printf("\n");
}


bool f_dash1(arb_ptr res, arb_ptr t0, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,t_error,diff_error;
  static arb_t sin_diff,neg_sin_diff,cos_diff,neg_cos_diff;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(step);
      arb_init(t);
      arb_init(t1);
      arb_init(t_delta);
      arb_mul_ui(step,Lc->one_over_A,Lu->stride,prec);
      arb_init(A);
      arb_div_ui(A,Lc->arb_A,Lu->stride,prec);
      arb_init(diff);
      arb_init(this_diff);
      arb_init(sin_diff);
      arb_init(neg_sin_diff);
      arb_init(cos_diff);
      arb_init(neg_cos_diff);
      arb_init(term);
      arb_init(t_error);
      arb_init(diff_error);
      arb_mul_2exp_si(t_error,Lc->one_over_A,-1);
      arb_set_d(diff_error,Lu->stride*2);
      arb_inv(diff_error,diff_error,prec);
    }

  int64_t n=nearest_n(diff,t0,Lc->arb_A,prec);
  if(n==BAD_64)
    return(false);
  //printf("Nearest n = %ld diff = ",n);arb_printd(diff,10);printf("\n");
  arb_mul(neg_sin_diff,diff,Lc->arb_A,prec);
  arb_sin_cos_pi(sin_diff,cos_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  arb_neg(neg_cos_diff,cos_diff);
  arb_mul_ui(t_delta,Lc->one_over_A,Lu->stride,prec);
  arb_mul_si(t1,Lc->one_over_A,2*n+1,prec);
  arb_mul_2exp_si(t1,t1,-1);
  arb_add_error(t1,t_error);
  arb_set(t,t1);
  int64_t nn=n+Lu->N*Lu->stride*2,nn1=nn; // offset into values

  //printf("nearest point is n=%lu\n",n);
  //printf("delta to t0 is ");arb_printd(diff,20);printf("\n");
  // do nearest point
  do_point_f_dash(res,Lu,Lc,Lf,Lfu,nn,t,diff,sin_diff,cos_diff,A,side,prec);
  //printf("nearest point contributed ");arb_printd(res,20);printf("\n");
  arb_set(this_diff,diff);
  // do Lu->N points to left of nearest
  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_sub(this_diff,this_diff,step,prec);
      arb_sub(t,t,t_delta,prec);
      nn-=Lu->stride;
      do_point_f_dash(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
      //printf("Res now ");arb_printd(res,20);printf("\n");
    }

  arb_set(t,t1); // point to middle again
  arb_set(this_diff,diff);
  nn=nn1;

  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_add(this_diff,this_diff,step,prec);
      arb_add(t,t,t_delta,prec);
      nn+=Lu->stride;
      do_point_f_dash(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
      //printf("Res now ");arb_printd(res,20);printf("\n");
    }
  // nothing rigorous about N-R
  //arb_add_error(res,Lu->upsampling_error);
  //printf("upsampled value is ");arb_printd(res,100);printf("\n");
  return(true);
}

bool f_dash1_div(arb_ptr res, arb_ptr t0, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,t_error,diff_error;
  static arb_t sin_diff,neg_sin_diff,cos_diff,neg_cos_diff;
  static bool init=false;
  if(!init)
    {
      printf("REMEMBER TO ADD ERROR TERM IN f_dash1_div.\n");
      init=true;
      arb_init(step);
      arb_init(t);
      arb_init(t1);
      arb_init(t_delta);
      arb_mul_ui(step,Lc->one_over_A,Lu->stride,prec);
      arb_init(A);
      arb_div_ui(A,Lc->arb_A,Lu->stride,prec);
      arb_init(diff);
      arb_init(this_diff);
      arb_init(sin_diff);
      arb_init(neg_sin_diff);
      arb_init(cos_diff);
      arb_init(neg_cos_diff);
      arb_init(term);
      arb_init(t_error);
      arb_init(diff_error);
      arb_mul_2exp_si(t_error,Lc->one_over_A,-1);
      arb_set_d(diff_error,Lu->stride*2);
      arb_inv(diff_error,diff_error,prec);
    }

  int64_t n=nearest_n(diff,t0,Lc->arb_A,prec);
  if(n==BAD_64)
    return(false);
  //printf("Nearest n = %ld diff = ",n);arb_printd(diff,10);printf("\n");
  arb_mul(neg_sin_diff,diff,Lc->arb_A,prec);
  arb_sin_cos_pi(sin_diff,cos_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  arb_neg(neg_cos_diff,cos_diff);
  arb_mul_ui(t_delta,Lc->one_over_A,Lu->stride,prec);
  arb_mul_si(t1,Lc->one_over_A,2*n+1,prec);
  arb_mul_2exp_si(t1,t1,-1);
  arb_add_error(t1,t_error);
  arb_set(t,t1);
  int64_t nn=n+Lu->N*Lu->stride,nn1=nn; // offset into values_div

  //printf("nearest point is n=%lu\n",n);
  //printf("delta to t0 is ");arb_printd(diff,20);printf("\n");
  // do nearest point
  do_point_f_dash_div(res,Lu,Lc,Lf,Lfu,nn,t,diff,sin_diff,cos_diff,A,side,prec);
  //printf("nearest point contributed ");arb_printd(res,20);printf("\n");
  arb_set(this_diff,diff);
  // do Lu->N points to left of nearest
  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_sub(this_diff,this_diff,step,prec);
      arb_sub(t,t,t_delta,prec);
      nn-=Lu->stride;
      do_point_f_dash_div(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
      //printf("Res now ");arb_printd(res,20);printf("\n");
    }

  arb_set(t,t1); // point to middle again
  arb_set(this_diff,diff);
  nn=nn1;

  for(uint64_t count=0;count<Lu->N;count++)
    {
      arb_add(this_diff,this_diff,step,prec);
      arb_add(t,t,t_delta,prec);
      nn+=Lu->stride;
      do_point_f_dash_div(term,Lu,Lc,Lf,Lfu,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,A,side,prec);
      //printf("Point at %ld contributed ",nn);arb_printd(term,10);printf("\n");
      arb_add(res,res,term,prec);
      //printf("Res now ");arb_printd(res,20);printf("\n");
    }
  //arb_add_error(res,Lu->upsampling_error);
  //printf("upsampled value is ");arb_printd(res,100);printf("\n");
  return(true);
}


#define N_NEWTON_ITS (5)
// using arb is overkill here as its not rigorous anyway
bool newton(arb_ptr res, arb_ptr t0, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static bool init=false;
  static arb_t f0,fd,t,tmp1;
  if(!init)
    {
      init=true;
      arb_init(f0);
      arb_init(fd);
      arb_init(t);
      arb_init(tmp1);
    }
  //printf("In Newton with t0=");arb_printd(t0,20);printf("\n");
  arb_set(t,t0);
  for(uint64_t n=0;n<N_NEWTON_ITS;n++)
    {      
      if(!upsample_stride(f0,t,Lu,Lc,Lf,Lfu,side,prec)) // f(t)
	return(false);
      //printf("f(");arb_printd(t,20);printf(")=");arb_printd(f0,20);printf("\n");
      if(!f_dash1(fd,t,Lu,Lc,Lf,Lfu,side,prec))
	return(false);
      //printf("f'(");arb_printd(t,20);printf(")=");arb_printd(fd,20);printf("\n");
      arb_div(tmp1,f0,fd,prec);
      arb_sub(t,t,tmp1,prec);
      //printf("New value of t is ");arb_printd(t,40);printf("\n");
    }
  arb_set(res,t);
  return(true);
}
