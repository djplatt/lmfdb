// return Q(s) is analytic conductor defined pg 387 col 2
void Q(acb_t res, acb_t s, L_family_t *Lf, uint64_t conductor, uint64_t prec)
{
  //printf("In logQ with s=");acb_printd(s,10);printf("\n");
  static bool init=false;
  static arb_t two_pi,tmp1,tmp2;
  static acb_t stmp1,stmp2,stmp3;
  if(!init)
    {
      init=true;
      arb_init(two_pi);arb_init(tmp1);arb_init(tmp2);
      acb_init(stmp1);acb_init(stmp2);acb_init(stmp3);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }
  acb_set_ui(stmp1,conductor);
  arb_set(acb_imagref(stmp2),acb_imagref(s));
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_set_d(tmp2,Lf->mus[j]);
      arb_add(acb_realref(stmp2),acb_realref(s),tmp2,prec);
      acb_mul(stmp3,stmp1,stmp2,prec);
      acb_div_arb(stmp1,stmp3,two_pi,prec);
    }
  //printf("Analytic conductor = ");acb_printd(stmp1,10);printf("\n");
  acb_set(res,stmp1);
}

// |Q(s)|
void absQ(arb_t res, acb_t s, L_family_t *Lf, uint64_t conductor, uint64_t prec)
{
  static bool init=false;
  static acb_t tmp1;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
    }
  Q(tmp1,s,Lf,conductor,prec);
  acb_abs(res,tmp1,prec);
}

// log |Q(s)|
void logQ(arb_t res, acb_t s, L_family_t *Lf, uint64_t conductor, uint64_t prec)
{
  //printf("In logQ with s=");acb_printd(s,10);printf("\n");
  static bool init=false;
  static arb_t tmp1;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
    }
  absQ(tmp1,s,Lf,conductor,prec);
  arb_log(res,tmp1,prec);
} /* logQ */

void abs_gamma_r(arb_t res, acb_t s, uint64_t prec)
{
  static bool init=false;
  static acb_t s_by_2,lg;
  static arb_t tmp,tmp0,log_pi;
  if(!init)
    {
      init=true;
      acb_init(s_by_2);
      acb_init(lg);
      arb_init(tmp);
      arb_init(tmp0);
      arb_const_pi(tmp,prec);
      arb_init(log_pi);
      arb_log(log_pi,tmp,prec);
    }

  //printf("in abs_gamma_r with s = ");acb_printd(s,10);printf("\n");
  acb_mul_2exp_si(s_by_2,s,-1); // s/2
  acb_lgamma(lg,s_by_2,prec);
  //printf("loggamma(s/2) = ");acb_printd(lg,10);printf("\n");

  arb_mul(tmp,log_pi,acb_realref(s_by_2),prec);

  arb_sub(tmp0,acb_realref(lg),tmp,prec);
  arb_exp(res,tmp0,prec);
  //printf("gamma_r returning ");arb_printd(res,prec);printf("\n");
}

// |gamma(s)| per top col 2 pg 387 without N^1/2(s-1/2)
void abs_gamma(arb_t res, acb_t s, L_family_t *Lf, uint64_t prec)
{
  acb_t tmp1,tmp2,tmp3;
  acb_init(tmp1);
  acb_init(tmp2);
  acb_init(tmp3);
  arb_t tmp;
  arb_init(tmp);

  arb_set_ui(res,1);
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_set_d(tmp,Lf->mus[j]);
      acb_add_arb(tmp1,s,tmp,prec);
      abs_gamma_r(tmp,tmp1,prec);
      arb_mul(res,res,tmp,prec);
    }
  //printf("abs_gamma returning ");arb_printd(res,10);printf("\n");
  acb_clear(tmp1);
  acb_clear(tmp2);
  acb_clear(tmp3);
  arb_clear(tmp);
}

bool init_ftwiddle_error(L_family_t *Lf, L_comp_t *Lc, uint64_t prec)
{
  acb_t s,s_plus_mu;
  arb_t tmp1,tmp2,tmp3,tmp4,E,two_pi,four_by_pi2,beta;
  acb_init(s);
  acb_init(s_plus_mu);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(tmp4);
  arb_init(E);
  arb_init(two_pi);
  arb_init(four_by_pi2);
  arb_init(beta);
  arb_const_pi(two_pi,prec);
  arb_inv(tmp1,two_pi,prec);
  arb_mul(four_by_pi2,tmp1,tmp1,prec);
  arb_mul_2exp_si(four_by_pi2,four_by_pi2,2); // 4/Pi^2
  arb_mul_2exp_si(two_pi,two_pi,1); // 2Pi

  arb_set_d(acb_realref(s),0.5);
  arb_set(acb_imagref(s),Lf->B); // 1/2+iB

  arb_set_d(tmp2,1.5);
  arb_zeta(tmp1,tmp2,prec);
  arb_pow_ui(tmp2,tmp1,Lf->r,prec);

  abs_gamma(tmp1,s,Lf,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  /*
  arb_set(acb_imagref(s_plus_mu),acb_imagref(s));
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_set_d(tmp1,Lf->mus[j]);
      arb_add(acb_realref(s_plus_mu),acb_realref(s),tmp1,prec);
      acb_abs(tmp2,s_plus_mu,prec);
      arb_div(tmp1,tmp2,two_pi,prec);
      arb_mul(E,E,tmp1,prec);
    }
  */
  absQ(tmp1,s,Lf,1,prec);
  arb_sqrt(tmp2,tmp1,prec);
  arb_mul(E,tmp2,tmp3,prec);

  // E without the N

  arb_mul_2exp_si(E,E,1);



  arb_mul_2exp_si(tmp1,two_pi,-3); //pi/4
  arb_mul_ui(beta,tmp1,Lf->r,prec); // pi r/4

  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_set_d(tmp1,0.5+Lf->mus[j]);
      arb_div(tmp2,tmp1,Lf->B,prec);
      arb_atan(tmp1,tmp2,prec);
      arb_sub(beta,beta,tmp1,prec);
    }

  arb_mul(tmp2,Lf->B,Lf->B,prec);
  arb_zero(tmp4);
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_set_d(tmp1,(0.5+Lf->mus[j])*(0.5+Lf->mus[j]));
      arb_sub(tmp3,tmp2,tmp1,prec);
      arb_inv(tmp1,tmp3,prec);
      arb_add(tmp4,tmp4,tmp1,prec);
    }
  arb_mul(tmp1,tmp4,four_by_pi2,prec);
  arb_sub(beta,beta,tmp1,prec);

  arb_mul(tmp1,beta,Lf->B,prec);
  arb_neg(tmp1,tmp1);
  arb_exp(tmp2,tmp1,prec);
  arb_sub_ui(tmp1,tmp2,1,prec);
  arb_neg(tmp1,tmp1);
  arb_div(Lf->ftwiddle_error,E,tmp1,prec);
  //printf("F twiddle error set to ");arb_printd(Lf->ftwiddle_error,10);
  //printf("\n");
  
  acb_clear(s);
  acb_clear(s_plus_mu);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
  arb_clear(tmp4);
  arb_clear(two_pi);
  arb_clear(E);
  arb_clear(four_by_pi2);

  return(true);
}

// multiply the error precomputed in L_family by sqrt{conductor} in Q(s)
void complete_ftwiddle_error(L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_sqrt_ui(tmp,Lfu->N,prec);
  arb_mul(Lfu->ftwiddle_error,Lf->ftwiddle_error,tmp,prec);
}
  

bool M_error(arb_t res, arb_t x, uint64_t M, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static arb_t C_bit,arb_M,tmp1,tmp2,tmp3,tmp4,u,X,half_log_N,pi,XM2r,Mc;
  if(!init)
    {
      init=true;
      arb_init(arb_M);
      arb_init(C_bit);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(u);
      arb_init(X);
      arb_init(half_log_N);
      arb_init(pi);
      arb_init(XM2r);
      arb_init(Mc);
      arb_init(tmp4);
      // pi = pi
      arb_const_pi(pi,prec);
      // C_bit = C2^(r/2)
      arb_sqrt_ui(tmp1,3,prec);
      arb_set_ui(tmp2,2);
      arb_set_d(tmp3,Lf->r/2.0);
      arb_pow(tmp4,tmp2,tmp3,prec);
      arb_mul(C_bit,tmp1,tmp4,prec);
    }
  // half_log_N = 1/2 log(N)
  arb_log_ui(half_log_N,Lfu->N,prec);
  arb_mul_2exp_si(half_log_N,half_log_N,-1);
  // u = x-1/2 log N
  arb_sub(u,x,half_log_N,prec);
  //printf("u=");arb_printd(u,10);printf("\n");
  // X = Pi r exp(2u/r)
  arb_div_ui(tmp1,u,Lf->r,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_exp(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,Lf->r,prec);
  arb_mul(X,tmp1,pi,prec);
  //printf("X=");arb_printd(X,10);printf("\n");
  // XM2r=XM^(2/r)
  arb_set_d(tmp1,2.0);
  arb_div_ui(tmp2,tmp1,Lf->r,prec);
  arb_set_ui(arb_M,M);
  arb_pow(tmp1,arb_M,tmp2,prec); // M^(2/r)
  arb_mul(XM2r,tmp1,X,prec); // XM^(2/r)
  arb_set_d(tmp1,Lf->r/2.0);
  arb_sub(tmp2,XM2r,tmp1,prec);
  if(!arb_is_positive(tmp2))
    {
      printf("Failed XM^(r/2)>r/2 test in Merror.\n");
      printf("X M^(2/r)=");arb_printd(XM2r,10);printf("\n");
      return(false);
    }
  arb_mul(tmp2,tmp1,Lf->c,prec);
  arb_sub_ui(tmp1,tmp2,1,prec);
  arb_sub(tmp2,XM2r,tmp1,prec);
  if(!arb_is_positive(tmp2))
    {
      printf("Failed XM^(r/2)>cr/2-1 test in Merror. Exiting.\n");
      printf("X M^(2/r)=");arb_printd(XM2r,10);printf("\n");
      return(false);
    }
  
  // Mc=M^c-1
  arb_sub_ui(tmp1,Lf->c,1,prec);
  arb_pow(Mc,arb_M,tmp1,prec);
  //printf("M^{c-1}=");arb_printd(Mc,10);printf("\n");
  arb_mul(tmp1,Lf->nu,u,prec);
  arb_sub(tmp3,tmp1,XM2r,prec); // nu*u-XM^(2/r)
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp1,Mc,prec); // M^{c-1} exp(nu*u-XM^(2/r))
  arb_add_ui(tmp1,XM2r,1,prec);
  arb_mul_2exp_si(tmp1,tmp1,1); // 2XM^(r/2)+2
  arb_mul_ui(tmp2,Lf->c,Lf->r,prec); // cr
  arb_sub(tmp4,tmp1,tmp2,prec);
  arb_mul_ui(tmp1,arb_M,Lf->r,prec);
  arb_div(tmp2,tmp1,tmp4,prec);
  arb_add_ui(tmp1,tmp2,1,prec);
  arb_mul(tmp2,tmp1,C_bit,prec);
  arb_mul(res,tmp2,tmp3,prec);
  //printf("before product=");arb_printd(res,10);printf("\n");
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_mul_ui(tmp2,Lf->nus[j],Lf->r,prec);
      arb_div(tmp1,tmp2,XM2r,prec);
      arb_add_ui(tmp2,tmp1,1,prec);
      arb_pow(tmp1,tmp2,Lf->nus[j],prec);
      arb_mul(res,res,tmp1,prec);
    }
  return(true);
}

bool F_hat_twiddle_error(arb_t res, arb_t x, L_family_t *Lf, L_func_t *Lfu, L_comp_t *Lc, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,u,X,half_log_N,pi,twoXbyr;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(u);
      arb_init(X);
      arb_init(half_log_N);
      arb_init(pi);
      arb_init(twoXbyr);
      // pi = pi
      arb_const_pi(pi,prec);
      // C_bit = Cr2^(r/2-1)
      // half_log_N = 1/2 log(N)
    }
  arb_log_ui(half_log_N,Lfu->N,prec);
  arb_mul_2exp_si(half_log_N,half_log_N,-1);
  // u = x-1/2 log N
  arb_sub(u,x,half_log_N,prec);
  //printf("u=");arb_printd(u,10);printf("\n");
  // X = Pi r exp(2u/r)
  arb_div_ui(tmp1,u,Lf->r,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_exp(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,Lf->r,prec);
  arb_mul(X,tmp1,pi,prec);
  // check X>r/2
  arb_set_d(tmp1,Lf->r/2.0);
  arb_sub(tmp2,X,tmp1,prec);
  if(!arb_is_positive(tmp2))
    {
      printf("Fhat twiddle error failed X>r/2 test. Exiting.\n");
      return(false);
    }
  //printf("X=");arb_printd(X,10);printf("\n");
  // 2X/r
  arb_div_ui(twoXbyr,X,Lf->r,prec);
  arb_mul_2exp_si(twoXbyr,twoXbyr,1);
  arb_zeta(tmp2,twoXbyr,prec);
  arb_pow_ui(tmp1,tmp2,Lf->r,prec);
  arb_set_d(tmp2,0.5);
  arb_sub(tmp3,tmp2,twoXbyr,prec);
  arb_mul(tmp2,tmp3,Lc->arb_A,prec);
  arb_mul(tmp3,tmp2,pi,prec);
  arb_mul_2exp_si(tmp3,tmp3,1); // 2piA(1/2-2X/r)
  arb_exp(tmp2,tmp3,prec);
  arb_sub_ui(tmp3,tmp2,1,prec);
  arb_neg(tmp3,tmp3);
  arb_div(tmp2,tmp1,tmp3,prec);
  arb_mul(tmp1,u,Lf->nu,prec);
  arb_sub(tmp3,tmp1,X,prec);
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  arb_sqrt_ui(tmp1,2,prec);
  arb_pow_ui(tmp2,tmp1,Lf->r,prec);
  arb_mul(res,tmp2,tmp3,prec);
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_mul_ui(tmp1,Lf->nus[j],Lf->r,prec);
      arb_div(tmp2,tmp1,X,prec);
      arb_add_ui(tmp1,tmp2,1,prec);
      arb_pow(tmp2,tmp1,Lf->nus[j],prec);
      arb_mul(res,res,tmp2,prec);
    }
  //printf("F_hat_twiddle_error returning ");arb_printd(res,10);printf("\n");
  return(true);
}
/*
// thanks Bober
void acb_reasonable_sqrt(acb_t out, acb_t in, slong prec) 
{
  if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) 
    {
      acb_neg(in, in);
      acb_sqrt(out, in, prec);
      acb_mul_onei(out, out);
      acb_neg(in, in);
    }
  else 
    acb_sqrt(out, in, prec);
}
*/
void acb_reasonable_arg(arb_t out, acb_t in, uint64_t prec)
{
  static bool init=false;
  static arb_t pi;
  static acb_t tmp;
  if(!init)
    {
      init=true;
      acb_init(tmp);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }

  if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) 
    {
      acb_neg(tmp, in);
      acb_arg(out, tmp, prec);
      arb_add(out,out,pi,prec);
    }
  else 
    acb_arg(out, in, prec);
}

// compute epsilon from x and y such that x*epsilon = conj(y*epsilon)  
void fix_epsilon(acb_t res, acb_t x, acb_t y,uint64_t prec)
{
  static arb_t th1,th2;
  static arb_t ep,two_pi;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
      arb_init(th1);
      arb_init(th2);
      arb_init(ep);
    }
  //printf("Doing fix_epsilon based on ");acb_printd(x,10);
  //printf(" ");acb_printd(y,10);
  acb_reasonable_arg(th1,x,prec);
  //printf("\narg(x)=");arb_printd(th1,10);printf("\n");
  if(arb_is_negative(th1))
    arb_add(th1,th1,two_pi,prec);
  //printf("arg(x) now =");arb_printd(th1,10);printf("\n");
  acb_reasonable_arg(th2,y,prec);
  //printf("arg(y)=");arb_printd(th2,10);printf("\n");
  if(arb_is_negative(th2))
    arb_add(th2,th2,two_pi,prec);

  if(arb_gt(th1,th2))
    {
      arb_sub(ep,th1,th2,prec);
      arb_mul_2exp_si(ep,ep,-1);
      arb_add(ep,ep,th2,prec);
    }
  else
    {
      arb_sub(ep,th2,th1,prec);
      arb_mul_2exp_si(ep,ep,-1);
      arb_add(ep,ep,th1,prec);
    }
  arb_neg(ep,ep);

  arb_sin_cos(acb_imagref(res),acb_realref(res),ep,prec);
}

// when this is called, we have convolved with all
// a_m from M0 to M and then done m=1->M0 simplistically
// on all n from 0 to hi_i
// thus n=0 has all m<=M
// n=n has all m<=max(floor(exp((hi_i-n+0.5)*2*Pi/B)*sqrt(N)),M0)

// given n, what is the last a_m have we summed into res[n]
// we need to add the error for sum m>=m+1
// returns 0 if no coefficients have been added in
// then use F_hat_twiddle bound
uint64_t inv_m(uint64_t hi_i, uint64_t n, uint64_t M0, double one_over_B, double dc)
{
  return(exp(((double) hi_i-(double) n+0.5)*2.0*M_PI*one_over_B)*dc);
}

bool do_pre_iFFT_errors(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, L_error_t *Le, uint64_t prec)
{
  static arb_t err,x,fhattwiddle,pi,two_pi_A,tmp,th,th1,th2;
  static acb_t ctmp1,ctmp2;
  static bool init=false;
  if(!init)
    {
      arb_init(th);
      arb_init(th2);
      arb_init(th2);
      arb_init(fhattwiddle);
      arb_init(err);arb_init(x);arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(two_pi_A);
      arb_mul(two_pi_A,pi,Lc->arb_A,prec);
      arb_mul_2exp_si(two_pi_A,two_pi_A,1);
      arb_init(tmp);
      acb_init(ctmp1);
      acb_init(ctmp2);
      init=true;
    }
  uint64_t i=0;
  arb_zero(x);
  for(i=0;i<=Lc->N/2;i++)
    {
      uint64_t M=inv_m(Lf->hi_i,i,Lc->M0,Lf->one_over_B,Lfu->dc);
      //printf("i=%lu => M=%lu ",i,M);
      if(M==0)
	break;
      if(M>Lfu->M)
	M=Lfu->M;
      if(!M_error(err,x,M+1,Lf,Lfu,prec))
	return(false);
      if(i==0)
	{
	  printf("M Error for n=%lu is ",i);
	  arb_printd(err,10);
	  printf("\n");
	}
      arb_add_error(acb_realref(Lc->res[i]),err);
      arb_add_error(acb_imagref(Lc->res[i]),err);
      if(i==1) // use same error for F_hat[-1] in case its needed for epsilon
	{
	  arb_add_error(acb_realref(Lc->res[Lc->N-1]),err);
	  arb_add_error(acb_imagref(Lc->res[Lc->N-1]),err);
	}
      arb_add(x,x,Lf->two_pi_by_B,prec);
    }

  // error sum k\neq 0 F_hat(x+2\pi k A)
  arb_sub(err,two_pi_A,x,prec);
  if(!F_hat_twiddle_error(fhattwiddle,err,Lf,Lfu,Lc,prec))
    return(false);
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  /*
  printf("F_hat_twiddle error at n=%lu = ",i);
  arb_printd(fhattwiddle,10);
  printf("\n");
  */

  arb_mul(err,Le->eq59,Lfu->sum_ans,prec);
  printf("Adding eq 5-9 error = ");arb_printd(err,10);printf("\n");

  for(uint64_t j=0;j<i;j++)
    {
      arb_add_error(acb_realref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_imagref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_realref(Lc->res[i]),err);
      arb_add_error(acb_imagref(Lc->res[i]),err);
    }

  // figure out epsilon
  
  if(acb_contains_zero(Lc->res[0])) // can't use F_hat[0]
    {
      // first set the error terms for F_hat[-1]
      arb_add_error(acb_realref(Lc->res[Lc->N-1]),fhattwiddle);
      arb_add_error(acb_imagref(Lc->res[Lc->N-1]),fhattwiddle);
      arb_add_error(acb_realref(Lc->res[Lc->N-1]),err);
      arb_add_error(acb_imagref(Lc->res[Lc->N-1]),err);
      printf("F_hat[1]=");acb_printd(Lc->res[1],30);printf("\n");
      printf("F_hat[-1]=");acb_printd(Lc->res[Lc->N-1],30);printf("\n");
      printf("Using F_hat[1] and F_hat[-1] to fix epsilon.\n");
      fix_epsilon(Lfu->epsilon,Lc->res[1],Lc->res[Lc->N-1],prec);
    }
  else // can use F_hat[0]
    {  
      printf("Using F_hat[0] to fix epsilon.\n");
      printf("F_hat[0]=");acb_printd(Lc->res[1],10);printf("\n");
      //printf("F_hat[1]=");acb_printd(Lc->res[1],10);printf("\n");
      //printf("F_hat[-1]=");acb_printd(Lc->res[Lc->N-1],10);printf("\n");
      if((arb_is_negative(acb_realref(Lc->res[0])))&&(arb_contains_zero(acb_imagref(Lc->res[0]))))
	{
	  acb_neg(ctmp1,Lc->res[0]);
	  acb_arg(th,ctmp1,prec);
	  arb_sin_cos(th1,th2,th,prec);
	  arb_set(acb_realref(Lfu->epsilon),th2);
	  arb_set(acb_imagref(Lfu->epsilon),th1);
	}
      else
	{
	  acb_arg(th,Lc->res[0],prec);
	  //acb_arg(th2,Lc->res[Lc->N/2],prec);
	  //arb_set(th,th1); //arb_intersect(th,th1,th2,prec);
	  arb_neg(th,th);
	  arb_sin_cos(th1,th2,th,prec);
	  arb_set(acb_realref(Lfu->epsilon),th2);
	  arb_set(acb_imagref(Lfu->epsilon),th1);
	  //fix_epsilon(ctmp1,Lc->res[1],Lc->res[Lc->N-1],prec);
	  //printf("Or we could have used ");acb_printd(ctmp1,10);printf("\n");
	}
    }
  //arb_one(acb_imagref(Lfu->epsilon));arb_zero(acb_realref(Lfu->epsilon));
  printf("epsilon set to ");
  acb_printd(Lfu->epsilon,10);
  printf("\n");
  acb_sqr(Lfu->epsilon_sqr,Lfu->epsilon,prec);
  for(uint64_t n=0;n<i;n++)
    acb_mul(Lc->res[n],Lc->res[n],Lfu->epsilon,prec);
  
  // the rest of the vector we just approximate with F_hat_twiddle error

  // error sum F_hat(x+2\pi k A)
  if(!F_hat_twiddle_error(fhattwiddle,x,Lf,Lfu,Lc,prec))
    return(false);
  // we have sum k\geq 0 F_hat(x+2\pi k A)
  // since x<\pi A we can just double this
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  printf("F_hat_twiddle error beyond n=%lu (x = ",i);
  arb_printd(x,10);
  printf(" ) = ");
  arb_printd(fhattwiddle,10);
  printf("\n");
  for(;i<=Lc->NN/2;i++)
    {
      acb_zero(Lc->res[i]);
      arb_add_error(acb_realref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_imagref(Lc->res[i]),fhattwiddle);
    }
  for(uint64_t n=Lc->NN/2+1;n<Lc->NN;n++)
    acb_conj(Lc->res[n],Lc->res[Lc->NN-n]);
  return(true);
}

