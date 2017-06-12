bool quick_test(L_upsample_t *Lu, int64_t side, int64_t last_n)
{
  if(arb_is_negative(Lu->values_div_start[side][0]))
    {
      for(uint64_t n=1;n<=last_n;n++)
	if(!arb_is_negative(Lu->values_div_start[side][n]))
	  return false;
      return true;
    }
  if(arb_is_positive(Lu->values_div_start[side][0]))
    {
      for(uint64_t n=1;n<=last_n;n++)
	if(!arb_is_positive(Lu->values_div_start[side][n]))
	  return false;
      return true;
    }
  return false;
}

bool proper_test(L_upsample_t *Lu, int64_t side, int64_t last_n, arb_t A, arb_t one_over_A, int64_t prec)
{
  static bool init=false;
  static arb_t fmax,fmin,tmp,tmp1,tmp2;

  if(!init)
    {
      init=true;
      arb_init(fmax);
      arb_init(fmin);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_set_d(tmp2,3.5);
      arb_mul(tmp2,tmp2,A,prec);
      arb_mul(tmp2,tmp2,A,prec);
      arb_mul(tmp2,tmp2,Lu->H,prec); // max possible f'(t)
      arb_mul(tmp2,tmp2,one_over_A,prec); // max possible change in f over one step
      printf("max delta in f set to ");arb_printd(tmp2,10);printf("\n");
    }

  arb_zero(fmax);
  for(uint64_t n=0;n<=last_n+2*Lu->N*Lu->stride;n++)
    {
      arb_abs(tmp,Lu->values_div[side][n]);
      arb_sub(tmp1,tmp,fmax,prec);
      if(!arb_is_negative(tmp1))
	arb_set(fmax,tmp);
    }

  arb_abs(fmin,Lu->values_div_start[side][0]);
  for(uint64_t n=1;n<=last_n;n++)
    {
      arb_abs(tmp,Lu->values_div_start[side][n]);
      arb_sub(tmp1,tmp,fmin,prec);
      if(!arb_is_positive(tmp1))
	arb_set(fmin,tmp);
    }
  printf("max |f| in proper_test = ");arb_printd(fmax,20);printf("\n");
  printf("min |f| in proper_test = ");arb_printd(fmin,20);printf("\n");
  return false;
}

void print_values_div(L_upsample_t *Lu, double A, int64_t last_n, uint64_t side)
{
  for(uint64_t n=0;n<=last_n;n++)
    {
      printf("%e ",((double) n-(1.0/DIV_UPSAMPLE_RATE*2))/A);
      arb_printd(Lu->values_div_start[side][n],20);
      printf("\n");
    }
}

//#define T_MINUS_RHO
// if defined, divide out by (s-rho)
// otherwise
// divide out by (1-t/rho)
void div_zero(L_upsample_t *Lu, L_comp_t *Lc, uint64_t n_zeros, uint64_t side, int64_t last_n, int64_t prec)
{
  static bool init=false;
  static arb_t tmp,t,tmp1;
  if(!init)
    {
      init=true;
      arb_init(t);
      arb_init(tmp);
      arb_init(tmp1);
    }
  for(uint64_t n=0;n<=last_n+2*Lu->N*Lu->stride;n++)
    {
      arb_mul_si(t,Lc->one_over_A,((n-Lu->N*Lu->stride)*2*DIV_UPSAMPLE_RATE)-1,prec); // compute t for this n
      arb_mul_2exp_si(t,t,-LOG_UPS-1);
      for(uint64_t z=0;z<n_zeros;z++)
	{
#ifdef T_MINUS_RHO
	  arb_sub(tmp1,t,Lc->zeros[side][z],prec);
#else
	  arb_div(tmp,t,Lc->zeros[side][z],prec);
	  arb_sub_ui(tmp1,tmp,1,prec);
	  arb_neg(tmp1,tmp1); // 1-t/rho
#endif
	  arb_div(Lu->values_div[side][n],Lu->values_div[side][n],tmp1,prec);
	}
    }
}

bool old_div_zeros(L_comp_t *Lc,L_family_t *Lf,L_upsample_t *Lu, L_func_t *Lfu,uint64_t side,uint64_t n_zeros, uint64_t prec)
{
  static bool init=false;
  static arb_t rho,tmp,tmp1,t,dt,ddt;
  static fmpz_t fz;
  if(!init)
    {
      init=true;
      arb_init(rho);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(t);
      arb_init(dt);
      arb_init(ddt);
      arb_div_ui(dt,Lc->one_over_A,DIV_UPSAMPLE_RATE,prec);
      arb_mul_2exp_si(ddt,dt,-1);

      fmpz_init(fz);
    }

  //printf("Not checking zero list rigorously.\n");return true;

  // where does the last zero to be output lie?
  arb_mul(tmp1,Lc->zeros[side][n_zeros-1],Lc->arb_A,prec);
  arb_ceil(tmp,tmp1,prec);
  if(arb_get_unique_fmpz(fz,tmp)==0) // spans no or more than one int
    {
      printf("Weird error in div_zeros. Skipping.\n");
      return(false);
    }
  int64_t last_n=fmpz_get_si(fz)+1; // the +1 is because we have shifted by
                                    // 1/(A*Lu->DIV_UPSAMPLE_RATE*2)
                                    // will fail if there is another zero
                                    // in the gap
  printf("Last n for values div = %ld.\n",last_n);
  fflush(stdout);
  div_zero(Lu,Lc,n_zeros,side,last_n,prec);
  printf("Divided by all known zeros.\n");
  if(quick_test(Lu,side,last_n))
    {
      printf("Quick test worked. Doing proper test now.\n");
    }
  else
    {
      printf("Don't like the look of zeros division. Skipping.\n");
      print_values_div(Lu,Lc->A,last_n,side);
      return false;
    }
  if(proper_test(Lu,side,last_n,Lc->arb_A,Lc->one_over_A,prec))
    printf("Completeness of zero list confirmed by dividing them out.\n");
  else
    {
      printf("Proper test of div_zeros failed. Skipping.\n");
      print_values_div(Lu,Lc->A,last_n,side);
      return false;
    }
  return true;
}

/*
bool new_new_div_zeros(L_comp_t *Lc,L_family_t *Lf,L_upsample_t *Lu, L_func_t *Lfu,uint64_t side,uint64_t n_zeros, uint64_t prec)
{
  static bool init=false;
  static arb_t rho,tmp,tmp1,t,dt,ddt;
  static fmpz_t fz;
  if(!init)
    {
      init=true;
      arb_init(rho);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(t);
      arb_init(dt);
      arb_init(ddt);
      arb_div_ui(dt,Lc->one_over_A,DIV_UPSAMPLE_RATE,prec);
      arb_mul_2exp_si(ddt,dt,-1);

      fmpz_init(fz);
    }

  // where does the last zero to be output lie?
  arb_mul(tmp1,Lc->zeros[side][n_zeros-1],Lc->arb_A,prec);
  arb_ceil(tmp,tmp1,prec);
  if(arb_get_unique_fmpz(fz,tmp)==0) // spans no or more than one int
    {
      printf("Weird error in div_zeros. Skipping.\n");
      return(false);
    }
  int64_t last_n=fmpz_get_si(fz)+1; // the +1 is because we have shifted by
                                    // 1/(A*Lu->DIV_UPSAMPLE_RATE*2)
                                    // will fail if there is another zero
                                    // in the gap
  printf("Last n for values div = %ld.\n",last_n);
  fflush(stdout);
  div_zero(Lu,Lc,n_zeros,side,last_n,prec);
  printf("Divided out known zeros.\n");
  fflush(stdout);
  arb_set(dt,Lc->one_over_A);
  arb_mul_2exp_si(t,dt,-LOG_UPS);
  //arb_add_error(t,t);
  for(uint64_t n=0,nn=Lu->N*Lu->stride;n<last_n;n++,nn++)
    {
      f_dash1_div(tmp,t,Lu,Lc,Lf,Lfu,side,prec);
      printf("t=");arb_printd(t,20);printf(" f(t)=");
      arb_printd(Lu->values_div[side][nn],20);printf(" f'(t)=");
      arb_printd(tmp,20);printf("\n");
      arb_add(t,t,Lc->one_over_A,prec);
    }
  return(false);
}

bool new_div_zeros(L_comp_t *Lc,L_family_t *Lf,L_upsample_t *Lu, L_func_t *Lfu,uint64_t side,uint64_t n_zeros, uint64_t prec)
{
  static bool init=false;
  static arb_t rho,tmp,tmp1,t,dt,ddt;
  static fmpz_t fz;
  if(!init)
    {
      init=true;
      arb_init(rho);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(t);
      arb_init(dt);
      arb_init(ddt);
      arb_div_ui(dt,Lc->one_over_A,DIV_UPSAMPLE_RATE,prec);
      arb_mul_2exp_si(ddt,dt,-1);

      fmpz_init(fz);
    }

  // where does the last zero to be output lie?
  arb_mul(tmp1,Lc->zeros[side][n_zeros-1],Lc->arb_A,prec);
  arb_ceil(tmp,tmp1,prec);
  if(arb_get_unique_fmpz(fz,tmp)==0) // spans no or more than one int
    {
      printf("Weird error in div_zeros. Skipping.\n");
      return(false);
    }
  int64_t last_n=fmpz_get_si(fz)+1; // the +1 is because we have shifted by
                                    // 1/(A*Lu->DIV_UPSAMPLE_RATE*2)
                                    // will fail if there is another zero
                                    // in the gap
  printf("Last n for values div = %ld.\n",last_n);
  fflush(stdout);
  div_zero(Lu,Lc,n_zeros,side,last_n,prec);
  printf("Divided out known zeros.\n");
  fflush(stdout);
  for(uint64_t n=0;n<DIV_UPSAMPLE_RATE-1;n++)
    {
      for(uint64_t nn=0;nn<Lc->NN;nn++)
	acb_zero(Lc->res[nn]);
      for(uint64_t nn=0;nn<=last_n+Lu->N*Lu->stride*2;nn++)
	arb_set(acb_realref(Lc->res[nn]),Lu->values_div[side][nn]);
      acb_convolve1(Lc->res,Lc->res,Lu->upsample_expsincs[n],Lc->NN,Lc->ww,prec);
      for(uint64_t nn=0;nn<last_n;nn++)
	{
	  printf("%lu ",nn);
	  arb_printd(acb_realref(Lc->res[nn+Lu->N*Lu->stride]),20);
	  printf("\n");
	}
      acb_t *div=Lc->res+Lu->N*Lu->stride;
      if(arb_contains_zero(acb_realref(div[0])))
	{
	  printf("f(0)/P(t) contains zero. Skipping.\n");
	  return(false);
	}
      if(arb_is_positive(acb_realref(div[0])))
	{
	  for(uint64_t nn=1;nn<last_n;nn++)
	    if(!arb_is_positive(acb_realref(div[nn])))
	      {
		printf("Sign change after dividing by zeros on sub-interval %lu at datum %lu. Skipping.\n",n,nn);
		return false;
	      }
	  continue;
	}
      for(uint64_t nn=1;nn<last_n;nn++)
	if(!arb_is_negative(acb_realref(div[nn])))
	  {
	    printf("Sign change after dividing by zeros on sub-interval %lu at datum %lu. Skipping.\n",n,nn);
	    return false;
	  }
    }
  printf("Completeness of list of zeros confirmed by dividing them out.\n");
  return true;
}
*/
