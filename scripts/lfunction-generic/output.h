// exp(Pi t/4)
void exp_pi_t_4(arb_t res, arb_t t, uint64_t prec)
{
  static bool init=false;
  static arb_t pi4,tmp;
  if(!init)
    {
      init=true;
      arb_init(pi4);
      arb_init(tmp);
      arb_const_pi(pi4,prec);
      arb_mul_2exp_si(pi4,pi4,-2); // Pi/4
    }
  arb_mul(tmp,t,pi4,prec);
  arb_exp(res,tmp,prec);
} /* exp_pi_t_4 */

// return f(t) * exp(-Pi*t/4) /|gamma(1/2+it)| 
// where gamma is per top col 2 pg 387 without N^1/2(s-1/2)
void normalise(arb_t res,arb_t ft,arb_t t,L_family_t *Lf,uint64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,tmp2;
  static acb_t s;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(s);
      arb_set_d(acb_realref(s),0.5);
    }
  exp_pi_t_4(tmp,t,prec);
  arb_set(acb_imagref(s),t);
  abs_gamma(tmp1,s,Lf,prec);
  arb_mul(tmp2,tmp,tmp1,prec);
  arb_div(res,ft,tmp2,prec);
} /* normalise */
  
// use this for zeros as we snap them to check for sign change
bool output_arb_no_snap(arb_t x, FILE *outfile, uint64_t OP_ACC)
{
  static bool init=false;
  static fmpz_t fz;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      fmpz_init(fz);
      arb_init(tmp);
    }

  arb_mul_2exp_si(tmp,x,OP_ACC);
  if(arb_get_unique_fmpz(fz,tmp)==0) // spans no or more than one int
    return(false);
  fmpz_fprint(outfile,fz);
  return(true);
}

// snap a point to Z/2^{OP_ACC+1} and write it
bool output_arb(arb_t x, FILE *outfile, uint64_t OP_ACC, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  arb_mul_2exp_si(tmp,x,OP_ACC+1);
  arb_add_ui(tmp1,tmp,1,prec);
  arb_mul_2exp_si(tmp1,tmp1,-1);
  arb_floor(tmp2,tmp1,prec);
  arb_get_mid_arb(tmp,tmp2);
  arb_mul_2exp_si(tmp,tmp,-OP_ACC);
  arb_sub(tmp1,tmp,x,prec);
  arb_abs(tmp1,tmp1); // should be less than 2^{-OP_ACC+1}
  //arb_printd(tmp1,20);printf("\n");
  arb_mul_2exp_si(tmp1,tmp1,OP_ACC+1);
  //arb_printd(tmp1,20);printf("\n");
  arb_sub_ui(tmp2,tmp1,1,prec);
  //arb_printd(tmp2,20);printf("\n");
  if(!arb_is_negative(tmp2))
    return(false);
  return(output_arb_no_snap(tmp,outfile,OP_ACC));
}

double output_zeros(arb_t *zeros, FILE *outfile,uint64_t prec, uint64_t OP_ACC)
{
  static bool init=false;
  static arf_t tmp;
  static arb_t tmp1;
  if(!init)
    {
      init=true;
      arf_init(tmp);
      arb_init(tmp1);
    }
  double zd;
  for(uint64_t z=0;z<OUTPUT_ZEROS;z++)
    {
      printf("Zero %lu: ",z);arb_printd(zeros[z],100);printf("\n");
      arb_get_mid_arb(tmp1,zeros[z]);
      arb_get_lbound_arf(tmp,tmp1,prec);
      zd=arf_get_d(tmp,ARF_RND_NEAR);
      fprintf(outfile,"%20.18e ",zd);
      if(!output_arb_no_snap(zeros[z],outfile,OP_ACC))
	{
	  arb_mul_2exp_si(zeros[z],zeros[z],OP_ACC+1);
	  printf("zero * 2^%lu is ",OP_ACC+1);arb_printd(zeros[z],100);printf("\n");
	  printf("Zero number %lu not precise enough. Skipping.\n",z);
	  return(0.0);
	}
      fprintf(outfile,"\n");
    }
  //printf("last zero is near %10.8e\n",zd);
  return(zd);
}

// print a double precision estimate of pt
void output_point(arb_t pt, FILE *outfile,uint64_t prec)
{
  static bool init=false;
  static arf_t tmp;
  static arb_t tmp1;
  if(!init)
    {
      init=true;
      arf_init(tmp);
      arb_init(tmp1);
    }

  arb_get_mid_arb(tmp1,pt);
  arb_get_lbound_arf(tmp,tmp1,prec);
  fprintf(outfile,"%20.18e\n",arf_get_d(tmp,ARF_RND_NEAR));
}

bool output_epsilon(acb_t eps, FILE *outfile, arb_t err, uint64_t op_acc, uint64_t prec)
{
  if(!output_arb(acb_realref(eps),outfile,op_acc,prec))
    {
      printf("Error outputting real part of epsilon. Skipping.\n");
      return(false);
    }
  fprintf(outfile," ");
  if(!output_arb(acb_imagref(eps),outfile,op_acc,prec))
    {
      printf("Error outputting imaginary part of epsilon. Skipping.\n");
      return(false);
    }
  fprintf(outfile,"\n");
  return(true);
}


bool write_data(const char *out_fname,L_family_t *L_family,L_func_t *L_func,
		L_comp_t *L_comp,L_upsample_t *L_upsample,uint64_t side, uint64_t prec)
{
  static bool init=false;
  static arb_t t,delta_t,tmp;
  if(!init)
    {
      init=true;
      arb_init(t);arb_init(delta_t);arb_init(tmp);
    }
      FILE *outfile=fopen(out_fname,"w");
      if(!outfile)
	{
	  printf("Failed to open %s for output. Skipping.\n",out_fname);
	  return false;
	}
      fprintf(outfile,"%lu\n",L_comp->OP_ACC+1);
      if(!output_epsilon(L_func->epsilon_sqr,outfile,L_comp->zero_prec,L_comp->OP_ACC,prec))
	{
	  fclose(outfile);
	  remove(out_fname);
	  return false;
	}
      arb_zero(t);
      normalise(tmp,L_upsample->values[side][L_upsample->N*L_upsample->stride*2],t,L_family,prec);
      printf("L(1/2)=");arb_printd(tmp,10);printf("\n");
      if(!output_arb(tmp,outfile,L_comp->OP_ACC,prec))
	{
	  printf("Error outputting L(1/2). Skipping.\n");
	  fclose(outfile);
	  remove(out_fname);
	  return false;
	}
      fprintf(outfile,"\n%lu\n",OUTPUT_ZEROS);
      double last_zero=output_zeros(L_comp->zeros[side],outfile,prec,L_comp->OP_ACC);
      if(last_zero==0.0)
	{
	  printf("Last zero was at 0.0. Skipping.\n");
	  fclose(outfile);
	  remove(out_fname);
	  return false;
	}

      
      double dt=1.0/L_comp->A; // spacing of points
      uint64_t step=1;
      while(dt*NUM_OUTPUT_POINTS*step<last_zero)
	step++;
      step+=step/2;
      if(step*NUM_OUTPUT_POINTS>L_comp->NN/OUTPUT_RATIO)
	step=L_comp->NN/(OUTPUT_RATIO*NUM_OUTPUT_POINTS);
      dt*=step;
      printf("Outputing every %lu'th point.\n",step);
      arb_mul_ui(delta_t,L_comp->one_over_A,step,prec);

      fprintf(outfile,"%10.8e\n",dt);
      fprintf(outfile,"%lu\n",NUM_OUTPUT_POINTS+1);
      for(uint64_t n=L_upsample->N*L_upsample->stride*2,nn=0;nn<=NUM_OUTPUT_POINTS;nn++,n+=step)
	{
	  normalise(tmp,L_upsample->values[side][n],t,L_family,prec);
	  output_point(tmp,outfile,prec);
	  arb_add(t,t,delta_t,prec);
	}
      fclose(outfile);

      return true;
}


