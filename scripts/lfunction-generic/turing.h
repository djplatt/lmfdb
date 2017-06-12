
// int im loggamma(sigma+It) t=a..b
bool im_int(arb_t res, arb_t sigma, arb_t a, arb_t b, uint64_t prec)
{
  //printf("Integrating with sig=");arb_printd(sigma,10);printf(" a=");arb_printd(a,10);printf(" and b=");arb_printd(b,10);printf("\n");
  static arb_t a2,b2,s2,b2s2,a2s2,tmp1,tmp2,tmp3,shalf;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(a2);
      arb_init(b2);
      arb_init(s2);
      arb_init(b2s2);
      arb_init(a2s2);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(shalf);
    }
  arb_set_d(tmp1,-0.5);
  arb_add(shalf,sigma,tmp1,prec); // sigma-1/2

  arb_mul(s2,sigma,sigma,prec);
  arb_mul(a2,a,a,prec);
  arb_mul(b2,b,b,prec);
  arb_add(tmp1,a2,s2,prec);
  arb_log(a2s2,tmp1,prec);
  arb_add(tmp1,b2,s2,prec);
  arb_log(b2s2,tmp1,prec);
  arb_sub(tmp1,s2,sigma,prec);
  arb_sub(tmp2,tmp1,b2,prec);
  arb_mul(tmp3,tmp2,b2s2,prec); // log(b^2+sig^2)(sig^2-b^2-sig)
  arb_sub(tmp2,tmp1,a2,prec);
  arb_mul(tmp1,tmp2,a2s2,prec); // log(a^2+sig^2)(sig^2-a^2-sig)
  arb_sub(tmp2,tmp3,tmp1,prec);
  arb_sub(tmp1,b2,a2,prec);
  arb_mul_ui(tmp3,tmp1,3,prec);
  arb_add(res,tmp2,tmp3,prec);
  arb_mul_2exp_si(res,res,-2);

  arb_div(tmp1,b,sigma,prec);
  arb_atan(tmp2,tmp1,prec);
  arb_mul(tmp1,tmp2,b,prec); // b atan(b/sig)
  arb_div(tmp2,a,sigma,prec);
  arb_atan(tmp3,tmp2,prec);
  arb_mul(tmp2,tmp3,a,prec); // a atan(a/sig)
  arb_sub(tmp3,tmp2,tmp1,prec); // a atan(a/sigma)-b atan(b/sigma)
  arb_mul(tmp2,tmp3,shalf,prec);

  arb_add(res,res,tmp2,prec);
  arb_neg(res,res);
  arb_div(tmp1,b,a,prec);
  arb_log(tmp2,tmp1,prec);
  arb_mul_2exp_si(tmp2,tmp2,-3);
  arb_add_error(res,tmp2);
  return(true);
} /* im_int */


// return upb for 2\int\limits_{t_0}^{t0+h} S(t) \dif t
// see Th 4.6
bool St_int(arb_t res, arb_t h, arb_t t0, L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static acb_t s;
  static arb_t Q1,Q2,l2_half,pi,rc_theta_etc,tmp,tmp1;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      acb_init(s);
      arb_init(pi);
      arb_set_d(acb_realref(s),1.5);
      arb_init(Q1);
      arb_init(Q2);
      arb_init(l2_half);
      arb_log_ui(Q1,2,prec);
      arb_set_d(Q2,-0.5);
      arb_add(l2_half,Q1,Q2,prec);
      arb_init(rc_theta_etc);
      arb_set_d(Q2,5.65056); // c_theta < 5.65055 in ARB Th 4.6
      arb_sqrt_ui(Q1,2,prec);
      arb_mul_ui(pi,Q1,Lf->X-5,prec);
      arb_inv(Q1,pi,prec); // 1.0/(sqrt(2)(X-5)
      arb_add(pi,Q1,Q2,prec); 
      arb_mul_ui(rc_theta_etc,pi,Lf->r,prec); // c_\theta r+r/(\sqrt{2}(X-5))
      arb_const_pi(pi,prec);
    }
  arb_set(acb_imagref(s),t0); // 3/2+it0
  logQ(Q2,s,Lf,Lfu->N,prec);
  arb_mul(tmp,Q2,l2_half,prec);
  arb_add(acb_imagref(s),acb_imagref(s),h,prec); // 3/2+it1
  logQ(Q1,s,Lf,Lfu->N,prec);
  arb_mul_2exp_si(Q1,Q1,-2);
  arb_add(tmp1,tmp,Q1,prec);
  arb_add(tmp,tmp1,rc_theta_etc,prec);
  arb_div(res,tmp,pi,prec);
  arb_mul_2exp_si(res,res,1);
  return(true);
} /* St_int */

void print_turing_region(L_comp_t *Lc, L_upsample_t *Ls,uint64_t side)
{
  printf("Turing region for side %lu\n",side);
  uint64_t nn=Ls->N*Ls->stride+Lc->NN/OUTPUT_RATIO;
  for(uint64_t n=0;n<=Lc->NN/TURING_RATIO;n++,nn++)
    {
      arb_printd(Ls->values[side][nn],10);
      printf("\n");
    }
  printf(".............End.\n");
}

sign_t sign(arb_t);
direction_t direction(arb_t,arb_t,uint64_t);
bool upsample_stride(arb_ptr, arb_ptr, L_upsample_t *, L_comp_t *, L_family_t *, L_func_t *, uint64_t , uint64_t );


bool stat_point_no_isolate(arb_t z1, arb_t z2, arb_t tt0, arb_t tt1, arb_t tt2, arb_t ff0, arb_t ff1, arb_t ff2, L_comp_t *Lc, L_upsample_t *Lu, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static bool init=false;
  static arb_t t0,t1,f0,f1,t2,f2,t01,f01,t12,f12;
  if(!init)
    {
      init=true;
      arb_init(t0);
      arb_init(t1);
      arb_init(t2);
      arb_init(f0);
      arb_init(f1);
      arb_init(f2);
      arb_init(t01);
      arb_init(t12);
      arb_init(f01);
      arb_init(f12);
    }

  arb_set(t0,tt0);
  arb_set(t1,tt1);
  arb_set(t2,tt2);
  arb_set(f0,ff0);
  arb_set(f1,ff1);
  arb_set(f2,ff2);
  sign_t s=sign(f0);

  while(true)
    {

      arb_printd(t0,20);printf(" ");arb_printd(f0,20);printf("\n");
      arb_printd(t1,20);printf(" ");arb_printd(f1,20);printf("\n");
      arb_printd(t2,20);printf(" ");arb_printd(f2,20);printf("\n\n");
      
      arb_add(t01,t0,t1,prec);
      arb_mul_2exp_si(t01,t01,-1);
      upsample_stride(f01,t01,Lu,Lc,Lf,Lfu,side,prec);
      sign_t s01=sign(f01);
      if(s01==UNK)
	{
	  printf("Indeterminate in stat_point.\n");
	  return(false);
	}
      if(s01!=s) // z1 is between t0 and t01, z2 t01 and t1
	{
	  arb_union(z1,t0,t01,prec);
	  arb_union(z2,t01,t1,prec);
	  return(true);
	}
      direction_t left=direction(f0,f01,prec);
      direction_t right=direction(f01,f1,prec);
      if((left==UNK)||(right==UNK))
	{
	  printf("Unknown direction in stat point.\n");
	  return(false);
	}
      if((left!=right)&&(((s==POS)&&(left==DOWN))||((s==NEG)&&(left==UP)))) // stat pt is between t0,t01,t1
	{
	  arb_set(t2,t1);
	  arb_set(f2,f1);
	  arb_set(t1,t01);
	  arb_set(f1,f01);
	  continue;
	}

      arb_add(t12,t1,t2,prec);
      arb_mul_2exp_si(t12,t12,-1);
      upsample_stride(f12,t12,Lu,Lc,Lf,Lfu,side,prec);
      sign_t s12=sign(f12);
      if(s12==UNK)
	{
	  printf("Indeterminate in stat_point.\n");
	  return(false);
	}
      if(s12!=s)
	{
	  arb_union(z1,t1,t12,prec);
	  arb_union(z2,t12,t2,prec);
	  return(true);
	}
      left=direction(f1,f12,prec);
      right=direction(f12,f2,prec);
      if((left==UNK)||(right==UNK))
	{
	  printf("Unknown direction in stat point.\n");
	  return(false);
	}
      if((left!=right)&&(((s==POS)&&(left==DOWN))||((s==NEG)&&(left==UP)))) // stat pt is t1,t12,t2
	{
	  arb_set(t0,t1);
	  arb_set(f0,f1);
	  arb_set(t1,t12);
	  arb_set(f1,f12);
	  continue;
	}
      else // stat pt is t01 t1 t12
	{
	  arb_set(t0,t01);
	  arb_set(f0,f01);
	  arb_set(t2,t12);
	  arb_set(f2,f12);
	  continue;
	}
    }
  return(true);
}


bool turing_int(arb_t res, arb_t h, L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp,ptr,err,t0,t1,t2,z1,z2;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(ptr);
      arb_init(err);
      arb_init(t0);
      arb_init(t1);
      arb_init(t2);
      arb_init(z1);
      arb_init(z2);
      arb_zero(err);
      arb_mul_2exp_si(tmp,Lc->one_over_A,-1);
      arb_add_error(err,tmp);
    }

  //print_turing_region(Lc,Ls,side);
  arb_zero(res);
  uint64_t n=0,nn=Ls->N*Ls->stride*2+Lc->NN/OUTPUT_RATIO;
  uint8_t this_sign,last_sign=sign(Ls->values[side][nn]);
  n++;nn++;
  while((last_sign==UNK)&&(n<Lc->NN/TURING_RATIO))
    {
      last_sign=sign(Ls->values[side][nn]);
      n++;nn++;
    }
  /*
    {
      printf("Indeterminate in turing_int. ");
      arb_printd(Ls->values[side][nn],10);
      printf("\n");
      return(false);
    }
  */
  if(last_sign!=UNK) // whole turing region was indeterminate!
    while(true)
      {
	this_sign=sign(Ls->values[side][nn]);
	while((this_sign==UNK)&&(n<Lc->NN/TURING_RATIO))
	  {
	    n++;nn++;
	    this_sign=sign(Ls->values[side][nn]);
	  }
	if(this_sign==UNK) // tail of Turing region was indeterminate
	  break;
	/*
	  if(this_sign==UNK)
	  {
	  printf("Indeterminate in turing_int.\n");
	  arb_printd(Ls->values[side][nn],10);
	  printf("\n");
	  return(false);
	  }
	*/
	if(this_sign!=last_sign)
	  {
	    //printf("Turing Zero found between %lu and %lu.\n",n-1,n);
	    arb_set_d(tmp,(double)n-0.5);
	    arb_mul(ptr,Lc->one_over_A,tmp,prec);
	    //div_zero(Ls,Lc,ptr,side,prec);
	    arb_add(tmp,ptr,err,prec);
	    arb_sub(ptr,h,tmp,prec);
	    arb_add(res,res,ptr,prec);
	    last_sign=this_sign;
	  }
	n++;nn++;
	if(n>Lc->NN/TURING_RATIO)
	  break;
      }
  //return(true); // don't do stat points bit
  printf("Int before stat points is ");arb_printd(res,10);printf("\n");
  nn=Ls->N*Ls->stride*2+Lc->NN/OUTPUT_RATIO;
  direction_t left=direction(Ls->values[side][nn],Ls->values[side][nn+1],prec);
  direction_t right=direction(Ls->values[side][nn+1],Ls->values[side][nn+2],prec);
  nn++;
  n=1;
  while(true)
    {
      if((left==UNK)||(right==UNK)) // data getting too noisy to use now
	{printf("Turing data too noisy for stat points.\n");return(true);}
      if(left!=right) // change in direction
	{
	  sign_t mid_sign=sign(Ls->values[side][nn]);
	  if((left==UP)&&(mid_sign==NEG))
	    {
	      if(!stat_point(z1,z2,nn-Ls->N*Ls->stride*2,Lc,Ls,Lf,Lfu,side,prec,false))
		{
		  printf("Failed to resolve stat point in Turing zone at:-\n");
		  printf("Continuing anyway.\n");
		}
	      else
		{
		  printf("Stat points found zeros at ");
		  arb_printd(z1,10); printf(" and ");
		  arb_printd(z2,10); printf("\n");
		  arb_sub(tmp,h,z1,prec);
		  arb_add(res,res,tmp,prec);
		  arb_sub(tmp,h,z2,prec);
		  arb_add(res,res,tmp,prec);
		  // now rebase z1,z2 to start from Turing region
		  arb_mul_ui(tmp,Lc->one_over_A,2*Lc->NN/OUTPUT_RATIO,prec);
		  arb_add(res,res,tmp,prec);
		}
	    }
	  else
	    {
	      if((left==DOWN)&&(mid_sign==POS))
		{
		  if(!stat_point(z1,z2,nn-Ls->N*Ls->stride*2,Lc,Ls,Lf,Lfu,side,prec,false))
		    {
		      printf("Failed to resolve stat point in Turing zone at:-\n");
		      printf("Continuing anyway.\n");
		    }
		  else
		    {
		      printf("Stat points found zeros at ");
		      arb_printd(z1,10); printf(" and ");
		      arb_printd(z2,10); printf("\n");
		      arb_sub(tmp,h,z1,prec);
		      arb_add(res,res,tmp,prec);
		      arb_sub(tmp,h,z2,prec);
		      arb_add(res,res,tmp,prec);
		      // now rebase z1,z2 to start from Turing region
		      arb_mul_ui(tmp,Lc->one_over_A,2*Lc->NN/OUTPUT_RATIO,prec);
		      arb_add(res,res,tmp,prec);
		    }
		}
	    }
	}
      nn++;n++;
      if(n==Lc->NN/TURING_RATIO)
	break;
      left=right;
      right=direction(Ls->values[side][nn],Ls->values[side][nn+1],prec);
    }
  return(true);
}


bool turing_count(arb_t res, L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static arb_t tint,tmp,tmp1,h,t0,sint,pi,rlogpi,t0hbit;
  if(!init)
    {
      init=true;
      arb_init(tint);
      arb_init(sint);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(h);
      arb_init(t0);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_div_ui(h,Lf->B,TURING_RATIO,prec);
      printf("Turing h set to ");arb_printd(h,10);printf("\n");
      arb_div_ui(t0,Lf->B,OUTPUT_RATIO,prec);
      printf("Turing t0 set to ");arb_printd(t0,10);printf("\n");
      arb_init(rlogpi);
      arb_log(tmp,pi,prec);
      arb_mul_ui(rlogpi,tmp,Lf->r,prec);
      arb_init(t0hbit);
      arb_mul(tmp,t0,h,prec);
      arb_mul_2exp_si(tmp,tmp,1);
      arb_mul(tmp1,h,h,prec);
      arb_add(sint,tmp1,tmp,prec);
      arb_div(t0hbit,sint,pi,prec); 
      arb_mul_2exp_si(t0hbit,t0hbit,-1); // (2t0h+h^2)/2Pi
    }
  if(!St_int(sint,h,t0,Lc,Lf,Ls,Lfu,prec))
    return(0);
  printf("St_int returned ");arb_printd(sint,10);printf("\n");

  if(!turing_int(tmp1,h,Lc,Lf,Ls,Lfu,0,prec))
    return(0);
  if(Lfu->self_dual_p)
    arb_mul_2exp_si(tint,tmp1,1);
  else
    {
      if(!turing_int(tmp,h,Lc,Lf,Ls,Lfu,1,prec))
	return(0);
      arb_add(tint,tmp1,tmp,prec);
    }

  printf("Turing int [*] = ");arb_printd(tint,20);printf("\n");

  arb_div(tmp1,Lf->imint,pi,prec);
  arb_add(tmp,tmp1,sint,prec);
  arb_sub(sint,tmp,tint,prec); //  

  arb_log_ui(tmp1,Lfu->N,prec);
  arb_sub(tmp,tmp1,rlogpi,prec);
  arb_mul(tmp1,tmp,t0hbit,prec);

  arb_add(tmp,tmp1,sint,prec);
  arb_div(res,tmp,h,prec);

  printf("Turing bound will try to return ");arb_printd(res,10);printf("\n");
  return(true);
}