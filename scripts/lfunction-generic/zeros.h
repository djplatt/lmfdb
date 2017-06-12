

void print_direction(direction_t d)
{
  printf("%s",(d==UP ? "+" : "-"));
}

#define dir_t uint8_t

direction_t direction(arb_t a, arb_t b, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_sub(tmp,a,b,prec);
  if(arb_contains_zero(tmp))
    return(UNK);
  if(arb_is_negative(tmp)) // b>a
    return(UP);
  return(DOWN); // b<a
}

sign_t sign(arb_t x)
{
  if(arb_contains_zero(x))
    return(UNK);
  if(arb_is_positive(x))
    return(POS);
  return(NEG);
}

// use simple binary chop
bool isolate_zero(arb_t res, arb_t tt0, arb_t tt1, arb_t ff0, arb_t ff1, int8_t s0, int8_t s1, L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t side, uint64_t prec, uint64_t op_acc)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,t0,t1,f0,f1;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(t0);
      arb_init(t1);
      arb_init(f0);
      arb_init(f1);

    }
  arb_set(t0,tt0);
  arb_set(t1,tt1);
  arb_set(f0,ff0);
  arb_set(f1,ff1);
  uint64_t iter=0;
  while(true)
    {
      arb_add(tmp1,t0,t1,prec);
      arb_mul_2exp_si(tmp1,tmp1,-1);

      // see if we are close enough yet
      arb_sub(tmp2,t1,t0,prec);
      arb_mul_2exp_si(tmp2,tmp2,op_acc+2); // +1 wasn't tight enough
      arb_sub_ui(tmp3,tmp2,1,prec);
      if(arb_is_negative(tmp3))
	{
	  printf("Zero resolved in binary chop.\n");
	  arb_set(res,tmp1);
	  return true;
	}

      iter++;
      /*      
      printf("t0 = ");arb_printd(t0,40);
      printf("\nt1 = ");arb_printd(t1,40);
      printf("\nf0 = ");arb_printd(f0,40);
      printf("\nf1 = ");arb_printd(f1,40);
      printf("\n");
      
      printf("\nupsampling at ");arb_printd(tmp1,40);
      */
      upsample_stride(tmp2,tmp1,Ls,Lc,Lf,Lfu,side,prec);
      /*
      printf("\ngot ");arb_printd(tmp2,40);printf("\n");

      printf("t0 = ");arb_printd(t0,40);
      printf("\nt1 = ");arb_printd(t1,40);
      printf("\nf0 = ");arb_printd(f0,40);
      printf("\nf1 = ");arb_printd(f1,40);
      printf("\n");
      */

      sign_t new_sign=sign(tmp2);
      if(new_sign==UNK)
	{
	  printf("Indeterminate sign in upsampling on iteration %lu. ",iter);
	  arb_printd(tmp2,10);
	  printf(" ");
	  arb_printd(tt0,10);
	  printf(" ");
	  arb_printd(ff0,10);
	  printf("\n");
	  return(false);
	}

      if(new_sign!=s0) // change is between t0 and tmp1
	{
	  arb_set(t1,tmp1);
	  arb_set(f1,tmp2);
	}
      else // change is between tmp1 and t1
	{
	  arb_set(t0,tmp1);
	  arb_set(f0,tmp2);
	}
    }
}
      
// pass it an exact x=N*2^{-OP_ACC}, N an integer
// check rigorously that there is a sign change between x\pm2^{-OP_ACC-1}
// if so, then zero is within 2^{-OP_ACC-1} of x
// this does Table Maker's dilemma
bool confirm_zero(arb_t x, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec,arb_t zero_prec, uint64_t side)
{
  static bool init=false;
  static arb_t fleft,fright,tmp;
  if(!init)
    {
      init=true;
      arb_init(fleft);
      arb_init(fright);
      arb_init(tmp);
    }
  //printf("Confirming zero at ");arb_printd(x,50);printf("\n");
  arb_sub(tmp,x,zero_prec,prec);
  if(arb_is_negative(tmp))
    {
      side^=1;
      arb_neg(tmp,tmp);
      //printf("computing F(");arb_printd(tmp,10);printf(") on side %lu\n",side);
      if(!upsample_stride(fleft, tmp, Lu, Lc, Lf, Lfu, side, prec))
	return(false);
      //printf("Left = ");arb_printd(fleft,20);printf("\n");
      side^=1;
    }
  else
    {
      //printf("computing F(");arb_printd(tmp,10);printf(") on side %lu\n",side);
      if(!upsample_stride(fleft, tmp, Lu, Lc, Lf, Lfu, side, prec))
	return(false);
      //printf("Left = ");arb_printd(fleft,20);printf("\n");
    }
  arb_add(tmp,x,zero_prec,prec);
  //printf("computing F(");arb_printd(tmp,10);printf(") on side %lu\n",side);
  if(!upsample_stride(fright, tmp, Lu, Lc, Lf, Lfu, side, prec))
    return(false);
  //printf("Right = ");arb_printd(fright,20);printf("\n");
  if((sign(fleft)&sign(fright))==0)
    return true;
  printf("Problem in confirm zero.\nf(t-delta)=");
  arb_printd(fleft,20);
  printf("\nf(t+delta=");
  arb_printd(fright,20);
  return false;

}

// called if we don't find enough zeros the obvious way
// locate maxima(minima) below(above) the line.
// go find the two zeros that lie therein.
uint64_t do_stat_points(L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  return(0);
  uint64_t count=0; // number of zeros found
  dir_t left,right;
  uint64_t mid_ptr=1;
  left=direction(Ls->values[side][0],Ls->values[side][1],prec);
  right=direction(Ls->values[side][1],Ls->values[side][2],prec);
  while(true)
    {
      if((left==UNK)||(right==UNK))
	{
	  printf("Unknown direction in stat points. Skipping.\n");
	  return(BAD_64);
	}
      if(left!=right) // change in direction
	{
	  sign_t mid_sign=sign(Ls->values[side][mid_ptr]);
	  if((left==UP)&&(mid_sign==NEG))
	    {
	      printf("Looks like a stationary point. Will seek zeros here.\n");
	      count+=2;
	    }
	  else
	    {
	      if((left==DOWN)&&(mid_sign==POS))
		{
		  printf("Looks like a stationary point. Will seek zeros here.\n");
		  count+=2;
		}
	    }
	}
      mid_ptr++;
      if(mid_ptr==Lc->NN/OUTPUT_RATIO)
	return(count);
      left=right;
      right=direction(Ls->values[side][mid_ptr],Ls->values[side][mid_ptr+1],prec);
    }
  return(BAD_64);
}

// snap rho to N*2^{-OP_ACC} where N is an integer
bool snap_point(arb_t res, arb_t rho, uint64_t prec, uint64_t op_acc)
{
  static bool init=false;
  static arb_t tmp;
  static mpfr_t a,b;
  static mpz_t z;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      mpfr_init2(a,prec);
      mpfr_init2(b,prec);
      mpz_init(z);
    }
  arb_mul_2exp_si(tmp,rho,op_acc);
  arb_get_interval_mpfr(a,b,tmp);
  mpfr_add(a,a,b,MPFR_RNDN);
  mpfr_div_ui(b,a,2,MPFR_RNDN);
  mpfr_get_z(z,b,MPFR_RNDN);
  mpfr_set_z(a,z,MPFR_RNDN);
  arb_set_interval_mpfr(tmp,a,a,prec);
  arb_mul_2exp_si(res,tmp,-op_acc);
  return(true);
}

uint64_t guess_rank(L_func_t *Lf, L_upsample_t *Ls, uint64_t side, uint64_t prec) 
{
  static bool init=false;
  static arb_t t;
  double ratio;

  if (!init) 
    {
      init=true;
      arb_init(t);
    }
  arb_div(t,Ls->values_off[side][2],Ls->values_off[side][1],prec);
  ratio = arf_get_d(arb_midref(t),ARF_RND_NEAR);
  printf("Ratio = %10.8e\n",ratio);
  printf("eps^2 = ");arb_printd(acb_realref(Lf->epsilon_sqr),20);printf("\n");
  if (arb_is_negative(acb_realref(Lf->epsilon_sqr)))
    {
      if (ratio <= 1) return 1;
      return 1+2*(int)floor(log(ratio/2)/log(4.0)+0.5);
    } 
  else 
    {
      if (ratio <= 1) return 0;
      return 2*(int)floor(log(ratio)/log(4.0)+0.5);
    }
}

int64_t find_zeros(L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, L_func_t *Lfu, uint64_t side, uint64_t prec)
{
  uint64_t op_acc=Lc->OP_ACC;
  int64_t count_factor=1;
  uint64_t count=0;
  static bool init=false;
  static arb_t tmp,tmp1,tmp2,tmp3,t0,zero_prec,z1,z2;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(t0);
      arb_init(zero_prec);
      arb_init(z1);
      arb_init(z2);
    }
  arb_set_ui(zero_prec,1);
  arb_mul_2exp_si(zero_prec,zero_prec,-op_acc-1);

  uint64_t n=0;
  sign_t this_sign,last_sign=sign(Ls->values_off[side][0]);
  direction_t this_dir,last_dir;
  last_dir=direction(Ls->values_off[side][0],Ls->values_off[side][1],prec);
  if(last_dir==UNK)
    {
      printf("Unknown direction at start of data.\n");
      return 0;
    }

  if((side==0)&&(Lf->ltype==GENUS3))
    {
      if(last_sign==UNK)
	{
	  Lfu->rank=guess_rank(Lfu,Ls,0,prec);
	  printf("Guessed analytic rank to be %lu\n",Lfu->rank);
	}
      else
	Lfu->rank=0;
    }

  // Genus2 and Elliptic curve data gives us the rank
  if(((Lf->ltype==GENUS3)||(Lf->ltype==GENUS2)||(Lf->ltype==ELLIPTIC))&&(Lfu->rank>0))
    {
      if(last_sign!=UNK)
	{
	  printf("Elliptic/Genus 2/Genus 3 curve of rank %lu but central value not 0.\n",Lfu->rank);
	  return 0;
	}
      count=Lfu->rank;
      if(side==0) count_factor=-1;
      for(uint64_t i=0;i<Lfu->rank;i++)
	arb_zero(Lc->zeros[side][i]);
      n++;
      last_sign=sign(Ls->values_off[side][n]);
    }
  else
    {
      if(last_sign==UNK) // we'll tolerate an UNK at the centre
	{
	  if(side==0)
	    {
	      printf("We strongly suspect a central zero.\n");
	      arb_zero(Lc->zeros[0][0]);
	      arb_zero(Lc->zeros[1][0]);
	      sign_t tmp_sign=sign(Ls->values_off[side][1]);
	      count=1;
	      count_factor=-1; // negate it to indicate we have a central zero
	      arb_zero(tmp);
	      if(!confirm_zero(tmp,Ls,Lc,Lf,Lfu,prec,zero_prec,0))
		{
		  printf("Looks like a double zero. Please check.\n");
		  count++;
		  arb_zero(Lc->zeros[0][1]);
		  arb_zero(Lc->zeros[1][1]);
		}
	      last_sign=tmp_sign;
	    }
	  else // just count it, side 0 checked it
	    {
	      if(arb_is_zero(Lc->zeros[0][1])) // it was a double zero
		count=2;
	      else
		count=1;
	      last_sign=sign(Ls->values_off[side][1]);
	    }
	  n=1;
	}
    }

  if(last_sign==UNK) // that was your last chance
    {
      printf("Indeterminate sign in output region at n=1.\n");
      return(0);
    }

  while(true)
    {
      n++;
      if(n>Lc->NN/OUTPUT_RATIO)
	return(count*count_factor);
      this_sign=sign(Ls->values_off[side][n]);
      if(this_sign==UNK)
	{
	  printf("Indeterminate sign in output region.\n");
	  return(0);
	}
      if(this_sign!=last_sign) // found a zero between n and n-1
	{
	  if(count<OUTPUT_ZEROS) // need to isolate it properly
	    {
	      arb_mul_ui(tmp1,Lc->one_over_A,n-1,prec);
	      arb_mul_ui(tmp2,Lc->one_over_A,n,prec);
	      //printf("zero found between ");arb_printd(tmp1,20);
	      //printf(" and ");arb_printd(tmp2,20);printf("\n");
	      arb_add(tmp3,tmp1,tmp2,prec);
	      arb_mul_2exp_si(tmp3,tmp3,-1);
	      if(!newton(tmp,tmp3,Ls,Lc,Lf,Lfu,side,prec)||!snap_point(tmp3,tmp,prec,op_acc)||!confirm_zero(tmp3,Ls,Lc,Lf,Lfu,prec,zero_prec,side))
		{
		  printf("Resorting to binary chop.\n");
		  if(!isolate_zero(tmp,tmp1,tmp2,Ls->values_off[side][n-1],Ls->values_off[side][n],last_sign,this_sign,Lc,Lf,Ls,Lfu,side,prec,op_acc))
		    return(0);
		  if(!snap_point(tmp3,tmp,prec,op_acc)) return 0;
		  if(!confirm_zero(tmp3,Ls,Lc,Lf,Lfu,prec,zero_prec,side))
		    {
		      printf("Failed to confirm zero near ");
		      arb_printd(tmp3,100);
		      printf("\n");
		      return(0);
		    }
		}
	      arb_set(Lc->zeros[side][count++],tmp3);
	    }
	  else
	    count++;
	  if((count==OUTPUT_ZEROS)&&(Lf->r>MAX_TURING_R))
	    return count*count_factor;
	  last_sign=this_sign;
	  last_dir=direction(Ls->values_off[side][n-1],Ls->values_off[side][n],prec);
	  if(last_dir==UNK)
	    {
	      printf("Unknown direction in find_zeros.\n");
	      return 0;
	    }
	  continue;
	}

      this_dir=direction(Ls->values_off[side][n-1],Ls->values_off[side][n],prec);
      if(this_dir==UNK)
	{
	  printf("Unknown direction in find_zeros.\n");
	  return 0;
	}

      if((this_dir&last_dir)==0) // change in direction
	{
	  if(((last_dir==UP)&&(this_sign==NEG))||((last_dir==DOWN)&&(this_sign==POS)))
	    {
	      printf("Stationary point detected.\n");
	      if(!stat_point(z1,z2,n-1,Lc,Ls,Lf,Lfu,side,prec,true))
		{
		  printf("Error resolving stat point.\n");
		  return false;
		}
	      if(count<OUTPUT_ZEROS)
		arb_set(Lc->zeros[side][count],z1);
	      if(count<OUTPUT_ZEROS-1)
		arb_set(Lc->zeros[side][count+1],z2);
	      count+=2;
	    }
			 
	  last_dir=this_dir;
	}
      last_dir=this_dir;
    }
}


// called with suspected stationary point between m-1, m and m+1
// if !isolate_p, just confirm the two zeros to minimal precision
bool stat_point(arb_t z1, arb_t z2, uint64_t m, L_comp_t *Lc, L_upsample_t *Lu, L_family_t *Lf, L_func_t *Lfu, uint64_t side, uint64_t prec, bool isolate_p)
{
  static bool init=false;
  static arb_t t0,t1,f0,f1,t2,f2,t01,f01,t12,f12,zero_prec,tmp;
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
      arb_init(zero_prec);
      arb_set_ui(zero_prec,1);
      arb_mul_2exp_si(zero_prec,zero_prec,-Lc->OP_ACC-1);
      arb_init(tmp);
    }

  arb_mul_ui(t0,Lc->one_over_A,m-1,prec);
  arb_mul_ui(t1,Lc->one_over_A,m,prec);
  arb_mul_ui(t2,Lc->one_over_A,m+1,prec);
  arb_set(f0,Lu->values_off[side][m-1]);
  arb_set(f1,Lu->values_off[side][m]);
  arb_set(f2,Lu->values_off[side][m+1]);
  sign_t s=sign(f0);

  while(true)
    {
      //printf("iterating with\n");
      //arb_printd(t0,20);printf(" ");arb_printd(f0,20);printf("\n");
      //arb_printd(t1,20);printf(" ");arb_printd(f1,20);printf("\n");
      //arb_printd(t2,20);printf(" ");arb_printd(f2,20);printf("\n");
      
      arb_add(t01,t0,t1,prec);
      arb_mul_2exp_si(t01,t01,-1);
      upsample_stride(f01,t01,Lu,Lc,Lf,Lfu,side,prec);
      //printf("left middle = ");arb_printd(f01,20);printf("\n");
      sign_t s01=sign(f01);
      if(s01==UNK)
	{
	  printf("Double zero suspected near ");arb_printd(t01,20);printf("\n");
	  return false;
	}
      if(s01!=s)
	{
	  if(!isolate_p)
	    {
	      arb_union(z1,t0,t01,prec);
	      arb_union(z2,t01,t1,prec);
	      return(true);
	    }
	  if(!isolate_zero(tmp,t0,t01,f0,f01,s,s01,Lc,Lf,Lu,Lfu,side,prec,Lc->OP_ACC))
	    {printf("Isolate zero failed.\n");return false;}
	  if(!snap_point(z1,tmp,prec,Lc->OP_ACC))
	    {printf("Snap point failed.\n");return false;}
	  if(!confirm_zero(z1,Lu,Lc,Lf,Lfu,prec,zero_prec,side))
	    {printf("Confirm zero failed.\n");return false;}
	  if(!isolate_zero(tmp,t01,t1,f01,f1,s01,s,Lc,Lf,Lu,Lfu,side,prec,Lc->OP_ACC))
	    {printf("Isolate zero failed.\n");return false;}
	  if(!snap_point(z2,tmp,prec,Lc->OP_ACC))
	    {printf("Snap point failed.\n");return false;}
	  if(!confirm_zero(z2,Lu,Lc,Lf,Lfu,prec,zero_prec,side))
	    {printf("Confirm zero failed.\n");return false;}
	  return(true);
	}
      direction_t left=direction(f0,f01,prec);
      direction_t right=direction(f01,f1,prec);
      if((left==UNK)||(right==UNK))
	{
	  printf("Double zero suspected near ");arb_printd(t1,20);
	  return false;
	}
      if(left!=right)
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
      //printf("right middle = ");arb_printd(f12,20);printf("\n");
      sign_t s12=sign(f12);
      if(s12==UNK)
	{
	  printf("Double zero suspected near ");arb_printd(t12,20);printf(" Skipping.\n");
	  return false;
	}
      if(s12!=s)
	{
	  if(!isolate_p)
	    {
	      arb_union(z1,t1,t12,prec);
	      arb_union(z2,t12,t2,prec);
	      return true;
	    }

	  if(!isolate_zero(tmp,t1,t12,f1,f12,s,s12,Lc,Lf,Lu,Lfu,side,prec,Lc->OP_ACC))
	    return(false);
	  if(!snap_point(z1,tmp,prec,Lc->OP_ACC))
	    return false;
	  if(!confirm_zero(z1,Lu,Lc,Lf,Lfu,prec,zero_prec,side))
	    return false;
	  if(!isolate_zero(tmp,t12,t2,f12,f2,s12,s,Lc,Lf,Lu,Lfu,side,prec,Lc->OP_ACC))
	    return(false);
	  if(!snap_point(z2,tmp,prec,Lc->OP_ACC))
	    return false;
	  if(!confirm_zero(z2,Lu,Lc,Lf,Lfu,prec,zero_prec,side))
	    return false;
	  return(true);
	}
      left=direction(f1,f12,prec);
      right=direction(f12,f2,prec);
      if((left==UNK)||(right==UNK))
	{
	  printf("Double zero suspected near ");arb_printd(t1,20);printf(" Skipping.\n");
	  return false;
	}
      if(left!=right)
	{
	  arb_set(t0,t1);
	  arb_set(f0,f1);
	  arb_set(t1,t12);
	  arb_set(f1,f12);
	  continue;
	}
      else
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
