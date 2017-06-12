#include "string.h"
#include "acb_poly.h"
#include "smalljac.h"

int64_t ell_qd_arb_prec;

// this will be called by smalljac for each prime p
int do_lpoly_elliptic_qd(smalljac_Qcurve_t C,unsigned long p,
		      int good,long aa[],int n,void *arg) 
{
  return 0;
}


// L->setup_string contains the line read from input
// need to set L->conductor
// L->curve
// L->rank
// some way of handling bad primes
// see elliptic_parse_line in g2_g3_ell.h for example
bool elliptic_qd_parse_line(L_func_t *L)
{
  return false;
}

bool do_lpolys_ell_qd(L_func_t *L,L_family_t *Lf,L_comp_t *Lc,int64_t prec)
{
  // initialise all the Dirichlet coefficients to 1
  for(uint64_t m=0;m<L->M;m++)
    acb_set_ui(L->ans[m],1);
  
  // global variable to tell the callback what precision to use
  ell_qd_arb_prec=prec;

  smalljac_Lpolys(L->curve,1,L->M,SMALLJAC_A1_ONLY_SQRT,do_lpoly_elliptic_qd,L);
  return false;
}



