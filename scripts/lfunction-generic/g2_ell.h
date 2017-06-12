#include "string.h"
#include "acb_poly.h"
#include "smalljac.h"
#include "../andy_elliptic/pari.c"

//**********************************************************************
// utility routines by Andy Booker
//**********************************************************************

/* invert the polynomial p to precision k */
static void inverse(int64_t *ans,int64_t *p,int64_t k) {
	int64_t i,j,c[32];

	c[0] = 1;
	for (j=1;j<=k;j++) c[j] = 0;
	for (i=0;i<=k;i++) {
		ans[i] = c[0];
		for (j=1;j<=k-i;j++)
			c[j-1] = c[j]-p[j]*ans[i];
	}
}

// floor(log(x)/log(p))
static inline int logp(uint64_t x,uint64_t p) {
	int k;
	for (k=0;x>=p;k++)
		x /= p;
	return k;
}

// structure to hold the bad primes data
#define MAX_DEGREE (6) // Genus 3 curves
typedef struct {
	long int p,f[MAX_DEGREE+1];
} bad_lfactors_t;

// there can't be more bad primes than this
bad_lfactors_t bad_lfactors[PRIMORIAL_MAX_COND+1];

// global so we don't have to pass prec all over the place
int64_t g2_ell_arb_prec;

// handle an Euler Polynomial of degree DEGREE
// polynomial coefficients are long ints
int do_lpoly_smalljac(smalljac_Qcurve_t C,unsigned long p,
	     int good,long aa[],int n,void *arg,int DEGREE) 
  {
    if(DEGREE>MAX_DEGREE)
      {
	printf("DEGREE exceeded MAX_DEGREE in do_lpoly. Exiting.\n");
	exit(0);
      }
	int64_t f[64],c[64],i,j;
	uint64_t q,k;

	L_func_t *L=(L_func_t *)arg;

	f[0] = 1; for (j=1;j<64;j++) f[j] = 0;
	if (good) 
	  {
	    switch(DEGREE)
	      {
	      case 2:
		f[1] = aa[0];
		f[2] = p;
		break;
	      case 4:
		f[1] = aa[0];
		f[2] = aa[1];
		f[3] = f[1]*p;
		f[4] = p*p;
		break;
	      case 6:
		if(p<=100)
		  printf("p=%d aa = {%d,%d,%d}\n",p,aa[0],aa[1],aa[2]);
		f[1] = aa[0];
		f[2] = aa[1];
		f[3] = aa[2];
		f[4] = p*aa[2];
		f[5] = p*p*aa[1];
		f[6] = p*p*p;
		break;
	      default: printf("bad DEGREE in do_lpoly_smalljac.\n");
		return(0);
	      }
	  }
	else 
	  {
	    for (i=0;bad_lfactors[i].p;i++)
	      if (bad_lfactors[i].p == p) {
		for (j=0;j<=DEGREE;j++)
		  f[j] = bad_lfactors[i].f[j];
		break;
	      }
	    if (!bad_lfactors[i].p) {
	      fprintf(stderr,"unexpected bad reduction at p=%d\n",p);
	      exit(1);
	    }
	  }

	// invert the polynomial to required precision
	k = logp(L->M,p);
	inverse(c,f,k);
	// compute Dirichlet coefficients with it
	use_inv_lpoly(L,p,c,g2_ell_arb_prec);
	return 1;
}

// these are the callback routines
// called by smalljac for each prime p with aa containing the lpoly
// arg is set to the L_func_t structure (contains the Dirichlet coefficients)

int do_lpoly_elliptic(smalljac_Qcurve_t C,unsigned long p,
		      int good,long aa[],int n,void *arg) 
{
  return(do_lpoly_smalljac(C,p,good,aa,n,arg,2));
}

int do_lpoly_g2(smalljac_Qcurve_t C,unsigned long p,
		int good,long aa[],int n,void *arg) 
{
  return(do_lpoly_smalljac(C,p,good,aa,n,arg,4));
}

int do_lpoly_g3(smalljac_Qcurve_t C,unsigned long p,
		int good,long aa[],int n,void *arg) 
{
  return(do_lpoly_smalljac(C,p,good,aa,n,arg,6));
}


extern "C" int i_poly_parse(long *,int,char *);

// for elliptic and genus 2 curves, the setup_string contains
// the necessary incantation to get smalljac to return Euler polys
// The following is a rip off of Andy Booker's code.
bool do_lpolys_ell_g2(L_func_t *L, L_family_t *Lf, L_comp_t *Lc, int64_t prec)
{
  char *s,PQ[BUFF_LEN],s1[BUFF_LEN],s2[BUFF_LEN],lpolys[BUFF_LEN],*pq;
  long disc,cond,p,last_cond=0,hash;

  int dummy,d,i,k,epsilon;
  smalljac_curve_t curve;
  switch(Lf->ltype)
    {
    case ELLIPTIC:
      if(sscanf(L->setup_string,"%s %s %s",s1,s2,PQ)!=3) return false;
      /*
      // setup smalljac we just need the string to pass
      if (!(s=strtok(L->setup_string," ")) ) return(false); // skip 
      if (!(s=strtok(NULL," ")) ) return(false); // skip
      if (!(s=strtok(NULL," ")) ) return(false);
      PQ=s; // here it is
      */
      printf("calling smalljac with %s\n",PQ);
      
      if ( !(curve=smalljac_curve_init(PQ,&dummy)) )
	return(false);
      printf("smalljac returned ok.\n");
      
      // handle the bad primes
      p = 2, k = 0;
      disc=L->N; // this is <= MAX_COND 
      for (;p<=disc;p++)
	if (disc % p == 0) {
	  do disc /= p; while (disc % p == 0);
	  bad_lfactors[k].p = p;
	  bad_lfactors[k].f[0] = 1;
	  sprintf(lpolys,"ellap(ellinit(%s),%ld)",PQ,p);
	  bad_lfactors[k].f[1] = -parilong(lpolys);
	  bad_lfactors[k].f[2] = 0;
	  k++;
	}
      bad_lfactors[k].p = 0;
      if (disc != 1) {
	fprintf(stderr,"something has gone wrong with bad L-factors\n");
	return false;
      }
      break;

    case GENUS2:
      printf("Doing genus 2 with str=|%s|\n",L->setup_string);
      if (!(s=strtok(L->setup_string," ")) || sscanf(s,"%ld",&disc) != 1) return false;
      if (!(s=strtok(NULL," ")) ) return false;
      if (strchr(s,'?') || sscanf(s,"%ld",&cond) != 1 || !cond) return false;
      if (!(s=strtok(NULL," ")) || sscanf(s,"%ld",&hash) != 1) return false;
      if (!(pq=strtok(NULL," "))) return false;
      if (!strtok(NULL," ")) return false; // skip disc sign
      if (!strtok(NULL," ")) return false; // skip Igusa-Clebsch invariants
      if (!(s=strtok(NULL," \n"))) return false;
      if (!strcmp(s,"1") || !strcmp(s,"-1")) {
	sscanf(s,"%d",&epsilon);
	if (!(s=strtok(NULL," ")) || !strcpy(lpolys,s)) return false;
	if (!(s=strtok(NULL,"\n"))) return false;
	
	printf("Calling smalljac_init with %s\n",pq);
	if ( !(curve=smalljac_curve_init(pq,&k)) )
	  return false;
	//printf("smalljac_init returned something.\n");
	
	p = 2, k = 0;
	for (s=strtok(lpolys,",");s;s=strtok(NULL,","))
	  {
	    //printf("handling bad l-factor %s\n",s);
	    for (;p<=disc;p++)
	      if (disc % p == 0) {
		do disc /= p; while (disc % p == 0);
		bad_lfactors[k].p = p;
		//printf("calling i_poly_parse with s=%s\n",s);
		d = i_poly_parse(bad_lfactors[k].f,4,s);
		for (i=d+1;i<=4;i++)
		  bad_lfactors[k].f[i] = 0;
		k++;
		break;
	      }
	  }
	bad_lfactors[k].p = 0;
	if (disc != 1) return false;
      }

    case GENUS3: /* looks like <disc>:<cond>:<curve> */
      printf("Doing genus 3 with str=|%s|\n",L->setup_string);
      if (!(s=strtok(L->setup_string,":")) || sscanf(s,"%ld",&disc) != 1) return false;
      if (!(s=strtok(NULL,":")) ) return false; // the conductor. ignored
      if (!(pq=strtok(NULL," "))) return false;
      printf("Calling smalljac_init with %s\n",pq);
      if ( !(curve=smalljac_curve_init(pq,&k)) )
	{
	  printf("smalljac_curve_init failed.\n");
      return false;
	}
      
      
      // this is a complete guess at the functionality required
      s[0]='1';
      s[1]=0;
      p = 2, k = 0;
      for (;p<=disc;p++)
	if (disc % p == 0) 
	  {
	    do disc /= p; while (disc % p == 0);
	    bad_lfactors[k].p = p;
	    d = i_poly_parse(bad_lfactors[k].f,6,s);
	    for (i=d+1;i<=6;i++)
	      bad_lfactors[k].f[i] = 0;
	    k++;
	  }
      bad_lfactors[k].p = 0;
      if (disc != 1) return false;
      break;
    default: printf("Bad ltype in do_lpolys_ell_g2_g3. Exiting.\n");
      exit(0);
    }

  // initialise all the Dirichlet coefficients to 1
  for(uint64_t m=0;m<L->M;m++)
    acb_set_ui(L->ans[m],1);

  g2_ell_arb_prec=prec;

  switch(Lf->ltype)
    {
    case ELLIPTIC: smalljac_Lpolys(curve,1,L->M,0,do_lpoly_elliptic,L);break;
    case GENUS2: smalljac_Lpolys(curve,1,L->M,0,do_lpoly_g2,L);break;
    case GENUS3: smalljac_Lpolys(curve,1,L->M,0,do_lpoly_g3,L);break;
    }
  /*
    for(uint64_t m=0;m<128;m++)
    {
    printf("a_%ld = ",m+1);arb_printd(acb_realref(L->ans[m]),20);printf("\n");
    }
  */
  smalljac_Qcurve_clear(curve);
  
  return(true);
}




