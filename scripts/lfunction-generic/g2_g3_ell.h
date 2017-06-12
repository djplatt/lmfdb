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
		f[3] = aa[0]*p;
		f[4] = p*p;
		break;
	      case 6:
		f[1] = aa[0];
		f[2] = aa[1];
		f[3] = aa[2];
		f[4] = p*aa[1];
		f[5] = p*p*aa[0];
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
	      printf("unexpected bad reduction at p=%d. Exiting.\n",p);
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

bool do_lpolys_ell_g2_g3(L_func_t *L, L_family_t *Lf, L_comp_t *Lc, int64_t prec)
{
  // initialise all the Dirichlet coefficients to 1
  for(uint64_t m=0;m<L->M;m++)
    acb_set_ui(L->ans[m],1);
  
  // global variable to tell the callback what precision to use
  g2_ell_arb_prec=prec;
  
  switch(Lf->ltype)
    {
    case ELLIPTIC: smalljac_Lpolys(L->curve,1,L->M,SMALLJAC_A1_ONLY_SQRT,do_lpoly_elliptic,L);break;
    case GENUS2: smalljac_Lpolys(L->curve,1,L->M,SMALLJAC_A1_ONLY_SQRT,do_lpoly_g2,L);break;
    case GENUS3: smalljac_Lpolys(L->curve,1,L->M,SMALLJAC_A1_ONLY_SQRT,do_lpoly_g3,L);break;
    }
  smalljac_Qcurve_clear(L->curve);
  return(true);
}

bool g2_parse_line(L_func_t *L)
{
  char *s,*pq,lpolys[BUFF_LEN];
  uint64_t disc;
  int p,k,epsilon,d,i;

  printf("Doing genus 2 with str=|%s|\n",L->setup_string);
  if (!(s=strtok(L->setup_string,":")) || sscanf(s,"%ld",&disc) != 1) return false;
  printf("Discriminant = %lu\n",disc);fflush(stdout);
  if (!(s=strtok(NULL,":")) || sscanf(s,"%lu",&L->N)!=1) return false;
  printf("Conductor = %lu\n",L->N);fflush(stdout);
  if (!(s=strtok(NULL,":")) || sscanf(s,"%ld",&L->hash) != 1) return false;
  printf("hash = %lu\n",L->hash);
  if (!(pq=strtok(NULL,":"))) return false;
  if (!strtok(NULL,":")) return false; // skip disc sign
  if (!strtok(NULL,":")) return false; // skip Igusa-Clebsch invariants
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":")) || !strcpy(lpolys,s)) return false;
  printf("lpolys = %s\n",lpolys);fflush(stdout);
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":"))) return false;
  if (!(s=strtok(NULL,":")) || sscanf(s,"%lu",&L->rank)!=1) return false;
  printf("Rank = %lu\n",L->rank);
  printf("Calling smalljac_init with %s\n",pq);
  if ( !(L->curve=smalljac_curve_init(pq,&k)) )
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
  return true;
}

bool g3_parse_line(L_func_t *L)
{
  char *s,*pq;
  int64_t disc;
  int k,p,d,i;
  printf("Doing genus 3 with str=|%s|\n",L->setup_string);fflush(stdout);
  if (!(s=strtok(L->setup_string,":")) || sscanf(s,"%ld",&disc) != 1) 
    {
      printf("Bad discriminant.\n");
      return false;
    }
  printf("Discriminant = %ld\n",disc);fflush(stdout);
  if (!(s=strtok(NULL,":")) || sscanf(s,"%ld",&L->N)!=1)
    {
      printf("Bad conductor.\n");
      return false;
    }
  printf("Conductor = %ld\n",L->N);fflush(stdout);

  if(!(pq=strtok(NULL,":"))) return false; // the hash
  if(sscanf(pq,"%lu",&L->hash)!=1) return false;
  printf("Hash = %lu\n",L->hash);
  if (!(pq=strtok(NULL,":"))) return false; // string desribing curveto smalljac
  printf("Calling smalljac_init with %s\n",pq);fflush(stdout);
  if ( !(L->curve=smalljac_curve_init(pq,&k)) )
    {
      printf("smalljac_curve_init failed.\n");
      return false;
    }
  if(!strtok(NULL,":")) return false; // not sure what this field is
  if (!(pq=strtok(NULL," "))) return false;
  printf("Bad prime data = %s\n",pq);fflush(stdout);
  if(pq[1]=='0')
    {
      printf("Don't know what to do with prime 2 for this curve.\n");
      return false;
    }
  p = 2, k = 0;
  for (s=strtok(pq+1,",");s;s=strtok(NULL,","))
    {
      for (;p<=disc;p++)
	if (disc % p == 0) {
	  do disc /= p; while (disc % p == 0);
	  bad_lfactors[k].p = p;
	  d = i_poly_parse(bad_lfactors[k].f,6,s);
	  for (i=d+1;i<=6;i++)
	    bad_lfactors[k].f[i] = 0;
	  k++;
	  break;
	}
    }
  bad_lfactors[k].p = 0;
  if (disc != 1) return false;
  return true;
}


bool elliptic_parse_line(L_func_t *L)
{
  char PQ[BUFF_LEN],s1[BUFF_LEN],s2[BUFF_LEN],lpolys[BUFF_LEN];
  int p,k;
  uint64_t disc,d1,d2,d3;

  if(sscanf(L->setup_string,"%lu %s %lu %s %lu %lu %lu",&L->N,s1,&d1,PQ,&L->rank,&d2,&d3)!=7) return false;
  printf("Conductor = %lu\nRank = %lu\n",L->N,L->rank);
  printf("calling smalljac with %s\n",PQ);
      
  if ( !(L->curve=smalljac_curve_init(PQ,&k)) )
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
  return true;
}
