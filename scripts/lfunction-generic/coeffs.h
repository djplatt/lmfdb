#include "string.h"
#include "slint.h" // thanks Bober/NTL
#include "acb_poly.h"

// callback routines
// given an inverted Euler polynomial
// use its coefficients to compute Dirichlet coefficients

// provide an inverted Euler polynomial for prime p as ints in c
// assumes there are enough coefficients in c to reach L->M
void use_inv_lpoly(L_func_t *L, uint64_t p, int64_t c[], uint64_t prec)
{
  // use inverted poly to populate Dirichlet coefficients
  uint64_t pn=p,pnn=p*p,j=1;
  while(pn<=L->M)
    {
      uint64_t ptr=pn,count=1;
      while(ptr<=L->M)
	{
	  if(count<p) // its not a higher prime power
	    {
	      acb_mul_si(L->ans[ptr-1],L->ans[ptr-1],c[j],prec);
	      count++;
	      ptr+=pn;
	    }
	  else // it is higher prime power, so skip it
	    {
	      ptr+=pn;
	      count=1;
	    }
	}
      pn*=p;
      pnn*=p;
      j++;
    }
}

// same thing, but the inverted polynomial is an arb_poly_t
void use_inv_lpoly(L_func_t *L, uint64_t p, arb_poly_t c, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  // use inverted poly to populate Dirichlet coefficients
  uint64_t pn=p,pnn=p*p,pow=1;
  while(pn<=L->M)
    {
      arb_poly_get_coeff_arb(tmp,c,pow);
      uint64_t ptr=pn,count=1;
      while(ptr<=L->M)
	{
	  if(count<p) // its not a higher prime power
	    {
	      acb_mul_arb(L->ans[ptr-1],L->ans[ptr-1],tmp,prec);
	      count++;
	      ptr+=pn;
	    }
	  else // it is higher prime power, so skip it
	    {
	      ptr+=pn;
	      count=1;
	    }
	}
      pn*=p;
      pnn*=p;
      pow++;
    }
}

// same thing, but the inverted polynomial is an acb_poly_t
void use_inv_lpoly(L_func_t *L, uint64_t p, acb_poly_t c, uint64_t prec)
{
  static bool init=false;
  static acb_t tmp;
  if(!init)
    {
      init=true;
      acb_init(tmp);
    }
  //if(p<10) {printf("p=%lu poly=",p);acb_poly_printd(c,20);printf("\n");}
  // use inverted poly to populate Dirichlet coefficients
  uint64_t pn=p,pnn=p*p,pow=1;
  while(pn<=L->M)
    {
      acb_poly_get_coeff_acb(tmp,c,pow);
      uint64_t ptr=pn,count=1;
      while(ptr<=L->M)
	{
	  if(count<p) // its not a higher prime power
	    {
	      acb_mul(L->ans[ptr-1],L->ans[ptr-1],tmp,prec);
	      count++;
	      ptr+=pn;
	    }
	  else // it is higher prime power, so skip it
	    {
	      ptr+=pn;
	      count=1;
	    }
	}
      pn*=p;
      pnn*=p;
      pow++;
    }
}

//************************************************************************
// Edit from here to add new L-function types
//***********************************************************************

// #define it here
#define MODFORM (1)
#define ARTIN (2)
#define ELLIPTIC (3)
#define GENUS2 (4)
#define GENUS3 (5)
#define ELLIPTIC_QD (6)

// Add a case here (and edit the default) to tell me how to print
// a description of your L-function
void print_ltype(L_family_t *lf)
{
  printf("Trying to process some ");
  switch(lf->ltype)
    {
    case MODFORM: printf("Modforms");break;
    case ARTIN: printf("Artin L-functions");break;
    case ELLIPTIC: printf("L-functions of elliptic curves");break;
    case GENUS2: printf("L-functions of genus 2 curves");break;
    case GENUS3: printf("L-functions of genus 3 curves");break;
    case ELLIPTIC_QD: printf("L-functions of elliptic curves over Q(sqrt(d))");
    default: printf("an unknown L-function type. Options are:-\n%d Modform\n%d Artin\n%d elliptic\n%d genus 2\n%d genus 3\nExiting.\n",MODFORM,ARTIN,ELLIPTIC,GENUS2,GENUS3);
      exit(0);
    }
  printf(".\n");
}

// Use L to be the L_func_t structure you are passed.
// add a prototype that will either
// a) put Dirichlet coefficients directly into L->ans[m]
//    (see modforms.h for an example of this method)
// b) generate inverted Euler polynomials, one per prime
//    in any order, to pass to use_inv_lpoly
//    I assume your poly has enough coefficients, i.e. p^(degree)>=L->M
//    (see artin.h or g2_g3_ell.h for examples of this method)
//
// when your handler is called, L_func->M will contain the number of
// Dirichlet coefficients I think I need. If you can't supply enough
// Euler polys to cover that, just reduce L_func->M accordingly
//
// If you are using the Euler polynomial version, set all the L_func->ans
// to 1+0i before calling use_inv_lpoly for the first prime
bool do_lpolys_ell_g2_g3(L_func_t *,L_family_t *,L_comp_t *,int64_t);
bool do_lpolys_artin(L_func_t *,L_family_t *,L_comp_t *,int64_t);
bool do_mf_dirichlet_coeffs(L_func_t *,L_family_t *,L_comp_t *,int64_t);
bool do_lpolys_ell_qd(L_func_t *,L_family_t *,L_comp_t *,int64_t);

// add a case for your new type here, calling the handler you prototyped
// above
bool do_lpolys(L_func_t *L, L_family_t *Lf, L_comp_t *Lc, int64_t prec)
{
  //printf("In do_lpolys with setup_string = |%s|.\n",L->setup_string);
  switch(Lf->ltype)
    {
    case MODFORM: return do_mf_dirichlet_coeffs(L,Lf,Lc,prec);
    case ELLIPTIC: 
    case GENUS2: 
    case GENUS3: return do_lpolys_ell_g2_g3(L,Lf,Lc,prec);
    case ARTIN: return do_lpolys_artin(L,Lf,Lc,prec);
    case ELLIPTIC_QD: return do_lpolys_ell_qd(L,Lf,Lc,prec);
    default: printf("Can't handle this type of L-function yet. Exiting.\n");
      exit(0);
    }
  return false;
}

// create a suitable name for the output file
// suffix will be "conj." for the conjugate L-function output
bool make_out_fname(char* fname, L_family_t *Lf, L_func_t *Lu, const char *suffix)
{
  char s1[BUFF_LEN],s2[BUFF_LEN],s3[BUFF_LEN],*s;
  uint64_t disc,cond,hash;
  switch(Lf->ltype)
    {
    case MODFORM:
      sprintf(fname,"%lu.%lu.%lu.%lu.%sout",Lu->N,Lu->wt,Lu->character,Lu->ver,suffix);return true;
    case ARTIN: sprintf(fname,"%s.%sout",Lu->setup_string,suffix);return true;
    case ELLIPTIC: 
      if(sscanf(Lu->setup_string,"%s %s %s",s1,s2,s3)!=3) return false;
      sprintf(fname,"%c%c/%s_%s_%s.%sout",s1[0],s1[1],s1,s2,s3,suffix);
      return true;
    case GENUS2:
    case GENUS3:
      sprintf(fname,"%lu_%lu.%sout",Lu->N,(long unsigned int)Lu->hash,suffix);
      return true;
    default: sprintf(fname,"%lu.%sout",Lu->N,suffix);return true;
    }
}

bool modform_parse_line(L_func_t *);
bool g2_parse_line(L_func_t *);
bool g3_parse_line(L_func_t *);
bool elliptic_parse_line(L_func_t *);
bool elliptic_qd_parse_line(L_func_t *);


// reads a line from infile
bool parse_line(FILE *infile, L_func_t *L, L_family_t *Lf)
{
  if(!fgets(L->setup_string,BUFF_LEN,infile))
    return false;
  switch(Lf->ltype)
    {
    case ELLIPTIC: return elliptic_parse_line(L);
    case MODFORM: return modform_parse_line(L);
    case GENUS2: return g2_parse_line(L);
    case GENUS3: return g3_parse_line(L);
    case ELLIPTIC_QD: return elliptic_qd_parse_line(L);
    default: printf("Unknown l function type.\n");return false;
    }
}

// I assume the Dirichlet coefficients are
// normalised for Lambda(s)=Lambda(wt-s)
// These coefficients will get normalised by n^((wt-1)/2)
// I'll need to divide by n^{-1/2} anyway for Booker's algorithm
// so this is at no extra charge
uint64_t get_weight(L_func_t *L, L_family_t *Lf)
{    
  switch(Lf->ltype)
    {
    case MODFORM: return L->wt;
    case ARTIN: return 1;
    case ELLIPTIC_QD:
    case ELLIPTIC:
    case GENUS2:
    case GENUS3: return 2;
    default: printf("Bad ltype in get_weight. Exiting.\n");exit(0);
    }
}

// is this type self dual (all the time)
// if not, set it to false
void self_dual_p(L_func_t *L, L_family_t *Lf)
{
  switch(Lf->ltype)
    {
    case ELLIPTIC_QD:
    case ELLIPTIC: 
    case GENUS2:
    case GENUS3: L->self_dual_p=true;break;
    default: L->self_dual_p=false;
    }
}

// put the code for your handler in a suitably named include file
// add it to the dependancies in makefile
#include "modform.h"
#include "artin.h"
#include "g2_g3_ell.h"
#include "ell_qd.h"
// you are done!

