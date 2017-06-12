#include "string.h"
#include "acb_poly.h"
#include "primesieve.h"

// when called, we have just read the '[' of a polynomial definition
bool make_artin_lpoly(acb_poly_t inv, char *s, uint64_t den, uint64_t k, uint64_t prec)
{
  static bool init=false;
  static acb_poly_t tmp1,tmp2;
  static acb_t ctmp;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      acb_poly_init(tmp1);
      acb_poly_init(tmp2);
      acb_poly_one(tmp2);
      arb_init(tmp);
      acb_init(ctmp);
    }
  //printf("In make artin lpoly with %s\n",s);
  acb_poly_one(tmp1);
  uint64_t num,offset=0;
  if(s[0]==']') // its the empty polynomial []
    {
      acb_poly_one(inv);
      return true;
    }
  while(true)
    {
      if(sscanf(s+offset,"%lu",&num)!=1) return false;
      arb_set_si(tmp,num*2);
      arb_div_ui(tmp,tmp,den,prec);
      arb_sin_cos_pi(acb_imagref(ctmp),acb_realref(ctmp),tmp,prec);
      acb_neg(ctmp,ctmp);
      acb_poly_set_coeff_acb(tmp2,1,ctmp);
      acb_poly_mul(tmp1,tmp1,tmp2,prec);
      offset++;
      while((s[offset]!=0)&&(s[offset]!=']')&&(s[offset]!=','))
	offset++;
      if(s[offset]==',') // there is another term
	offset++;
      else // there isn't
	break;
    }
  // invert into tmp and return
  acb_poly_inv_series(inv,tmp1,k,prec);
  return true;
}

// for Artin L-functions, L->setup_string names a file containing 
// the Euler product info.
//
//format is
// <conductor>   we should already know this, so we just double check
// [u1, u2, ..., ur] the gamma_r factors. again we already know, so just check
//                   we see the right number. Don't bother checking the mus
// []
// <d>
// [poly1, poly2, poly3, ...., polym] the polynomial definitions
// which poly to use for p=2
// which poly to use for p=3
// ....
// a polynomial looks like 
// [] the polynomial 1
// [n_1, n_2, n_3 ... , n_k] where 1<=k<=r
// polys are defined as prod(1-e(n_i/d)T)
//
// since the polynomials are fixed and few in number, we invert them once and
// store them.
#define MAX_POLYS (32) // maximum number of defining polynomials allowed
bool do_lpolys_artin(L_func_t *L, L_family_t *Lf, L_comp_t *Lc, int64_t prec)
{
  static bool init=false;
  static acb_poly_t invs[MAX_POLYS]; // to hold the inverted polynomials
  static acb_poly_t tmp1,tmp2;
  static primesieve_iterator it;
  if(!init)
    {
      init=true;
      primesieve_init(&it);
      acb_poly_init(tmp1);
      acb_poly_init(tmp2);
      for(uint64_t i=0;i<MAX_POLYS;i++)
	acb_poly_init(invs[i]);
    }
  FILE *infile;
  if(!(infile=fopen(L->setup_string,"r")))
    {
      printf("Failed to open Artin Euler Product file %s. Skipping.\n",L->setup_string);
      return(false);
    }
  uint64_t cond;
  char buff[BUFF_LEN],*s;
  if(!fgets(buff,BUFF_LEN,infile)) return false;
  if(sscanf(buff,"%lu",&cond)!=1) return false;
  if(cond!=L->N)
    {
      printf("conductor %lu in file %s doesn't match.\n",cond,L->setup_string);
      return false;
    }
  if(!fgets(buff,BUFF_LEN,infile)) return false; // [mu, mu, mu...]
  if(buff[0]!='[') return false;
  uint64_t count=1;
  if(!(s=strtok(buff+1,",]"))) return false;
  while((s=strtok(NULL,",]\n"))) count++;
  if(count!=Lf->r) {printf("# Gamma factors in file %s don't match those in g-data.\n",L->setup_string);return false;}
  if(!fgets(buff,BUFF_LEN,infile)) return false; // should be []
  if(buff[0]!='['||buff[1]!=']'||buff[2]!='\n') return false;
  if(!fgets(buff,BUFF_LEN,infile)) return false; // the denominator
  uint64_t den;
  if(sscanf(buff,"%lu\n",&den)!=1) return false;

  uint64_t k=1;
  for(uint64_t M=L->M;M;M<<=1,k++); // log_2(M) max number of coeffs needed
                                    // in inverted polynomials

  if(!fgets(buff,BUFF_LEN,infile)) return false; // the Euler factors (nums)
  if(buff[0]!='[') return false; // syntax error
  uint64_t l=1; // buff[l]='['
  uint64_t poly_ptr=0;
  while(true)
    {
      if(buff[l]!='[') return false; // syntax error
      l++;
      if(!make_artin_lpoly(invs[poly_ptr++],buff+l,den,k,prec)) return false;
      while((buff[l]!=0)&&(buff[l]!=']')) l++; // scan to ]
      if(buff[l]==0) return false;
      l++;
      if(buff[l]==']') break; // the final closing ] ignore anything which follows
      if(buff[l]==',') l++; // another poly coming
      if(buff[l]==' ') l++;
    }

  // initialise all the Dirichlet coefficients to 1
  for(uint64_t m=0;m<L->M;m++)
    acb_one(L->ans[m]);


  uint64_t p=2,polyn;
  primesieve_skipto(&it,2,L->M);
  while(p<=L->M)
    {
      if(fscanf(infile,"%lu\n",&polyn)!=1) // ran out of polynomials
	{
	  uint64_t np=primesieve_next_prime(&it);
	  if(np<L->M)
	    {
	      printf("Ran out of polynomials with p=%lu. M set to %lu.\n",p,np-1);
	      L->M=np-1;
	    }
	  break;
	}
      // callback to populate the Dirichlet characters using the inverted
      // polynomial for this prime
      use_inv_lpoly(L,p,invs[polyn-1],prec);
      p=primesieve_next_prime(&it);
    }      
  fclose(infile);
  printf("Finished processing Euler factors.\n");
  fflush(stdout);
  return true;
}
