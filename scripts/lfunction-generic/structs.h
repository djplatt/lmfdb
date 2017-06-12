#include <inttypes.h>
#include "acb.h"
#include "smalljac.h"


// structure describing L function family
// see p 387 col 2 of Booker 2006
// and conditions for Lemma 5.2, 5,5
typedef struct{
  uint64_t r; // number of Gamma_R 0 means error
  double *mus; // the r mu_j's
  acb_t mu; // -1/2+1/r(1+sum mus[j])
  //int64_t m; // number of poles, all on Re(s)=1 -1 means error
  //arb_t *lambdas; // the imaginary parts of the poles
  //acb_t *residues; // residues at those poles
  arb_t *nus; // nu[j]=(Re mu[j]-{0.5,1})/2
  arb_t nu;
  arb_t alpha; // >= 1/r
  arb_t c; // Re mu +1/2+alpha
  //arb_t C;
  arb_t c_dash; // max(cr/2-1,0)
  double one_over_B;
  arb_t B;
  arb_t two_pi_by_B; // spacing of the Gs in the file
  int64_t low_i; // lowest value of i for which G^(k)(2 pi i/max_B) is present
                 // = round(log(1/sqrt(N))*max_B/2/pi)
  int64_t hi_i; // last i present
  uint64_t max_K; // max G^(k) provided
  arb_t **Gs; // the values of G^(k)(i/B)/k! read to be used
  uint64_t G_len; // number of G values to use =N/2-offset+1
  uint64_t X; // see eq 4-10
  arb_t imint;
  arb_t ftwiddle_error; // error of Lemma 5.7 (t= \pm B) over N
  int ltype; // what type of l-function is is
  // without the conductor factor from Q(s) which varies by L func
} L_family_t;

// structure with data for this particular L function
typedef struct{
  uint64_t N; // conductor 0 means error
  arb_t one_over_root_N;
  double dc;
  uint64_t allocated_M; // how much space have we allocated for ans?
  uint64_t M; // how many Dirichlet coefficients a_n do we have?
  bool self_dual_p;
  bool even_p;
  uint64_t character;
  int64_t conj_char; // might get -1 if bad character number
  acb_t *ans; // the coefficients a_n/sqrt(n)
  arb_t sum_ans; // sum over |a_n/sqrt(n)|
  acb_t epsilon; // the square root of the root number
  acb_t epsilon_sqr; // the square root of the root number
  // we'll compute epsilons from the output.
  arb_t ftwiddle_error; // error of Lemma 5.7 (t= \pm B)
  uint64_t wt;
  uint64_t ver;
  uint64_t hash;
  uint64_t rank;
  char setup_string[BUFF_LEN];
  char *sj_str1; // smalljac strings
  char *sj_str2;
  smalljac_curve_t curve;
} L_func_t;

// structure with parameters for the computation
typedef struct{
  uint64_t N; // length of convolutions we are going for a positive power of 2 please
  uint64_t NN; // length of large iFFT we are going to do
  int64_t offset; // lowest value of i for which G^(k)(i/B) is needed
  acb_t *w; // e(i/N) for i=0..N/2 the DFT
  acb_t *ww; // e(i/NN)
  double A; // NN/B
  uint64_t M0;
  arb_t arb_A;
  arf_t arf_A;
  arb_t one_over_A;
  arf_t arf_one_over_A;
  uint64_t K; // truncation bound for Taylor approx 5-9 and 5-10
  double eta; // in (0,1) ARB does (-1,1)
  arb_t delta; // pi/2*(1-eta) (|eta| in ARB but see above)
  acb_t *G,*kres,*res,**skm;
  arb_t *zeros[2];
  int64_t OP_ACC;
  arb_t zero_prec; // 0+/-2^{-OP_ACC-1}
} L_comp_t;

typedef struct{
  arb_t lem54;
  arb_t lem56;
  arb_t lem57;
  arb_t eq59;
} L_error_t;

/*
void print_family(L_family_t *L)
{
  printf("# Gamma factors=%lu\n",L->r);
  for(uint64_t i=0;i<L->r;i++)
    printf("   mu_%lu=%10.8e\n",i,L->mus[i]);
  for(uint64_t i=0;i<L->r;i++)
    {
      printf("   nu_%lu=",i);
      arb_printd(L->nus[i],10);
      printf("\n");
    }
  printf("mu=");acb_printd(L->mu,10);
  printf("c=");arb_printd(L->c,10);
  printf("\nc'=");arb_printd(L->c_dash,10);printf("\n");
}
*/

typedef struct{
  arb_t H;
  arb_t inv_2H2; // -1.0/2H^2
  arb_t upsampling_error;
  arb_t div_upsampling_error;
  uint64_t N; // number of samples to use each side of t
  arb_t *exps; // exp(-(n/A)^2/(2H^2)) n= 0 .. N-1
  uint64_t stride; // use every stride'th sample
  arb_t *values[2]; // normalised values of F(t)
  arb_t *values_off[2]; // offset so [0] is F(0)
  arb_t *values_div[2]; // normalised values of F(t) for poly division
  acb_t *expsincs; // fft's values of exp()sinc()
  acb_t *upsample_expsincs[DIV_UPSAMPLE_RATE-1];
  arb_t *values_start[2]; // where is F(0)
  arb_t *values_div_start[2]; // where is F(delta)
  uint64_t no_values;
  uint64_t no_values_div;
} L_upsample_t;

bool stat_point(arb_t, arb_t, uint64_t, L_comp_t *, L_upsample_t *, L_family_t *, L_func_t *, uint64_t, uint64_t, bool);
