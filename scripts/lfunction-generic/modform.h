// code to read Dirichlet coefficients from Bober generated files
// NB We do NOT go via Euler factors, we populate L->ans directly


// structures and I/O as per Bober's email 
typedef struct{
    uint32_t version; // == 3229261 (ascii string "MF1\0")
    uint32_t level;
    uint32_t weight;
    uint32_t chi;
    uint32_t orbit;
    uint32_t j;
    int32_t prec;
    int32_t exponent;
    uint32_t ncoeffs;
    char reserved[92];
} mfheader; // For a grand total of 128 bytes, since that is a nice round number.

int write_mfheader(FILE * outfile, mfheader * header) {
    if(!fwrite((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if(!fwrite((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fwrite((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fwrite((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fwrite((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fwrite((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fwrite((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fwrite((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fwrite((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fwrite((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

int read_mfheader(FILE * outfile, mfheader * header) {
    if(!fread((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if(header->version != 3229261) return 0;
    if(!fread((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fread((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fread((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fread((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fread((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fread((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fread((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fread((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fread((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

#define MF_TRIVIAL (1)

// read a single Dirichlet coefficient
bool mf_get_coeff(FILE *infile,acb_t an ,uint32_t chi, int32_t exponent, arb_t err,int64_t prec)
{
  static bool init=false;
  static fmpz_t z;
  if(!init)
    {
      init=true;
      fmpz_init(z);
    }
  //printf("chi=%u exponent=%u\n",chi,exponent);
  if(!fmpz_inp_raw(z,infile)) // the real part
    {
      printf("Error reading coefficient.\n");
      return false;
    }
  //fmpz_print(z);printf("\n");
  //exit(0);
  arb_set_fmpz(acb_realref(an),z);
  arb_mul_2exp_si(acb_realref(an),acb_realref(an),exponent);
  arb_add_error(acb_realref(an),err);
  if(chi==MF_TRIVIAL) // trivial character
    {
      arb_zero(acb_imagref(an));
      return true;
    }
  if(!fmpz_inp_raw(z,infile)) // the imaginary part
    {
      printf("Error reading coefficient.\n");
      return false;
    }
  arb_set_fmpz(acb_imagref(an),z);
  arb_mul_2exp_si(acb_imagref(an),acb_imagref(an),exponent);
  arb_add_error(acb_imagref(an),err);

  return true;
}

#define MF_EXACT (2147483647)
bool do_mf_dirichlet_coeffs(L_func_t *L, L_family_t *Lf, L_comp_t *Lc, int64_t prec)
{ 
  static bool init=false;
  static arb_t err;
  if(!init)
    {
      init=true;
      arb_init(err);
    }

  mfheader mf;
  char buff[BUFF_LEN];
  if(!sprintf(buff,"/scratch/jb12407/mf/%lu/%lu/%lu.%lu.%lu.%lu",L->N,L->wt,L->N,L->wt,L->character,L->ver))
    {
      printf("Failed to create Modform filename.\n");
      return false;
    }
  FILE *infile=fopen(buff,"r");
  if(!infile)
    {
      printf("Failed to open file %s.");
      return false;
    }
  if(!read_mfheader(infile,&mf))
    {
      printf("Failed to read mf header.");
      return false;
    }
  if(mf.ncoeffs<L->M)
    {
      printf("File only contains %lu coefficients, reducing M to that.\n");
      L->M=mf.ncoeffs;
    }

  if(mf.prec==MF_EXACT)
    arb_zero(err);
  else
    {
      arb_one(err);
      arb_mul_2exp_si(err,err,mf.prec);
    }
  printf("Modform coefficient error set to ");arb_printd(err,20);printf("\n");
  for(uint64_t m=0;m<L->M;m++)
    {
      if(!mf_get_coeff(infile,L->ans[m],mf.chi,mf.exponent,err,prec))
	{
	  if(m==0) // failed at first fence
	    {
	      printf("Failed to read any coefficients.\n");
	      return false;
	    }
	  printf("Failed to read coefficient at %lu. Reducing M.\n",m);
	  L->M=m;
	  return true;
	}
    }

  return true;
}

bool modform_parse_line(L_func_t *L)
{
  if(sscanf(L->setup_string,"%lu %lu %lu %lu\n",&L->N,&L->wt,&L->character,&L->ver)!=4)
    {
      printf("Error parsing Modform descriptor. Skipping.\n");
      return false;
    }
  // get rid of the spurious newline
  L->setup_string[strlen(L->setup_string)-1]=0;
  // this isn't used anywhere, but it could be
  printf("Character Number=%lu\n",L->character);
  L->conj_char=InvMod(L->character,L->N);
  if(L->conj_char==-1) // it had no inverse
    {
      printf("Bad character number %lu mod %lu. Skipping.\n",L->character,L->N);
      return false;
    }
  return true;
}
  
