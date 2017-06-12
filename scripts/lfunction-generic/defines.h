#define OUTPUT_RATIO (8) // we will analyse this portion of B
#define TURING_RATIO (8) // (16) we will use this portion for Turing
#define NUM_OUTPUT_POINTS (256) // output this + 1 Z values for plotting
#define OUTPUT_ZEROS (10) // number of zeros to rigorously isolate and output
#define MAX_TURING_R (4) // (3) largest degree our Turing code can cope with

#define POS (1)
#define NEG (2)
#define UNK (3)  
#define UP (1)
#define DOWN (2)
#define BOTH (3)
#define sign_t uint8_t
#define direction_t uint8_t

#define MAX_ZEROS (256)
#define BAD_64 (1LL<<62)

//#define HI_PREC

#define MAX_COND (1LL<<41)
#define PRIMORIAL_MAX_COND (12)

#define BUFF_LEN (1024)

#define DIV_UPSAMPLE_RATE (8)
#define LOG_UPS (3) // log 2 of the above

#define MAX_M_ERROR (1e-50) // the largest M error we can tolerate at x=0

