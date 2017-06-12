#include "smalljac.h"

int main()
{
  smalljac_curve_t curve;
  int k;

  curve=smalljac_curve_init("[x^5+x^3,x^4+x^2+1]",&k);
  smalljac_Qcurve_clear(curve);
  curve=smalljac_curve_init("[-x^7+3*x^6-4*x^5-x^4+8*x^3-14*x^2+9*x-4,x^2+x+1]",&k);
  smalljac_Qcurve_clear(curve);

  return 0;
}
