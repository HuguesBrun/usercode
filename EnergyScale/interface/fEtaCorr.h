#ifndef INC_FETACORR
#define INC_FETACORR
namespace fcorr{
  //
  //Method to correct energy loss due to the leackage of a shower from the crystals
  //
  inline double f5x5( double iEta ) {
    if ( iEta < 40.2198 ) return 1;
    return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
  }

  inline double f3x3( double iEta ) {
    return 1 + 1.15581e-6*iEta*iEta - 1.35535e-10*iEta*iEta*iEta*iEta;
  }

  inline double fCorrEta( double iEta) {
    // old corrections DO NOT USE!

    Double_t fcorreta = 0;
    double p0 = 56.82;
    if ( fabs(iEta) < p0 ) {
      fcorreta = 0.99879;
    }
    if ( fabs(iEta) > p0 && fabs(iEta) <85.0 ) {
      double p1 = 1.006;
      double p2 = -2.227E-6;
      double p3 = 7.592E-11;
      fcorreta = p1 + p2 * fabs(iEta)*fabs(iEta) + p3 * fabs(iEta)*fabs(iEta)*fabs(iEta)*fabs(iEta);
    }
    if ( fcorreta == 0 )
      {
	std::cout << "Something is not right with Correction Function of Eta!!!" << std::endl;
	//return -100.000;
      }
    return fcorreta;
  }
}
#endif //INC_FETACORR
