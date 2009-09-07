#ifndef INC_FCORR
#define INC_FCORR

namespace fcorr{
  inline double electron_br1(double brLinear, double e) {
    // These corrections are tuned on CMSSW 18x dynamicHybridSuperClusters!
    // YM: 02/05/2008
    //
    // first parabola (for brLinear < threshold)
    // p0*x^2 + p1*x + p2
    // second parabola (for brLinear >= threshold)
    // ax^2 + bx + c, make y and y' the same in threshold
    // y = p0*threshold^2 + p1*threshold + p2
    // yprime = p1 + 2*p0*threshold
    // a = p3
    // b = yprime - 2*a*threshold
    // c = y - a*threshold^2 - b*threshold
    // final result is multiplied by cos(p5/br + p6*br + p7)

    // make NO correction if brLinear is invalid!
    if ( brLinear == 0 ) return e;
    //

    // this is for 21X
    if ( brLinear < 1.1 ) brLinear = 1.1;
    if ( brLinear > 8 ) brLinear = 8.0;

    double p0 = -0.05185;
    double p1 = 0.1354;
    double p2 = 0.9165;
    double p3 = -0.0005626;
    double p4 = 1.385;    
    /*double p0 = -0.04382;
    double p1 = 0.1169;
    double p2 = 0.9267;
    double p3 = -0.0009413;
    double p4 = 1.419;    */


    double threshold = p4;
    
    double y = p0*threshold*threshold + p1*threshold + p2;
    double yprime = 2*p0*threshold + p1;
    double a = p3;
    double b = yprime - 2*a*threshold;
    double c = y - a*threshold*threshold - b*threshold;

    double fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return e/fCorr;

  }

  inline double electron_br1_complete(double et, double eta) {
    double fCorr = 0;

    // hybrid SC 
    // fBrem
    double c0 = 1.002;
    double c1 = -0.7424;
    double c2 = 0;
    double c3 = -0.00181;


    // fEtEta
    double c4 = 0.5558;
    double c5 = 2.375; 
    double c6 = 0.1869;
 
    // final fitting
    double c7 = 1.081;  // curve point in eta distribution
    double c8 = 7.6;     // sharpness of the curve

    double p0 = c0 + c1/(et + c2);
    double p1 = c4/(et + c5) + c6/(et*et);

    fCorr = p0 + p1*atan(c8*(c7 - fabs(eta))) + c3*fabs(eta);

    return et/fCorr;
  }

}
#endif //INC_FCORR
