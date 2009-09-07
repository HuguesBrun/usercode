#ifndef INC_FCORR_EE
#define INC_FCORR_EE

namespace fcorr_ee{
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
    // make a flat correction if brLinear is too big (>9)

    // FM with preshower
    if ( brLinear > 6.5 ) brLinear = 6.5;
    if ( brLinear < 0.9 ) brLinear = 0.9;

    // ============= Fixed Matrix With Preshower SC
    double p0 = -0.1214;   //-0.04163;
    double p1 = 0.2362;    //0.08551;
    double p2 = 0.8847;     //0.95048;
    double p3 = -0.00193; //-0.002308;
    double p4 = 1.057;      //1.077;

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
    double fCorr = 0.;

    double c0 = 2.213;
    double c1 = -17.29;
    double c2 = -0.599;

    double c3 = 8.874;
    double c4 = 0.09632; 
    double c5 = -1.457;

    double c6 = -0.7584;
    double c7 = 10.29;


    double p0 = c0 + c1/sqrt(et);
    double p1 = c2 + c3/sqrt(et);
    double p2 = c4 + c5/sqrt(et);
    double p3 = c6 + c7/sqrt(et);

    fCorr = p0 + p1*fabs(eta) + p2*eta*eta + p3/fabs(eta);

    return et/fCorr;
  }

}
#endif //INC_FCORR_EE


