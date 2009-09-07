#ifndef INC_FINVMASS
#define INC_FINVMASS
//Calculate invariant Mass
namespace finvm{
  inline float fInvMass(float Energy, float PX, float PY, float PZ) {
    if ( Energy * Energy - PX * PX - PY * PY - PZ * PZ > 0 ) {
      float InvMass = sqrt(Energy * Energy - PX * PX - PY * PY - PZ * PZ );
      return InvMass;
    } else {
      std::cout << "Negative value of the invariant mass in the code" << std::endl;
      return -sqrt(fabs(Energy*Energy - PX*PX - PY*PY - PZ*PZ));
    }
  }
}
#endif
