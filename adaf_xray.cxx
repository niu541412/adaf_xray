//#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <xsTypes.h>
//#include <XSFunctions/funcWrappers.h>
#include <funcWrappers.h>
#include "isisWrappers.h"

extern "C" 
void adafp(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr,
	       const string& initString)
{
   /*
     c Parameters:
     c    ==================================================================
     c    param(1) = slope of differential EM, alpha
     c    param(2) = minmum temperature, tmin (keV)
     c    param(3) = maximum temperature, tmax (keV)
     c    param(4) = Metal abundances (He fixed at cosmic), only for apec components due to intrisic bremsstrahlung ignoring it.
     c    param(5) = Redshift, z
         ==================================================================
     
    Shu Niu, Apr, 2020, merge adaf_apec and adaf_brem into one model which has a very hot inner region.
    Because the maximum temperature of apec is about 64(or 86)keV, the inner and hotter area is caculated by bremsstrahlung.
    Tint is the temparature at the interface between apec or bremsstrahlung. 
    Scale is the normalization which makes the connection at the boundary at Tint smoothly, i.e. differential emission
    measure is continuous.   
    Note that, bremsstrahlung doesn't consider metal abundance.  


  */

  const Real Tint = 64.0;
  
  RealArray aparams(5);
  for (size_t i=0; i<5; i++) aparams[i] = params[i];
  aparams[2]= Tint;
  
  RealArray bparams(4);
  bparams[0] = params[0];
  bparams[1] = Tint;
  bparams[2] = params[2];
  bparams[3] = params[4];
  
  /*printf("x address: %p\n", (void *)&params);
  printf("a address: %p\n", (void *)&aparams); 
  printf("b address: %p\n", (void *)&bparams); 
  FunctionUtility::xsWrite("ZZZZZ "+ std::to_string(params.size()) +" @@@",2);
  FunctionUtility::xsWrite("AAAAA "+ std::to_string(aparams.size()) +" @@@",2);
  FunctionUtility::xsWrite("BBBBB "+ std::to_string(bparams.size()) +" @@@",2);
*/
 
  const Real ab = FunctionUtility::getAbundance(2);
  //scaled norm of adaf_brem, which makes the connection at the boundary at Tint smoothly.
  Real scale = 3.02/10*(1+ab)*pow(Tint/params[1], params[0]);

  RealArray tempflux(0.0,energyArray.size()-1);
  RealArray tempfluxErr(0.0,energyArray.size()-1);
  
  adaf(energyArray, aparams, spectrumNumber, tempflux, tempfluxErr, initString);
  CXX_badaf(energyArray, bparams, spectrumNumber, flux, fluxErr, initString);

 // for (size_t i=0; i<flux.size(); i++) flux[i] = tempflux[i]+flux[i]*scale;
 // for (size_t i=0; i<fluxErr.size(); i++) fluxErr[i] = tempfluxErr[i]+fluxErr[i]*scale;
  
  flux = tempflux+flux*scale;
  fluxErr = tempfluxErr+fluxErr*scale;

  return;
}
