//#ifndef FUNCWRAPPERS_H
//#define FUNCWRAPPERS_H

#include    <funcType.h>
extern "C"
{
   /*void C_vadaf(const double* energy, int nFlux, const double* params,
             int spectrumNumber, double* flux, double* fluxError,
             const char* initStr);
   void vadaf(const RealArray& energyArray, const RealArray& params,
               int spectrumNumber, RealArray& flux, RealArray& fluxErr,
               const string& initString);
  
   void C_adaf(const double* energy, int nFlux, const double* params,
             int spectrumNumber, double* flux, double* fluxError,
             const char* initStr); 
   void adaf(const RealArray& energyArray, const RealArray& params,
               int spectrumNumber, RealArray& flux, RealArray& fluxErr,
               const string& initString);
   */
   xsccCall     C_vadaf;
   XSCCall      vadaf;
   xsccCall     C_adaf;
   XSCCall      adaf;
   
   xsf77Call    badaf_;
   xsccCall     C_badaf;
   XSCCall      CXX_badaf;
   
   XSCCall      adafp;
   xsccCall     C_adafp;
/*   void C_badaf(const double* energy, int nFlux, const double* params,
                int spectrumNumber, double* flux, double* fluxError,
                const char* initStr);
   void CXX_badaf(const RealArray& energyArray, const RealArray& params,
                int spectrumNumber, RealArray& flux, RealArray& fluxErr,
                const string& initString);
  
   void C_adafp(const double* energy, int nFlux, const double* params,
             int spectrumNumber, double* flux, double* fluxError,
             const char* initStr);
   void adafp(const RealArray& energyArray, const RealArray& params,
               int spectrumNumber, RealArray& flux, RealArray& fluxErr,
               const string& initString);
*/
}

//#endif
//#endif
