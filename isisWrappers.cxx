#include <XSFunctions/Utilities/XSCall.h>
#include <xsTypes.h>
#include <cfortran.h>
#include "isisWrappers.h"

// Shu Niu, Jul, 2024, fix bug introduced since Xspec 12.13.0.

void cppModelWrapper(const double *energy, int nFlux, const double *params,
                     int spectrumNumber, double *flux, double *fluxError, const char *initStr,
                     int nPar, void (*cppFunc)(const RealArray &, const RealArray &, int, RealArray &, RealArray &, const string &));
void fcppModelWrapper(const float *energy, int nFlux, const float *params,
                      int spectrumNumber, float *flux, float *fluxError,
                      int nPar, void (*cppFunc)(const RealArray &, const RealArray &, int, RealArray &, RealArray &, const string &));
void fModelWrapper(const double *energy, int nFlux, const double *params,
                   int nPar, int spectrumNumber, double *flux, double *fluxError,
                   xsf77Call fFunc);

// adaf_apec,"C_" wrapper seems not necessary.
void f_adaf(const float *energy, int nFlux, const float *params,
            int spectrumNumber, float *flux, float *fluxError);

FCALLSCSUB6(f_adaf, ADAF, adaf, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
FCALLSCSUB7(C_adaf, DADAF, dadaf, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, STRING)

void f_adaf(const float *energy, int nFlux, const float *params,
            int spectrumNumber, float *flux, float *fluxError)
{
        const size_t nPar = 5;
        fcppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                         nPar, adaf);
}

void C_adaf(const double *energy, int nFlux, const double *params,
            int spectrumNumber, double *flux, double *fluxError, const char *initStr)
{
        const size_t nPar = 5;
        cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                        initStr, nPar, adaf);
}

// adaf_vapec
void f_vadaf(const float *energy, int nFlux, const float *params,
             int spectrumNumber, float *flux, float *fluxError);

FCALLSCSUB6(f_vadaf, VADAF, vadaf, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
FCALLSCSUB7(C_vadaf, DVADAF, dvadaf, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, STRING)

void f_vadaf(const float *energy, int nFlux, const float *params,
             int spectrumNumber, float *flux, float *fluxError)
{
        const size_t nPar = 18;
        fcppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                         nPar, vadaf);
}

void C_vadaf(const double *energy, int nFlux, const double *params,
             int spectrumNumber, double *flux, double *fluxError, const char *initStr)
{
        const size_t nPar = 18;
        cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                        initStr, nPar, vadaf);
}

// adaf_brem
void C_badaf(const double *energy, int nFlux, const double *params,
             int spectrumNumber, double *flux, double *fluxError, const char *initStr)
{
        const size_t nPar = 4;
        fModelWrapper(energy, nFlux, params, nPar, spectrumNumber, flux, fluxError,
                      badaf_);
}

void CXX_badaf(const RealArray &energyArray, const RealArray &params,
               int spectrumNumber, RealArray &flux, RealArray &fluxErr, const string &initString)
{
        XSCallBase *funcObj = new XSCall<xsf77Call>(badaf_);
        (*funcObj)(energyArray, params, spectrumNumber, flux, fluxErr, initString);
        delete funcObj;
}

// adaf_xray
void f_adafx(const float *energy, int nFlux, const float *params,
             int spectrumNumber, float *flux, float *fluxError);

FCALLSCSUB6(f_adafx, ADAFX, adafx, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
FCALLSCSUB7(C_adafx, DADAFX, dadafx, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, STRING)

void f_adafx(const float *energy, int nFlux, const float *params,
             int spectrumNumber, float *flux, float *fluxError)
{
        const size_t nPar = 5;
        fcppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                         nPar, adafx);
}

void C_adafx(const double *energy, int nFlux, const double *params,
             int spectrumNumber, double *flux, double *fluxError, const char *initStr)
{
        const size_t nPar = 5;
        cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                        initStr, nPar, adafx);
}

// adaf_vxray
void f_vadafx(const float *energy, int nFlux, const float *params,
              int spectrumNumber, float *flux, float *fluxError);

FCALLSCSUB6(f_vadafx, VADAFX, vadafx, FLOATV, INT, FLOATV, INT, FLOATV, FLOATV)
FCALLSCSUB7(C_vadafx, DVADAFX, dvadafx, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, DOUBLEV, STRING)

void f_vadafx(const float *energy, int nFlux, const float *params,
              int spectrumNumber, float *flux, float *fluxError)
{
        const size_t nPar = 18;
        fcppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                         nPar, vadafx);
}

void C_vadafx(const double *energy, int nFlux, const double *params,
              int spectrumNumber, double *flux, double *fluxError, const char *initStr)
{
        const size_t nPar = 18;
        cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError,
                        initStr, nPar, vadafx);
}
