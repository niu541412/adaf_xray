      SUBROUTINE badaf(ear, ne, param, ifl, photar, photer)

         INTEGER ne, ifl
         REAL ear(0:ne), param(4), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C simple thermal bremsstrahlung, based on GSFC routine
C FBG, as derived by GSFC for low temperature (<20 keV) plasmas
C the number of photons in the bin is estimated by a simple 2
C point approximation to the integral. Includes variable He/H.
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 1
C      1	T(keV)	Plasma temperature in keV
C      2	n(He)/n(H)  Helium to Hydrogen ratio
c    param(1) = slope of differential EM, alpha
c    param(2) = minmum temperature, tmin (keV)
c    param(3) = maximum temperature, tmax (keV)
c    He abundance is fixed at cosmic, intrisic bremsstrahlung ignores metal, so no emission lines.
c    param(4) = Redshift, z
C intrinsic energy range:
C      Emine=epsilon, Emax=infinity
C algorithm:
C      Uses the function Gaunt, which returns the photon spectrum,
C      based on the GSFC routine FBG
C---
C
C Shu Niu, Apr, 2020, add new parameter "Redshift", now the He abundance is defined by
C                     the intrisic cosmic abundance in Xspec, not the fixed value 0.085.
C 28 July 1988 - kaa
C 16 Dec 1994 - kaa         switched to Gaunt routine and included
C                           1/sqrt(T) factor to make the normalization
C                           independent of T. Normalization is now :
C                             K = (3.02e-15/4/pi/R^2) Int n_e n_H dV
C
C---
c     INCLUDE 'xspec.inc'

         REAL t, ab, zfac, elow, alow, ehi, ahi, tfac
         REAL alpha, tmin, tmax, logtemp, k, ltmin, ltmax
         INTEGER ie, i, nt

         REAL gaunt
         EXTERNAL gaunt


c suppress a warning message from the compiler
         ie = ifl

c this model does not calculate errors
         DO ie = 1, ne
            photer(ie) = 0.0
         ENDDO

         zfac = 1.0+param(4)

         k = 8.6171e-8
         alpha = param(1)
         tmin = param(2)
         tmax = param(3)
         ltmin = log10(tmin/k)
         ltmax = log10(tmax/k)
         logstep=30
         nt = INT((ltmax - ltmin)*logstep + 1)

c      print *,tmin,tmax,ltmin,ltmax,nt

         DO i = nt, 1, -1

            logtemp = ltmax - 1.0*(i-1)/logstep
            t = (10**logtemp)*k
c         print *,t
            tfac = 1./SQRT(t)

            ab = FGABNZ(2)

            elow = ear(0)*zfac

            alow = 0.
            IF ( (t.GT.0) .AND. (elow.LT.50*t) .AND. (elow.GT.0.) )
     &       alow = (gaunt(elow, t, 1.) + 4*ab*gaunt(elow, t, 2.))
     &             *tfac*EXP(-elow/t)/elow

            DO ie = 1, ne

               ehi = ear(ie)*zfac

               ahi = 0.
               IF ( (t.GT.0) .AND. (ehi.LT.50*t) .AND. (ehi.GT.0.) )
     &         ahi = (gaunt(ehi, t, 1.) + 4*ab*gaunt(ehi, t, 2.))
     &                *tfac*EXP(-ehi/t)/ehi

               photar(ie) = photar(ie) + ((t/tmin)**alpha * 0.1)
     &            * (0.5*(alow+ahi)*(ehi-elow))

               elow = ehi
               alow = ahi

            ENDDO
         ENDDO

c Shift the energies back and correct the model

         DO ie = 0, ne
            ear(ie) = ear(ie)/zfac
         ENDDO

         DO ie = 1, ne
            photar(ie) = photar(ie)/zfac
         ENDDO

         RETURN
      END
