#pragma once

#include "castor/matrix.hpp"

#ifndef CXX17_COMPLIANT

extern "C"
{
    #include "gsl/gsl_sf_bessel.h"
    #include "gsl/gsl_sf_legendre.h"
}

#else
#include <cmath>
#endif

namespace castor
{

/// Choose the type of field computation:\n
/// - infinity: compute the far-field amplitude\n
/// - boundary: compute the scattered field on the surface of the scatterer\n
/// - domain: compute the scattered field in the propagation domain
enum radtyp{infinity,boundary,domain};

/// Choose the type of boundary condition:\n
///  - dirichlet: homogeneous Dirichlet condition
///  - neumann: homogeneous Neumann condition
enum bndc{dirichlet,neumann};


// [AnalyticalScattering]
/// Class for the analytical computation of scattered waves when the incident
/// field is a plane wave propagating along the -e_z axis. The obstacle is
/// assumed to be a non-penetrable sphere with various boundary conditions.\n
/// There are three types of computation:\n
///  - field on the surface,\n
///  - field in the propagation domain,\n
///  - farfield amplitude (amplitude of the 0th-order term of the asymptotic 
/// expansion at infinity)
template<typename T>
class AnalyticalScattering
{
public:
    /// Computes the exact scattered acoustic field when the obstacle is a sphere. The arguments are:\n
    ///  - \e Xobs : observation nodes. If \e bc == farfield, the \e Xobs are the directions in which the farfield amplitude will be computed.
    ///  - \e rho : radius of the sphere.
    ///  - \e k : wavenumber.
    ///  - \e rt : type of computation (\e boundary, \e domain, \e infinity)
    ///  - \e bc : type of boundary condition (\e Dirichlet, \e Neumann)
    matrix<std::complex<T>> sphereHelmholtz(matrix<T> Xobs, T rho, T k, bndc bc, radtyp rt)
    {
        std::size_t nobs = size(Xobs,1);
        std::size_t n;
        double ka = static_cast<double>(std::abs(k));
        matrix<double> X(nobs,1), Y(nobs,1), Z(nobs,1);
        matrix<double> theta, phi, r;
        matrix<std::complex<double>> uu(nobs,1);
        
        // convert to spherical coordinates
        for(std::size_t ix=0; ix<nobs; ++ix)
        {
            X(ix) = Xobs(ix,0);
            Y(ix) = Xobs(ix,1);
            Z(ix) = Xobs(ix,2);
        }
        std::tie(phi,theta,r) = cart2sph(X,Y,Z);
        for(std::size_t ix=0; ix<nobs; ++ix) theta(ix) = M_PI/2. - theta(ix);

        switch(bc)
        {
            case dirichlet:
                uu = sphereHelmholtzDir(theta,r,rho,ka,rt);
                break;
            case neumann:
                uu = sphereHelmholtzNeu(theta,r,rho,ka,rt);
                break;
            default:
                error(__FILE__, __LINE__, __FUNCTION__,"Boundary condition not valid.");
                break;
        }
        if(k > 0)
        {
            for(std::size_t ix=0; ix<nobs; ++ix) uu(ix) = std::conj(uu(ix));
        }
        return cast<std::complex<T>>(uu);
    }

private:
    // Spherical bessel functions
    inline double sphBessel(std::size_t n, double x)
    {
#ifndef CXX17_COMPLIANT
        return std::sqrt(M_PI/(2*x))*gsl_sf_bessel_Jnu(n+0.5,x);
#else
        return std::sqrt(M_PI/(2*x))*std::cyl_bessel_j(n+0.5,x);
#endif
    }

    inline double sphNeumann(std::size_t n, double x)
    {
#ifndef CXX17_COMPLIANT
        return std::sqrt(M_PI/(2*x))*gsl_sf_bessel_Ynu(n+0.5,x);
#else
        return std::sqrt(M_PI/(2*x))*std::cyl_neumann(n+0.5,x);
#endif
    }

    // Functions for Dirichlet conditions
    inline std::complex<double> alphaDir(std::size_t n, double x)
    {
        double sphj = sphBessel(n,x);
        return -sphj/std::complex<double>(sphj,-sphNeumann(n,x));
    }

    matrix<std::complex<double>> sphereHelmholtzDir(matrix<double>const& theta, matrix<double>const& r, double rho, double k, radtyp rt)
    {
        std::size_t nobs = length(theta);
        std::size_t n=0;
        double maxR = max(r);
        double kr = k*rho;
        double Pn;
        std::complex<double> im = std::sqrt(std::complex<double>(-1.,0.));
        matrix<std::complex<double>> add(nobs,1);
        // output
        matrix<std::complex<double>> u(1,nobs);

        if(rt == boundary || rt == domain)
        {
            while(std::abs(alphaDir(n,k*maxR)) > 1e-12)
            {
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
#ifndef CXX17_COMPLIANT
                    Pn     = gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
                    Pn     = std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
                    u(ix) += std::pow(im,n)*(2.*n+1.)*Pn*alphaDir(n,kr)*std::complex<double>(sphBessel(n,k*r(ix)),-sphNeumann(n,k*r(ix)));
                }
                ++n;
            }
            if(rt == domain)
            {
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
                    if(r(ix) <= rho) u(ix) = 0;
                }
            }
        }
        else if(rt == infinity)
        {
#ifndef CXX17_COMPLIANT
            for(std::size_t ix=0; ix<nobs; ++ix) add(ix) = alphaDir(n,kr)*gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
            for(std::size_t ix=0; ix<nobs; ++ix) add(ix) = alphaDir(n,kr)*std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
            while(std::abs(norm(add,"inf")) > 1e-12)
            {
                ++n;
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
                    u(ix) += add(ix);
#ifndef CXX17_COMPLIANT
                    add(ix) = std::pow(-1,n)*(2.*n+1.)*alphaDir(n,kr)*gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
                    add(ix) = std::pow(-1,n)*(2.*n+1.)*alphaDir(n,kr)*std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
                }
            }
            u *= im/k;
        }
        // the end
        return u;
    }

    // Functions for Neumann conditions
    // derivative of spherical bessel function
    inline double djsph(std::size_t n, double x)
    {
        double sphj_nm1;
        if(n>=1) sphj_nm1 = sphBessel(n-1,x);
        else
        {
#ifndef CXX17_COMPLIANT
            sphj_nm1 = -std::sqrt(M_PI/(2*x))*gsl_sf_bessel_Ynu(0.5,x);
#else
            sphj_nm1 = -std::sqrt(M_PI/(2*x))*std::cyl_neumann(0.5,x);
#endif
        }
        return (n*sphj_nm1 - (n+1)*sphBessel(n+1,x))/(2.*n+1.);
    }

    // derivative of spherical Hankel function
    inline std::complex<double> dhsph(std::size_t n, double x)
    {
        return -(std::complex<double>(0,1./(x*x)) - djsph(n,x)*std::complex<double>(sphBessel(n,x),-sphNeumann(n,x)))/sphBessel(n,x);
    }

    inline std::complex<double> alphaNeu(std::size_t n, double x)
    {
        return -djsph(n,x)/dhsph(n,x);
    }

    matrix<std::complex<double>> sphereHelmholtzNeu(matrix<double>const& theta, matrix<double>const& r, double rho, double k, radtyp rt)
    {
        std::size_t nobs = length(theta);
        std::size_t n=0;
        double maxR = max(r);
        double kr = k*rho;
        double Pn;
        std::complex<double> im = std::sqrt(std::complex<double>(-1.,0.));
        matrix<std::complex<double>> add(nobs,1);
        // output
        matrix<std::complex<double>> u(1,nobs);

        if(rt == boundary || rt == domain)
        {
            while(std::abs(alphaNeu(n,k*maxR)) > 1e-12)
            {
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
#ifndef CXX17_COMPLIANT
                    Pn     = gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
                    Pn     = std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
                    u(ix) += std::pow(im,n)*(2.*n+1.)*Pn*alphaNeu(n,kr)*std::complex<double>(sphBessel(n,k*r(ix)),-sphNeumann(n,k*r(ix)));
                }
                ++n;
            }
            if(rt == domain)
            {
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
                    if(r(ix) <= rho) u(ix) = 0;
                }
            }
        }
        else if(rt == infinity)
        {
#ifndef CXX17_COMPLIANT
            for(std::size_t ix=0; ix<nobs; ++ix) add(ix) = alphaNeu(n,kr)*gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
            for(std::size_t ix=0; ix<nobs; ++ix) add(ix) = alphaNeu(n,kr)*std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
            while(std::abs(norm(add,"inf")) > 1e-12)
            {
                ++n;
                for(std::size_t ix=0; ix<nobs; ++ix)
                {
                    u(ix) += add(ix);
#ifndef CXX17_COMPLIANT
                    add(ix) = std::pow(-1,n)*(2.*n+1.)*alphaNeu(n,kr)*gsl_sf_legendre_Plm(n,0,std::cos(theta(ix)));
#else
                    add(ix) = std::pow(-1,n)*(2.*n+1.)*alphaNeu(n,kr)*std::assoc_legendre(n,0,std::cos(theta(ix)));
#endif
                }
            }
            u *= im/k;
        }
        // the end
        return u;
    }

};

// end of namespace
}