#ifndef _CONVERSION_H
#define _CONVERSION_H


using function = void(const double &, const double &, const double &, double &);
using derfunction = void (const double &, const double &, const double &, double DW[]);

void GaussianKernel(const double &x, const double &y, const double &h, double &W);
void GradientGaussianKernel(const double &x, const double &y, const double &h, double DW[]);
double Density(const double &x, const double &y, const double &m, const double&h, function Kernel);
double TaitStateEquation(const double &rho0, const double &P0, const double &cs, const double &gamma, const double &rho);
void LennardJonesPotential(const double pos[], const double &r0, const double &D, const double &p1, const double &p2, double fl[]);
double ArtificialViscosity(const double vel[], const double pos[], const double &cs, const double &rhoij, const double &hij, const double &alpha, const double &eps);
double TransitionFunction(const double &t, const double &t_trans);
void RelativeDistance(const double pos[][2], const int &N, const int &i, double relpos[][2]);
void CreateFluidParticles(const double &dx, const double &x, const double &y, const double &P0, const double &rho0, double pos[][2], double vel[][2], double P[], double rho[]);
void CreateBoundaryParticles(const double &dx, const double &H, const double &W, double pos[][2]);
void WallAcceleration(const int &N, const int &Nw, const double pospar[][2], const double poswall[][2], const double &r0, const double &D, const double &p1, const double &p2, double fl[][2]);
void Acceleration(double rho[], double P[], double Drho[], double acc[][2], const double pos[][2], const double vel[][2], const double fl[][2],
                  const double &h, const int &N, const double&m, const double &rho0, const double &P0, const double &cs, const double &gamma, const double &alpha, const double &eps, 
                  function Kernel, derfunction GradientKernel);

void CubicKernel(const double &x, const double &y, const double &h, double &W);
void GradientCubicKernel(const double &x, const double &y, const double &h, double DW[]);
int sgn(const double &x);

#endif