#include <iostream>
#include <cmath>

using function = void(const double &, const double &, const double &, double &);
using derfunction = void (const double &, const double &, const double &, double DW[]);

int sgn(const double &x){
    int s{(x > 0) - (x < 0)};
    return s;
}

void GaussianKernel(const double &x, const double &y, const double &h, double &W){
    /*---------------------------------------------------------------------------------------------
    Truncated Gaussian smoothing function
    ---------------------------------------------------------------------------------------------*/
    const double H{h*h};
    const double x2{x*x+y*y};
    W = 0.0;
    if (x2/H <= 9){
        W = 1.0/(H*M_PI)*exp(-x2/H);
    }
}

void CubicKernel(const double &x, const double &y, const double &h, double &W){
    /*---------------------------------------------------------------------------------------------
    Truncated Gaussian smoothing function
    ---------------------------------------------------------------------------------------------*/
    const double a{1/h};
    const double q{sqrt(x*x+y*y)/h};
    W = 0.0;
    if (q < 1.){
        W = a*(2./3 -q*q+0.5*pow(q,3));
    }
    if (q <= 2.){
        W = a*1./6.0 *pow(2-q,3);
    }
}

void GradientGaussianKernel(const double &x, const double &y, const double &h, double DW[]){
    /*---------------------------------------------------------------------------------------------
    Gradient of the Truncated Gaussian smoothing function
    ---------------------------------------------------------------------------------------------*/
    const double x2{x*x+y*y};
    double W{0.0};
    double H{h*h};
    if (x2/H <= 9){
        W = -2.0/(pow(H,2)*M_PI)*exp(-x2/H);
    }
    DW[0] = x*W;
    DW[1] = y*W;
}

void GradientCubicKernel(const double &x, const double &y, const double &h, double DW[]){
    /*---------------------------------------------------------------------------------------------
    Gradient of the Truncated Gaussian smoothing function
    ---------------------------------------------------------------------------------------------*/
    const double a{1/h};
    const double q{sqrt(x*x+y*y)/h};
    double W{0.0};
    if (q < 1.){
        W = sgn(q)*a*(-2*q+1.5*q*q);
    }
    if (q <= 2.){
        W = -sgn(q)*a* 0.5 *(2-q)*(2-q);
    }
    DW[0] = x*W;
    DW[1] = y*W;
}

double Density(const double &x, const double &y, const double &m, const double&h, function Kernel) {
    /*---------------------------------------------------------------------------------------------
    SPH approximation for density equation. Particles have equal mass m.
    ---------------------------------------------------------------------------------------------*/
    double W{0.0};
    Kernel(x, y, h, W);
    return m*W;
}

double TaitStateEquation(const double &rho0, const double &P0, const double &cs, const double &gamma, const double &rho) {
    /*---------------------------------------------------------------------------------------------
    Tait State equation: Relate liquid density to hydrostatic pressure.
    -----------------------------------------------------------------------------------------------
    Arguments:
    rho0    :   Density at pressure P0
    P0      :   Reference pressure
    cs      :   Speed of sound
    gamma   :   Material parameter
    rho     :   Density
    -----------------------------------------------------------------------------------------------
    Returns the Pressure for a given density rho. 
    ---------------------------------------------------------------------------------------------*/
    double P {cs*cs*rho0/gamma *(pow(rho/rho0, gamma)-1)+P0};
    //std::cout << "P:\t" << cs*cs*rho0/gamma << "\t" << pow(rho/rho0, gamma)<< "\t " << P <<std::endl;
    return P;
}

void LennardJonesPotential(const double pos[], const double &r0, const double &D, const double &p1, const double &p2, double fl[]){
    double MagPos = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    double Fl{0.0};
    if (MagPos<=r0) {
        Fl = D*(pow(r0/MagPos, p1)-pow(r0/MagPos, p2))/pow(MagPos, 2);
    }
    fl[0] = pos[0]*Fl;
    fl[1] = pos[1]*Fl;
}

double ArtificialViscosity(const double vel[], const double pos[], const double &cs, const double &rhoij, const double &hij, const double &alpha, const double &eps){
        double vis{(alpha*cs*hij/rhoij)*(vel[0]*pos[0]+vel[1]*pos[1])/((pow(pos[0], 2)+pow(pos[1], 2))+eps*hij*hij)};
        return std::max(-vis, 0.0);
}

double TransitionFunction(const double &t, const double &t_trans){
    double zeta{1.0};
    if (t< t_trans){
        zeta = 0.5*(sin((-0.5+t/t_trans)*M_PI)+1);
    }
    return zeta;
}

void RelativeDistance(const double pos[][2], const int &N, const int &i, double relpos[][2]){
    for (int jj=0; jj<(N-i-1); jj++){
        relpos[jj][0] = pos[i][0]-pos[i+1+jj][0];
        relpos[jj][1] = pos[i][1]-pos[i+1+jj][1];
    }
}

void CreateFluidParticles(const double &dx, const double &x, const double &y, const double &P0, const double &rho0, double pos[][2], double vel[][2], double P[], double rho[]){
    int Ny{int(y/dx)};
    int Nx{int(x/dx)};
    for (int ii=0; ii<Ny; ii++){
        for (int jj=0; jj<Nx; jj++){
            pos[ii*Nx+jj][0] = dx*(jj+1);
            pos[ii*Nx+jj][1] = dx*(ii+1);
            vel[ii*Nx+jj][0] = 0;
            vel[ii*Nx+jj][1] = 0;
            P[ii*Nx+jj] = P0;
            rho[ii*Nx+jj] = rho0;
        }
    }
}

void CreateBoundaryParticles(const double &dx, const double &H, const double &W, double pos[][2]){
    int Ny{int(H/dx)};
    int Nx{int(W/dx)};
    for (int ii=0; ii<(Nx+1); ii++){
        pos[ii][0] = ii*dx;
        pos[ii][1] = 0.0;
    }
    for (int jj=1; jj<(Ny+1); jj++){
        pos[Nx+jj][0]=0;
        pos[Nx+jj][1]=jj*dx;
        pos[Nx+Ny+jj][0]= W;
        pos[Nx+Ny+jj][1]=pos[Nx+jj][1];
    }
}

void WallAcceleration(const int &N, const int &Nw, const double pospar[][2], const double poswall[][2], const double &r0, const double &D, const double &p1, const double &p2, double fl[][2]){
    for (int ii=0; ii<N; ii++){
        for (int jj=0; jj<Nw; jj++){
            double dr[2]{pospar[ii][0] - poswall[jj][0], pospar[ii][1] - poswall[jj][1]};
            double f_lj[2]{0.0, 0.0};
            LennardJonesPotential(dr, r0, D, p1, p2, f_lj);
            fl[ii][0] += f_lj[0];
            fl[ii][1] += f_lj[1];
        }
    }   
}

void Acceleration(double rho[], double P[], double Drho[], double acc[][2], const double pos[][2], const double vel[][2], const double fl[][2],
                  const double &h, const int &N, const double&m, const double &rho0, const double &P0, const double &cs, const double &gamma, const double &alpha, const double &eps, 
                  function Kernel, derfunction GradientKernel) {
    /*---------------------------------------------------------------------------------------------
    Hydrodynamics SPH equations to evolve the system.
    ---------------------------------------------------------------------------------------------*/
    double rho_0{Density(0, 0, m, h, Kernel)};
    for (int ii=0; ii < N; ii++){
        //rho[ii] = rho_0;
        acc[ii][0] = fl[ii][0];
        acc[ii][1] = fl[ii][1];
        Drho[ii] = 0.;
    }
    // Density
    /*for (int ii=0; ii<(N-1); ii++){
        double dr[N-ii-1][2]{0.};
        RelativeDistance(pos, N, ii, dr);
        for (int jj=ii+1; jj<N; jj++){
            if ((pow(dr[jj-ii-1][0], 2)+pow(dr[jj-ii-1][1],2))/(h*h) <= 9){
                double rho_ij{Density(dr[jj-ii-1][0], dr[jj-ii-1][1], m, h, Kernel)};
                rho[ii] += rho_ij;
                rho[jj] += rho_ij;
            }
        }
    }*/
    // Pressure
    for (int ii=0; ii<N; ii++) {
        P[ii] = TaitStateEquation(rho0, P0, cs, gamma, rho[ii]);
    }
    //Acceleration
    for (int ii=0; ii<(N-1); ii++) {
        double dr[N-ii-1][2]{0.};
        double dv[N-ii-1][2]{0.};
        RelativeDistance(pos, N, ii, dr);
        RelativeDistance(vel, N, ii, dv);
        for (int jj=ii+1; jj<N; jj++){
            if ((pow(dr[jj-ii-1][0], 2)+pow(dr[jj-ii-1][1],2))/(h*h) <= 9){
                double DW[2]{0.0, 0.0};
                GradientKernel(dr[jj-ii-1][0], dr[jj-ii-1][1], h,  DW);
                double Drho_ij{m*(dv[jj-ii-1][0]*DW[0]+dv[jj-ii-1][1]*DW[1])};
                Drho[ii] += rho[ii]/rho[jj] * Drho_ij;
                Drho[jj] += rho[jj]/rho[ii] * Drho_ij;
                double rhoij{0.5*(rho[ii]+rho[jj])};
                double ArtVis{ArtificialViscosity(dv[jj-ii-1], dr[jj-ii-1], cs, rhoij, h, alpha, eps)};
                double acc_ij[2]{0.0, 0.0};
                acc_ij[0] = m*(P[ii]/pow(rho[ii],2) + P[jj]/pow(rho[jj],2)+ArtVis);
                acc_ij[1] = acc_ij[0]*DW[1];
                acc_ij[0] *= DW[0];
                acc[ii][0] -= acc_ij[0];
                acc[ii][1] -= acc_ij[1];
                acc[jj][0] += acc_ij[0];
                acc[jj][1] += acc_ij[1];
            }
        }
    }
}