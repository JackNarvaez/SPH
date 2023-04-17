#include <iostream>
#include <stdlib.h>     // atof
#include <cmath>
#include <fstream>  // std::ifstream; // std::ofstream
#include <sstream>  // std::istringstream
#include <chrono>
#include "SPHfunctions.h"

int main(int argc, char **argv)
{
    auto start{std::chrono::steady_clock::now()};
    std::cout.precision(10);

    // Parameters
    const double dx = 0.05;
    const double X = 1;
    const double Y = 2;
    const double H = 3;
    const double W = 4;
    const int N = int(X*Y/(dx*dx));
    const int Nw = int(2*H/dx + W/dx);
    const double rho0 = 1000;
    const double P0 = 0;
    const double m = rho0*dx*dx;
    const double alpha=0.01;
    const double eps = 0.01;
    const double g = 9.81;
    const double cs = sqrt(200*g*Y);
    const double gamma = 7;
    const double p1 = 4;
    const double p2 = 2;
    const double D = 7.0/(4*M_PI*H*H);
    const double t0 = 0;
    const int T = 50000;
    const double h = 2*dx;
    const double dt=0.00002;
    const double tf = dt*T;
    const double t_trans{5*dt};

    double pos[N][2]{0.0};
    double vel[N][2]{0.0};
    double P[N]{0.0};
    double rho[N]{0.0};
    double posWall[Nw][2]{0.};
    std::ofstream Position, Boundary, Velocity, Press, Density, Accel;
    Position.open("./Files/Position.txt");
    Boundary.open("./Files/Boundary.txt");
    Velocity.open("./Files/Velocity.txt");
    Press.open("./Files/Press.txt");
    Density.open("./Files/Density.txt");
    //Create Particles 

    CreateFluidParticles(dx, X, Y, P0, rho0, pos, vel, P, rho);
    CreateBoundaryParticles(dx, H, W, posWall);

    // Save initial setting
    for (int ii = 0; ii<N; ii++){
        Position << pos[ii][0] << " " << pos[ii][1] << " ";
        Velocity << vel[ii][0] << " " << vel[ii][1] << " ";
        Press << P[ii] << " ";
        Density << rho[ii] << " ";
    }
    Position << "\n";
    Velocity << "\n";
    Press << "\n";
    Density << "\n";

    //Evolution
    for (int tt=1; tt < T; tt++){
        double fl[N][2]{0.0};
        double Drho[N]{0.0};
        double acc[N][2]{0.0};

        //Body forces
        for (int ii=0; ii<N; ii++){
            fl[ii][0] = 0;
            fl[ii][1] = -g*TransitionFunction(dt*tt, t_trans);
        }
        WallAcceleration(N, Nw, pos, posWall, 2*dx, D, p1, p2, fl);
        Acceleration(rho, P, Drho, acc, pos, vel, fl, h, N, m, rho0, P0, cs, gamma, alpha, eps, GaussianKernel, GradientGaussianKernel);

        double Temppos[N][2]{0.0};
        double Tempvel[N][2]{0.0};
        for (int ii = 0; ii<N; ii++){
            Tempvel[ii][0] = vel[ii][0] + 0.5*dt*acc[ii][0];
            Tempvel[ii][1] = vel[ii][1] + 0.5*dt*acc[ii][1];
            Temppos[ii][0] = pos[ii][0] + 0.5*dt*Tempvel[ii][0];
            Temppos[ii][1] = pos[ii][1] + 0.5*dt*Tempvel[ii][1];
            pos[ii][0] = Temppos[ii][0] + 0.5 *dt * Tempvel[ii][0];
            pos[ii][1] = Temppos[ii][1] + 0.5 *dt * Tempvel[ii][1];

            rho[ii] += dt*Drho[ii];
            
            vel[ii][0] = Tempvel[ii][0] + 0.5*dt * acc[ii][0];
            vel[ii][1] = Tempvel[ii][1] + 0.5*dt * acc[ii][1];
            
            Position << pos[ii][0] << " " << pos[ii][1] << " ";
            Velocity << vel[ii][0] << " " << vel[ii][1] << " ";
            Accel << acc[ii][0] << " " << acc[ii][1] << " ";
            Press << P[ii] << " ";
            Density << rho[ii] << " ";
        }
        if (tt%10==0){
            Position << "\n";
            Velocity << "\n";
            Press << "\n";
            Density << "\n";
        }
    }
    for (int ii=0; ii<Nw; ii++){
        Boundary << posWall[ii][0] << " " << posWall[ii][1] << " ";
    }
    
    Position.close();
    Boundary.close();
    Velocity.close();
    Press.close();
    Density.close();

    auto end{std::chrono::steady_clock::now()};
    std::cout << "Elapsed time in ms: "
        << std::chrono::duration_cast<std::chrono::seconds>(end-start).count()
        << "\n";
    return 0;
    return 0;
}