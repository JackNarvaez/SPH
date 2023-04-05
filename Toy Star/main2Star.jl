#using Plots
using LaTeXStrings
using Random
include("ToyStar.jl")
include("Integrator.jl")
include("Plots.jl")

# Simulation parameters
seed = 1234 # Random seed
t0 = 0      # Initial time
tEnd = 20   # Final time
dt = 0.04   # Timestep
d  = 2      # dimensions
k  = 0.1    # Pressure constant
n  = 1      # Polytropic index
nu = 1      # Viscocity coefficient

# Star 1
M1 = 4      # Star mass
R1 = 0.75   # Star radius
N1 = 400    # Number of Particles
m1 = M1/N1  # Particle mass

# Star 2
M2 = 1      # star mass
R2 = 0.5      # star radius
N2 = 100   # Number of Particles
m2 = M2/N2  # Particle mass


h = 0.04/sqrt(0.5*(N1+N2)/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, 0.5*(M1+M2), 0.5*(R1+R2))
m = 0.5*(m1+m2)

rho = zeros(T, N1+N2)
P = zeros(T, N1+N2)
vel = zeros(T, N1+N2, 2)
pos = zeros(T, N1+N2, 2)
a = zeros(T, N1+N2, 2)

pos[1, 1:N1, :] = Init_Dis(N1, R1, -1, -1, seed)
pos[1, N1+1:end, :] = Init_Dis(N2, R2, 1, 1, seed);

Leap_Frog(T, dt, Acceleration, pos, vel, a, N1+N2, k, n, lmbda, nu, m, h, Gaussian_Kernel, Gradient_Gaussian_Kernel)

record = 2

pmin = minimum(rho) 
pmax = maximum(rho)

animation(T, record, dt, pos, rho, "Density", pmin, pmax, 1.2(R1+R2), "2Star")