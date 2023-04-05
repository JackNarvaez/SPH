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
M  = 2      # Star mass
R  = 0.75   # Star radius
k  = 0.1    # Pressure constant
n  = 1      # Polytropic index
nu = 1      # Viscocity coefficient
N  = 100   # Number of Particles
m  = M/N    # Particle mass

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, M, R)
rho = zeros(T, N) # Density

P = zeros(T, N) # Pressure
vel = zeros(T, N, 2) # Velocity
pos = zeros(T, N, 2) # Position
a = zeros(T, N, 2) # Acceleration

pos[1, :, :] = Init_Dis(N, R, 0, 0, seed)

Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Gaussian_Kernel, Gradient_Gaussian_Kernel)

record = 2

pmin = minimum(rho) 
pmax = maximum(rho)

animation(T, record, dt, pos, rho, "Density", pmin, pmax, R, "1Star")