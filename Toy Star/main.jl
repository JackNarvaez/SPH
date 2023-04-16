using DelimitedFiles
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
nu = 1      # Viscosity coefficient
N  = 100   # Number of Particles

m  = M/N    # Particle mass

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, M, R)
rho = zeros(N) # Density

P = zeros(N) # Pressure
vel = zeros(N, 2) # Velocity
pos = zeros(N, 2) # Position
a = zeros(N, 2) # Acceleration

pos = Init_Dis(N, R, 0, 0, seed)

#Files
ioPos = open("./Files/1StarPos.txt", "w");
ioRho = open("./Files/1StarRho.txt", "w");
writedlm(ioPos, [pos[:, 1]])
writedlm(ioPos, [pos[:, 2]])
writedlm(ioRho, [rho])

Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Gaussian_Kernel, Gradient_Gaussian_Kernel, ioPos, ioRho)

close(ioPos)
close(ioRho)