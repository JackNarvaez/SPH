#using DelimitedFiles
include("ToyStar.jl")
include("Integrator.jl")
include("Plots.jl")

t0 = 0      # Initial time
tEnd = 20   # Final time
dt = 0.04   # Timestep
N  = 100   # Number of Particles
R  = 0.75   # Star radius

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
record = 2

ioPos = open("./Files/1StarPos.txt", "r");
ioRho = open("./Files/1StarRho.txt", "r");

pmin = 0.0 
pmax = 3.0

animation(T, record, dt, ioPos, ioRho, "Density", pmin, pmax, R, "1Star")