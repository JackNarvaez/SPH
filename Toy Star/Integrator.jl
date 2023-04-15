function Euler_Cromer(T, dt, pos, vel, a, N, k, n, lmbda, nu, m, h, Acceleration, Kernel, Gradient_Kernel)
    #=---------------------------------------------------------------
    Euler Cromer method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    ---------------------------------------------------------------=#
    for i in 1:T-1
        pos[i+1, :, :] = pos[i, :, :] .+ dt.*vel[i, :, :]
        rho[i+1, :], P[i+1, :], a[i+1, :, :] = Acceleration(pos[i+1, :, :], vel[i, :, :], N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        vel[i+1, :, :] = vel[i, :, :] .+ dt.*a[i+1, :, :]
        
    end
end;

function Verlet_Pos(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    #=---------------------------------------------------------------
    Position Verlet method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    ---------------------------------------------------------------=#
    for i in 1:T-1
        Temp_Pos = pos[i, :, :] .+ 0.5*dt.*vel[i, :, :] 
        rho[i, :], P[i, :], a[i, :, :] = Acceleration(Temp_Pos, vel[i, :, :], N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        vel[i+1, :, :] = vel[i, :, :] .+ dt.*a[i, :, :]
        pos[i+1, :, :] = Temp_Pos .+ 0.5*dt.*vel[i+1, :, :]
    end
end;

function Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    #=---------------------------------------------------------------
    Leap Frog method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    ---------------------------------------------------------------=#
    rho[1, :], P[1, :], a[1, :, :] = Acceleration(pos[1, :, :], vel[1, :, :], N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    Temp_Vel1 = vel[1, :, :]
    for i in 1:T-1
        Temp_Vel2 = Temp_Vel1 .+ dt.*a[i, :, :]
        pos[i+1, :, :] = pos[i, :, :] .+ dt.*Temp_Vel2
        vel[i+1, :, :] = 0.5.*(Temp_Vel1 .+ Temp_Vel2)
        Temp_Vel1 = Temp_Vel2
        rho[i+1, :], P[i+1, :], a[i+1, :, :] = Acceleration(pos[i+1, :, :], vel[i+1, :, :], N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    end
end;