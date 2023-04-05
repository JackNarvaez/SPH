using Random

function Init_Dis(N, R, x0, y0, seed =1234)
    #=---------------------------------------------------------------
    Random distribution of particles inside a circle of ratio R that
    is centered at (x0, y0)
    -----------------------------------------------------------------
    Arguments:
    N:  Number of particles
    R:  Radius
    (x0, y0) = Center of circle
    -----------------------------------------------------------
    Return:
    [x y]: Random positions of N particles.
    ---------------------------------------------------------------=#
    Random.seed!(1234)
    r = R*sqrt.(rand(Float64, N))
    theta = 2*pi*rand(Float64, N)
    x = r.*cos.(theta) .+ x0
    y = r.*sin.(theta) .+ y0
    return [x y]
end;

function Gaussian_Kernel(x, y, h)
    #=---------------------------------------------------------------
    Gaussian kernel in 2D 
    ---------------------------------------------------------------=#
    r2 = x^2 + y^2
    H = h^2
    W = 1.0/(H*pi)*exp(-r2/H)
    return W
end;

function Gradient_Gaussian_Kernel(x, y, h)
    #=---------------------------------------------------------------
    Gradient of the Gaussian kernel in 2D
    ---------------------------------------------------------------=#
    r2 = x^2 + y^2
    n = -2/(h^4*pi)*exp(-r2/h^2)
    return [n*x n*y]
end;

function Density(dx, dy, m, h, Kernel)
    #=---------------------------------------------------------------
    Mass density of particles of mass m
    ---------------------------------------------------------------=#
    rho = m*Kernel(dx, dy, h)
    return rho
end;

function Pressure(rho, k, n)
    #=---------------------------------------------------------------
    Polytropic equation of state
    -----------------------------------------------------------------
    Arguments:
    rho:  Density
    k:    Pressure constant
    n:    Polytropic index
    -----------------------------------------------------------------
    Return:
    P:    Pressure
    ---------------------------------------------------------------=#
    P = k * rho^(1+1.0/n)
    return P
end;

function Coeff_static_grav_potential(k, n, M, R, d=2)
    #=---------------------------------------------------------------
    Calculate the coefficient of static gravity potential
    -----------------------------------------------------------------
    Arguments:
    k:    Pressure constant
    n:    Polytropic index
    M:    Star mass
    R:    Star Radius
    d:    Dimensions
    -----------------------------------------------------------------
    Return:
    lmbda: coefficient of static gravity potential in 2D
    ---------------------------------------------------------------=#
    lmbda = 2*k/(pi^(1/n)) * (M*(1+n)/R^2)^(1+1/n) / M
    return lmbda
end;

function Acceleration(pos, vel, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    #=---------------------------------------------------------------
    SPH approximation of acceleration
    -----------------------------------------------------------------
    Arguments:
    pos:     Position of particles. Size: N x 2
    vel:     Velocity of particles. Size: N x 2
    N:       Number of particles.
    k:       Pressure constant
    n:       Polytropic index
    lmbda:   Coefficient of static gravity potential
    nu:      Viscocity coefficient
    m:       Mass particle
    h:       Smoothing kernel length
    -----------------------------------------------------------
    Return:
    rho:     Density at each particle
    P:       Pressure
    a:       Acceleration
    ---------------------------------------------------------=#
    
    rho = ones(N)*Density(0, 0, m, h, Kernel)
    P = zeros(N)
    a = -lmbda.*pos .- nu.*vel
    for i in 1:N
        dr = pos[i, :]' .- pos[i+1:end, :]
        for j = (i+1):N
            if (dr[j-i, 1]^2 + dr[j-i, 2]^2)/h^2 <= 4
                rho_ij = Density(dr[j-i, 1], dr[j-i, 2], m, h, Kernel)
                rho[i] += rho_ij
                rho[j] += rho_ij
            end
        end
    end
    for i in 1:N
        P[i] = Pressure(rho[i], k, n)
    end
    for i in 1:N
        dr = pos[i, :]' .- pos[i+1:end, :]
        for j in i+1:N
            if (dr[j-i, 1]^2 + dr[j-i, 2]^2)/h^2 <= 4
                acc_ij = (m*(P[i]/rho[i]^2 + P[j]/rho[j]^2)).*Gradient_Kernel(dr[j-i, 1], dr[j-i, 2], h)'
                a[i, :] -= acc_ij
                a[j, :] += acc_ij
            end
        end
    end
    return rho, P, a
end;