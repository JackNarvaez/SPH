function Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    for i in 1:T-1
        rho[i, :], P[i, :], a[i, :, :] = Acceleration(pos[i, :, :], vel[i, :, :], N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        vel[i+1, :, :] = vel[i, :, :] .+ 0.5*dt.*a[i, :, :]
        pos[i+1, :, :] = pos[i, :, :] .+ dt.*vel[i+1, :, :]
    end
end;