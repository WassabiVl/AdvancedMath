def heatEquation(nx, ny):
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            uxx = (u0[i + 1, j] - 2 * u0[i, j] + u0[i - 1, j]) / dx2
            uyy = (u0[i, j + 1] - 2 * u0[i, j] + u0[i, j - 1]) / dy2
            u[i, j] = u0[i, j] + dt * D * (uxx + uyy)
