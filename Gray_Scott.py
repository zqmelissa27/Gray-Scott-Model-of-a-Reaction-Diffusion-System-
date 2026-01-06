import os
import numpy as np
import matplotlib.pyplot as plt
import imageio

def gray_scott_model(Du, Dv, f, k, Lx, Ly, mx, my, NT, dt, St, output_dir, cmap):
    N = mx  #mx = my
    V = np.ones((N, N))
    U = np.zeros((N, N))
    dU = np.zeros((N, N))
    dV = np.zeros((N, N))

    def initial_state(U, V):
        U.fill(1.0)
        V.fill(0.0)
        for i in range(N // 3, 2 * N // 3):
            for j in range(N // 3, 2 * N // 3):
                U[i, j] = 0.5 * (1 + np.random.uniform(-1, 1))
                V[i, j] = 0.25 * (1 + np.random.uniform(-1, 1))
        return U, V

    def timestep(U, V, dU, dV, Du, Dv, f, k):
        for i in range(N):
            for j in range(N):
                u = U[i, j]
                v = V[i, j]
                lapU = (U[(i-1)%N,j] + U[(i+1)%N,j] + U[i,(j-1)%N] + U[i,(j+1)%N] - 4*u)
                lapV = (V[(i-1)%N,j] + V[(i+1)%N,j] + V[i,(j-1)%N] + V[i,(j+1)%N] - 4*v)
                dU[i, j] = Du * lapU - u * v**2 + f * (1 - u)
                dV[i, j] = Dv * lapV + u * v**2 - (f + k) * v

        U += dU * dt
        V += dV * dt
        return U, V


    U, V = initial_state(U, V)

    #Create output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    images = []

    for t in range(NT):
        U, V = timestep(U, V, dU, dV, Du, Dv, f, k)

        if t % int(St / dt) == 0:
            plt.clf()
            plt.imshow(1 - U, cmap=cmap)
            plt.title(f'Simulación Gray-Scott\nIteración: {t}, F = {f}, K = {k}')
            plt.axis('off')
            image_filename = f"{output_dir}/frame_{t:04d}.png"
            plt.savefig(image_filename)
            images.append(imageio.imread(image_filename))

    
    gif_filename = os.path.join(output_dir, 'gray_scott_simulation.gif')
    imageio.mimsave(gif_filename, images, format='GIF', duration=2)

    print(f"Simulación completada. GIF guardado como {gif_filename}")

simulaciones = [
   # {"Du": 0.16, "Dv": 0.08, "f": 0.035, "k": 0.060, "cmap": 'plasma'},
  #  {"Du": 0.21, "Dv": 0.105, "f": 0.034, "k": 0.056, "cmap": 'inferno'},
    #{"Du": 0.21, "Dv": 0.105, "f": 0.039, "k": 0.058, "cmap": 'magma'},
   # {"Du": 0.21, "Dv": 0.105, "f": 0.026, "k": 0.051, "cmap": 'viridis'},
   # {"Du": 0.16, "Dv": 0.08, "f": 0.060, "k": 0.062, "cmap": 'Spectral'},
    {"Du": 0.21, "Dv": 0.105, "f": 0.029, "k": 0.057, "cmap": 'twilight'},
   # {"Du": 0.14, "Dv": 0.06, "f": 0.035, "k": 0.065, "cmap": 'coolwarm'},
   # {"Du": 0.21, "Dv": 0.105, "f": 0.025, "k": 0.06, "cmap": 'terrain'},
  #  {"Du": 0.21, "Dv": 0.105, "f": 0.030, "k": 0.060, "cmap": 'Spectral'},
]


Lx = 100.0  # Longitud en el eje x
Ly = 100.0  # Longitud en el eje y
mx = 350 # Número de celdas en el eje x
my = 350  # Número de celdas en el eje y
NT = 5050  # Número de iteraciones
dt = 1.0  # Paso de tiempo
St = 50  # Guardar cada St s

output_base_dir = 'gray_scott_simulations'
for i, params in enumerate(simulaciones):
    output_dir = f'{output_base_dir}/simulation_{i+1}'
    print(f"Ejecutando simulación {i+1} con parámetros: Du={params['Du']}, Dv={params['Dv']}, f={params['f']}, k={params['k']}, cmap={params['cmap']}")

    gray_scott_model(params['Du'], params['Dv'], params['f'], params['k'], Lx, Ly, mx, my, NT, dt, St, output_dir, params['cmap'])
