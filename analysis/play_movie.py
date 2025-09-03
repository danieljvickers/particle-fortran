import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from tqdm import tqdm

plt.style.use('dark_background')

data_path = "/home/dan/Documents/repos/particle-fortran/build/data"

AU_in_meters = 149597870700
marker_size_scale = 0.05


def main():
    # get files
    file_nums = [int(f.split('.')[0]) for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))]
    file_nums.sort()
    files = [os.path.join(data_path, f"{num}.bin") for num in file_nums]

    particles = []
    norm_factor = 0
    for file in files:
        file_data = np.fromfile(file, dtype='float64')
        x_vals = file_data[::5]
        y_vals = file_data[1::5]
        px_vals = file_data[2::5]
        py_vals = file_data[3::5]
        r_vals = file_data[4::5]
        if not norm_factor:
            norm_factor = np.mean(r_vals)

        while len(x_vals) > len(particles):
            particles.append({'x':[], 'y':[], 'r': []})

        for i in range(len(x_vals)):
            particles[i]['x'].append(x_vals[i] / AU_in_meters)
            particles[i]['y'].append(y_vals[i] / AU_in_meters)
            particles[i]['r'].append(marker_size_scale * (r_vals[i] / norm_factor)**3)

            if len(particles[i]['x']) < len(particles[0]['x']):
                print(f"Found {len(particles[i]['x'])} time steps at particle {i}")
        

    def animation_function(i):
        plot_data = []
        radii = []
        for j in range(len(particles)):
            if len(particles[j]['x']) > i:
                plot_data.append((particles[j]['x'][i], particles[j]['y'][i]))
                radii.append(particles[j]['r'][i])
        scatter.set_offsets(plot_data)
        scatter.set_sizes(radii)
    
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    plt.xlabel('x (AU)', fontsize=18)
    plt.ylabel('y (AU)', fontsize=18)
    plt.xlim([-2, 2])
    plt.ylim([-2, 2])

    scatter = ax.scatter([], [], c='white', s=0.5)
    plt.scatter((0), (0), c='orange', s=200)

    ani = animation.FuncAnimation(fig, animation_function, repeat=False, frames=len(files))
    writer = animation.FFMpegWriter(fps=30)
    with tqdm(total=len(files), desc='Saving video') as progress_bar:
        ani.save('out.mp4', writer=writer, progress_callback=lambda i, n: progress_bar.update(1))


if __name__ == '__main__':
    main()
