import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from tqdm import tqdm

plt.style.use('dark_background')

data_path = "/home/dan/Documents/repos/particle-fortran/build/data"

AU_in_meters = 149597870700
marker_size_scale = 4000


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
            particles.append({'x':[], 'y':[], 's': []})

        for i in range(len(x_vals)):
            particles[i]['x'].append(x_vals[i] / AU_in_meters)
            particles[i]['y'].append(y_vals[i] / AU_in_meters)
            particles[i]['s'].append(marker_size_scale * (np.pi * (r_vals[i] / AU_in_meters)**2))

            if len(particles[i]['x']) < len(particles[0]['x']):
                print(f"Found {len(particles[i]['x'])} time steps at particle {i}")
        

    def animation_function(i):
        ax.cla()
        plot_x = []
        plot_y = []
        area = []
        for j in range(len(particles)):
            if len(particles[j]['x']) > i:
                plot_x.append(particles[j]['x'][i])
                plot_y.append(particles[j]['y'][i])
                area.append(particles[j]['s'][i])
        plot_z = np.zeros_like(plot_x)
        ax.scatter(plot_x, plot_y, plot_z, c='white', s=area)
        ax.scatter([0], [0], [0], c='orange', s=200)

        ax.set_xlabel('x (AU)', fontsize=18)
        ax.set_ylabel('y (AU)', fontsize=18)
        ax.set_zlabel('z (AU)', fontsize=18)
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        ax.set_zlim([-0.4, 0.4])
        ax.grid(False)
        ax.view_init(elev=30, azim=45)
        ax.xaxis.pane.set_facecolor("black")
        ax.yaxis.pane.set_facecolor("black")
        ax.zaxis.pane.set_facecolor("black")

        ax.xaxis.pane.set_edgecolor("black")
        ax.yaxis.pane.set_edgecolor("black")
        ax.zaxis.pane.set_edgecolor("black")
        
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect((5, 5, 1))

    scatter = ax.scatter([], [], [], c='white', s=0.5)
    plt.scatter((0), (0), c='orange', s=200)

    ani = animation.FuncAnimation(fig, animation_function, repeat=False, frames=len(files))
    writer = animation.FFMpegWriter(fps=60)
    with tqdm(total=len(files), desc='Saving video') as progress_bar:
        ani.save('out.mp4', writer=writer, progress_callback=lambda i, n: progress_bar.update(1))


if __name__ == '__main__':
    main()
