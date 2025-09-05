import numpy as np
import matplotlib.pyplot as plt
import os

data_path = "/home/dan/Documents/repos/particle-fortran/build/data"
plot_particles = 8

colors = [
    'blue',
    'red',
    'green',
    'cyan',
    'purple',
    'darkgoldenrod',
    'pink',
    'black',
    'darkgreen',
    'darkslatgrey',
    'saddlebrown'
]

AU_in_meters = 149597870700


def main():
    # get files
    file_nums = [int(f.split('.')[0]) for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))]
    file_nums.sort()
    files = [os.path.join(data_path, f"{num}.bin") for num in file_nums]

    particles = []
    for file in files:
        file_data = np.fromfile(file, dtype='float64')
        x_vals = file_data[::5]
        y_vals = file_data[1::5]
        px_vals = file_data[2::5]
        py_vals = file_data[3::5]
        r_vals = file_data[4::5]

        while len(x_vals) > len(particles):
            particles.append({'x':[], 'y':[]})

        for i in range(len(x_vals)):
            particles[i]['x'].append(x_vals[i] / AU_in_meters)
            particles[i]['y'].append(y_vals[i] / AU_in_meters)
    
    plt.figure(figsize=(8, 8))
    for particle_index in range(plot_particles):
        plt.plot(particles[particle_index]['x'], particles[particle_index]['y'], color=colors[particle_index])
        plt.scatter(particles[particle_index]['x'][-1], particles[particle_index]['y'][-1], c=colors[particle_index])
    plt.scatter((0), (0), c='orange', s=200)

    plt.xlabel('x (AU)', fontsize=18)
    plt.ylabel('y (AU)', fontsize=18)

    # plt.show()
    plt.savefig('out.png')


if __name__ == '__main__':
    main()
