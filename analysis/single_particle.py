import numpy as np
import matplotlib.pyplot as plt
import os

data_path = "/home/dan/Documents/repos/particle-fortran/build/data"
particle_index = 0


def main():
    # get files
    file_nums = [int(f.split('.')[0]) for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))]
    file_nums.sort()
    files = [os.path.join(data_path, f"{num}.bin") for num in file_nums]

    particles = []
    for file in files:
        file_data = np.fromfile(file, dtype='float64')
        x_vals = file_data[::2]
        y_vals = file_data[1::2]

        while len(x_vals) > len(particles):
            particles.append({'x':[], 'y':[]})

        for i in range(len(x_vals)):
            particles[i]['x'].append(x_vals[i])
            particles[i]['y'].append(y_vals[i])
    
    plt.plot(particles[particle_index]['x'], particles[particle_index]['y'])
    plt.scatter((0), (0), c='orange', s=200)
    plt.scatter(particles[particle_index]['x'][-1], particles[particle_index]['y'][-1], c='black')
    plt.show()


if __name__ == '__main__':
    main()
