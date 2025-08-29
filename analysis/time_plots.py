import numpy as np
import matplotlib.pyplot as plt
import os

cases = [256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144]
cpu_time = [
    5.0000000000000002E-005,
    3.0000000000000003E-004,
    6.0000000000000006E-004,
    1.7500000000000003E-003,
    5.7499999999999999E-003,
    2.1350000000000001E-002,
    8.3699999999999997E-002,
    0.33424999999999999,
    1.3162499999999999,
    5.2058499999999999,
    20.434049999999999
]

def main():
    plt.figure(figsize=(8, 8))
    plt.loglog(cases, cpu_time, linewidth=3)
    plt.xlabel("# particles", fontsize=18)
    plt.ylabel("time (s)", fontsize=18)
    plt.show()

if __name__ == '__main__':
    main()