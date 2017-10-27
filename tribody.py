"""
Deprecated
"""


from ctypes import *
import numpy as np
import time
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

#
# class TribodySolver:
#     def __init__(self):
#         self.solver()
#



fig = plt.figure()
ax = p3.Axes3D(fig)
ax.autoscale()
solver = CDLL('./solver_c.so')
solver.iter_step.argtypes = [
    POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_int),
    POINTER(c_double), POINTER(c_double), POINTER(c_double)
]
solver.init.argtypes = [
    POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_int),
    POINTER(c_double), POINTER(c_double), POINTER(c_double)
]
p_num = 4
dots = []
p_m = np.empty(p_num, dtype=c_double, order='F')
p_x = np.empty((p_num, 3), dtype=c_double, order='F')
p_v = np.empty((p_num, 3), dtype=c_double, order='F')

# sample data
p_m[:] = np.array([1.0, 1.0, 1.0, 0.0])
p_x[:] = np.array([[0.0, 3.0**0.5, 0.0],
                   [-1.0, 0.0, 0.0],
                   [1.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0]])
p_v[:] = np.array([[-(0.5)**0.5, 0.0, 0.0],
                   [(0.5)**0.5*0.5, -(1.5)**0.5*0.5, 0.0],
                   [(0.5)**0.5*0.5, (1.5)**0.5*0.5, 0.0],
                   [0.0, 0.0, 0.0]])
h = c_double(1.0 / 60.0)
t = c_double(0.0)
delta_t = c_double(1.0 / 60.0)
err = c_double(0.0)
tol = c_double(1.0E-16)
cur_time = time.process_time()


def fort_iter_init():
    global dots
    solver.init(h, t, delta_t, err, tol, c_int(p_num),
                p_m.ctypes.data_as(POINTER(c_double)),
                p_x.ctypes.data_as(POINTER(c_double)),
                p_v.ctypes.data_as(POINTER(c_double)))
    temp = p_x.T
    dots, = ax.plot(temp[0], temp[1], temp[2], 'bo')
    return dots,


#fort_iter_init()


def fort_iter_run():
    solver.iter_step(h, t, delta_t, err, tol, c_int(p_num),
                     p_m.ctypes.data_as(POINTER(c_double)),
                     p_x.ctypes.data_as(POINTER(c_double)),
                     p_v.ctypes.data_as(POINTER(c_double)))
    # for a_dot in p_x:
    #     print(*a_dot, sep=' ', end='    ')
    # print('')
    return p_x.T


def plot_iter(i):
    global dots
    dot_list = fort_iter_run()
    dots.set_data(dot_list[0], dot_list[1])
    dots.set_3d_properties(dot_list[2])
    return dots,


ani = animation.FuncAnimation(fig, plot_iter, interval=1000 / 60, blit=True, init_func=fort_iter_init)
plt.show()

