from ctypes import *
import numpy as np
import time


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
p_num = 3
p_m = np.empty(p_num, dtype=c_double)
p_x = np.empty((p_num, 3), dtype=c_double)
p_v = np.empty((p_num, 3), dtype=c_double)

# sample data
p_m[:] = np.array([1.0, 1.0, 1.0])
p_x[:] = np.array([[0.0, 3.0**0.5, 0.0],
                   [-1.0, 0.0, 0.0],
                   [1.0, 0.0, 0.0]])
p_v[:] = np.array([[-(0.5)**0.5, 0.0, 0.0],
                   [(0.5)**0.5*0.5, -(1.5)**0.5*0.5, 0.0],
                   [(0.5)**0.5*0.5, (1.5)**0.5*0.5, 0.0]])
h = c_double(1.0 / 60.0)
t = c_double(0.0)
delta_t = c_double(1.0 / 60.0)
err = c_double(0.0)
tol = c_double(1.0E-16)
cur_time = time.process_time()
solver.init(h, t, delta_t, err, tol, c_int(p_num),
            p_m.ctypes.data_as(POINTER(c_double)),
            p_x.ctypes.data_as(POINTER(c_double)),
            p_v.ctypes.data_as(POINTER(c_double)))
for i in range(0, int(60.0 * 200 * 1 / 60.0 / h.value)):
    solver.iter_step(h, t, delta_t, err, tol, c_int(p_num),
                     p_m.ctypes.data_as(POINTER(c_double)),
                     p_x.ctypes.data_as(POINTER(c_double)),
                     p_v.ctypes.data_as(POINTER(c_double)))
print(time.process_time() - cur_time)
