import numpy as np
import ctypes
import multiprocessing
import threading
import copy
import time
import sys
from PyQt5.QtWidgets import *
import MainWin_UI

CDLL = ctypes.CDLL
POINTER = ctypes.POINTER
c_double = ctypes.c_double
c_int = ctypes.c_int
Pipe = multiprocessing.Pipe
Thread = threading.Thread
Process = multiprocessing.Process


class GravityState:
    def __init__(self, m=None, x=None, v=None, debug=False):
        self.debug = debug
        if m is None or x is None or v is None:
            if self.debug:
                m = [1.0, 1.0, 1.0, 0.0]
                x = [[0.0, 3.0 ** 0.5, 0.0],
                     [-1.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0],
                     [0.0, 0.0, 0.0]]
                v = [[-0.5 ** 0.5, 0.0, 0.0],
                     [0.5 ** 0.5 * 0.5, -1.5 ** 0.5 * 0.5, 0.0],
                     [0.5 ** 0.5 * 0.5, 1.5 ** 0.5 * 0.5, 0.0],
                     [0.0, 0.0, 0.0]]
        p_num = c_int(len(m))
        p_m = np.empty(p_num.value, dtype=c_double, order='F')
        p_x = np.empty((p_num.value, 3), dtype=c_double, order='F')
        p_v = np.empty((p_num.value, 3), dtype=c_double, order='F')
        p_m[:] = np.array(m)
        p_x[:] = np.array(x)
        p_v[:] = np.array(v)
        self.p_m = p_m
        self.p_x = p_x
        self.p_v = p_v
        self.p_num = p_num

    def set_point(self, mp=None, xp=None, vp=None, pos=-1):
        if pos == -1:
            self.p_m = np.append(self.p_m, [mp], axis=0)
            self.p_x = np.append(self.p_x, [xp], axis=0)
            self.p_v = np.append(self.p_v, [vp], axis=0)
            self.p_num = c_int(self.p_num.value + 1)
        else:
            if mp is not None:
                self.p_m[pos] = mp
            if xp is not None:
                self.p_x[pos] = xp
            if vp is not None:
                self.p_v[pos] = vp

    def get_point(self, pos):
        return {'m': float(self.p_m[pos]),
                'x': self.p_x[pos].tolist(),
                'v': self.p_v[pos].tolist()}


class FortranSolver:
    def __init__(self):
        self.s = CDLL('./solver_c_new.so')
        self.s.iter_step.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_int),
            POINTER(c_double), POINTER(c_double)
        ]
        self.s.init.argtypes = [
            POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_int),
            POINTER(c_double), POINTER(c_double), POINTER(c_double)
        ]

    def set_param(self, state=None, step_len=1.0 / 60.0, t=0.0, tolerance=1.0E-16):
        if state is None:
            self.state = GravityState(debug=True)
        else:
            self.state = state
        self.t = c_double(t)
        self.delta_t = c_double(step_len)
        self.err = c_double(0.0)
        self.energy = c_double(0.0)
        self.tol = c_double(tolerance)
        self.s.init(self.t, self.delta_t, self.tol, self.state.p_num,
                    self.state.p_m.ctypes.data_as(POINTER(c_double)),
                    self.state.p_x.ctypes.data_as(POINTER(c_double)),
                    self.state.p_v.ctypes.data_as(POINTER(c_double)))

    def step_forward(self):
        self.s.iter_step(self.t, self.err, self.energy, self.state.p_num,
                         self.state.p_x.ctypes.data_as(POINTER(c_double)),
                         self.state.p_v.ctypes.data_as(POINTER(c_double))
                         )
        return {'s': {'m': self.state.p_m.tolist(),
                      'x': self.state.p_x.tolist(),
                      'v': self.state.p_v.tolist()},
                't': self.t.value,
                'err': self.err.value,
                'E': self.energy.value}


def calc_worker(state, tol, fps, pipe):
    fort_worker = FortranSolver()
    fort_worker.set_param(state=state, step_len=1.0 / fps, tolerance=tol)
    stop = False
    pause = False
    while not stop:
        if pipe.poll():
            req = pipe.recv()
            if req == 'STOP':
                stop = True
            elif req == 'PAUSE':
                pause = True
            elif req == 'CONTINUE':
                pause = False
        if not pause:
            pipe.send(fort_worker.step_forward())
        else:
            time.sleep(0.1)


if __name__ == '__main__':
    pa, pb = Pipe()
    t = Process(target=calc_worker, args=(GravityState(debug=True), 1.0E-6, 60.0, pb))
    t.start()
    for i in range(int(60 * 10.28)):
        res = pa.recv()
        print(res['t'], res['E'], res['err'], end=' ')
        for xi in res['s']['x']:
            print(*xi, end=' ')
        for vi in res['s']['v']:
            print(*vi, end=' ')
        print()
    pa.send('STOP')
