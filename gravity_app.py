import numpy as np
import ctypes
import multiprocessing
import threading
import copy
import time
import sys
import json
from PyQt5.QtWidgets import *
import MainWin_UI, ParticleAddForm_UI, ParticlePresetForm_UI

CDLL = ctypes.CDLL
POINTER = ctypes.POINTER
c_double = ctypes.c_double
c_int = ctypes.c_int
Pipe = multiprocessing.Pipe
Thread = threading.Thread


class GravityState:
    def __init__(self, m=None, x=None, v=None, tag=None, tol=1E-3, debug=False):
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
                tag = ['A', 'B', 'C', 'D']
        p_num = c_int(len(m))
        p_m = np.empty(p_num.value, dtype=c_double, order='F')
        p_x = np.empty((p_num.value, 3), dtype=c_double, order='F')
        p_v = np.empty((p_num.value, 3), dtype=c_double, order='F')
        p_tag = [i for i in tag]
        p_m[:] = np.array(m)
        p_x[:] = np.array(x)
        p_v[:] = np.array(v)
        self.p_m = p_m
        self.p_x = p_x
        self.p_v = p_v
        self.p_tag = p_tag
        self.p_num = p_num
        self.tol = c_double(tol)

    def set_point(self, mp=None, xp=None, vp=None, tp=None, pos=-1):
        if pos == -1:
            self.p_m = np.append(self.p_m, [mp], axis=0)
            self.p_x = np.append(self.p_x, [xp], axis=0)
            self.p_v = np.append(self.p_v, [vp], axis=0)
            self.p_tag.append(tp)
            self.p_num = c_int(self.p_num.value + 1)
        else:
            if mp is not None:
                self.p_m[pos] = mp
            if xp is not None:
                self.p_x[pos] = xp
            if vp is not None:
                self.p_v[pos] = vp
            if tp is not None:
                self.p_tag[pos] = tp

    def deled_point(self, pos=-1):
        if pos == -1:
            pass
        else:
            m = []
            x = []
            v = []
            tag = []
            for i in range(self.p_num.value):
                if pos == i:
                    pass
                else:
                    m.append(float(self.p_m[i]))
                    x.append(self.p_x[i].tolist())
                    v.append(self.p_v[i].tolist())
                    tag.append(self.p_tag[i])
            return GravityState(m=m, x=x, v=v, tag=tag, debug=self.debug)

    def get_point(self):
        for i in range(self.p_num.value):
            yield {'m': float(self.p_m[i]),
                   'x': self.p_x[i].tolist(),
                   'v': self.p_v[i].tolist(),
                   'tag': self.p_tag[i]}


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

    def set_param(self, state=None, step_len=1.0 / 60.0, t=0.0):
        if state is None:
            self.state = GravityState(debug=True)
        else:
            self.state = state
        self.h = c_double(step_len)
        self.t = c_double(t)
        self.delta_t = self.h
        self.err = c_double(0.0)
        self.energy = c_double(0.0)
        self.tol = c_double(self.state.tol)
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


class PresetParticleForm(QWidget, ParticlePresetForm_UI.Ui_ParticlePresetForm):
    def __init__(self):
        super(PresetParticleForm, self).__init__()
        self.setupUi(self)


class AddParticleForm(QWidget, ParticleAddForm_UI.Ui_ParticleAddForm):
    def __init__(self):
        super(AddParticleForm, self).__init__()
        self.setupUi(self)


class MainWin(QWidget, MainWin_UI.Ui_MainWin):
    def __init__(self):
        super(MainWin, self).__init__()
        self.setupUi(self)
        self.gstate = GravityState(debug=True)
        self.refresh_mainwin()
        self.show()
        self.add = AddParticleForm()
        self.preset = PresetParticleForm()
        self.add_button.clicked.connect(self.add.show)
        self.add.ok_button.clicked.connect(self.try_add_particle)
        self.del_button.clicked.connect(self.try_del_particle)
        self.preset_json = None
        self.refresh_preset()
        self.preset_button.clicked.connect(self.open_preset)
        self.preset.preset_list.itemClicked.connect(self.show_preset_description)
        self.preset.use_button.clicked.connect(self.load_preset)

    def refresh_preset(self):
        self.preset.preset_list.clear()
        try:
            self.preset_json = json.load(open('presets.json', encoding='utf-8'))
        except FileNotFoundError:
            self.error_hint.setText('No presets found!')
            return
        except:
            self.error_hint.setText('Errors in \'presets.json\'!')
            return
        for i in self.preset_json:
            self.preset.preset_list.addItem(i['name'])

    def open_preset(self):
        self.refresh_preset()
        self.preset.show()

    def show_preset_description(self):
        if self.preset.preset_list.selectedItems():
            self.preset.detail_label.setText(
                self.preset_json[self.preset.preset_list.row(self.preset.preset_list.selectedItems()[0])][
                    'description'])

    def load_preset(self):
        if self.preset.preset_list.selectedItems():
            preset = self.preset_json[self.preset.preset_list.row(self.preset.preset_list.selectedItems()[0])]
            self.gstate = GravityState(m=preset['m'], x=preset['x'], v=preset['v'], tag=preset['tag'])
            self.preset.hide()
            self.refresh_mainwin()

    def try_del_particle(self):
        if self.part_list.selectedItems():
            if self.gstate.p_num.value == 1:
                self.error_hint.setText('At least have one particle!')
                return
            else:
                self.error_hint.setText('')
                self.gstate = self.gstate.deled_point(pos=self.part_list.row(self.part_list.selectedItems()[0]))
                self.part_list.takeItem(self.part_list.row(self.part_list.selectedItems()[0]))

    def refresh_mainwin(self):
        self.part_list.clear()
        self.error_hint.setText('')
        for i in self.gstate.get_point():
            self.part_list.addItem('%s;M:%s;X:[%s,%s,%s];V:[%s,%s,%s]' % (i['tag'], i['m'], *i['x'], *i['v']))

    def try_add_particle(self):
        try:
            m = float(self.add.mass_edit.text())
            x = [float(self.add.posx_edit.text()),
                 float(self.add.posy_edit.text()),
                 float(self.add.posz_edit.text())]
            v = [float(self.add.velx_edit.text()),
                 float(self.add.vely_edit.text()),
                 float(self.add.velz_edit.text())]
            tag = self.add.tag_edit.text()
        except:
            self.add.error_hint.setText('Invalid input')
            return
        self.add.error_hint.setText('')
        self.gstate.set_point(mp=m, xp=x, vp=v, tp=tag)
        self.add.hide()
        self.refresh_mainwin()


def calc_worker(state, fps, pipe):
    fort_worker = FortranSolver()
    fort_worker.set_param(state=state, step_len=1.0 / fps)
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
    app = QApplication(sys.argv)
    main_win = MainWin()
    sys.exit(app.exec_())
