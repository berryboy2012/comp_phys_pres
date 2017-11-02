import numpy as np
import ctypes
import multiprocessing
import threading
import copy
import time
import sys
import json
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from PyQt5.QtWidgets import *
import MainWin_UI, ParticleAddForm_UI, ParticlePresetForm_UI, PlotForm_UI

CDLL = ctypes.CDLL
POINTER = ctypes.POINTER
c_double = ctypes.c_double
c_int = ctypes.c_int
Pipe = multiprocessing.Pipe
Thread = threading.Thread
Process = multiprocessing.Process


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
        self.s.iter_step_alter.argtypes = [
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
        self.tol = self.state.tol
        self.s.init(self.t, self.delta_t, self.tol, self.state.p_num,
                    self.state.p_m.ctypes.data_as(POINTER(c_double)),
                    self.state.p_x.ctypes.data_as(POINTER(c_double)),
                    self.state.p_v.ctypes.data_as(POINTER(c_double)))

    def step_forward(self):
        #print('pingping')
        self.s.iter_step_alter(self.t, self.err, self.energy, self.state.p_num,
                         self.state.p_x.ctypes.data_as(POINTER(c_double)),
                         self.state.p_v.ctypes.data_as(POINTER(c_double))
                         )
        #print('pongpong')
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


class EditParticleForm(AddParticleForm):
    def __init__(self):
        super(EditParticleForm, self).__init__()
        self.setWindowTitle('Edit a particle')
        self.edit_pos = 0


class MyMplCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.axes = fig.add_subplot(111, projection='3d')
        #self.daxe = fig.add_subplot(122, projection='3d')
        self.axes.mouse_init()
        self.compute_initial_figure()
        self.setParent(parent)

    def compute_initial_figure(self):
        pass


class AnimationWidget(QWidget):
    def __init__(self):
        super(AnimationWidget, self).__init__()

        vbox = QVBoxLayout()
        self.canvas = MyMplCanvas(self, width=12, height=8, dpi=100)
        vbox.addWidget(self.canvas)

        hboxL = QVBoxLayout()
        self.ehint_label = QLabel(self)
        self.ehint_label.setText('Energy:')
        self.errhint_label = QLabel(self)
        self.errhint_label.setText('Relative error:')
        self.fpshint_label = QLabel(self)
        self.fpshint_label.setText('Frame time:')
        hboxL.addWidget(self.ehint_label)
        hboxL.addWidget(self.errhint_label)
        hboxL.addWidget(self.fpshint_label)
        hboxR = QVBoxLayout()
        self.error_label = QLabel(self)
        self.energy_label = QLabel(self)
        self.fps_label = QLabel(self)
        hboxR.addWidget(self.energy_label)
        hboxR.addWidget(self.error_label)
        hboxR.addWidget(self.fps_label)
        hbox = QHBoxLayout()
        hbox.addLayout(hboxL)
        hbox.addLayout(hboxR)
        vbox.addLayout(hbox)
        self.setLayout(vbox)


class MainWin(QWidget, MainWin_UI.Ui_MainWin):
    def __init__(self):
        super(MainWin, self).__init__()
        self.setupUi(self)

        self.gstate = GravityState(debug=True)

        self.refresh_mainwin()
        self.show()

        self.add = AddParticleForm()
        self.preset = PresetParticleForm()
        self.editor = EditParticleForm()

        self.add_button.clicked.connect(self.add.show)
        self.add.ok_button.clicked.connect(self.try_add_particle)

        self.del_button.clicked.connect(self.try_del_particle)

        self.edit_button.clicked.connect(self.edit_open)
        self.editor.ok_button.clicked.connect(self.try_edit_particle)

        self.preset_json = None
        self.refresh_preset()
        self.preset_button.clicked.connect(self.open_preset)
        self.preset.preset_list.itemClicked.connect(self.show_preset_description)
        self.preset.use_button.clicked.connect(self.load_preset)

        self.plot = AnimationWidget()
        dots = self.__state2dots()
        self.dots = self.plot.canvas.axes.scatter(*dots[0:2], zs=dots[2])
        self.pa, self.pb = Pipe()
        self.worker = Process(target=calc_worker, args=(self.gstate, 60.0, self.pb))
        #self.phase_graph = self.plot.canvas.daxe.plot()

        self.playing = False
        self.play_button.clicked.connect(self.toggle_play)

    def edit_open(self):
        if not self.part_list.selectedItems():
            return
        plist = []
        for i in self.gstate.get_point():
            plist.append(i)
        p = plist[self.part_list.row(self.part_list.selectedItems()[0])]
        self.editor.mass_edit.setText(str(p['m']))
        self.editor.posx_edit.setText(str(p['x'][0]))
        self.editor.posy_edit.setText(str(p['x'][1]))
        self.editor.posz_edit.setText(str(p['x'][2]))
        self.editor.velx_edit.setText(str(p['v'][0]))
        self.editor.vely_edit.setText(str(p['v'][1]))
        self.editor.velz_edit.setText(str(p['v'][2]))
        self.editor.tag_edit.setText(p['tag'])
        self.editor.edit_pos = self.part_list.row(self.part_list.selectedItems()[0])
        self.editor.show()

    def try_edit_particle(self):
        try:
            m = float(self.editor.mass_edit.text())
            x = [float(self.editor.posx_edit.text()),
                 float(self.editor.posy_edit.text()),
                 float(self.editor.posz_edit.text())]
            v = [float(self.editor.velx_edit.text()),
                 float(self.editor.vely_edit.text()),
                 float(self.editor.velz_edit.text())]
            tag = self.editor.tag_edit.text()
        except:
            self.editor.error_hint.setText('Invalid input')
            return
        self.editor.error_hint.setText('')
        self.gstate.set_point(mp=m, xp=x, vp=v, tp=tag, pos=self.editor.edit_pos)
        self.editor.hide()
        self.refresh_mainwin()

    def __state2dots(self):
        x = []
        y = []
        z = []
        for i in self.gstate.get_point():
            x.append(i['x'][0])
            y.append(i['x'][1])
            z.append(i['x'][2])
        return x, y, z

    def toggle_play(self):
        if self.playing:
            self.playing = False
            return self.stop_playback()
        else:
            self.playing = True
            return self.start_playback()

    def stop_playback(self):
        res = self.pa.recv()
        self.pa.send('STOP')
        self.ani._stop()
        self.gstate = GravityState(m=res['s']['m'], x=res['s']['x'], v=res['s']['v'], tag=self.gstate.p_tag
                                   , tol=self.gstate.tol.value)
        self.part_list.setEnabled(True)
        self.refresh_mainwin()

    def start_playback(self):
        self.pa, self.pb = Pipe()
        self.worker = Process(target=calc_worker, args=(self.gstate, 60.0, self.pb))
        dots = self.__state2dots()
        self.dots, = self.plot.canvas.axes.plot(*dots, 'bo')
        self.part_list.setEnabled(False)
        self.worker.start()
        self.plot.canvas.axes.clear()
        self.ani = animation.FuncAnimation(self.plot.canvas.figure, self.update_plot, blit=True, interval=1000.0 / 60.0)
        self.plot.canvas.axes.relim()
        self.plot.canvas.axes.autoscale_view()
        self.plot.show()

    def update_plot(self, i):
        ct = time.time()
        res = self.pa.recv()
        # self.plot.canvas.axes.can_zoom()
        # self.plot.canvas.axes.clear()
        # self.dots, = self.plot.canvas.axes.plot(*list(zip(*res['s']['x'])), 'bo')
        self.dots.set_data(*list(zip(*res['s']['x']))[:2])
        self.dots.set_3d_properties(list(zip(*res['s']['x']))[2])
        # self.plot.canvas.axes.relim()
        self.plot.error_label.setText(str(res['err']))
        self.plot.energy_label.setText(str(res['E']))
        self.plot.fps_label.setText(str(time.time() - ct))
        return self.dots,

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
        ct = time.time()
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
            print(ct - time.time())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_win = MainWin()
    sys.exit(app.exec_())
