from PyQt5.QtWidgets import *
import test_UI
import sys


app_args = sys.argv
app = QApplication(app_args)
m = QDialog()
ui = test_UI.Ui_Dialog()
ui.setupUi(m)
m.show()
sys.exit(app.exec_())
