import sys
from PyQt5 import QtWidgets

from model import Model
from presenter import Presenter
from view import View


def main():
    app = QtWidgets.QApplication(sys.argv)
    view = View()
    model = Model()
    presenter = Presenter(view, model)
    view.connect(presenter)
    view.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
