from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

from dialog import Ui_Dialog


# dendrogram dialog
class Dendrogram(QtWidgets.QDialog):
    def __init__(self, parent, linkageTable, labels, colorThreshold):
        super().__init__(parent)
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # return value
        self.threshold = None

        # create dendrogram canvas
        self.figureDendrogram = plt.figure()
        self.axesDendrogram = self.figureDendrogram.add_subplot(111)
        self.canvasDendrogram = FigureCanvas(self.figureDendrogram)
        # self.toolbarDendrogram = NavigationToolbar(self.canvasDendrogram, self)
        self.ui.frameDendrogram.layout().addWidget(self.canvasDendrogram)
        # self.ui.frameDendrogram.layout().addWidget(self.toolbarDendrogram)
        self.axesDendrogram.cla()

        # draw dendrogram
        hierarchy.dendrogram(linkageTable, ax=self.axesDendrogram, labels=labels, color_threshold=colorThreshold)
        self.axesDendrogram.set_xlabel("Motif Length")
        self.axesDendrogram.set_ylabel("Cluster Distance")
        self.line = self.axesDendrogram.axhline(color='red')
        self.line.set_ydata(colorThreshold)
        self.figureDendrogram.tight_layout()
        self.canvasDendrogram.draw()

        # connect moving line to cursor event
        self.canvasDendrogram.mpl_connect("motion_notify_event", self.mouseMove)

        # connect button press
        self.canvasDendrogram.mpl_connect("button_press_event", self.mousePress)

    def mouseMove(self, event):
        if event.inaxes:
            self.line.set_ydata(event.ydata)
            self.canvasDendrogram.draw()

    def mousePress(self, event):
        if event.inaxes:
            self.threshold = event.ydata
            self.accept()

    def getSelectedThreshold(self):
        return self.threshold

    def accept(self):
        plt.close(self.figureDendrogram)
        super().accept()

    def reject(self):
        plt.close(self.figureDendrogram)
        super().reject()
