import subprocess
from pathlib import Path
from PyQt5.QtCore import QObject, pyqtSignal, QThread
from scipy import stats

from utils import loadTS, getMotifPoints, getRanges
from model import Model
from view import View


class Presenter:
    def __init__(self, view: View, model: Model):
        self.view = view
        self.model = model

        self.runThread = None
        self.worker = None

    def openFile(self):
        fileName = self.view.openFileNameDialog()
        self.loadData(fileName)

    def loadData(self, filePath):
        if not self.model.running and filePath is not None:  # TODO feedback if running; check if it is a file
            self.model.tsPath = filePath
            self.model.ts = loadTS(filePath)  # TODO catch exceptions
            self.model.clusteredMotifs = []
            self.model.motifs = []

            self.view.plotTs(self.model.ts, Path(filePath).name)
            self.view.setSliderRanges(10, len(self.model.ts) // 2)
            self.view.clearMotifPlot()
            self.view.clearMotifTable()

            self.view.setProgress(0, "ready")
            self.view.setRunnable(True)

    def setRunning(self, running):
        self.view.setRunning(running)
        self.model.running = running

    def clearMotifs(self):
        self.view.showRanges([])
        self.view.clearMotifPlot()
        self.view.clearMotifTable()

    def run(self):
        if not self.model.running:
            # clear
            self.clearMotifs()
            self.model.clusteredMotifs = []
            self.model.motifs = []
            # load settings
            self.model.settingsCurrentRun = self.view.getSettings()
            minLength, maxLength, correlation = self.model.settingsCurrentRun
            # run
            self.setRunning(True)
            self.view.setProgress(0, "started")
            # create thread
            self.runThread = QThread()
            self.worker = Worker(self.model.tsPath, minLength, maxLength, correlation)
            self.worker.moveToThread(self.runThread)
            # setup
            self.runThread.started.connect(self.worker.run)
            self.worker.finished.connect(self.runThread.quit)
            self.worker.finished.connect(self.worker.deleteLater)
            self.runThread.finished.connect(self.runThread.deleteLater)
            # connect signals and slots
            self.worker.progress[int].connect(lambda value: self.view.setProgress(value=value))
            self.worker.progress[str].connect(lambda text: self.view.setProgress(text=text))
            self.worker.progress[int, str].connect(lambda value, text: self.view.setProgress(value, text))
            self.worker.motif.connect(lambda motif, window: self.addMotif(motif, window))
            self.worker.cluster.connect(lambda motif, window: self.addCluster(motif, window))
            self.runThread.finished.connect(self.finishRun)
            # start thread
            self.runThread.start()
        else:
            # cancel
            self.setRunning(False)
            self.view.setProgress(text="cancelled")
            self.runThread.join()  # TODO stop

    def addMotif(self, motif, window):
        self.model.motifs.append((motif, window))
        if not self.model.showClusters:
            self.view.addMotifRow(window, len(motif))

    def addCluster(self, motif, window):
        self.model.clusteredMotifs.append((motif, window))
        if self.model.showClusters:
            self.view.addMotifRow(window, len(motif))

    def finishRun(self):
        self.view.setProgress(100, "finished")
        self.setRunning(False)

    def showMotif(self):
        motifId = self.view.getSelectedMotif()
        if motifId is not None:
            if self.model.showClusters:
                motif, window = self.model.clusteredMotifs[motifId]
            else:
                motif, window = self.model.motifs[motifId]

            ranges = getRanges(getMotifPoints(motif, window))
            self.view.showRanges(ranges)

            subsequences = []
            mean = [0] * window
            for i in motif:
                subsequence = stats.zscore(self.model.ts[i:i + window])
                subsequences.append(subsequence)
                mean += subsequence
            mean /= len(motif)
            self.view.plotMotif(mean, subsequences)

    def save(self):
        pass  # TODO

    def switchMotifTable(self):
        # clear
        self.clearMotifs()

        self.model.showClusters = not self.model.showClusters
        if self.model.showClusters:
            self.view.setClusterButtonText("Show All")
            for entry in self.model.clusteredMotifs:
                motif, window = entry
                self.view.addMotifRow(window, len(motif))
        else:
            self.view.setClusterButtonText("Filter")
            for entry in self.model.motifs:
                motif, window = entry
                self.view.addMotifRow(window, len(motif))


class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal([int], [str], [int, str])
    motif = pyqtSignal(list, int)
    cluster = pyqtSignal(list, int)

    def __init__(self, tsPath, windowMin, windowMax, correlation):
        super().__init__()
        self.tsPath = tsPath
        self.windowMin = windowMin
        self.windowMax = windowMax
        self.correlation = correlation

    def run(self):
        runs = self.windowMax - self.windowMin + 1
        run = 0
        currentCluster = []
        # get motifs with VLCM
        args = ["../vlcm/build/VLCM", self.tsPath, str(self.windowMin), str(self.windowMax),
                str(self.correlation)]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, text=True)
        for line in process.stdout:
            splits = line.split(";")
            w = int(splits[0])
            motif = sorted([int(i) for i in splits[1].split(",")])
            if len(motif) > 1:
                self.motif.emit(motif, w)
                if len(currentCluster) > 0 and motifSimilarity(motif, w, currentCluster[-1][0], w - 1) < 0.8:  # TODO
                    # cluster done
                    maxMotif = (0, -1)  # cover, i
                    for i in range(len(currentCluster)):
                        cover = len(currentCluster[i][0]) * currentCluster[i][1]
                        if cover >= maxMotif[0]:
                            maxMotif = (cover, i)
                    self.cluster.emit(*currentCluster[maxMotif[1]])
                    currentCluster = []
                else:
                    currentCluster.append((motif, w))
            run += 1
            self.progress[int].emit(int(100.0 * run / runs))
        return_code = process.wait()
        if return_code:  # TODO handle error
            raise subprocess.CalledProcessError(return_code, args)
        self.finished.emit()


def motifSimilarity(m1, w1, m2, w2):
    pointerM1 = 0
    pointerM2 = 0

    onlyM1 = 0
    onlyM2 = 0
    shared = 0

    currM1range = [m1[0], m1[0] + w1]
    currM2range = [m2[0], m2[0] + w2]

    while pointerM1 < len(m1) or pointerM2 < len(m2):
        if pointerM1 >= len(m1):
            onlyM2 += currM2range[1] - currM2range[0]
            pointerM2 += 1
            if pointerM2 >= len(m2):
                break
            currM2range = [m2[pointerM2], m2[pointerM2] + w2]
            continue
        if pointerM2 >= len(m2):
            onlyM1 += currM1range[1] - currM1range[0]
            pointerM1 += 1
            if pointerM1 >= len(m1):
                break
            currM1range = [m1[pointerM1], m1[pointerM1] + w1]
            continue
        if currM1range[1] <= currM2range[0]:
            onlyM1 += currM1range[1] - currM1range[0]
            pointerM1 += 1
            if pointerM1 < len(m1):
                currM1range = [m1[pointerM1], m1[pointerM1] + w1]
            continue
        if currM2range[1] <= currM1range[0]:
            onlyM2 += currM2range[1] - currM2range[0]
            pointerM2 += 1
            if pointerM2 < len(m2):
                currM2range = [m2[pointerM2], m2[pointerM2] + w2]
            continue
        if currM1range[0] < currM2range[0]:
            onlyM1 += currM2range[0] - currM1range[0]
            currM1range[0] = currM2range[0]
        elif currM2range[0] < currM1range[0]:
            onlyM2 += currM1range[0] - currM2range[0]
            currM2range[0] = currM1range[0]
        if currM1range[1] < currM2range[1]:
            shared += currM1range[1] - currM1range[0]
            currM2range[0] = currM1range[1]
            pointerM1 += 1
            if pointerM1 < len(m1):
                currM1range = [m1[pointerM1], m1[pointerM1] + w1]
        elif currM2range[1] < currM1range[1]:
            shared += currM2range[1] - currM2range[0]
            currM1range[0] = currM2range[1]
            pointerM2 += 1
            if pointerM2 < len(m2):
                currM2range = [m2[pointerM2], m2[pointerM2] + w2]
        else:
            shared += currM1range[1] - currM1range[0]
            pointerM1 += 1
            if pointerM1 < len(m1):
                currM1range = [m1[pointerM1], m1[pointerM1] + w1]
            pointerM2 += 1
            if pointerM2 < len(m2):
                currM2range = [m2[pointerM2], m2[pointerM2] + w2]
    p = shared / (shared + onlyM1)
    r = shared / (shared + onlyM2)
    if p == 0 or r == 0:
        f1 = 0
    else:
        f1 = 2 / (1 / p + 1 / r)
    return f1
