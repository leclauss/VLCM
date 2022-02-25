import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from PyQt5.QtCore import QObject, pyqtSignal, QThread
from scipy import stats

from utils import loadTS, getMotifPoints, getRanges, motifDistance, getMaxCoverMotif, maxMotifDistanceHigherThreshold, \
    cover
from model import Model
from view import View


class Presenter:
    def __init__(self, view: View, model: Model):
        self.view = view
        self.model = model

        self.runThread = None
        self.worker = None

    def openFile(self):
        fileName = self.view.openFileDialog()
        self.loadData(fileName)

    def reset(self):
        self.model.motifs = []
        self.model.prefilteredMotifs = []
        self.model.showClustered = False
        self.view.setDendrogramButtonEnabled(False)
        self.model.clusteredMotifs = []
        self.model.clustering = None
        self.model.currentClusterThreshold = None

    def clearMotifs(self):
        self.view.showRanges([])
        self.view.clearMotifPlot()
        self.view.clearMotifTable()

    def loadData(self, filePath):
        if filePath is not None:
            if self.model.running:
                self.view.showWarningMessage("Already Running",
                                             "The program is already running. Please cancel the current run before loading new data.")
            else:
                try:
                    ts = loadTS(filePath)
                except Exception as e:
                    self.view.showWarningMessage("Error", "There was a problem loading the file:\n" + str(e))
                    return

                self.reset()
                self.model.ts = ts
                self.model.tsPath = filePath

                self.view.plotTs(self.model.ts, Path(filePath).name)
                self.view.setSliderRanges(10, len(self.model.ts) // 2)
                self.clearMotifs()

                self.view.setProgress(0, "ready")
                self.view.setError(False)
                self.view.setRunnable(True)

    def setRunning(self, running):
        self.view.setRunning(running)
        self.model.running = running

    def run(self):
        if not self.model.running:
            # clear
            self.clearMotifs()
            self.reset()
            # load settings
            self.model.settingsCurrentRun = self.view.getSettings()
            minLength, maxLength, correlation, prefilter = self.model.settingsCurrentRun
            # run
            self.setRunning(True)
            self.view.setProgress(0, "started")
            self.view.setError(False)
            # create thread
            self.runThread = QThread()
            self.worker = Worker(self.model.tsPath, minLength, maxLength, correlation, prefilter)
            self.worker.moveToThread(self.runThread)
            # setup
            self.runThread.started.connect(self.worker.run)
            self.worker.finished.connect(self.runThread.quit)
            self.worker.finished.connect(self.worker.deleteLater)
            self.runThread.finished.connect(self.runThread.deleteLater)
            # connect signals and slots
            self.worker.error.connect(self.error)
            self.worker.progress[int].connect(lambda value: self.view.setProgress(value=value))
            self.worker.progress[str].connect(lambda text: self.view.setProgress(text=text))
            self.worker.progress[int, str].connect(lambda value, text: self.view.setProgress(value, text))
            self.worker.motif.connect(lambda motif, window: self.addMotif(motif, window))
            self.worker.prefiltered.connect(lambda motif, window: self.addPrefiltered(motif, window))
            self.worker.clustering.connect(lambda clustering: self.processClustering(clustering))
            self.runThread.finished.connect(self.finishRun)
            # start thread
            self.runThread.start()
        else:
            # cancel
            self.worker.kill()

    def error(self):
        self.view.setError(True)
        self.view.setProgress(text="cancelled")

    def addMotif(self, motif, window):
        self.model.motifs.append((motif, window))
        if not self.model.showFiltered:
            self.view.addMotifRow(window, len(motif))

    def addPrefiltered(self, motif, window):
        self.model.prefilteredMotifs.append((motif, window))
        if self.model.showFiltered:
            self.view.addMotifRow(window, len(motif))

    def processClustering(self, clustering):
        self.model.clustering = clustering
        self.model.showClustered = True
        clusterThreshold = 0.25  # default cluster threshold
        self.model.currentClusterThreshold = clusterThreshold
        self.updateClusteredMotifs(clusterThreshold)
        self.view.setDendrogramButtonEnabled(True)

    def updateClusteredMotifs(self, cutoff):
        self.model.clusteredMotifs = self.getMotifsFromClustering(cutoff)
        if self.model.showFiltered:
            self.updateMotifTable()

    def getMotifsFromClustering(self, cutoff):
        # reconstruct clusters
        prefilteredMotifs = self.model.prefilteredMotifs.copy()
        clustering = self.model.clustering

        prefilteredNumber = len(prefilteredMotifs)
        clusterMap = {}
        for i in range(clustering.shape[0]):
            if clustering[i, 2] > cutoff:
                # next clusters to merge have too high distance
                break
            # merge two clusters
            clusterId = prefilteredNumber + i

            # get representatives of both clusters
            index1 = int(clustering[i, 0])
            if index1 >= prefilteredNumber:
                index1 = clusterMap[index1]
            index2 = int(clustering[i, 1])
            if index2 >= prefilteredNumber:
                index2 = clusterMap[index2]
            motif1 = prefilteredMotifs[index1]
            motif2 = prefilteredMotifs[index2]

            # only keep motif with highest cover
            if cover(motif1) > cover(motif2):
                mergeMotif = index1
                prefilteredMotifs[index2] = None
            else:
                mergeMotif = index2
                prefilteredMotifs[index1] = None
            clusterMap[clusterId] = mergeMotif

        resultMotifs = []
        for motif in prefilteredMotifs:
            if motif is not None:
                resultMotifs.append(motif)
        return resultMotifs

    def finishRun(self):
        self.setRunning(False)

    def getCurrentMotifs(self):
        if self.model.showFiltered:
            if self.model.showClustered:
                return self.model.clusteredMotifs
            else:
                return self.model.prefilteredMotifs
        else:
            return self.model.motifs

    def showMotif(self):
        motifId = self.view.getSelectedMotif()
        if motifId is not None:
            motif, window = self.getCurrentMotifs()[motifId]

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
        motifs = self.getCurrentMotifs()
        if len(motifs) == 0:
            self.view.showWarningMessage("No Motifs", "There are no motifs to save.")
            return
        fileName = self.view.writeFileDialog()
        if fileName is not None:
            minLength, maxLength, correlation, prefilter = self.model.settingsCurrentRun
            with open(fileName, "w") as file:
                file.write("time series;" + str(self.model.tsPath) +
                           "\ntime series length;" + str(len(self.model.ts)) +
                           "\nlength min;" + str(minLength) +
                           "\nlength max;" + str(maxLength) +
                           "\ncorrelation;" + str(correlation) +
                           "\n\nlength;size;occurrences\n")
                for entry in motifs:
                    motif, window = entry
                    file.write(str(window) + ";" + str(len(motif)) + ";" + str(motif) + "\n")

    def switchMotifTable(self):
        self.model.showFiltered = not self.model.showFiltered
        self.updateMotifTable()

    def updateMotifTable(self):
        # clear
        self.clearMotifs()

        if self.model.showFiltered:
            self.view.setClusterButtonText("Show All")
        else:
            self.view.setClusterButtonText("Show Filtered")

        motifs = self.getCurrentMotifs()
        for entry in motifs:
            motif, window = entry
            self.view.addMotifRow(window, len(motif))

    def showDendrogram(self):
        labels = []
        for motif in self.model.prefilteredMotifs:
            labels.append(motif[1])
        threshold = self.view.showDendrogramWindow(self.model.clustering, labels, self.model.currentClusterThreshold)
        if threshold is not None:
            self.model.currentClusterThreshold = threshold
            self.updateClusteredMotifs(threshold)


class Worker(QObject):
    finished = pyqtSignal()
    error = pyqtSignal()
    progress = pyqtSignal([int], [str], [int, str])
    motif = pyqtSignal(list, int)
    prefiltered = pyqtSignal(list, int)
    clustering = pyqtSignal(object)

    def __init__(self, tsPath, windowMin, windowMax, correlation, prefilter):
        super().__init__()
        self.tsPath = tsPath
        self.windowMin = windowMin
        self.windowMax = windowMax
        self.correlation = correlation
        self.prefilter = prefilter
        self.process = None

    def run(self):
        try:
            # get motifs with VLCM
            args = ["../vlcm/build/VLCM", self.tsPath, str(self.windowMin), str(self.windowMax), str(self.correlation)]
            self.process = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)

            totalMotifs = self.windowMax - self.windowMin + 1
            calculatedMotifs = 0

            currentMotifs = []
            prefilteredMotifs = []

            # process motif as soon as it is output
            for line in self.process.stdout:
                splits = line.split(";")
                w = int(splits[0])
                m = sorted([int(i) for i in splits[1].split(",")])
                calculatedMotifs += 1
                self.progress[int].emit(int(math.ceil(95.0 * calculatedMotifs / totalMotifs)))

                if len(m) > 1:
                    self.motif.emit(m, w)

                    # pre-clustering
                    motif = [m, w]
                    if len(currentMotifs) > 0 and maxMotifDistanceHigherThreshold(motif, currentMotifs, self.prefilter):
                        maxMotif = getMaxCoverMotif(currentMotifs)
                        self.prefiltered.emit(*maxMotif)
                        prefilteredMotifs.append(maxMotif)
                        currentMotifs = [motif]
                    else:
                        currentMotifs.append(motif)

            if len(currentMotifs) > 0:
                maxMotif = getMaxCoverMotif(currentMotifs)
                self.prefiltered.emit(*maxMotif)
                prefilteredMotifs.append(maxMotif)

            return_code = self.process.wait()
            if return_code:
                # error / killed
                self.error.emit()
                return

            if len(prefilteredMotifs) >= 2:
                # agglomerative clustering
                motifDistanceMatrix = pdist(np.array(prefilteredMotifs, dtype=object), metric=motifDistance)
                clustering = linkage(motifDistanceMatrix, 'average')
                self.clustering.emit(clustering)

            self.progress[int, str].emit(100, "finished")
        finally:
            self.finished.emit()

    def kill(self):
        if self.process is not None:
            try:
                self.process.kill()
            except ProcessLookupError:
                pass
