class Model:
    def __init__(self):
        self.ts = None
        self.tsPath = None
        self.running = False
        self.showFiltered = True
        self.showClustered = False
        self.motifs = []
        self.prefilteredMotifs = []
        self.clusteredMotifs = []
        self.settingsCurrentRun = None
        self.clustering = None
        self.currentClusterThreshold = None
