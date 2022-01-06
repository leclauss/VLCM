from multiprocessing import Process, Queue


class ParallelRun:
    def __init__(self, threads, outputFunction=None, log=lambda *_, **__: None):
        self.log = log
        self.outputFunction = outputFunction

        self.argsQueue = Queue()
        self.workers = [Process(target=self.__work, args=(i,)) for i in range(threads)]
        if outputFunction is not None:
            self.resultQueue = Queue()
            self.outputWorker = Process(target=self.__outputWorker)
            self.outputWorker.start()
        for worker in self.workers:
            worker.start()

    def __work(self, threadId):
        while True:
            work = self.argsQueue.get()
            if work is None:
                self.argsQueue.put(None)
                break
            function, args = work
            functionName = function.__name__
            self.log("Thread", threadId, "running", functionName, *args)
            result = function(*args)
            if self.outputFunction is not None:
                self.resultQueue.put((functionName, args, result))
        self.log("Thread", threadId, "done")

    def __outputWorker(self):
        while True:
            results = self.resultQueue.get()
            if results is None:
                break
            self.outputFunction(*results)

    def run(self, function, args):
        self.argsQueue.put((function, args))

    def join(self):
        self.argsQueue.put(None)
        for worker in self.workers:
            worker.join()
        if self.outputFunction is not None:
            self.resultQueue.put(None)
            self.outputWorker.join()
