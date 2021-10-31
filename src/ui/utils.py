def loadTS(tsPath):
    with open(tsPath) as tsFile:
        ts = [float(line) for line in tsFile]
    return ts


def getMotifPoints(motif, window):
    motifPoints = set()
    for index in motif:
        motifPoints = motifPoints.union(set(range(index, index + window)))
    return motifPoints


def getRanges(points):
    sortedList = sorted(points)
    if len(sortedList) == 0:
        return []
    startPoint = sortedList[0]
    lastPoint = sortedList[0]
    ranges = []
    for point in sortedList[1:]:
        if point > lastPoint + 1:
            ranges.append((startPoint, lastPoint))
            startPoint = point
        lastPoint = point
    ranges.append((startPoint, lastPoint))
    return ranges
