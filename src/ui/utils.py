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


def cover(motif):
    return len(motif[0] * motif[1])


def getMaxCoverMotif(motifs):
    maxMotif = None
    maxCover = 0.0
    for motif in motifs:
        currentCover = cover(motif)
        if currentCover >= maxCover:
            maxCover = currentCover
            maxMotif = motif
    return maxMotif


def motifSimilarity(t1, t2):
    m1, w1 = t1
    m2, w2 = t2

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


def motifDistance(t1, t2):
    return 1 - motifSimilarity(t1, t2)


def maxMotifDistanceHigherThreshold(motif, motifList, threshold):
    for otherMotif in motifList:
        if motifDistance(motif, otherMotif) > threshold:
            return True
    return False
