import ast
import re
import resource
import shutil
import subprocess
import sys
from pathlib import Path

from parallelRun import ParallelRun

SCRIPT_PATH = Path(__file__).parent.absolute()
PROJECT_PATH = SCRIPT_PATH.parent
EXE_PATHS = {
    "VLCM": PROJECT_PATH / "src/vlcm/build/VLCM",
    "HIME": PROJECT_PATH / "algorithms/HIME/bin/HIME_release.jar",
    "GV": PROJECT_PATH / "algorithms/grammarviz2_src/target/grammarviz2-1.0.1-SNAPSHOT-jar-with-dependencies.jar",
    "VALMOD": PROJECT_PATH / "algorithms/VALMOD_src/VALMOD"
}
UCR_PATH = PROJECT_PATH / "data/UCR_Anomaly_FullData"

RUN_TIMEOUT = 3600  # timeout per run in seconds
MAX_MEMORY = 8589934592  # max memory for each run in bytes


def checkPrerequisites():
    print("Checking prerequisites ...")
    failed = False
    for alg, path in EXE_PATHS.items():
        present = path.is_file()
        print(alg, "executable in", str(path), "...", present)
        if not present:
            failed = True
    present = UCR_PATH.is_dir() and len(list(UCR_PATH.glob("*"))) == 250
    print("UCR_Anomaly dataset in", str(UCR_PATH), "...", present)
    if not present:
        failed = True
    if failed:
        print("Error: missing executables or datasets\nPlease read the instructions in the README file.")
        return False
    else:
        return True


def run(args, timeout=RUN_TIMEOUT, maxMemory=MAX_MEMORY, output=True, cwd=None):
    # limit memory
    def setLimits():
        resource.setrlimit(resource.RLIMIT_AS, (maxMemory, maxMemory))

    preFunc = None
    if maxMemory is not None:
        if args[0] == "java":
            args.insert(1, "-Xmx" + str(maxMemory))
        else:
            preFunc = setLimits

    # limit time
    if timeout is not None:
        timeoutArgs = ["timeout", "--kill-after=15", str(timeout)]
    else:
        timeoutArgs = []

    # run command
    procInfoId = "PROC_INFO"
    allArgs = ["/usr/bin/time", "-f" + procInfoId + ":%e,%M,%P"] + timeoutArgs + [str(arg) for arg in args]
    outFile = subprocess.PIPE if output else subprocess.DEVNULL
    proc = subprocess.Popen(allArgs, preexec_fn=preFunc, stdout=outFile, stderr=subprocess.PIPE,
                            universal_newlines=True, cwd=cwd)
    stdout, stderr = proc.communicate()

    # parse time output to get elapsed time, max memory consumption and cpu share
    info = [None, None, None]
    for line in stderr.splitlines():
        if line.startswith(procInfoId):
            try:
                infoSplits = line.split(":")[1].split(",")
                time = float(infoSplits[0])
                memory = int(infoSplits[1])
                cpu = int(infoSplits[2].split("%")[0])
                info = [time, memory, cpu]
            except:
                pass
    return proc.returncode, info, stdout


def VLCM(tsPath, minLength, maxLength, correlation, naive=False):
    args = [EXE_PATHS["VLCM"], tsPath, minLength, maxLength, correlation]
    if naive:
        args.append("n")
    returnCode, info, output = run(args)
    motifs = []
    for line in output.splitlines():
        try:
            splits = line.split(";")
            window = int(splits[0])
            motif = sorted([int(i) for i in splits[1].split(",")])
            motifs.append((window, motif))
        except:
            pass
    return returnCode == 0, info, motifs


def VLCMNaive(tsPath, minLength, maxLength, correlation):
    return VLCM(tsPath, minLength, maxLength, correlation, naive=True)


def VALMOD_Base(tsPath, minLength, maxLength, tsLength):
    runId = str(Path(tsPath).stem) + "_" + str(minLength) + "_" + str(maxLength) \
            + "_" + str(abs(hash(tsPath)) % (10 ** 8))
    runName = "tmp_" + runId
    runPath = SCRIPT_PATH / runName
    try:
        runPath.mkdir()
        args = [EXE_PATHS["VALMOD"], tsLength, minLength, maxLength, 0, 50, runName, tsPath]
        returnCode, info, _ = run(args, output=False, cwd=runPath)
        motifs = []
        patternFirst = re.compile("with length (\d+):.*Q:(\d+) , D:(\d+)")
        pattern = re.compile("Size (\d+) is over.*q:(\d+), d:(\d+)")
        try:
            with open(runPath / "Logrun.log") as logFile:
                for line in logFile:
                    values = patternFirst.findall(line) + pattern.findall(line)
                    try:
                        window = int(values[0][0])
                        motif = sorted([int(values[0][1]), int(values[0][2])])
                        motifs.append((window, motif))
                    except:
                        pass
        except:
            pass
        return returnCode == 1, info, motifs
    finally:
        shutil.rmtree(runPath)


def VALMOD(tsPath, tsLength, minLength, maxLength):
    return VALMOD_Base(tsPath, minLength, maxLength, tsLength)


def VALMOD_LogSteps(tsPath, tsLength, minLength, maxLength):
    mergedResults = (True, (0, 0, 0), [])
    currentMinLength = minLength
    currentMaxLength = minLength * 2
    while currentMinLength <= maxLength:
        runMaxLength = min(currentMaxLength, maxLength)
        result = VALMOD_Base(tsPath, currentMinLength, runMaxLength, tsLength)
        mergedResults = mergeResults([mergedResults, result])
        if not result[0] or len(result[2]) == 0:
            break
        if mergedResults[1][0] > RUN_TIMEOUT:
            mergedResults = (False, mergedResults[1], mergedResults[2])
            break
        currentMinLength = currentMaxLength + 1
        currentMaxLength = currentMinLength * 2
    return mergedResults


def HIME(tsPath, minLength):
    args = ["java", "-Xmx8g", "-jar", EXE_PATHS["HIME"], tsPath, minLength]
    returnCode, info, output = run(args)
    motifs = []
    for line in output.splitlines():
        if line.startswith("Motif:"):
            try:
                splits = line.split(" ")
                window = int(splits[5])
                motif = sorted([int(splits[1]), int(splits[3])])
                motifs.append((window, motif))
            except:
                pass
    return returnCode == 0, info, motifs


def extended(function, tsPath, minLength, maxLength):
    mergedResults = (True, (0, 0, 0), [])
    maxMotif = minLength - 1
    while maxMotif < maxLength:
        result = function(tsPath, maxMotif + 1)
        mergedResults = mergeResults([mergedResults, result])
        if not result[0] or len(result[2]) == 0:
            break
        if mergedResults[1][0] > RUN_TIMEOUT:
            mergedResults = (False, mergedResults[1], mergedResults[2])
            break
        maxMotif = max(result[2], key=lambda item: item[0])[0]
    return mergedResults


def HIME_Extended(tsPath, minLength, maxLength):
    return extended(HIME, tsPath, minLength, maxLength)


def GrammarViz(tsPath, minLength):
    args = ["java", "-cp", EXE_PATHS["GV"], "net.seninp.grammarviz.cli.TS2SequiturGrammar",
            "-o", "/dev/stdout", "--strategy", "EXACT", "-d", tsPath, "-w", minLength]
    returnCode, info, output = run(args)
    motifs = []
    currentLines = []
    for line in output.splitlines():
        if line.startswith("mean length"):
            try:
                window = int(line.split("mean length ")[1])
                motif = sorted(ast.literal_eval(currentLines[2].split("subsequences starts: ")[1]))
                motifs.append((window, motif))
            except:
                pass
            currentLines = []
        else:
            currentLines.append(line)
    return returnCode == 0, info, motifs


def GrammarViz_Extended(tsPath, minLength, maxLength):
    return extended(GrammarViz, tsPath, minLength, maxLength)


def mergeResults(results):
    success = True
    timeSum = 0
    maxMemory = 0
    cpuSum = 0
    allMotifs = []
    try:
        for result in results:
            success = success and result[0]
            timeSum += result[1][0]
            maxMemory = max(maxMemory, result[1][1])
            cpuSum += result[1][0] * result[1][2]
            allMotifs += result[2]
        return success, [timeSum, maxMemory, cpuSum / timeSum], allMotifs
    except:
        return False, [None, None, None], []


def getTsLength(tsPath):
    length = 0
    with open(str(tsPath)) as tsFile:
        for line in tsFile:
            try:
                float(line)
                length += 1
            except:
                pass
    return length


def main(argv):
    if not len(argv) == 5:
        print("Usage: python3 benchmark.py <THREADS> <START-ID> <END-ID> <OUTPUT-FILE>")
        return 1
    if not checkPrerequisites():
        return 1

    threads = int(argv[1])
    startId = int(argv[2])
    endId = int(argv[3])
    outFilePath = argv[4]

    def output(functionName, args, result):
        motifs = result[2] if result[0] else None
        with open(outFilePath, "a") as outFile:
            outFile.write(";".join(
                [functionName, ",".join([str(arg) for arg in args]), str(result[0]),
                 ";".join([str(res) for res in result[1]]), str(motifs)]) + "\n")

    pool = ParallelRun(threads, log=print, outputFunction=output)

    for filePath in UCR_PATH.glob("*"):
        if filePath.is_file():
            try:
                tsId = int(str(filePath.name).split("_")[0])
                if startId <= tsId <= endId:
                    ts = str(filePath)
                    tsLength = getTsLength(ts)
                    minLength = 32
                    maxLength = tsLength // 2

                    pool.run(HIME, (ts, minLength))
                    pool.run(HIME_Extended, (ts, minLength, maxLength))

                    pool.run(GrammarViz, (ts, minLength))
                    pool.run(GrammarViz_Extended, (ts, minLength, maxLength))

                    pool.run(VLCMNaive, (ts, minLength, maxLength, 0.99))

                    for correlation in [0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7]:
                        pool.run(VLCM, (ts, minLength, maxLength, correlation))

                    pool.run(VALMOD, (ts, tsLength, minLength, maxLength))
                    pool.run(VALMOD_LogSteps, (ts, tsLength, minLength, maxLength))
            except:
                print("Unknown file in UCR benchmark:", filePath)

    pool.join()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
