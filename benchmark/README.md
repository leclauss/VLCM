# Benchmark

This benchmark compares VLCM against other state-of-the-art algorithms.

### Prerequisites

Create the folder _algorithms_ in the project root and download the repositories of the other algorithms:

```
mkdir algorithms
cd algorithms
git clone https://github.com/flash121123/HIME
git clone https://github.com/GrammarViz2/grammarviz2_src
```
You also need to download the VALMOD source code from [here](https://helios2.mi.parisdescartes.fr/~mlinardi/VALMOD.html) and unzip it into the same folder.
Then, follow the respective instructions to compile the algorithms if necessary.

Download the benchmark from [here](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/UCR_TimeSeriesAnomalyDatasets2021.zip). Extract the folder _UCR_Anomaly_FullData_ and put it into the _data_ directory.

### Run

Use the following command to run the full benchmark. The script checks whether the prerequisites were done correctly.

```
python3 benchmark.py <threads> 1 250 <output-path>
```
