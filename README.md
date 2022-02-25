# VLCM

Variable-Length CliqueMotif (VLCM) finds frequently occurring patterns, so-called motifs, in a time series. It is the first exact algorithm that discovers latent motifs of variable length.

### Compile
```
cd src/vlcm
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

### Run
```
./VLCM <time-series-path> <window-min> <window-max> <correlation> [v]
```

## GUI

### Install Requirements

```
cd src/ui
pip3 install -r requirements.txt
```
Alternatively, use a virtul environment:
```
python3 -m venv vlcm_venv/
source vlcm_venv/bin/activate
cd src/ui
python3 -m pip install -r requirements.txt
```

### Run
```
python3 app.py
```

## References

[1] Linardi, Michele, et al. "[Matrix Profile X: VALMOD - Scalable Discovery of Variable-Length Motifs in Data Series.](https://dl.acm.org/doi/pdf/10.1145/3183713.3183744)" *Proceedings of the 2018 International Conference on Management of Data*. 2018. Source: https://helios2.mi.parisdescartes.fr/~mlinardi/VALMOD.html

[2] Zimmerman, Zachary, et al. "[Matrix Profile XIV: Scaling Time Series Motif Discovery with GPUs to Break a Quintillion Pairwise Comparisons a Day and Beyond.](https://www.cs.ucr.edu/~eamonn/SCAMP-camera-ready-final1.pdf)" *Proceedings of the ACM Symposium on Cloud Computing*. 2019. Source: https://github.com/zpzim/SCAMP

[3] Jiang, Hua, Chu-Min Li, and Felip Manyà. "[Combining Efficient Preprocessing and Incremental MaxSAT Reasoning for MaxClique in Large Graphs.](http://www.mis.u-picardie.fr/~cli/ecai2016PublishedVersion.pdf)" *Proceedings of the twenty-second European conference on artificial intelligence*. 2016. Source: https://home.mis.u-picardie.fr/~cli/EnglishPage.html


