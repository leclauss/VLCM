# VLCM

Compile:
cd src/vlcm
mkdir build
cd build
cmake ..
cmake --build . --config Release

Run:
./VLCM <time-series-path> <window-min> <window-max> <correlation> [v]


## GUI

install requirements

Run:
cd src/ui
python3 app.py
