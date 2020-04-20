# The Mackey-Glass Anomaly Benchmark

Image here

This repository contains the Mackey-Glass anomaly benchmark (MGAB), which is composed of synthetic Mackey-Glass time series with non-trivial anomalies. Mackey-Glass time series are known to exhibit chaotic behavior under certain conditions. MGAB contains 10 MG time series of length 10<sup>5</sup>. Into each time series 10 anomalies are inserted with a procedure as described below.
In contrast to other synthetic benchmarks,  it is very hard for the human eye to distinguish the introduced anomalies from the normal (chaotic) behavior. 
An excerpt of a time series caontaining 3 anomalies is shown in the graph below. The location of the anomalies are revealed in the last graph of this page.

## The Benchmark Files
The labeled data for the time series 1-10 can be found in the CSV files [1-10].csv. Each file contains 4 columns: 
1. `time`: Time in seconds in the range 0 to 10<sup>5</sup> -1.
2. `value`: Value in the range [0.26, 0.66].
3. `is_anomaly`: Binary values (0/1) indicate which data points are considered as normal (0) or anomalous (1). For each anomaly, a range of 400 points, the so called anomaly window, is flagged. Detections within the anomaly window are to be considered as correct.
4. `is_ignored`: Since some algorithms might require a "warm-up" phase when processing the time series, the is_ignored column indicates for the initial 256 time steps that false detections can be ignored. It was ensured that no anomalies were placed at the beginning of any time series.


## Anomaly Insertion Process
