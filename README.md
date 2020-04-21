![Anomalous Mackey-Glass Time Series][header]


# The Mackey-Glass Anomaly Benchmark

This repository contains the Mackey-Glass anomaly benchmark (MGAB), which is composed of synthetic Mackey-Glass time series with non-trivial anomalies. Mackey-Glass time series are known to exhibit chaotic behavior under certain conditions. MGAB contains 10 MG time series of length 10<sup>5</sup>. Into each time series 10 anomalies are inserted with a procedure as described below.
In contrast to other synthetic benchmarks,  it is very hard for the human eye to distinguish the introduced anomalies from the normal (chaotic) behavior. 
An excerpt of a time series caontaining 3 anomalies is shown in the graph above. The location of the anomalies are revealed in the last plot of this page.

## The Benchmark Files
The labeled data for the time series 1-10 can be found in the CSV files `[1-10].csv`. Each file contains a table with 4 columns: 
1. `time`: Time in seconds in the range 0 to 10<sup>5</sup> -1.
2. `value`: Value of x(t) in the range [0.26, 1.66].
3. `is_anomaly`: Binary values (0/1) indicate which data points are considered as normal (0) or anomalous (1). For each anomaly, a range of 400 points, the so called anomaly window, is flagged. Detections within the anomaly window are to be considered as correct.
4. `is_ignored`: Since some algorithms might require a "warm-up" phase when processing the time series, the is_ignored column indicates for the initial 256 time steps that false detections can be ignored. It was ensured that no anomalies were placed at the beginning of any time series.

## Time Series Generation
We use the following Mackey-Glass equation (a non-linear time delay differential equation, DDE) to generate our time series:

![\frac{dx}{dt} = \beta \cdot \frac{x(t-\tau)}{1+x(t-\tau)^n} - \gamma x(t)](https://render.githubusercontent.com/render/math?math=%5Cfrac%7Bdx%7D%7Bdt%7D%20%3D%20%5Cbeta%20%5Ccdot%20%5Cfrac%7Bx(t-%5Ctau)%7D%7B1%2Bx(t-%5Ctau)%5En%7D%20-%20%5Cgamma%20x(t))

The parameters are real numbers which we set to τ=18, n=10, β=0.25, γ=0.1. Additionally, a constant history parameter is required which is set to h=0.9. A sufficiently long time series is generated using the [JiTCDDE](https://github.com/neurophysik/jitcdde) solver ( with an integration stepsize of one) which is then divided into 10 new time series. The time delay embedding of such a time series is illustrated in Fig. 1.

![Pseudo-code of the anomaly insertion procedure for Mackey-Glass time series.][timedelay]<br>**Figure 1**: Time  delay  embedding  of  the  Mackey-Glass  attractor.

## Anomaly Insertion Process
The main idea of the anomaly insertion process is to randomly remove segments from each time series in a way that this will be hardly visible later. To do so, we try to find two points (with a minimal and maximal distance) in a random segment of the time series so that the values of these two points as well as their derivatives closely match. Then, we remove the segment between those two points and "stich" the remaining parts together again. The exact procedure is as follows:

1. For the time series sequence ![x(t)](https://render.githubusercontent.com/render/math?math=x(t)) estimate the first 3 derivatives ![dx/dt](https://render.githubusercontent.com/render/math?math=dx%2Fdt), ![d^2x/dx^2](https://render.githubusercontent.com/render/math?math=d%5E2x%2Fdx%5E2) and ![d^3x/dx^3](https://render.githubusercontent.com/render/math?math=d%5E3x%2Fdx%5E3) by numerical differentiation of ![x(t)](https://render.githubusercontent.com/render/math?math=x(t)). Then, stack the original time series ![x(t)](https://render.githubusercontent.com/render/math?math=x(t)) and the three derivatives in a four-dimensional time series ![\mathbf x(t)](https://render.githubusercontent.com/render/math?math=%5Cmathbf%20x(t)).
2. Randomly select a position ![t_i](https://render.githubusercontent.com/render/math?math=t_i) in ![\mathbf x(t)](https://render.githubusercontent.com/render/math?math=%5Cmathbf%20x(t)). This will be the first split point
3. Starting at ![t_i'=t_i+m](https://render.githubusercontent.com/render/math?math=t_i'%3Dt_i%2Bm), with ![m=100](https://render.githubusercontent.com/render/math?math=m%3D100), for all ![k\in K$, $K=\{0,1,\ldots,100\}](https://render.githubusercontent.com/render/math?math=k%5Cin%20K%24%2C%20%24K%3D%5C%7B0%2C1%2C%5Cldots%2C100%5C%7D), compare ![\mathbf x(t)](https://render.githubusercontent.com/render/math?math=%5Cmathbf%20x(t)) to ![\mathbf x(t_i'+k)](https://render.githubusercontent.com/render/math?math=%5Cmathbf%20x(t_i'%2Bk)) and compute the euclidean norm ![d(k)=||\mathbf x(t_i) - \mathbf x(t_i'+k)||](https://render.githubusercontent.com/render/math?math=d(k)%3D%7C%7C%5Cmathbf%20x(t_i)%20-%20%5Cmathbf%20x(t_i'%2Bk)%7C%7C).
4. The index ![k](https://render.githubusercontent.com/render/math?math=k) which minimizes the distance ![d(k)](https://render.githubusercontent.com/render/math?math=d(k)) will give us the second split point ![t_j=t_i'+k_m](https://render.githubusercontent.com/render/math?math=t_j%3Dt_i'%2Bk_m), where ![k_m=\argmin_{k\in K} d(k)](https://render.githubusercontent.com/render/math?math=k_m%3D%5Cargmin_%7Bk%5Cin%20K%7D%20d(k))
5. Construct a new manipulated time series ![x_m(t)](https://render.githubusercontent.com/render/math?math=x_m(t)), which is ![x(t)](https://render.githubusercontent.com/render/math?math=x(t)), for ![t\le t_i](https://render.githubusercontent.com/render/math?math=t%5Cle%20t_i) and ![x(t+m+k_m)](https://render.githubusercontent.com/render/math?math=x(t%2Bm%2Bk_m)) for ![t > t_i](https://render.githubusercontent.com/render/math?math=t%20%3E%20t_i).

The procedure is also summarized in Algorithm 1 below.
![Pseudo-code of the anomaly insertion procedure for Mackey-Glass time series.][algorithm1]

An example of how such an anomaly, which was generated using the described procedure, could look like is illustrated in Fig. 2. For the human eye it would be almost impossible to spot the anomaly.

![Example Anomaly in a MG time series 1][mgexample1] ![Example Anomaly in a MG time series 1][mgexample2]<br>
**Figure 2**: Top: Example for the creation of a Mackey-Glass time series with a temporal anomaly. The original time series (dashed line) is manipulated in such a way, that a segment is removed and the two remaining ends are joined together. In this example, the interval [21562,21703] is removed from the original curve. The resulting manipulated time series (solid line) has a smooth point of connection, but significantly differs from the original. Bottom: Zoomed-In. The red shaded area indicates the position where the anomaly was inserted.

## Final Notes
For the 10 time series of this benchmark, in total 100 anomalies were inserted (10 anomalies per time series). In the last step, in order to increase the complexity of the anomaly detection task slightly, we add noise drawn from a random uniform distribution with the range [-0.01, 0.01] to each point of all time series. 



![Anomalous Mackey-Glass Time Series with revealed Anomalies][footer]
**Figure 3**: This graph shows the same section of a Mackey-Glass time series as the first graph on this page, but now reveals the location of the anomalies in the time series. The anomalies are at t<sub>1</sub> = 40388, t<sub>2</sub>=40917 and t<sub>3</sub>=41550. The positions are indicated by the black crosses in the plot.

[header]: https://github.com/MarkusThill/MGAB/blob/master/img/mg_example_anomalies_hidden.png?raw=true "Anomalous Mackey Glass Time Series"
[footer]: https://github.com/MarkusThill/MGAB/blob/master/img/mg_example_anomalies_revealed.png?raw=true "Anomalous Mackey Glass Time Series with revealed Anomalies"
[algorithm1]: https://github.com/MarkusThill/MGAB/blob/master/img/pseudocode.png?raw=true "Pseudo-code of the anomaly insertion procedure for Mackey-Glass time series."
[timedelay]: https://github.com/MarkusThill/MGAB/blob/master/img/time-delay-embedding.png?raw=true "Time  delay  embedding  of  the  Mackey-Glass  attractor."
[mgexample1]: https://github.com/MarkusThill/MGAB/blob/master/img/mg-anomaly-example.png?raw=true "Example Anomaly in a MG time series 1"
[mgexample2]: https://github.com/MarkusThill/MGAB/blob/master/img/mg-anomaly-example-zoom.png?raw=true "Example Anomaly in a MG time series 2"
