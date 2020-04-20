---
header-includes:
  - \usepackage[ruled,vlined,linesnumbered]{algorithm2e}
---

# The Mackey-Glass Anomaly Benchmark

## Anomaly Insertion Process
\begin{algorithm}[tbp]
\caption{Pseudo-code of the anomaly insertion procedure for Mackey-Glass time series.
}       
\label{alg:create-anomaly}
   \begin{algorithmic}[1]
		 \State \textbf{Define parameters}: 
		 \State \indent $m \leftarrow100$: Minimum distance between both split points $t_i$ and $t_j$
		 \State \indent $K \leftarrow 100$: Length of the comparison window 
		 \State
		 \Function{insertAnomaly}{$x(t)$, $m$, $K$} \Comment Time series $x(t)$ of length $T$
		    %\State $\vec x(t) = \Big(x(t), \frac{d}{dt}x(t),  \frac{d^2}{dt^2}x(t),  \frac{d^3}{dt^3}x(t) \Big)$\Comment Stack derivatives in $\vec x(t) \in \mathbb{R}^4$
		    \State $\vec x: \mathbb{N} \rightarrow \mathbb{R}^4, t \mapsto \Big(x(t); \frac{d}{dt}x(t);  \frac{d^2}{dt^2}x(t);  \frac{d^3}{dt^3}x(t) \Big)$ \Comment{Stack $x(t)$ and its derivatives}
		    \State $t_i \leftarrow$ Random value from range $[0,T]$ \Comment First split point
		    %\State $t_i' \leftarrow t_i+m$ \Comment There should be a min. distance between both split points %, where $t_i$ is randomly drawn from range $[0,T]$
		    \State $d(k) \leftarrow ||\vec x(t_i) - \vec x(t_i+m+k)||$ \Comment euclidean distance between $\vec x(t_i)$ \& $\vec x(t_i'+k)$
		    \State $\mathbb{K} \leftarrow \{1,\ldots, K\}$
		    \State $k_{min} \leftarrow \argmin\limits_{k\in \mathbb{K}} d(k)$ \Comment gives us 2\textsuperscript{nd} split point $t_j=t_i+m+k_{min}$
		    \State $x_m(t) \leftarrow \begin{cases} 
                                x(t), & t < t_i \\
                                x(t+m+k_{min}), & t > t_i \\
                            \end{cases}$ \Comment Manipulated time series with anomaly
		 %\For  {$i \in \{1\ldots k\}$} 
		 %	\State $x^2+x$
		  %\Call{predict}{dd}
		 %\EndFor
		 \State \Return $x_m(t)$
		 \EndFunction 
	\end{algorithmic}
\end{algorithm}
