# Monocapa de grafeno

El grafeno es un alótropo de carbono, y consiste de una única capa de C $(sp^2)$ empacados en una estructura  bidimensional.

## Geometría de red

### Espacio Real

La monocapa de grafeno está conformada por varios carbonos arreglados en una estructura tipo panal. Se considera el sistema de coordenadas procurando que el arreglo en zig-zag esté alineado al eje x.

![AQUI](SLGreal.png)

Cada celda unitaria contiene 2 átomos de carbono pertenecientes a diferentes subredes, A y B. El arreglo de varias celdas unitarias forma una red hexagonal de Bravais $\{\textbf{R}\}$, con posiciones:

\begin{equation}
\textbf{R} = n_1 \textbf{a}_1 + n_2 \textbf{a}_2, \quad n_1, n_2 \epsilon \mathrm{\mathbb{Z}}	
\end{equation}

donde los vectores base $a_1$ y $a_2$ se definen como:

\begin{equation}
\textbf{a}_1 =a \left( \frac{1}{2},\frac{\sqrt3}{2}\right) 
\hspace{5mm},
\hspace{5mm}
\textbf{a}_2 =a \left(-\frac{1}{2},\frac{\sqrt3}{2}\right)
\end{equation}

y $a \simeq 2.46 \mathring A $ es el parámetro de red, ligado a la distancia carbono-carbono, $d$, por $a = \sqrt 3 d$.

### Espacio Recíproco

Dada la red en el espacio real, podemos definir el conjunto de puntos $\{\textbf{G}\}$ tales que $ e^{i \{\textbf{R}\} \cdot \{\textbf{G}\}} = 1$. Estos puntos forman la red recíproca y pueden escribirse en términos de sus vectores base, de la siguiente manera:

\begin{equation}
\textbf{G} = m_1 \textbf{b}_1 + m_2 \textbf{b}_2 \hspace{8mm} m_1, m_2  \hspace{1mm} \epsilon \hspace{1mm} \mathbb{Z}	
\end{equation}

donde los vectores $\textbf{b}_1$ y $\textbf{b}_2$ obedecen por definición:
\begin{equation}
\textbf{a}_i \cdot \textbf{b}_j = 2 \pi \delta_{i,j} 
\end{equation}

Por lo tanto, en este caso $\textbf{b}_1$ y $\textbf{b}_2$ se espresan como:

\begin{equation}
\textbf{b}_1 = \frac{4\pi}{3d} \left(\frac{\sqrt3}{2},\frac{1}{2}\right) 
\hspace{5mm},\hspace{5mm} 
\textbf{b}_2 =\frac{4\pi}{3d} \left(\frac{-\sqrt3}{2},\frac{1}{2}\right)
\end{equation}

![](SLGrec.png)

Ya que los estados de Bloch no se modifican por traslaciones en el espacio de momentos, toda la información del sistema puede obtenerse limitando el enfoque a la celda unitaria en el espacio recíproco (1ª zona de Brillouin). Los puntos de interés para este trabajo son aquellos de alta simetría, los cuales son:

\begin{equation}
\textbf{K}^{\pm}= \pm \frac{4\pi}{3 a} \left(1,0\right)
\hspace{5mm},\hspace{5mm} 
\textbf{M}= \frac{\pi}{a}  \left( 1,\frac{-1}{\sqrt{3}}       \right)
\hspace{5mm},\hspace{5mm} 
\boldsymbol{\Gamma} = \left(0,0 \right)
\end{equation}