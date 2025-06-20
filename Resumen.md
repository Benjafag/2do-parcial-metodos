# Unidad 3 - Sistemas de ecuaciones lineales

## 🔍 Objetivo

Resolver sistemas del tipo:

$$
A \cdot x = b
$$

donde $A \in \mathbb{R}^{n \times n}$, $x$ es el vector incógnita y $b$ el vector de constantes. Se usan métodos numéricos para evitar errores acumulativos y reducir el costo computacional.

## 🔸 Tipos de Matrices

**DENSAS:** Son aquellas que poseen pocos elementos nulos y son de orden bajo\
**RALAS (SPARSE):** Son aquellas que poseen muchos ceros y son de orden alto

## 🔸 Tipos de Métodos

### 1. **Métodos Directos**

Son aquellos que nos conducirían a la solución exacta luego de un número finito de operaciones elementales, si no hubiera errores de redondeo

### 🟢 **Eliminación de Gauss**

#### 🔍 Objetivo

Transforma el sistema en forma escalonada superior usando operaciones elementales a la matriz ampliada (modifica b) y luego se resuelve con **sustitución hacia atrás**.

**Errores comunes:**

- El error de redondeo se acumula en cada paso.
- Matrices mal condicionadas pueden amplificar errores.

### 💡 Ventajas

- Método **general**: se puede aplicar a cualquier sistema no singular.
- Base de otros métodos como LU, Thomas (caso tridiagonal), y muchos algoritmos numéricos.
- Se puede mejorar con:

  - **Pivoteo parcial** (mayor estabilidad).
  - **Escalamiento** (evita errores por magnitudes muy diferentes).

### ⚠️ Requisitos

- La matriz debe ser **no singular**.
- Sin pivoteo, puede fallar o generar grandes errores numéricos si hay ceros o valores muy pequeños en la diagonal.

### 💻 Complejidad

- Eliminación (reducción a triangular superior):$ O\left( \frac{2}{3} n^3 \right)$
- Sustitución hacia atrás:$ O(n^2) $
- **Total:** $O(n^3)$

---

### 🟢 **Factorización LU**

#### 🔍 Objetivo

Descompone $A = L \cdot U$, donde $A$ debe ser **no singular**, y para LU sin pivoteo, debe ser posible evitar ceros en la diagonal de $U$:

- $L$: matriz triangular inferior con unos en la diagonal
- $U$: matriz triangular superior
- Luego resuelve $L \cdot y = b$, y $U \cdot x = y$

#### 💡 Ventajas

- Útil para resolver muchos sistemas con la misma matriz $A$ y distintos vectores $b$.
- Se puede calcular $|A| = |L.U| = |U|$
- Cálculo de $A^{-1}$ usando como b las columnas canónicas
- Se puede combinar con **pivoteo parcial escalado** para mayor estabilidad

### 💻 Complejidad

- **Factorización LU:** $ O\left( \frac{2}{3} n^3 \right) $
- **Sustituciones (hacia adelante y hacia atrás):** $ O(n^2) $
- **Total por un sistema:** $ O(n^3) $
- **Total si ya tenés la factorización (solo resolver para otro $b$):** $ O(n^2) $

### **Pivoteo parcial**

- Consiste en, en cada paso de la eliminación de Gauss, buscar en la columna actual **entre las filas restantes** (desde la fila pivote hacia abajo) el elemento con mayor valor absoluto con el objetivo de evitar divisiones por números muy pequeños (o $0$) de tal forma que los multiplicadores $|m_{ik}| \leq 1$, evitando aumentar errores de redondeo.
- **Vector $P$ (permutaciones):**
  Registra el orden actual de las filas luego de los intercambios durante la eliminación. Inicialmente, $P = [1, 2, ..., n]^t$. Cada vez que se intercambian filas $i$ y $k$, se intercambian también $P_i \leftrightarrow P_k$. Útil para mantener la correspondencia entre filas originales y filas actuales, sin mover físicamente toda la matriz (optimización).

### **Escalamiento Implicito**

- Wilkinson propone que una matriz debe equilibrarse antes de aplicar una algoritmo de solucion. Se basa en **normalizar la comparación de valores en cada fila según la magnitud relativa en su fila** para elegir mejor el pivote. Por lo que para cada fila $i$, se calcula su factor de escala en $S$:

- **Vector $S$ (escalas):**
  Guarda para cada fila $i$ el valor máximo absoluto de sus coeficientes en la matriz original $A$:

$$
S_i = \max_{1 \leq j \leq n} |a_{ij}|
$$

- Al elegir el pivote, no se mira solo el valor absoluto del elemento de la columna $k$, sino el cociente: $ ck =\max{1 \leq j \leq n} \frac{|a\_{ij}^k| }{S_i} $ haciendo una comparación justa y relativa

### ⚙️ Funcionamiento resumido del ppe:

1. Se calcula $S$ antes de comenzar la eliminación y $P$ se inicializa desde 1 hasya n.
2. En cada paso $k$, para elegir fila pivote $p$, se evalúa: $ c*i = \frac{|a*{ik}|}{S_i}$

3. Se intercambian filas $k$ y $p$ tanto en la matriz como en $P$.
4. Se continúa con la eliminación habitual.
5. Al resolver queda: $PLU x = b \implies L y = P b \implies U x = y$. Es decir, debemos aplicar la permutacion a b.

#### 📝 Beneficios

- Permite hacer pivoteo con una comparación justa considerando la escala de cada fila.
- Mejora estabilidad y evita errores numéricos.
- El vector $P$ facilita reconstruir la solución o aplicar permutaciones posteriores sin perder datos.

---

### 🟢 Metodo de Gauss Jordan

Es similar al método de Gauss, la diferencia es que se diagonaliza la matriz

Al finalizar el algoritmo tenemos $𝑥 = b^n$

Desventaja: Costo aumenta en 50%

### 🟢 Descomposición de Cholesky

#### 🔍 Objetivo

Resolver $A x = b$ cuando $A$ es **simétrica y definida positiva**.

#### 📐 Descomposición:

$$
A = L \cdot L^T
$$

- $L$: matriz triangular inferior con coeficientes reales, incluso en su diagonal.
- Luego se resuelven:

$$
L y = b \quad \text{(sustitución hacia adelante)}
$$

$$
L^T x = y \quad \text{(sustitución hacia atrás)}
$$

### 🧮 Fórmulas por componentes

Para construir $L$, se recorren filas y columnas de forma incremental. Los elementos se calculan así:

#### 🟩 Diagonal ($i = j$):

$$
\ell_{ii} = \sqrt{ a_{ii} - \sum_{k=1}^{i-1} \ell_{ik}^2 }
$$

#### 🟦 Debajo de la diagonal ($i > j$):

$$
\ell_{ji} = \frac{1}{\ell_{ii}} \left( a_{ji} - \sum_{k=1}^{i-1} \ell_{jk} \ell_{ik} \right)
$$

### 💡 Ventajas

- Es **más eficiente** y **más estable** que la LU tradicional si se cumplen las condiciones.
- Requiere casi la **mitad del trabajo** que LU.
- No necesita pivoteo
- Ideal para matrices simétricas dispersas.

### 💻 Complejidad

- Descomposición: $O\left( \frac{1}{3} n^3 \right)$ (mitad que LU)
- Sustituciones: $O(n^2)$

---

### 🟢 Método de Thomas (tridiagonal)

#### 🔍 Objetivo

Resolver sistemas donde $A$ es **tridiagonal**, sólo tiene elementos distintos de cero en la diagonal principal y las diagonales adyacentes, es decir:

$$
\begin{bmatrix}f_1 & g_1 & 0   & \cdots & 0 \\e_2 & f_2 & g_2 & \cdots & 0 \\0   & e_3 & f_3 & \ddots & \vdots \\\vdots & \ddots & \ddots & \ddots & g_{n-1} \\0 & \cdots & 0 & e_n & f_n\end{bmatrix} .
\begin{bmatrix}x_1 \\x_2 \\x_3 \\\vdots \\x_{n-1} \\x_n\end{bmatrix} =
\begin{bmatrix}r_1 \\r_2 \\r_3 \\\vdots \\r_{n-1} \\r_n\end{bmatrix}
$$

### 🧠 Algoritmo

#### Descomposición

Para $k = 2$ hasta $n$: $e_k = \frac{e_k}{f_{k-1}}$; $f_k = f_k - e_k \cdot g_{k-1} $

#### Sustitución hacia adelante

Para $k = 2$ hasta $n$: $r*k = r_k - e_k \cdot r*{k-1} $

#### Sustitución hacia atrás

$x_n = \frac{r_n}{f_n}$ \
Para $k = n-1$ hasta $1, -1$ : $x_k = \frac{r_k - g_k \cdot x_{k+1}}{f_k} $

Es una forma optimizada de **eliminación de Gauss** que aprovecha la estructura tridiagonal para:

1. **Eliminar los elementos debajo de la diagonal** de forma eficiente
2. **Resolver con sustitución hacia atrás**

### 💡 Ventajas

- No requiere almacenar toda la matriz → solo 3 vectores
- Precisión superior a métodos genéricos al reducir operaciones

### 💻 Complejidad

- Fase de eliminación: $O(n)$
- Fase de sustitución: $O(n)$
- **Total:** $O(n)$

---

### 📌 Comparativa final de métodos directos

| Método       | Requisitos principales                      | Complejidad                                    | Ventajas clave                                                           |
| ------------ | ------------------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------ |
| **Gauss**    | No singular y $a_{ii} \neq 0$ (sin pivoteo) | $O(n^3)$                                       | General, robusto, base de muchos otros métodos                           |
| **LU**       | No singular y $a_{ii} \neq 0$ (sin pivoteo) | $O(n^3)$ (una vez) <br> $O(n^2)$ (por sistema) | Reutilizable para varios $b$, más eficiente que Gauss en ese caso        |
| **Cholesky** | **Simétrica definida positiva**             | $O\left(\frac{1}{3} n^3 \right)$               | Más rápido y estable que LU si aplica; usa menos operaciones             |
| **Thomas**   | **Tridiagonal**                             | $O(n)$                                         | Extremadamente eficiente; ideal para sistemas grandes con esa estructura |

---

### 2. Metodos iterativos

## 🔍 Objetivo

Resolver sistemas lineales $A x = b$ mediante **aproximaciones sucesivas**, comenzando con una estimación inicial $x^{(0)}$, que en principio convergen a la solucion x:

$$
x^{(k+1)} = B(x^{(k)}) + C
$$

B se llama matriz de iteración, es una generalización del método de punto fijo.

Se aplican principalmente cuando el sistema es de **gran tamaño** o **disperso (sparse)**.

### 💡 Ventajas respecto a los metodos directos

- Número de operaciones
- Posiciones de memoria
- Errores de redondeo

## 🔸 Tipos de Métodos Iterativos

- **Jacobi**
- **Gauss-Seidel**
- **SOR (Successive Over-Relaxation)**

---

## 🟢 Método de Jacobi

### 🔧 Ecuaciones base

Dado el sistema:

$$
A x = b \quad \Rightarrow \quad a_{ii} x_i^{(k+1)} = b_i - \sum_{j \neq i} a_{ij} x_j^{(k)}
$$

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)
$$

> Cada componente se calcula **usando solo valores de la iteración anterior** $\rightarrow$ Aproximaciones simultaneas.

### 💡 Ventajas

- Muy sencillo de implementar.
- Fácil de paralelizar (cada ecuación es independiente de las demás en cada iteración).

### ⚠️ Requisitos

- La matriz no debe tener ceros en la diagonal.

### 💻 Complejidad

- Cada iteración: $O(n^2)$
- Número de iteraciones depende de la convergencia (condición del sistema).

---

## 🟢 Método de Gauss-Seidel

### 🔧 Ecuación de iteración

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
$$

> Se usa **inmediatamente** cada nuevo valor calculado $x_j^{(k+1)}$ para actualizar el siguiente.

### 💡 Ventajas

- Mejora notablemente la velocidad de convergencia respecto a Jacobi.
- Menor consumo de memoria (no requiere vector auxiliar).
- Adecuado para sistemas **grandes y dispersos**.

### ⚠️ Requisitos

- No siempre converge si A tiene mala condición.

### 💻 Complejidad

- Cada iteración: $O(n^2)$
- Menos iteraciones que Jacobi (aunque más difíciles de paralelizar).

### 📉 Observaciones (Chapra)

- En la práctica, se usa como base para aceleración (como en SOR).
- Más susceptible a errores si hay malas estimaciones iniciales.

---

## 🟢 Método SOR (Relajación Sucesiva)

Este método usa un factor de ponderación para mejorar el valor calculado.

### 🔧 Ecuación iterativa

$$
x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
$$

- $\omega \in (0, 2)$: parámetro de relajación

  - $\omega = 1$ → Gauss-Seidel
  - $\omega < 1$ → subrelajación
  - $\omega > 1$ → **sobre-relajación** (acelera convergencia si está bien elegido)

### 💡 Ventajas

- Puede acelerar la convergencia de Gauss-Seidel de forma significativa.
- Muy útil en sistemas grandes y bien condicionados.

### ⚠️ Requisitos

- Requiere prueba empírica o análisis para encontrar el $\omega$ óptimo.

### 💻 Complejidad

- Cada iteración: $O(n^2)$
- Menos iteraciones que Gauss-Seidel si $\omega$ es adecuado

### 📉 Observaciones (Chapra)

- En muchos problemas prácticos, un $\omega \in [1.1, 1.5]$ acelera fuertemente la convergencia.
- Chapra sugiere experimentar con distintos valores y observar el número de iteraciones requeridas.

## 📌 Comparativa final de métodos iterativos

| Método           | Requisitos para convergencia                      | Convergencia | Paralelización | Velocidad relativa | Notas clave                          |
| ---------------- | ------------------------------------------------- | ------------ | -------------- | ------------------ | ------------------------------------ |
| **Jacobi**       | Diagonal dominante                                | Lenta        | Muy fácil      | 🟡 Lenta           | Simple, pero puede diverger          |
| **Gauss-Seidel** | Diagonal dominante, simétrica y definida positiva | Rápida       | Difícil        | 🟢 Media           | Usa valores actualizados al instante |
| **SOR**          | Idem + $\omega$ adecuado                          | Muy rápida   | Difícil        | 🟢🟢 Rápida        | Acelera G-S si $\omega$ bien elegido |

---

## Número de condición

$K(A) = ||A|| \cdot ||A^{-1}||$

- Es un medida cuantitativa del grado de mal condicionamiento de la matriz de coeficientes. Mide cuan cerca está una matriz de ser singular
- Se usa para calcular como afectan los errores relativos en A y/o b el cálculo de x.
- Si A y b tienen t cifras significativas y κ(A) es de un orden 10 s entonces la precisión del resultado será 10 s-t
- Se puede demostrar que:

  $$\frac {||\delta x||} {||x||} \le K(A) \frac {||\delta b||}{||b||} $$
  $$\frac {||\delta x||} {||x + \delta x||} \le K(A) \frac {||\delta A||}{||A||} $$

Wilkinson estudió el efecto del redondeo en el método de Eliminación de Gauss, considerando la triangulación con pivoteo y la solución de los dos sistemas triangulares, concluyendo que es un proceso muy estable, considerando que la matriz A no sea mal condicionada.

Una forma de chequear esto es controlando los elementos
de U, si crecen mucho es una señal de mala condición de
la misma

¿Cómo afecta el numero de condicion a los tipos de métodos?

- En **métodos directos**, el mal condicionamiento afecta la **exactitud** de la solución.
- En **métodos iterativos**, el mal condicionamiento afecta la **eficiencia y convergencia** del proceso.

## Comparación entre metodos directos e iterativos

|                                 | **Métodos Directos**           | **Métodos Iterativos**                                               |
| ------------------------------- | ------------------------------ | -------------------------------------------------------------------- |
| **Tiempo de ejecución**         | $\mathcal{O}(n^3)$             | $O(n^2 \times iteraciones)$                                          |
| **Almacenamiento**              | $n \times n$ (matriz completa) | $n$ (solo vectores y diagonales)                                     |
| **Errores de redondeo**         | Grandes                        | Despreciables (menos acumulativos)                                   |
| **Tiempo de ejecución (total)** | Finito (se conoce a priori)    | Indeterminado (depende de la convergencia)                           |
| **Tareas adicionales (TI)**     | "Barato"                       | "Caro" (más iteraciones o ajustes)                                   |
| **Aplicaciones típicas**        | Problemas generales            | Problemas específicos (matrices sparse, diagonales dominantes, etc.) |

---

# Unidad 4 Interpolación

## 🔍 Objetivo

Dado un conjunto de puntos conocidos $(x_i, f(x_i))$, encontrar una función que pase exactamente por ellos. La **interpolación** permite:

- Estimar valores intermedios de una función.
- Aproximar funciones complejas.
- Base para derivación, integración y resolución de ecuaciones.

## 🔸 Tipos de Interpolación

- **Interpolación Polinómica** (global)
- **Interpolación por tramos (Spline)**
- **Interpolación lineal y cuadrática simple**
- **Interpolación de Newton / Lagrange**

---

## 🟢 Interpolación Polinómica

### 🔧 Forma general

Dado $n+1$ puntos, existe un único polinomio de grado ≤ $n$ que los interpola:

$$
P_n(x) = a_0 + a_1 x + a_2 x^2 + \dots + a_n x^n
$$

El sistema se puede construir y resolver usando:

- **Forma de Vandermonde**
- **Forma de Lagrange**
- **Forma de Newton**

### ⚠️ Problemas (Chapra y Burden)

- Para muchos puntos ($n$ grande), el polinomio oscila fuertemente (**Fenómeno de Runge**).
- Poca estabilidad numérica si los puntos están muy cerca.
- Mejor usar interpolación por tramos o nodos Chebyshev.

---

## 🟢 Forma de Lagrange

### 🔧 Fórmula

$$
P_n(x) = \sum_{i=0}^{n} f(x_i) \cdot L_i(x)
$$

donde:

$$
L_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}
$$

### 💡 Ventajas

- No requiere resolver sistemas.
- Forma explícita del polinomio.

### ⚠️ Desventajas

- Requiere recalcular todo si se añade un nuevo punto.
- No se reutiliza cálculo.

---

## 🟢 Forma de Newton (Diferencias Divididas)

### 🔧 Forma general

$$
P_n(x) = a_0 + a_1(x - x_0) + a_2(x - x_0)(x - x_1) + \dots
$$

Los coeficientes $a_k$ se calculan mediante **diferencias divididas**:

$$
f[x_i, x_{i+1}] = \frac{f(x_{i+1}) - f(x_i)}{x_{i+1} - x_i}
$$

$$
f[x_i, x_{i+1}, x_{i+2}] = \frac{f[x_{i+1}, x_{i+2}] - f[x_i, x_{i+1}]}{x_{i+2} - x_i}
$$

### 💡 Ventajas

- Reutilizable si se agregan puntos.
- Útil para tabulación incremental.

---

## 🟢 Interpolación Lineal y Cuadrática

### 🔧 Interpolación lineal

Usa dos puntos para construir una recta:

$$
f(x) \approx f(x_0) + \frac{f(x_1) - f(x_0)}{x_1 - x_0}(x - x_0)
$$

### 🔧 Interpolación cuadrática

Usa tres puntos para un polinomio de grado 2:

$$
P_2(x) = a_0 + a_1(x - x_0) + a_2(x - x_0)(x - x_1)
$$

Usando diferencias divididas.

---

## 🟢 Splines (Interpolación por Tramos)

### 🔧 Objetivo

Construir polinomios de bajo grado (generalmente cúbicos) en cada intervalo $[x_i, x_{i+1}]$, garantizando **suavidad**:

- Continúa en primera y segunda derivada
- Evita oscilaciones del polinomio global

### 📐 Spline cúbico natural

- Condiciones de borde: $S''(x_0) = S''(x_n) = 0$
- Se resuelve un sistema tridiagonal → **Método de Thomas**

### 💡 Ventajas

- Alta precisión y suavidad
- Muy usado en gráficos, ingeniería y simulaciones

---

## 📉 Errores en la Interpolación

### 🔺 Error en interpolación polinómica (forma de Newton):

$$
f(x) - P_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^{n}(x - x_i)
$$

- Depende de derivada de orden $n+1$
- Crece con el número de nodos si no son bien distribuidos

### 🧠 Recomendación (Burden/Chapra)

- Usar **splines o polinomios de bajo grado por tramos** para alta precisión
- Evitar polinomios globales de grado alto

---

## 📌 Comparativa de Métodos de Interpolación

| Método            | Tipo       | Reutilizable | Precisión | Estabilidad | Observaciones                                 |
| ----------------- | ---------- | ------------ | --------- | ----------- | --------------------------------------------- |
| **Lagrange**      | Global     | ❌           | Alta      | Baja        | No se puede agregar puntos fácilmente         |
| **Newton DD**     | Global     | ✅           | Alta      | Media       | Buen desempeño incremental                    |
| **Lineal**        | Por tramos | ✅           | Media     | Alta        | Muy simple, baja continuidad                  |
| **Cuadrático**    | Por tramos | ✅           | Mejor     | Alta        | Aumenta suavidad y precisión                  |
| **Spline cúbico** | Por tramos | ✅           | Muy alta  | Muy alta    | Suave en derivadas, requiere resolver sistema |

---

# Unidad 5 – Integración Numérica

## 🔍 Objetivo

Aproximar integrales definidas de la forma:

$$
I = \int_a^b f(x)\,dx
$$

cuando:

- $f(x)$ no tiene antiderivada elemental
- Solo se conoce $f(x)$ en puntos discretos
- Se desea una solución aproximada con control del error

---

## 🔸 Tipos de Métodos

1. **Reglas de Newton-Cotes (puntos equiespaciados):**

   - Regla del rectángulo
   - Regla del trapecio
   - Regla de Simpson (1/3 y 3/8)

2. **Reglas compuestas:** aplican las anteriores por subintervalos

3. **Cuadratura Gaussiana:** nodos y pesos óptimos

---

## 🟢 Regla del Trapecio

### 🔧 Aproximación

Se interpola una recta entre $(a, f(a))$ y $(b, f(b))$:

$$
\int_a^b f(x)\,dx \approx \frac{b-a}{2} [f(a) + f(b)]
$$

### 💡 Ventajas

- Simple de implementar
- Útil para funciones lineales o suavemente curvadas

### ⚠️ Error

$$
E_T = -\frac{(b-a)^3}{12} f''(\xi)
$$

- Depende de la segunda derivada de $f$
- Se puede reducir dividiendo en subintervalos

## 💻 Regla del Trapecio Compuesta

$$
\int_a^b f(x)\,dx \approx \frac{h}{2} \left( f(x_0) + 2\sum_{i=1}^{n-1} f(x_i) + f(x_n) \right)
$$

con $h = \frac{b-a}{n}$

$$
E_T^{(comp)} = -\frac{(b-a)h^2}{12} f''(\xi)
$$

---

## 🟢 Regla de Simpson 1/3

### 🔧 Aproximación

Se interpola un polinomio cuadrático entre 3 puntos:

$$
\int_a^b f(x)\,dx \approx \frac{b-a}{6} [f(a) + 4f\left( \frac{a+b}{2} \right) + f(b)]
$$

## 💻 Regla de Simpson 1/3 Compuesta

Requiere **n par** subintervalos:

$$
\int_a^b f(x)\,dx \approx \frac{h}{3} \left[ f(x_0) + 4\sum_{\text{impares}} f(x_i) + 2\sum_{\text{pares}} f(x_i) + f(x_n) \right]
$$

### ⚠️ Error

$$
E_S = -\frac{(b-a)^5}{180} f^{(4)}(\xi)
$$

- Mucho más preciso que el trapecio si $f$ es suave
- Usa segunda y cuarta derivada

---

## 🟢 Cuadratura de Gauss

### 🔧 Objetivo

Aproximar:

$$
\int_{-1}^{1} f(x)\,dx \approx \sum_{i=1}^{n} w_i f(x_i)
$$

donde:

- $x_i$: raíces del polinomio de Legendre de grado $n$
- $w_i$: pesos asociados

> Puede integrarse con exactitud polinomios de grado $2n-1$ con solo $n$ puntos.

### 💡 Ventajas

- Mucho más precisa que Newton-Cotes con pocos puntos
- No requiere que los nodos estén equiespaciados
- Chapra y Burden recomiendan para integrales complicadas

### 💻 Cambio de intervalo

Para transformar $[a,b]$ a $[-1,1]$:

$$
\int_a^b f(x)\,dx = \frac{b-a}{2} \int_{-1}^{1} f\left( \frac{b-a}{2}x + \frac{a+b}{2} \right)\,dx
$$

### 📘 Tablas (para Gauss-Legendre)

Para $n = 2$:

$$
x_1 = -\frac{1}{\sqrt{3}}, \quad x_2 = \frac{1}{\sqrt{3}} \quad\text{y}\quad w_1 = w_2 = 1
$$

Para $n = 4$, los nodos y pesos son fracciones con raíces cuadradas (ver tabla en apunte anterior).

---

## 📉 Comparativa de Métodos

| Método             | Orden  | Nodos equiespaciados | Precisión relativa | Observaciones         |
| ------------------ | ------ | -------------------- | ------------------ | --------------------- |
| Trapecio simple    | 2      | ✅                   | Baja               | Sencillo              |
| Trapecio compuesto | 2      | ✅                   | Media              | Mejora con n          |
| Simpson 1/3        | 4      | ✅ (n par)           | Alta               | Muy usado             |
| Simpson 3/8        | 4      | ✅ (n múltiplo de 3) | Similar a 1/3      | Poco más complejo     |
| Gauss-Legendre     | $2n-1$ | ❌                   | Muy alta           | Nodos y pesos óptimos |

---

## 📌 Observaciones finales (Chapra/Burden)

- A mayor grado del polinomio, mayor el riesgo de oscilaciones → usar con cuidado
- Métodos compuestos (Simpson/Trapecio) son preferibles a globales
- Cuadratura de Gauss es ideal cuando se busca **alta precisión con pocos puntos**
- Siempre tener en cuenta el comportamiento de las derivadas al estimar el error
