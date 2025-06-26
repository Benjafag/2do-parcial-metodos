## 📋 Preguntas de repaso – Teórico integrador

### Conceptuales generales

1. ¿Por qué el error del método del trapecio involucra la segunda derivada de la función? ¿Y por qué el de Simpson involucra la cuarta?
2. ¿Qué significa que una fórmula de integración sea interpolatoria? ¿Qué métodos vistos se basan en esta idea?
3. ¿Cuál es la ventaja principal de la cuadratura de Gauss frente a los métodos de Newton-Cotes?
4. ¿Por qué Newton-Cotes no puede alcanzar una precisión de grado \(2n - 1\) como sí lo hace Gauss-Legendre?
5. Si una función es muy oscilante en un intervalo, ¿qué métodos de integración conviene usar y cuáles podrían fallar? Justificá.
6. ¿Qué significa que el error sea de orden \(O(h^4)\)? ¿Eso garantiza que el método será más preciso que uno con error \(O(h^2)\)?
7. ¿Qué condiciones deben cumplirse para que una fórmula de cuadratura tenga grado de precisión \(2n - 1\)?
8. ¿Por qué las raíces de los polinomios ortogonales son ideales como nodos en cuadratura de Gauss?
9. ¿Se puede usar cuadratura de Gauss si solo tenés datos tabulados de la función? ¿Por qué sí o por qué no?
10. ¿Cómo mejora la extrapolación de Richardson una estimación de integral? ¿Qué rol cumple el orden del error en ese proceso?
11. ¿Qué mejora aporta el método de Romberg frente al método del trapecio compuesto aplicado varias veces?
12. ¿Por qué la regla de Simpson 1/3 puede integrar exactamente polinomios de grado 3, si interpola con un polinomio de grado 2?
13. ¿Qué implica que el término cúbico se anule al integrar el polinomio interpolante?
14. Si tenés tres valores: \(f(0), f(0.5), f(1)\), ¿qué regla de integración podés usar? ¿Qué orden de error tiene?

---

### Verdadero o falso – Justificá

1. Dada \(f(x) = 0\), la iteración generada por el método de la secante siempre converge a la raíz.
2. Dados los puntos \((x_0,y_0), (x_1,y_1), ..., (x_n,y_n)\), existe un único polinomio de grado \(n\) que pasa por esos puntos.
3. En el pivoteo parcial solo se intercambian filas si el pivote es menor que algún otro elemento.
4. Para determinar si un sistema lineal está mal condicionado alcanza con analizar el determinante de la matriz.
5. Todo número con representación decimal finita tiene representación binaria finita.
6. El método de extrapolación de Richardson se puede usar para mejorar resultados de cualquier fórmula de Newton-Cotes.
7. El método de Newton-Raphson siempre tiene orden de convergencia cuadrático.
8. La diferencia de dos números aproximadamente iguales debe evitarse en algoritmos numéricos.
9. Cuantos más nodos interpolatorios se consideren, mejor aproximación se logra con el polinomio interpolante.
10. Un método inestable siempre corresponde a un problema mal condicionado.
11. El principal error en las fórmulas de cuadratura es el de truncamiento.
12. El error en Newton-Cotes puede minimizarse haciendo pequeños subintervalos con integración compuesta.
13. El error del polinomio de Newton es menor que el de Lagrange.
14. La principal ventaja de los métodos de Gauss es que permiten usar menos puntos para igual o mejor resultado.
15. Newton-Cotes es más robusto que cuadratura de Gauss.
16. Si se resuelve un sistema con eliminación de Gauss y se halla solución exacta, entonces el sistema era bien condicionado.
17. Si se tiene una función del tipo \(f(x) = t(t-1)(t-2)...(t-n)\), su integral exacta puede obtenerse con una fórmula de cuadratura.
18. Las fórmulas de cuadratura de Gauss son más eficientes que las de Newton-Cotes.

---

### De parciales reales (fragmentos de imágenes)

1. ¿Cuál es el tipo de error predominante en las fórmulas de Newton-Cotes? ¿Cómo se puede reducir?
2. Caracterizá a los métodos de un paso para resolver EDO con valor inicial.
3. ¿Qué relación hay entre error local y error global en métodos para EDO? ¿Cómo se relacionan gráficamente?
4. Hablar sobre las fórmulas de cuadratura vistas, comparándolas brevemente.
5. Verdadero o falso (justificá):
   - A. Si se tiene un número grande de puntos conocidos, siempre conviene usar splines.
   - B. Dados los puntos \((x_0,y_0), ..., (x_n,y_n)\), existe un único polinomio de grado \(n\) que pasa por esos puntos.
   - C. Cuantos más nodos interpolatorios se consideren, mejor aproximación se logra.
   - D. La principal ventaja de Gauss frente a Newton-Cotes es que permite integrar funciones con discontinuidades.
   - E. El error predominante en cuadratura es el truncamiento.
   - F. Para usar Runge-Kutta de orden alto, hay que usar otro método para calcular los primeros pasos.

