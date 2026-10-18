# Modelo UMZ-SIRS: análisis numérico de un sistema acoplado de desinformación y epidemias

**Universidad de Ingeniería y Tecnología (UTEC)**  
**Curso:** CC2104 – Métodos Numéricos  
**Semestre:** 2026-2  
**Grupo:** 8

## Descripción del proyecto

Este repositorio contiene el desarrollo de un proyecto de investigación numérica basado en el modelo **UMZ-SIRS**, que combina la propagación de información y desinformación con la dinámica de una enfermedad infecciosa.

El sistema está compuesto por seis variables de estado:

- **U:** población no informada.
- **M:** población que difunde desinformación.
- **Z:** población escéptica.
- **S:** población susceptible.
- **I:** población infectada.
- **R:** población recuperada.

El objetivo es estudiar el comportamiento del modelo mediante diferentes métodos numéricos, analizar sus propiedades de estabilidad y evaluar cómo los parámetros asociados a la información modifican la dinámica epidemiológica.

## Métodos numéricos implementados

Se desarrollaron y compararon los siguientes métodos:

1. Euler explícito.
2. Runge-Kutta de cuarto orden (RK4).
3. Runge-Kutta-Fehlberg (RKF45) con paso adaptativo.
4. Euler implícito o Euler hacia atrás.
5. Newton-Raphson para resolver los sistemas no lineales del método implícito.
6. Gauss-Seidel para resolver los sistemas lineales internos.

También se utiliza `ode15s` de MATLAB como solución de referencia en los experimentos de rigidez.

## Estructura del repositorio

```text
CC2104-UMZ-SIRS-Grupo8/
├── matlab/
│   ├── MN_Grupo8_01_Modelo_UMZ_SIRS.m
│   ├── MN_Grupo8_02_Euler_RK4.m
│   ├── MN_Grupo8_03_RKF45.m
│   ├── MN_Grupo8_04_Euler_Implicito.m
│   ├── MN_Grupo8_05_Comparacion_Metodos.m
│   ├── MN_Grupo8_06_Rigidez_Estabilidad.m
│   ├── MN_Grupo8_07_Escenarios_Epidemiologicos.m
│   └── MN_Grupo8_08_Analisis_Parametrico.m
├── figuras/
├── resultados/
├── paper/
│   ├── main.tex
│   ├── 01_introduccion.tex
│   ├── 02_metodologia.tex
│   ├── 03_resultados.tex
│   ├── 04_discusion.tex
│   ├── 05_referencias.tex
│   └── imagenes/
├── entregables/
│   ├── avances/
│   └── entrega1/
└── README.md
```

Los programas de MATLAB también disponen de versiones Live Script (`.mlx`) para facilitar su ejecución y revisión.

## Experimentos desarrollados

| Programa | Descripción |
|---|---|
| 01 | Formulación e implementación del modelo UMZ-SIRS |
| 02 | Integración mediante Euler explícito y RK4 |
| 03 | Implementación del método RKF45 adaptativo |
| 04 | Euler implícito con Newton-Raphson y Gauss-Seidel |
| 05 | Comparación de métodos numéricos |
| 06 | Análisis de rigidez y estabilidad para distintos valores de ε |
| 07 | Escenarios epidemiológicos asociados a desinformación y escepticismo |
| 08 | Análisis de sensibilidad de los parámetros del modelo |

## Resultados y visualizaciones

La carpeta `resultados/` contiene los archivos CSV generados en los experimentos, mientras que `figuras/` almacena las visualizaciones utilizadas en el informe.

Se analizan principalmente:

- La precisión y convergencia de los métodos numéricos.
- El comportamiento del sistema ante diferentes escalas temporales.
- La estabilidad de los métodos explícitos e implícitos.
- Los efectos de la desinformación y el escepticismo sobre los contagios.
- La sensibilidad del modelo ante cambios en sus parámetros.

## Ejecución

**Requisitos:** MATLAB con las funciones numéricas utilizadas en los programas.

1. Clonar el repositorio:

   ```bash
   git clone https://github.com/jmtimana/CC2104-UMZ-SIRS-Grupo8.git
   ```

2. Abrir MATLAB y establecer la carpeta `matlab/` como directorio de trabajo.

3. Abrir el programa `.m` o `.mlx` correspondiente al experimento.

4. Ejecutar el programa y revisar los resultados numéricos y las gráficas.

5. Consultar los archivos CSV en `resultados/` y las figuras exportadas en `figuras/`.

**Nota:** algunos programas guardan sus resultados mediante rutas relativas al directorio de ejecución. Si se ejecutan desde `matlab/`, es posible que los archivos generados aparezcan dentro de esa carpeta y deban trasladarse a `resultados/` o `figuras/`.

## Documentación y entregables

- **`paper/`**: código fuente editable en LaTeX del artículo de investigación.
- **`entregables/avances/`**: informes de los avances del proyecto.
- **`entregables/entrega1/`**: primera entrega oficial del curso.
- **`figuras/`**: imágenes generadas mediante MATLAB.
- **`resultados/`**: datos numéricos exportados en formato CSV.

## Estado del proyecto

El repositorio contiene los experimentos numéricos desarrollados hasta el tercer avance, junto con el artículo en LaTeX y los entregables correspondientes.

El trabajo continuará con la revisión, validación y preparación de la entrega final.

## Equipo

**Grupo 8 – Métodos Numéricos, UTEC, 2026-2.**