# cuaderno de bitácora

Registro en el que dejo por escrito ideas, planes futuros, líneas de investigación seguidas, por seguir y demás.

## slay

dos métodos distintos para optimizar órbitas

:::mermaid
graph TD
    A(Analítico) --> B(Secciones
    cónicas)  
    B --> C(Optimización) --> D
    D((Comparación))
    E(Numérico) --> F(RK-4) --> G(Optimización) --> D
:::

## mejor así

se usa el método analítico para enfocar con el numérico

:::mermaid
graph TD
    A(**Analítico**) --> B(*Secciones cónicas*
    se aproxima a elípticas
    en radios de influencia
    para poder obtener
    expresiones analíticas)  
    B --> C(*Optimización*
    es rápida, se analiza todo
    el espacio y se extraen
    ventanas óptimas) --> F
    C --> H(*Afinamiento analítico*
    se usa únicamente el
    modelo analítico para
    encontrar trayectorias
    óptimas) --> D
    E(**Numérico**) -->
    F(*RK-4 ~dt*
    más lento pero fiel a la
    realidad, se analizan las
    ventanas obtenidas) --> G(*Afinamiento numérico*
    se mejora la precisión,
    afinando parámetros y
    dejando margen de)
    G --> D(*Análisis de resultados*)
:::

[2026-01-10 13:24] Me acabo de dar cuenta de que el parámetro de impacto hay que recalcularlo porque ahora se está usando, erróneamente, la fase de Venus. Planea calcularlo en la función de simulación.

[2026-01-11 13:17] La parte analítica está acabada. Queda explicar el código de la simulación numérica y los resultados.
