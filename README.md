# Proyecto final de Genómica Funcional
## Reproducción de resultados

### Autoras:
* Mayra Lizeth Bárcenas Guevara (**GitHub:** https://github.com/mayra-bioinfo; **Página personal:** )
* Mariana Gómez Becerra (**GitHub:** https://github.com/Marianagomez-uaq; **Página personal:** https://marianagomez-uaq.github.io/mi-sitio-quarto/)

### Fecha de entrega:
Miércoles 20 de mayo de 2026

## Pregunta biológica

Las comunidades sintéticas (SynComs) son parte de un campo muy reciente de la biología, que integra ecología, análisis microbiano y metagenómica (Shayanthan et al., 2022). Las SynComs se construyen mediante co-cultivo de microorganismos de distintos grupos taxonómicos, en condiciones que simulan un microbioma; busca reducir la cantidad de microorganismos de la comunidad, sin perder los efectos positivos que tiene la interacción microbioma-hospedante. Estas comunidades han sido diseñadas y aplicadas para mejorar el crecimiento de plantas en invernadero y en campo, o para tratar desórdenes gastrointestinales en humanos (Van Leeuwen et al., 2023).

Las interacciones de los microorganismos en comunidades sintéticas generan propiedades emergentes que no son predecibles basándose en análisis de cada cepa aislada, por lo que estas interacciones forman redes complejas (Gastélum et al., 2025).

Dada la complejidad de estas redes, es importante generar metodologías que permitan estudiarlas y aprovechar de mejor manera el potencial de las interacciones. Esto es lo que se propuso Delgado‐Baquerizo (2022) en su artículo: Simplifying the complexity of the soil microbiome to guide the development of next-generation SynComs.

Delgado‐Baquerizo desarrolló una metodología basada en datos crudos de secuenciación de ARN, que tras la asignación taxonómica, realiza un filtrado basado en distintos parámetros, y finalmente clusteriza las comunidades mediante análisis con redes de correlación; con esta metodología, se identifican los microorganismos con mayor relevancia ecológica, por su función y presencia nativa en distintos tipos de suelo. El autor plantea que aunque esta metodología se aplicó en muestras tomadas de suelo, puede extrapolarse a otros tipos de microbioma, para facilitar el diseño de SynComs en distintos campos de estudio.

En este proyecto, nos planteamos reproducir los resultados que obtuvo el autor, con los siguientes objetivos:
* Verificar la reproducibilidad de la metodología propuesta
* Identificar diferencias entre los resultados originales y los obtenidos en este proyecto
* Realizar análisis adicionales en la red generada, para comprender mejor las propiedades topológicas de la misma

## Enfoque metodológico

Reproduciremos la metodología con programas distintos a los utilizados por el autor, reemplazandolos por los vistos durante el semestre. Se espera que los resultados sean similares, demostrando que el éxito de la metodología no depende del uso de un programa específico, y puede ser reproducido por investigadores alrededor del mundo con sus programas de preferencia.

Los análisis serán realizados utilizando paquetes de R (R Core Team, 2025):
* DADA2: Análisis de datos de secuenciación y asignación taxonómica.
* igraph: Generación y análisis de la red
* networkD3: Visualización de la red filtrada (solo si la máquina lo permite)

## Datos a utilizar



## Bibliografía

* R Core Team (2025). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
* Delgado‐Baquerizo, M. (2022). Simplifying the complexity of the soil microbiome to guide the development of next‐generation SynComs. Journal of Sustainable Agriculture and Environment, 1(1), 9–15. https://doi.org/10.1002/sae2.12012
* Gastélum, G., Gómez-Gil, B., Olmedo-Álvarez, G., & Rocha, J. (2025). Harnessing emergent properties of microbial consortia for Agriculture: Assembly of the Xilonen SynCom. Biofilm, 9, 100284. https://doi.org/10.1016/j.bioflm.2025.100284
* Shayanthan, A., Ordoñez, P. a. C., & Oresnik, I. J. (2022). The Role of Synthetic Microbial Communities (SynCom) in Sustainable Agriculture. Frontiers in Agronomy, 4. https://doi.org/10.3389/fagro.2022.896307
* Van Leeuwen, P. T., Brul, S., Zhang, J., & Wortel, M. T. (2023). Synthetic microbial communities (SynComs) of the human gut: design, assembly, and applications. FEMS Microbiology Reviews, 47(2). https://doi.org/10.1093/femsre/fuad012
