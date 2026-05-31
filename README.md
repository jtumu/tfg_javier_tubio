# tfg_javier_tubio
En este repositorio se encuentra todo el material utilizado para el TFG titulado "Aplicación de modelos de lenguaje de proteínas para la expansión funcional de genomas".

El archivo titulado "Análisis_EColi_csv.R" es un _script_ de R preparado para analizar los archivos .csv que FANTASIA Lite genera tras recibir los proteomas de E. Coli. El _script_ "Análisis_SCerevisiae_csv.R" es el mismo código que el de E. Coli pero las figuras, ontologías y tablas están adaptadas para S. Cerevisiae (el análisis se hizo igual para los dos organismos).

El archivo llamado "aleatorizador_copia.py" es el _script_ de Python utilizado para aleatorizar el orden de los aminoácidos de las proteínas otorgadas.

En la carpeta "Samples.zip" están los archivos .csv que FANTASIA Lite devuelve y que se analizaron en dichos _scripts_: los archivos .tsv son los datos originales disponibles en UniProt, los archivos "results_control..." son los archivos del proteoma original pasado por FANTASIA Lite y los archivos "..._random..." son los archivos de cada aleatorización pasada por FANTASIA Lite (están nombrados de la siguiente forma: organismo _ random _ semilla utilizada para aleatorizar).
