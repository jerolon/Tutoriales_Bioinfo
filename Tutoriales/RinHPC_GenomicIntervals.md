# Uso de R en la terminal del HPC del Laboratorio Nacional de Investigación Avanzada

## Configuración de R
Mi consejo es utilizar mobaxterm, pues trae la configuración X11 de cajón. Además, su editor de texto es bastante bueno para editar scripts remotos.

Lo primero que requerimos es cargar el módulo de R:

``` bash
module avail
```
Te mostrará todos los módulos disponibles en el clúster, incluyendo múltiples versiones de R. Utilizemos la versión más actual instalada en este momento 4.4.1

``` bash
module load r/4.4.1
```

Si usas R a menudo o lo usarás mucho aunque sea por un periodo corto, esto da mucha pereza. Por ello, añadiremos el comando anterior a nuestro archivo .bashrc. Este es un script oculto que se ejecuta en cada login para configurar ciertos valores del bash.

El archivo se encuentra en la "casa" de cada quien `/home/usuario`. Para ir allí, si no estás ya, escribe `cd ~` o simplemente `cd`. Una vez en nuestra casa, muestra los archivos oculos con `ls -a` (-a de all).
Aparecerán algunos archivos que comienzan con un punto como `.Xauthority, .bash_history, .bashrc`.

Cargamos varias librerias del tidyverse para hacernos mas facil la vida
con los dataframes

``` r
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(tximport)
library(stats)
library("RColorBrewer")
library(ggplot2)
library(ggfortify)
```
