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

Utiliza vi para editar el archivo .bashrc: `vi .bashrc`. Utiliza las flechas de arriba y abajo del teclado para ir al final del archivo. Ahora, presiona la letra <kbd>i</kbd> en el teclado para entrar al modo **insertar** que es donde se puede editar, en la parte de abajo de tu pantalla debe aparecer la leyenda "-- INSERT --".

Ahora, simplemente presiona <kbd>Enter</kbd> para ir a la última línea del archivo y escribe o copia el comando `module load r/4.4.1`. Para salvar y salir, hay que presionar la tecla <kbd>Esc</kbd> que nos saca del modo de inserción de texto (la leyenda de INSERT desaparecerá). Acto seguido, escribe <kbd>:</kbd> luego <kbd>w</kbd> luego <kbd>q</kbd> y luego <kbd>Enter</kbd>. Mientras hagas esto, irán apareciendo los caracteres :wq indicando que escribiremos el archivo (write) y saldremos del editor (quit). Al presionar <kbd>Enter</kbd>, saldremos del editor vi y volveremos a la terminal de linux.

Listo, ahora, entremos a la partición interactiva como nos indica Luis: 

``` bash
srun -p interactive --pty bash
```
El prompt de tu terminal debe cambiar de algo como `[username@login01 ~]` a `[username@node12 ~]` indicando que hemos pasado al nodo de cómputo 12 que debe tener recursos suficientes de RAM y procesamiento para utilizar R con datasets grandes.

Sin hacer más nada, el comando `module list` debería mostrar "r/4.4.1". Ahora, ve a la carpeta donde trabajarás tus datos o créala. Esta **NO debe estar en tu home donde pusimos .bashrc porque hay almacenamiento limitado allí**. Vamos a `cd /mnt/data/nombredetuPI/username/carpetadondetrabajaras`.

## Preparando el directorio.
Lo mejor es, una vez dentro del directorio del proyecto, mantener ordenadas las cosas pero eso queda a gusto de cada quien. Una manera es crear varios subdirectorios con el comando 
``` bash
mkdir mkdir {scripts,graficas,outputs,data}
```
Para mantener organizados los scripts, las graficas, los archivos de salida y los datos de entrada, respectivamente. Revisa este [tutorial](https://jerolon.github.io/Tutoriales/bash_training.html) si necesitas ayuda con bash.
El archivo del genoma y la anotación están disponibles para todos, para guardar el espacio en disco, haremos un softlink a éstos en lugar de que cada quien haga una copia para cada proyecto.

``` bash
#Genoma secuencia
ln -s /mnt/data/alfredvar/30-Genoma/Deroceras_laeve_genome_GCA_051403575.fasta genoma.fasta
#Anotacion en GFF
ln -s /mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff anotacion.gff
```

Esto forma links con color azul claro en tu directorio, de manera que puedes acceder a los archivos facilmente sin que estos ocupen espacio extra. Mira el resultado con `ls` y luego con `ls -lh`.

Ya casi estamos listos. El último paso es crear un script de R vacío donde iremos copiando y llenando los comandos del análisis.

``` bash
touch analisis_intervalos.R
```

Esto crea el archivo. Para mí, lo más fácil es navegar al directorio actual en el panel lateral de MobaXterm y hacer doble click en el archivo `analisis_intervalos.R` lo cual lo abre en una ventana nueva en el programa MobaTextEditor. Esto es conveniente porque todos los cambios a este archivo se guardan en el directorio remoto dentro del clúster.

Ahora sí, simplemente escribe el comando `R` en la terminal para entrar a la terminal de R.


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
