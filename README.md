# Final Project: Mushroom Fairy Rings
### Aron Smith-Donovan, Macalester College
***
<br>

### Description:
This repository contains two C++ files for modeling the growth of a mycelium network in a given area; one file generates the model with sequential processing (`fungi-seq.cpp`), and the other generates the model with parallel processing (`fungi-omp.cpp`) using the OpenMP API. The graphic model outputs directly to the terminal feed.

"Fairy rings" are a naturally-occuring ring or arc of mushrooms connected by underground mycelia. The term 'fungi' refers generally to multicellular, spore-producing organisms, and mushrooms are the fruiting body of certain types of fungi. Mushroom-producing fungi, however, have composed of much more than the mushrooms themselves: thin, branching tubules called hyphae grow underground in search of nutrients, and are capable of branching out and connecting with other hyphae. Collectively, a network of hyphae is called the mycelium (pl. mycelia). After developing from a spore, the mycelium develops and grows radially outward, sometimes sprouting mushrooms before it depletes the soil of nutrients and continues further outward, creating the ring shape. Mycelium networks can connect with others to create even larger compund ring structures.

In this model, the fungal growth described has been simulated with a structured grid approach on a cellular automata; a 2D grid is created and initialized with some cells being empty and other having spores. Each cell can be one of 11 states:
 * EMPTY = empty ground containing no spore or hyphae
 * SPORE = contains at least one spore
 * YOUNG = young hyphae that cannot form mushrooms yet 
 * MATURING = maturing hyphae that cannot form mushrooms yet
 * MUSHROOM = older hyphae with mushrooms
 * OLDER = older hyphae with no mushrooms
 * DECAYING = decaying hyphae with exhausted nutrients
 * DEAD = newly dead hyphae with exhausted nutrients
 * DEADER = hyphae that have been dead for a while
 * DEPLETED = area whose nutrients have previously been depleted by fungal growth
 * INERT = inert area where plants cannot grow

Spores develop into young hyphae capable of spreading to empty area around them before aging, potentially sprouting mushrooms, and dying, leaving the cell depleted. Depleted cells have a low chance of becoming empty again, and an even lower chance of becoming a spore, to acknowledge the lack of nutrients in that site (nutrient use itself has not been included in the model). The inert state has been implemented, but there is not yet a way to have a cell become inert.

This project was inspired by a project description from *Introduction to Computational Science: Modeling and Simulation for the Sciences* (Shiflet and ShifletPrinceton University Press 2014) and completed as the course project for COMP445: Parallel and Distributed Processing for the spring 2021 semester by Macalester College undergraduate student Aron Smith-Donovan under the guidance of Prof. Libby Shoop. The code in conjunction with the written report meet the requirements for a capstone project for undergraduates pursuing a Bachelor's in Computer Science.

Special thanks to Prof. Shoop, the C documentation, and the Stack Overflow forums.
<br>

***
<br>

### File structure:

```
  project-aronsmithdonovan\
      .gitignore
      README.md
      Makefile
      fungi-seq.cpp
      fungi-omp.cpp
      seq_time.h
      report\
         report.pdf
         fungi-state-diagram.png
         pseudocode.md
```
<br>

***
<br>

<blockquote>

### Running the simulation:

   <blockquote>

   **Option 1: sequential processing**<br>
   * set flags in Makefile
      * with DEBUG flag enabled, disable COLOR flag in Makefile for numerical output
      * with DEBUG flag enabled, enable COLOR flag in Makefile for color-coded output
      * disable DEBUG and COLOR flag for just the runtime as output
   * navigate to the main directory in the terminal
   * execute `$ make seq.fungi`
   * execute `$ ./seq.fungi -r R -c C -s S` where `R` is the number of rows, `C` is the number of columns, and `S` is the number of time steps
   * `R`, `C`, and `S` must all be positive nonzero integers (an error will be thrown at runtime if the arguments supplied do not meet this criteria)

   </blockquote>
   <br>
   <blockquote>

   **Option 2: parallel processing**<br>
   * set flags in Makefile
      * with DEBUG flag enabled, disable COLOR flag in Makefile for numerical output
      * with DEBUG flag enabled, enable COLOR flag in Makefile for color-coded output
      * disable DEBUG and COLOR flag for just the runtime as output
   * navigate to the main directory in the terminal
   * execute `$ make omp.fungi`
   * execute `$ ./omp.fungi -r R -c C -s S -t T` where `R` is the number of rows, `C` is the number of columns, `S` is the number of time steps, and `T` is the number of threads
   * `R`, `C`, `S`, and `T` must all be positive nonzero integers (an error will be thrown at runtime if the arguments supplied do not meet this criteria)

   </blockquote>
   <br>
</blockquote>
<!-- end of file -->