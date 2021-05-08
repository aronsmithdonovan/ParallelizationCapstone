# Pseudocode for Cell States + Probability Values


```
// CELL STATES:

0 = EMPTY      // empty ground containing no spore or hyphae
1 = SPORE      // contains at least one spore
2 = YOUNG      // young hyphae that cannot form mushrooms yet 
3 = MATURING   // maturing hyphae that cannot form mushrooms yet
4 = MUSHROOMS  // older hyphae with mushrooms
5 = OLDER      // older hyphae with no mushrooms
6 = DECAYING   // decaying hyphae with exhausted nutrients
7 = DEAD       // newly dead hyphae with exhausted nutrients
8 = DEADER     // hyphae that have been dead for a while
9 = DEPLETED   // area whose nutrients have previously been depleted by fungal growth
10 = INERT     // inert area where plants cannot grow
```

```
// PROBABILITY VALUES

probSpore = 0.001                  // probability that a site initially is SPORE
probSporeToYoung = 0.25            // probability that a SPORE will become YOUNG at the next time step
probSpread = 0.6                   // probability that a EMPTY with a neighbor that is YOUNG will become YOUNG at the next time step
probMaturingToMushroom = 0.7       // probability that a MATURING will become MUSHROOMS at the next time step (otherwise it becomes OLDER)
probDepletedToSpore = 0.0001       // probability that a DEPLETED will become SPORE at the next time step
probDepletedToEmpty = 0.5          // probability that a DEPLETED will become EMPTY at the next time step
```

```
// GROWTH MODEL:

    // initial
    for each cell {
        prob = random number between 0 and 1;

        if prob <= probSpore {
            cell is SPORE
        } else {
            cell is EMPTY
        }
    }

    // EMPTY --> YOUNG
    if cell is EMPTY {
        check cell's neighbors

        if cell has at least one neighbor that is YOUNG {
            prob = random number between 0 and 1;

            if prob <= probSpread {
                cell becomes YOUNG
            } else {
                cell stays EMPTY
            }
        }
    }

    // SPORE --> YOUNG
    if cell is SPORE {
        prob = random number between 0 and 1;

        if prob <= probSporeToYoung {
            cell becomes YOUNG
        } else {
            cell stays SPORE
        }
    }


    // YOUNG --> MATURING
    if cell is YOUNG {
        cell becomes MATURING
    }

    // MATURING --> MUSHROOMS | OLDER
    if cell is MATURING {
        prob = random number between 0 and 1;

        if prob <= probMaturingToMushrooms {
            cell becomes MUSHROOMS
        } else {
            cell becomes OLDER
        }
    }

    // MUSHROOMS | OLDER --> DECAYING
    if cell is MUSHROOMS or OLDER {
        cell becomes DECAYING
    }

    // DECAYING --> DEAD
    if cell is DECAYING {
        cell becomes DEAD
    }

    // DEAD --> DEADER
    if cell is DEAD {
        cell becomes DEADER
    }

    // DEADER --> DEPLETED
    if cell is DEADER {
        cell becomes DEPLETED
    }

    // DEPLETED --> SPORE | EMPTY
    if cell is DEPLETED {
        prob = random number between 0 and 1;

        if prob <= probDepletedToSpore {
            cell becomes SPORE
        } else if prob <= probDepletedToEmpty {
            cell becomes EMPTY
        } else {
            cell stays DEPLETED
        }
    }

    // INERT --> INERT
    if cell is INERT {
        cell stays INERT
    }
```