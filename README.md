# The 01-TB-KP

This repository contains the instances and the source code to solve the **0-1 Time-bomb Knapsack Problem**.
The problem is introduced in the following paper:

```bib
@online{MPS2020,
    title={{The 0--1 Time-bomb Knapsack Problem}},
    author={Monaci, Michele and Pike-Burke Ciara, and Santini, Alberto},
    year={2020},
    version={1}
}
```

## Instances

Instances are in folder `data/generated-instances`.
Folder `data` also contains original instances for the deterministic 0-1 Knapsack Problem and an instance generator (in Ruby).

## Solver

Source code for the implementation of algorithms solving the 01-TB-KP is in folder `src`.
You can build the source code using CMake and the provided `CMakeLists.txt` file.

## License

The code is released under the GNU Public License version 3 (see file `LICENSE`).