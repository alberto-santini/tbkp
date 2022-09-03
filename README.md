# The 01-TB-KP

This repository contains the instances and the source code to solve the **0-1 Time-bomb Knapsack Problem**.
The problem is introduced in the following paper ([preprint](https://santini.in/files/papers/monaci-pike-burke-santini-2021.pdf), [journal version](https://www.sciencedirect.com/science/article/pii/S0305054822001253)):

```bib
@article{MPS2020,
    title={Exact algorithms for the 0-1 time-bomb knapsack problem},
    author={Monaci, Michele and Pike-Burke, Ciara and Santini, Alberto},
    journal={{Computers \& Operations Research}},
    volume=145,
    doi={10.1016/j.cor.2022.105848},
    year=2022
}
```

You can cite this repository via Zenodo:

```bib
@misc{tbkp_github,
    title={Code and instances for the 01-TB-KP},
    author={Santini, Alberto},
    date={2020-11-02},
    year={2020},
    howpublished={Github repository},
    doi={10.5281/zenodo.4193102},
    url={https://github.com/alberto-santini/tbkp/}
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
