# CM5-Charge Manual

## Table of Contents

- [About CM5-Charge](#about-cm5-charge)
- [Install and Usage](#install-and-usage)
- [Code Structure](#code-structure)
- [Requirements](#requirements)

## About CM5-Charge

A repo for calculating the CM5-charge, including the `chargemol`, `cm5pac` and `CM5.sh` (which is use to automatically obtain the CM5-charge under a VASP job directory) files in [package](package) directory.

## Install and Usage

You can following the steps for install and use this tool.

1. typed the following command in the root directory (with the `CMakeLists.txt`):

```bash
cmake . -B build
```

then the `Makefile` will generated in the `build` directory.

2. compile the `chargemol`, `cm5pac`, `CM5` in `build` directory and install them in `bin` directory use the following command (if you want to see the compile process, use `make verbose=1`):

```bash
cd build && make && make install
```

> **Note**
> The `atomic_densities` will also be installed in `lib` directory.

3. (optional) typed the following command to perform a embedded test (**Note**: make sure you are in `build` directory)

```bash
make test
```

4. add the `bin` directory in `PATH` env

```bash
export PATH=$PATH:/the_directory_to_CM5_Charge/bin
```

then the `CM5` command can be used to calculate the CM5 charge!

## Code Structure

- [chargemol](package/chargemol): source files of chargemol
- [cm5pac](package/cm5pac): source files of cm5pac
- [CM5.sh](package/CM5.sh): automated script

## Requirements

- cmake [version >= 3.1]
- Fortran Compiler (*e.g.*, gfortran)

Copyright © 2022-2023 `Hui Zhou` All rights reserved.
