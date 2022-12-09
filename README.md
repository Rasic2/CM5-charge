# CM5-Charge Manual

## Table of Contents

- [About CM5-Charge](#about-cm5-charge)
- [Install and Usage](#install-and-usage)
- [Code Structure](#code-structure)

## About CM5-Charge

A repo for calculating the CM5-charge, including the `chargemol`, `cm5pac` and `CM5.sh` (which is use to automatically obtain the CM5-charge under a VASP job directory) files.

## Install and Usage

You can following the steps for install and use this tool.

1. compile the chargemol in [chargemol/src](chargemol/src/) directory (just typed the `make` command, then the chargemol_parallel will be generated)

2. compile the cm5pac in [cm5pac](cm5pac) use the following command:

```bash
gfortran -o cm5pac.exe cm5pac.f
```

3. change the `chargemol` and `cm5pac` variables in CM5.sh to **real executable path**

## Code Structure

- [chargemol](chargemol): source files of chargemol
- [cm5pac](cm5pac): source files of cm5pac
- [CM5.sh](CM5.sh): automated script

Copyright © 2022 `Hui Zhou` All rights reserved.
