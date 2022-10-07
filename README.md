[![IRSN](https://circleci.com/gh/IRSN/SalinoPhD.svg?style=shield)](https://circleci.com/gh/IRSN/SalinoPhD)

# Introduction

The production realized during this PhD thesis is centralized on this repository. This principle of openness is intended to allow, in a pragmatic and effective way:
* the verifiable and complete reproducibility of this publicly funded research,
* the transparency necessary for a rigorous peer review,
* the broadest dissemination of the ideas developed, with the aim of enhancing nuclear safety.

# Downloading

In order to retrieve all the necessary components, it is recommended to recover the material with:
```
git clone --recurse-submodules https://github.com/IRSN/SalinoPhD.git
```

# Monte-Carlo Serpent2 datasets

Five real reactor datasets are made available in [Serpent/FullCore](Serpent/FullCore). They may indeed be reused, with due citation. The [measurements](Measurements) results, found in (quoted) public data, are also available.

# Deterministic Version5 datasets execution

The [Version5](https://www.polymtl.ca/merlin/version5.htm) (Dragon5 and Donjon5) datasets can be executed with its [beta revision 1923](https://www.polymtl.ca/merlin/development.htm), as made available by its developers (École Polytechnique de Montréal). For the sake of convenience, this latter revision is also included in this repository. It can be compiled and then executed with the following commands:
```
cd SalinoPhD/Drakkar/Version5/Donjon/src
make
cd ../../../data
../runVersion5.sh Fessenheim-Bugey.x2m
```

This example is a simulation of the first start-up tests of [Fessenheim-1, Fessenheim-2 and Bugey-2](https://inis.iaea.org/collection/NCLCollectionStore/_Public/18/076/18076909.pdf). It is based on JEFF-3.3 nuclear data evaluation, with 295 energy groups. Other Draglibs can be downloaded [here](https://www.polymtl.ca/merlin/libraries.htm). From any ENDF-6 files, it is also possible to produce your own Draglibs using Python scripts in PyNjoy2016 (see more [there](https://github.com/IRSN/PyNjoy2016)).

![Radial view: first start-up cores of Almaraz-2, Bugey-2, Fessenheim-1 and 2](https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Radial_Almaraz_Fessenheim_Bugey.png "Radial view: first start-up cores of Almaraz-2, Bugey-2, Fessenheim-1 and 2")

More images can be downloaded [here](https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_images.tgz).

# Software requirements

For completeness, the precise version that was used successfully is indicated in parentheses. However, it is expected that many other versions would be compatible.

* Git (2.15.0)
* Bash (4.2.46(2))
* Make (3.16.2)
* CMake (3.16.2)
* GCC & GFortran (8.2.0)
* Python (3.6.9 ; see also the packages requirements [here](Plots/requirements.txt), [there](Serpent/D2S/requirements.txt) and [over there](https://github.com/luca-fiorito-11/sandy/blob/f944a451fa3d151a0b8a95583e398a4520289923/requirements.txt))
* Serpent (2.1.32)
* [serpentTools](https://github.com/CORE-GATECH-GROUP/serpent-tools) (0.9.3)
* Slurm (17.11.13-2). This software is only really necessary for massively parallel computations, i.e. nuclear data sampling (also called Total Monte-Carlo) either with deterministic or Monte-Carlo codes. Best estimate computations (without nuclear data sampling) either with Drakkar or Serpent will not benefit much from this software.

# Data release

The data produced during this PhD thesis are openly available on [this page](https://github.com/IRSN/SalinoPhD/releases/tag/Data-Release), along with numerous details ensuring their reusability.

## Disclaimer

The user of these files is solely responsible for their adequacy to his needs, the precautions to be taken, the qualification of his personnel and the use he makes of the results he obtains.
IRSN gives no guarantee concerning the results resulting from the use of these files as to their accuracy, reliability, up-to-date status or otherwise. Users use these files under their own responsibility and waive the right to take any action against IRSN regarding the consequences of their use of these files, particularly in the event of infringement proceedings brought against them by third parties.

L'utilisateur de ces fichiers est seul responsable de leur adéquation à ses besoins, des précautions à prendre, de la qualification de son personnel et de l'usage qu'il fait des résultats qu'il obtient.
L'IRSN ne donne aucune garantie concernant les résultats découlant de l'emploi de ces fichiers quant à leur exactitude, fiabilité, actualité ou autre. L'utilisateur exploite sous sa propre responsabilité ces fichiers, renonce à exercer à l'encontre de l'IRSN tout recours relatif aux conséquences de l'utilisation par elle-même des fichiers, et notamment en cas d'action en contrefaçon émanant de tiers à son encontre.

