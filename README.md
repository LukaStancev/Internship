[![IRSN](https://circleci.com/gh/IRSN/SalinoPhD.svg?style=shield)](https://circleci.com/gh/IRSN/SalinoPhD)

## Introduction

The production realized during this PhD thesis is centralized on this repository. This principle of openness is intended to allow, in a pragmatic and effective way:
* the verifiable and complete reproducibility of this publicly funded research,
* the transparency necessary for a rigorous peer review,
* the broadest dissemination of the ideas developed, with the aim of enhancing nuclear safety.

## Version5 datasets execution

The [Version5](https://www.polymtl.ca/merlin/version5.htm) datasets can be executed with its [beta revision 1923](https://www.polymtl.ca/merlin/development.htm), as made available by its developers (École Polytechnique de Montréal). For the sake of convenience, this latter revision is also included in this repository. It can be downloaded, compiled and then executed with the following commands:
```
git clone https://github.com/IRSN/SalinoPhD.git
cd SalinoPhD/Drakkar/Version5/Donjon/src
make
cd ../../../data
../runVersion5.sh Tihange.x2m
```

This latter example is a simulation of the first start-up tests of [Tihange-1](https://inis.iaea.org/collection/NCLCollectionStore/_Public/11/511/11511367.pdf). It is based on JEFF-3.1.1 evaluation, with 172 energy groups (so-called XMAS energy mesh). Other Draglibs can be downloaded from: https://www.polymtl.ca/merlin/libraries.htm

From any ENDF files, it is also possible to produce your own Draglibs using Python scripts in PyNjoy2016 (see more [here](https://github.com/IRSN/PyNjoy2016)).

The power distribution can be accessed through that [kind of command](https://github.com/IRSN/SalinoPhD/blob/1abc854045630af1af45fc0e682fb4aee5cea29e/Drakkar/Reference/Diff.sh#L22) (for a quick peek) or plotted with Matplotlib through PyGan (a Python interface for Version5 ; an example [here](https://github.com/IRSN/SalinoPhD/blob/master/Plots/2Dpow.py)).

## Disclaimer

The user of these files is solely responsible for their adequacy to his needs, the precautions to be taken, the qualification of his personnel and the use he makes of the results he obtains.
IRSN gives no guarantee concerning the results resulting from the use of these files as to their accuracy, reliability, up-to-date status or otherwise. Users use these files under their own responsibility and waive the right to take any action against IRSN regarding the consequences of their use of these files, particularly in the event of infringement proceedings brought against them by third parties.

L'utilisateur de ces fichiers est seul responsable de leur adéquation à ses besoins, des précautions à prendre, de la qualification de son personnel et de l'usage qu'il fait des résultats qu'il obtient.
L'IRSN ne donne aucune garantie concernant les résultats découlant de l'emploi de ces fichiers quant à leur exactitude, fiabilité, actualité ou autre. L'utilisateur exploite sous sa propre responsabilité ces fichiers, renonce à exercer à l'encontre de l'IRSN tout recours relatif aux conséquences de l'utilisation par elle-même des fichiers, et notamment en cas d'action en contrefaçon émanant de tiers à son encontre.
