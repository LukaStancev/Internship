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

<details><summary>Reveal/hide details on the data release by clicking here</summary>

<!---
https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-collapsed-sections#creating-a-collapsed-section
-->

## Initial data in ENDF-6 format

TENDL-2019 (random files), JEFF-3.3, JEFF-3.2, JEFF-3.1.2 and JEFF-3.1.1 nuclear data can be downloaded and properly arranged with:

```
cd ~/SalinoPhD/ENDF
./DownloadENDF.sh
```

## Introductory remarks

From this initial data in ENDF-6 format, all of the following (below) data should be reproducible. To speed up and improve the reproducibility and understanding of this PhD thesis, the intermediate and final files are also provided below.

To access any of the data below, it is recommended to download and unpack directly in `~/SalinoPhD` directory:

```
cd ~/SalinoPhD
wget [...] # See below for specific links
tar -xvzf *.tgz
```

If done as such, the unpacking will place these data in the directory where they must be for the various tools to find them.

In case `tendl.web.psi.ch` website would come to disappear sooner than this repository, the random data retrieved from TENDL-2019 has been duplicated in this archive:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_TENDL-2019.tgz
```

## Intermediate data, also in ENDF-6 format

300 ENDF-6 files per isotope randomly sampled by SANDY within the variance-covariance matrices of JEFF-3.3 can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_H1-B-O16.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_U235_part1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_U235_part2.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_U238_part1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_U238_part2.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ENDF_SANDY_U238_part3.tgz
```

## Intermediate data : random ACE files and random Draglib files

For continuous energy Monte-Carlo codes like Serpent, the random ACE files (produced from the ENDF-6 files of SANDY or TENDL-2019) can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_B10-B11.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_H1_H2O.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_H1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_O16.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U235_part1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U235_part2.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U238_part1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U238_part2.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U238_part3.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_SANDY_U238_part4.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_AgIn.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_Cd.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_Fe54-Fe56.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_Ni58-Cr52.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_Zr90-Zr91.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/ACE_TENDL-2019_Zr92-Zr94-Zr96.tgz
```

The corresponding xsdata file (required by Serpent) and xsdir file (required by MCNP) can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/xsdata_xsdir.tgz
```

For the Dragon5 deterministic code, the random Draglib (produced consistently from the same ENDF-6 files, as above) can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_H.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_U23x.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_Zr-FeNiCr-AgInCd-B-O16.tgz
```

The following Draglib are the same ones but in ASCII format. It is easier to read for humans or for plotting, but they require binarization (i.e. the latter Draglib files) before they can be used in Dragon5.

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_H1_H2O.tgz.xaa
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_H1_H2O.tgz.xab
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_H1.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_O16-B10-B11.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_U235.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_SANDY_U238.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Draglib_ASCII_TENDL-2019_Zr-FeNiCr-AgInCd.tgz
```

Moreover, best estimate nuclear data were consistently used to produce [these ACE files](https://github.com/IRSN/PyNjoy2016/releases/tag/ACE_JEFF-3.x) and [those Draglib files](https://github.com/IRSN/PyNjoy2016/releases/tag/JEFF-3.x) (in binary format).

### A side note

One (single) archive had to be splitted in two. Before uncompressing it, it must be reunited with:

```
cat Draglib_ASCII_SANDY_H1_H2O.tgz.xaa Draglib_ASCII_SANDY_H1_H2O.tgz.xab > Draglib_ASCII_SANDY_H1_H2O.tgz
```

## Intermediate data : Serpent2 datasets

The [D2S](Serpent/D2S/D2S.py) tool has been fed with the following Version5 objects:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/D2S_input_Almaraz.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/D2S_input_Fessenheim-Bugey.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/D2S_input_Tihange.tgz
```

The resulting datasets from D2S are described above, in the [*Monte-Carlo Serpent2 datasets* section](https://github.com/IRSN/SalinoPhD#monte-carlo-serpent2-datasets).

## Final data : Drakkar output

Deterministic Drakkar outputs for JEFF-3.3 best estimate nuclear data can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Drakkar_Output_BestEstimate.tgz
```

This latter archive includes Tihange-1 case (assemblies in an infinite lattice, full cores), Almaraz-2 case and Bugey-2, Fessenheim-1 and 2 cases.

Deterministic Drakkar outputs for randomly sampled nuclear data (also called "Total Monte-Carlo") can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Drakkar_Output_TIH_TMC.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Drakkar_Output_FSH-BUG_TMC.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Drakkar_Output_FSH-BUG_TMC_interactions.tgz
```

These three archives include, respectively :
* Tihange-1 case, by varying the nuclear data of each isotope one by one,
* Bugey-2, Fessenheim-1 and 2 cases, by varying the nuclear data of each isotope one by one,
* Bugey-2, Fessenheim-1 and 2 cases, with simultaneous variations of all nuclear data (all isotopes), thus taking into account the interactions between the uncertainties of the different isotopes.

## Final data : Serpent2 output

Monte-Carlo Serpent2 outputs for JEFF-3.3 best estimate nuclear data can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_Tihange-Kinf-Assemblies.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_TIH-ARO_BatchHist.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_TIH-D-bank_BatchHist.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_TIH-CD-banks_BatchHist.tgz
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_Almaraz-Fessenheim-Bugey.tgz
```

These five archives include, respectively :
* assemblies in an infinite lattice for Tihange-1 case,
* full core for Tihange-1 case, all rods out,
* full core for Tihange-1 case, with D bank inserted,
* full core for Tihange-1 case, with C and D banks inserted,
* full core for Almaraz-2, Bugey-2, Fessenheim-1 and 2 cases.

Monte-Carlo Serpent2 outputs for randomly sampled nuclear data (also called "Total Monte-Carlo") can be downloaded with:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Serpent_Output_TIH_TMC.tgz
```

Only the isotopes determined as major had their nuclear data sampled, i.e.: <sup>238</sup>U, <sup>1</sup>H, <sup>56</sup>Fe, <sup>16</sup>O, <sup>235</sup>U, <sup>90</sup>Zr and <sup>107</sup>Ag.

From all the files mentioned above (from the initial ENDF-6 files up to the final output), the plots can be generated again.

## Plots

Only a careful selection of the figures produced are shown in the thesis and its appendices. However, the figures are exhaustively available in this archive:

```
wget https://github.com/IRSN/SalinoPhD/releases/download/Data-Release/Plots.tgz
```

## MD5 checksums

If an archive fails to unpack, it may be useful to check its MD5 checksum, available on [this page](https://github.com/IRSN/SalinoPhD/releases/tag/Data-Release).

</details>

## Disclaimer

The user of these files is solely responsible for their adequacy to his needs, the precautions to be taken, the qualification of his personnel and the use he makes of the results he obtains.
IRSN gives no guarantee concerning the results resulting from the use of these files as to their accuracy, reliability, up-to-date status or otherwise. Users use these files under their own responsibility and waive the right to take any action against IRSN regarding the consequences of their use of these files, particularly in the event of infringement proceedings brought against them by third parties.

L'utilisateur de ces fichiers est seul responsable de leur adéquation à ses besoins, des précautions à prendre, de la qualification de son personnel et de l'usage qu'il fait des résultats qu'il obtient.
L'IRSN ne donne aucune garantie concernant les résultats découlant de l'emploi de ces fichiers quant à leur exactitude, fiabilité, actualité ou autre. L'utilisateur exploite sous sa propre responsabilité ces fichiers, renonce à exercer à l'encontre de l'IRSN tout recours relatif aux conséquences de l'utilisation par elle-même des fichiers, et notamment en cas d'action en contrefaçon émanant de tiers à son encontre.

