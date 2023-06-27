\page EABFMOD Extended-System Adaptive Biasing Force 

<!--
description: Methods for performing eABF or DRR method to calculate PMF along CVs
authors: Haochuan Chen, Haohao Fu
reference: \cite Chen2018 \cite Lelievre2007 \cite Lesage2016 \cite Fu2016
-->

\section Overview

This module contains the eABF/DRR method to do free energy calculation or enhance sampling along CVs.

\section Installation

This module is not installed by default and depends on the boost serialization module. Please make sure the boost serialization library is compiled and installed in your system before trying to compile this module. Add '\-\-enable-modules=drr \-\-enable-boost_serialization' to your './configure' command when building PLUMED to enable these features.

\section Usage

Please read \ref drr_tool and \ref DRR for more information.

\section Module Contents

- \subpage EABFMODBias
- \subpage EABFMODCLTools

\page EABFMODBias Biases Documentation

The following list contains descriptions of biases developed for the eABF module. They can be used in combination with other biases outside of the eABF module.

@EABFMOD_BIAS@

\page EABFMODCLTools Command Line Tools

The following list contains the command line tools available in the eABF module.

@EABFMOD_TOOLS@
