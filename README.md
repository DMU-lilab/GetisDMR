GetisDMR: DMR Detection 
========================================================

GetisDMR: Detection differentially methylated region in whole genome bisulfite sequencing (WGBS) data.

Contact Information
-------------------
Fushun Chen
chen_fushun@126.com

Installation
* Before attempting to compile GetisDMR please make sure that GNU Scientific Library (http://www.gnu.org/software/gsl/), the statistical software R and boost library are installed on your system.
* You need to install Rcpp, RInside, Rarmardillo.

Compiling (Example)

g++ GetisDMR_p_value.cpp FindFileHome.cpp FindDMR.cpp SplitString.cpp compare2Prop.cpp Getis_Ord.cpp RealDataMain.cpp one_vs_one.cpp BetaBinomialReg.cpp multi_vs_multi.cpp regression.cpp -I /usr/lib64/R/include/ -I ~/armadillo-5.100.2/include/ -I /usr/lib64/R/library/Rcpp/include/ -I /usr/lib64/R/library/RInside/include/ -I /usr/local/include/gsl/ -L /usr/local/lib/ -lgsl -lgslcblas -lm -std=c++11 -o3 -pipe -g -Wall /usr/lib64/libboost_program_options.so -lboost_system -lboost_filesystem -L /usr/lib64/R/lib -lR -lblas -llapack -L /usr/lib64/R/library/RInside/lib/ -lRInside -Wl,-rpath,/usr/lib64/R/library/RInside/lib -o GetisDMR

Usage
-----

Basic Command Line Parameters
Options:
*  -h [ --help ]            Print help messages
*  -c [ --Comparison ] arg  1: 1 vs 1 (no biological replicates) 
                           2: with biological replicates but no covariates 
                           3: biological replicates and covariates
*  --input1 arg             The file including the locations of files under 
                           treatment condition 1
*  --input2 arg             The file including the locations of files under 
                           treatment condition 2
*  -t [ --tolerance ] arg   tolerance
*  -l [ --length ] arg      length
*  --cutoff arg             cutoff
*  -o [ --outfiledir ] arg  The output folder: Default is the current folder
*  -s [ --sorted ] arg      The mtbr is sorted according to physical location
*  --CLINK_CPPFLAG arg      RcppArmadillo link
*  --MydistSoftware arg     MydistSoftware directory
*  --NewDistCppSoftware arg NewDistCppSoftware directory
*  --FindDMRfile arg        FindDMRSoftware.r directory
*  --sorted arg             If the mtbr is sorted according to physical 
                           position, currently not support un-sorted
*  --MINTOT arg             Minimum of the total number of reads required for 
                           each sample, default 5
*  --cov1 arg               The covariates file of treatment condition 1
*  --cov2 arg               The covariates file of treatment condition 2
*  -p [ --P_value ] arg     Calculate p-value or not

License
-------
Copyright (C) 2015 University of Auckland and Yalu Wen

    Authors: Yalu Wen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
