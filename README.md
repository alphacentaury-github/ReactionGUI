# ReactionGUI
GUI for easy and intuitive reaction calculation

Written by Y.-H. Song(yhsong@ibs.re.kr), I.J. Shin(geniean@ibs.re.kr)

2021-03-03. added GUI for fusion(CCFULL) and knockout(MOMDIS) reaction

2020-11-06. first version
    
## Requirements
Recommend to install Anaconda( https://www.anaconda.com/products/individual ) for python environment 

1. Python 3.8 or higher
2. PyQt5 ( install with command "pip install pyqt5")
3. matplotlib ( install with command "pip install matplotlib")
4. numpy ( install with command "pip install numpy")
5. scipy ( install with command "pip install scipy")

---followings are included and thus no separate install/download is necessary for WINDOWS. ---
---one have to prepare executable files for other OS.                                      ---

6. executable file of FRESCO and SFRESCO. 
  ( from http://www.fresco.org.uk/ ) 
7. executable file of OMGET and required data files. 
  ( from https://www-nds.iaea.org/RIPL-3/ (in the OPTICAL tab) ) 
8. density profile data files from HFB-14 calculation. 
  ( from https://www-nds.iaea.org/RIPL-3/ ( in the MASSES tab) )
9. executable file of CCFULL.
  ( from http://www2.yukawa.kyoto-u.ac.jp/~kouichi.hagino/ccfull.html ) 
10. executable file of MOMDIS.   
  ( from https://www.sciencedirect.com/science/article/pii/S0010465506001937?via%3Dihub )
  
## Usage
For OSX and linux : one have to set the executable file path using the 'Configure' menu.     

There are two ways to run the code 

1. type 'python main.py' in the (anaconda) terminal. 
2. modify 'main.bat' file for the correct path and run batch file in Windows. 

Input examples are available in example folder. (Use 'Open' in the menu.) 

## Limitation
ReactionGUI can compute (with a interface to Fresco )

1. elasctic scattering using optical potential 
2. inelastic scattering using Distorted Wave Born Approximation in rotor model
3. transfer reaction using Distorted Wave Born Approximation in cluster model
4. fusion reaction in complex optical potential or incoming wave boundary condition in coupled channel method.
5. one nucleon knockout reaction in eikonal approximation.

Capability of doing coupled cluster calculation, CDCC, radiative capture reaction
are planned but not yet available. 

## Known Problems
1. Symptom : "DLL load failed while importing xxx" 

   Solution: update anaconda(or python packages) to newest version 
   
2. Symptom : "No experimental data for fitting"

   Solution : One have to click "exp data" button before fitting.

3. Symptom: Global optical potential is not available.

   Solution: Check the path of executable of omget in the Configure menu
      
## Notice 
    Copyright (c) 2020, All rights reserved.
    
    ReactionGUI is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    ReactionGUI is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.   
