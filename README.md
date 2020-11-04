# ReactionGUI
GUI for easy and intuitive reaction calculation

Written by Y.-H. Song(yhsong@ibs.re.kr), I.J. Shin(geniean@ibs.re.kr)

## Disclamer 
    Copyright (c) 2020, All rights reserved.
    
    ReactionGUI is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    ReactionGUI is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
## Requirements
Recommend to install Anaconda( https://www.anaconda.com/products/individual ) for python environment 

1. Python 3.8 (The code is compiled with cpython 3.8. Thus, there may be a problem with lower version of python.)  
2. PyQt5 ( install with command "pip install pyqt5")
3. matplotlib ( install with command "pip install matplotlib")
4. numpy ( install with command "pip install numpy")
5. scipy ( install with command "pip install scipy")
---followings are included and thus no separate install/download is not necessary. ---
6. executable file of FRESCO and SFRESCO. 
  ( precompiled executable files are included.  
    Or Download/Compile from http://www.fresco.org.uk/ 
    and place it in the same folder of main.py ) 
7. executable file of OMGET and required data files. 
  ( precompiled executable files and data are included.
    Or Download/Compile from https://www-nds.iaea.org/RIPL-3/ (in the OPTICAL tab)
    and place it at the folder "omget_RIPL3/" ) 
8. density profile data files from HFB-14 calculation. 
  ( Located at folder "DoubleFolding/density-hfb14/". 
    Or Download from https://www-nds.iaea.org/RIPL-3/ ( in the MASSES tab)
   )
   
## Usage
There are two ways to run the code 

1. type 'python main.py' in the terminal.( within python environment) 
2. modify 'main.bat' file for the correct path and run batch file. 
3. For OSX and linux : set the executable file path using the 'Configure' menu.     

## Known Errors
1. Symptom : DLL load failed while importing xxx 

   Solution: update anaconda(or python packages) to newest version 
2. Symptom : No experimental data for fitting

   Solution : One have to click "exp data" button before fitting. 
   
   
   
