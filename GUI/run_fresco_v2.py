"""
!***********************************************************************
!     
!    Copyright (c) 2020, ReactionGUI GUI for reaction calculation
!                        Produced at Rare Isotope Science Project 
!                        Written by Young-Ho Song, yhsong@ibs.re.kr 
!                               and Ik-Jae Shin,  geniean@ibs.re.kr 
!                        All rights reserved.
!          
!    This file is part of ReactionGUI.
!
!    ReactionGUI is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!     
!    ReactionGUI is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!     
!    You should have received a copy of the GNU General Public License
!    along with FRESCO. If not, see <http://www.gnu.org/licenses/>.
!     
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!     
!    Our Preamble Notice
!
!      A. This work was produced at the Rare Isotope Science Project 
!         of Institute for Basic Science, 
!         funded by Ministry of Science and ICT and 
!         by National Research Foundation of Korea (2013M7A1A1075764).
!      B. Neither the Rare Isotope Science Project nor any of their employees, 
!         makes any warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!        
!***********************************************************************
"""
import os
import shutil
import sys
import numpy as np
import numpy.linalg as npla

from subprocess import (call,Popen)
import myutil
from myutil import read_fresco_res ,clean_comm

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]

#-----------------------------------------------------------------
class Fresco_Head:
    """
    Head part of the Fresco input

    Follows general names of FRESCO namelist format

    Property
    .data = Dictionary of Fresco Head information

    .set_item(keyward arguments) = add/modify Dictionary items
    .remove_item(key) = delete item from dictionary
    .write()          = return Fresco input text
    """
    # all possible items
    # if item name starts 'i,j,l,m,n' -> integer
    #                     'c'    -> complex
    #                     others -> float
    # but not always...
    float_variables = ['HCM','RMATCH','RINTP','HNL','RNL','CENTRE',
                       'RASYM','ACCRCY','SWITCH','AJSWITCH',
                       'JTMIN','JTMAX','ABSEND',
                       'JBORD',
                       'PP',
                       'THMIN','THMAX','THINC',
                       'CUTL','CUTR','CUTC',
                       'IPS',
                       'ERANGE','DK','SMALLCHAN','SMALLCOUP',
                       'ELAB',
                       'UNITMASS','FINEC'
                       ]
    int_variables = ['JSETS','NEARFA',
                     'JUMP',
                     'KQMAX','KOORDS',
                     'IT0','ITER','IBLOCK' ,'PADE',
                     'NNU','MAXL','MINL','MTMIN','EPC',
                     'INH', 'PLANE',
                     'INITWF',
                     'CHANS','LISTCC','TRENEG','CDETR',
                     'NLPL','WAVES','LAMPL','VEFF','KFUS',
                     'WDISK', 'BPM',
                     'PEL','EXL','LAB','LIN','LEX',
                     'NLAB',
                     'PSET','JSET',
                     'SMATS','XSTABL' ]
    str_variables = ['DRY','RELA','ISO',
                     'FATAL','NOSOL','PSIREN','TMP','MASFIL' ]
    def __init__(self, title):
        self.data = {'TITLE': title}
        # set default values
        self.data['HCM'] = 0.01
        self.data['RMATCH'] = 20.0
        self.data['JTMIN'] = 0.0
        self.data['JTMAX'] = 100.0
        self.data['THMIN'] = 0.0
        self.data['THMAX'] = 180.0
        self.data['THINC'] = 1.0
        self.data['CHANS'] = 1
        self.data['SMATS'] = 2
        self.data['XSTABL'] =1
        self.data['RNL'] =10.0
        self.data['ITER'] =1
        self.data['IBLOCK']=0

    def set_items(self,**attributes):
        for item,value in attributes.items():
            if (item in (self.str_variables
                         +self.int_variables
                         +self.float_variables
                         +['TITLE']) ): # title is an exception
                self.data[item]=value
            else :
                print('invalid input item name %s\n'%item)

    def remove_item(self,item):
        if item in self.data.keys():
            del self.data[item]

    def write(self,):
        """
        Write Fresco input according to heads information
        """
        text = self.data['TITLE']+'\n'
        text = text+"NAMELIST\n"
        text = text+'&FRESCO \n'
        keys = list(self.data.keys())
        keys.remove('TITLE') # title is already used
        count=0
        for item in keys:
            if item in self.str_variables:
                text = text+' {}=\'{}\''.format(item,self.data[item])
            else:
                text = text+' {}={} '.format(item,self.data[item])
            count = count+1
            if count==5: # maximum number of items in one line
                text = text+'\n'
                count = 0

        return text + '   /\n'

#----------------------------------------------------------------------------
class Fresco_Partitions():
    """
    Partition part of Fresco input
    only one partition with many states
    """
    #str_variables=['NAMEP','NAMET','PWF']
    str_variables=['NAMEP','NAMET'] # PWF=T or PWF=F without \'
    float_variables=['MASSP','MASST','ZP','ZT','QVAL',
                     'JP','JT','EP','ET','TP','TT',
                     'KKP','KKT' ]
    int_variables=['NEX','COPYP','COPYT','BANDP','BANDT',
                   'CPOT'  ]
    partition_variables = ['NAMEP','MASSP','ZP','NAMET','MASST','ZT',
                          'NEX', 'QVAL', 'PWF' ]
    state_variables = ['JP','JT','EP','ET','COPYP','COPYT',
                        'BANDP','BANDT','KKP','KKT',
                        'TP','TT','CPOT']
                     #'PWF' is not clear how to use
    def __init__(self,NAMEP='proj',MASSP=0.0,ZP=0.0,
                 NAMET='targ',MASST=0.0,ZT=0.0
                 ,QVAL=0.0,NEX=0):
        self.data = [ ]
        self.data.append( {'NAMEP': NAMEP,'MASSP':MASSP,'ZP':ZP,
                          'NAMET': NAMET,'MASST':MASST,'ZT':ZT,
                          'NEX': NEX, 'QVAL':QVAL,
                          'STATES': [] } )
                       # states are treated as a list with index
    def add_partition(self,NAMEP='proj',MASSP=0.0,ZP=0.0,
                 NAMET='targ',MASST=0.0,ZT=0.0
                 ,QVAL=0.0,NEX=0,**attributes):
        """
        adding new partition
        Use Keyward arguments of Fresco input in Capital letters
        """
        temp_dictionary={'NAMEP': NAMEP,'MASSP':MASSP,'ZP':ZP,
                          'NAMET': NAMET,'MASST':MASST,'ZT':ZT,
                          'NEX': NEX, 'QVAL':QVAL,
                          'STATES': [] }
        for item, value in attributes.items():
            if item in self.partition_variables:
                temp_dictionary[item]= value
            else:
                print('item: %s is not valid'%(item))
        # rename partition name
        try:
            if temp_dictionary['NAMEP']=='proj':
                temp_dictionary['NAMEP']='{}{}'.format(
                round(temp_dictionary['MASSP']),
                element_names[temp_dictionary['ZP']]
                )
            if temp_dictionary['NAMET']=='targ':
                temp_dictionary['NAMET']='{}{}'.format(
                round(temp_dictionary['MASST']),
                element_names[temp_dictionary['ZT']]
                )
        except:
            pass
        self.data.append(temp_dictionary) # need to check ?

    def change_partition(self,partition_index,**attributes):
        """
        change partition info of existing partition
        """
        for item, value in attributes.items():
            if item in self.partition_variables:
                self.data[partition_index][item] = value
            else:
                print('item: %s is not valid'%(item))
        # rename partition name
        if self.data[partition_index]['NAMEP']=='proj':
            self.data[partition_index]['NAMEP']='{}{}'.format(
                round(self.data[partition_index]['MASSP']),
                element_names[self.data[partition_index]['ZP']])
        if self.data[partition_index]['NAMET']=='targ':
            self.data[partition_index]['NAMET']='{}{}'.format(
                round(self.data[partition_index]['MASST']),
                element_names[self.data[partition_index]['ZT']])

    def remove_partition(self,partition_index):
        """
        remove partition
        """
        if partition_index < len(self.data):
            del self.data[partition_index]
            
    def clear_partition(self,partition_index):
        """
        clear partition part information (This removes existing states)
        but, does not remove parition index 
        """
        self.data[partition_index]={} 

    def add_state(self,partition_index,**attributes):
        """
        add new state into the partition
        """
        #set default values
        temp_dictionary ={}
        for item, value in attributes.items():
            if item in self.state_variables:
                temp_dictionary[item]= value
        self.data[partition_index]['STATES'].append(temp_dictionary)
        self.data[partition_index]['NEX'] += 1

    def remove_state(self,partition_index,state_index):
        """
        remove existing state from partition
        """
        del self.data[partition_index]['STATES'][state_index]
        self.data[partition_index]['NEX'] -= 1

    def change_state(self,partition_index,state_index,**attributes):
        """
        change existing state in partition
        """
        for item, value in attributes.items():
            if item in self.state_variables:
                self.data[partition_index]['STATES'][state_index][item] = value
                
    def clear_state(self,partition_index,state_index):
        """
        clear all previous state information 
        """            
        self.data[partition_index]['STATES'][state_index]={} 

    def write(self,):
        """
        print Fresco input for one partition part
        """
        text = ""
        for partition in self.data :
            text = text + '&PARTITION '
            keys = list( partition.keys() )
            keys.remove('STATES') # states will be written later
            # write partition line
            for key in keys:
                if key in self.str_variables:
                    text = text +" {}=\'{}\'".format(key,partition[key])
                elif key in self.float_variables:
                    text = text +' {}={:.5f}'.format(key,partition[key])
                else:
                    text = text +' {}={}'.format(key,partition[key])
                    
            text = text +'/\n'
            # write states(excitation pairs)
            for state in partition['STATES']:
                text = text+'  &STATES '
                for key in state.keys():
                    if key in self.str_variables:
                        text = text +" {} = \'{}\'".format(key,state[key] )
                    elif key in self.float_variables:
                        text = text +' {}={:.5f}'.format(key,state[key])    
                    else:
                        text = text +' {} = {}'.format(key,state[key] )
                text = text+' / \n'

        return text +'&PARTITION /\n'

#-----------------------------------------------------------------------------
class Fresco_Potentials:
    """
    Note that in FRESCO the order of potential terms is important!
    Thus, it is stored as a list.

    The part of potential is composed of potentials and terms.

    potentials is collection of terms which have the same kp index.
    each terms are collection of functions for given kp.

    In other words,
    each line of potential parts of Fresco input are 'terms'.
    And 'terms' belongs to a certain 'potentials' of kp index.

    data = [ pot1=[term1,term2] ,  pot2=[term1,term2...] , ...  ]

    each term is a dictionary.
    Keywords follows the FRESCO convention in Captical letters except RC.
    Use P1,P2,P3,P4 instead of (AT,AP,RC,AC)

    Unlike other parts of Fresco input,
    .add_potential adds a Coulomb term.
    Thus, to input Coulomb potential one either have to
    (1) add_potential with Coulomb parameters or
    (2) use change_terms
    """
    str_variables =[]
    int_variables =['KP','TYPE','SHAPE','ITT']
    float_variables =['AT','AP','RC','AC',
                      'P1','P2','P3','P4','P5','P6','P7',
                      'DEFP','DEFT','MNET','MNEP']
    def __init__(self):
        """
        initialize potential with empty list.

        externalPotential= {(potential_index,term_index):
                             {'R': array, 'V': complex_array }, .. }

        """
        self.data = []
        self.externalPotential = {}

    def add_potential(self,TYPE=0,SHAPE=0,P1=0,P2=0,P3=0,P4=0,P5=0,P6=0,P7=0):
        """
        if (pot_type,shape)==(0,0),(AT,AP,R0C,AC) <- p1,p2,p3,p4
        """
        self.data.append([{'TYPE': TYPE,'SHAPE':SHAPE,
                          'P1':P1,'P2':P2,'P3':P3,'P4':P4,'P5':P5,
                          'P6':P6,'P7':P7 }]  )  # add list not dictionary

    def remove_potential(self,pot_index):
        del self.data[pot_index]
    # no need of change_potential.change_term will do

    def add_term(self,pot_index,TYPE=0,SHAPE=0,P1=0,P2=0,P3=0,
                 P4=0,P5=0,P6=0,P7=0):
        """
        add a new term to a potential
        
        What if the term use external potential?? 
        """
        self.data[pot_index].append({'TYPE': TYPE,'SHAPE':SHAPE,
                          'P1':P1,'P2':P2,'P3':P3,'P4':P4,'P5':P5,
                          'P6':P6,'P7':P7 }) #add dictionary
        return

    def change_term(self, pot_index,term_index,**attributes ) :
        """
        change/set parameter values of existing terms
        """
        for item, value in attributes.items():
            if (item in self.int_variables
                        +self.float_variables+self.str_variables):
                self.data[pot_index][term_index][item] = value
            else:
                print('item: %s is not valid'%(item))
        return

    def remove_term(self,pot_index, term_index):
        del self.data[pot_index][term_index]

    def set_external_potential(self,pot_index,term_index,r_points,pot_array):
        """
        bind external potential to a term

        Parameters
        ----------
        pot_index : integer

        term_index : integer

        r_points : 1-d array

        pot_array : 1-d array
            if term shape=7,8 it is a real
            if term shape=9 it sould be complex

        Returns
        -------
        None.

        """
        self.externalPotential[( pot_index,term_index)]= {'R': r_points,
                                   'V': pot_array }
        return

    def write(self,):
        """
        write a text input for FRESCO potential

        """
        text=""
        KP = 0
        for (pindx, potential) in enumerate(self.data) :
            KP = KP +1
            for (tindx,term) in enumerate(potential) :
                text = text + '&POT KP={}'.format(KP)
                if (term['TYPE']==0 and term['SHAPE']==0):  #Coulomb term case
                    text = text +' TYPE=%i SHAPE=%i AP=%.3f AT=%.3f RC=%.3f / \n'%(
                                    term['TYPE'],term['SHAPE']
                                    ,term['P1'],term['P2'],term['P3'])
                elif (term['SHAPE'] in [7,8,9]) : # external file input
                    fmt=' TYPE=%i SHAPE=%i P1=%.3f P2=%.3f P3=%.3f  / \n'
                    text = text +fmt%(term['TYPE'],term['SHAPE']
                                    ,term['P1'],term['P2'],term['P3'] )
                    # need to write external potential to 'fort.4' files
                    dat = np.column_stack(
                        (self.externalPotential[(pindx,tindx)]['R'],
                         self.externalPotential[(pindx,tindx)]['V']) )
                    radial_points = self.externalPotential[(pindx,tindx)]['R']
                    x_range= (radial_points[0],radial_points[-1]+1,
                              radial_points[1]-radial_points[0])
                    if term['SHAPE'] in [7,8]:
                        num_type='R'
                    elif term['SHAPE']==9 :
                        num_type='C'
                    myutil.fresco_input_formfactor(dat,x_range,
                            'fort.4',fill_value='extrapolate',num_type= num_type,
                            comment='test folding potential',
                            write='a')
                else :
                    fmt=' TYPE=%i SHAPE=%i P1=%.3f P2=%.3f P3=%.3f P4=%.3f P5=%.3f P6=%.3f P7=%.3f / \n'
                    text = text +fmt%(term['TYPE'],term['SHAPE']
                                    ,term['P1'],term['P2'],term['P3'],
                                    term['P4'],term['P5'],term['P6'],term['P7'])
        return text +'&POT /\n'

#----------------------------------------------------------------------------
class Fresco_Overlaps():
    """
    Overlap part of Fresco input

    Because 'in' variable cannot be used ,
    here all variables are in Capital letters.

    """
    str_variables =['CH1','KEEP']
    int_variables =['KN','KN1','KN2','IC1','IC2','IN','KIND',
                    'NN','L','LMAX','IA','IB',
                    'KBPOT','KRPOT','ISC','IPC','NFL','NAM',
                    'NK','NLAG']
    float_variables =['SN','J','JN','BE','AMPL','ER','RSMIN','RSMAX',
                      'RSALPHA','PHASE','DM','E']
    def __init__(self,):
        self.data = []
        self.num_form_factor =0

    def add_overlap(self,**attributes):
        """
        For the moment, it assumes that kn2 = 0 case only.
        """
        temp_dictionary ={}
        for item, value in attributes.items():
            if ( (item in self.str_variables)
                or (item in self.int_variables)
                or (item in self.float_variables)):
                temp_dictionary[item]= value
        self.data.append(temp_dictionary)

    def change_overlap(self,overlap_index,**attributes):
        for item, value in attributes.items():
            if ( (item in self.str_variables)
                or (item in self.int_variables)
                or (item in self.float_variables)):
                self.data[overlap_index][item]= value

    def remove_overlap(self,overlap_index):
        del self.data[overlap_index]

    def write(self,):
        text = ""
        for overlap in self.data :
            text = text + '&OVERLAP '
            keys = overlap.keys()
            for key in keys:
                if key in self.str_variables:
                    text = text + " {}=\'{}\'".format( key, overlap[key] )
                elif key in self.float_variables:
                    text = text +' {}={:.5f}'.format(key,overlap[key])    
                else :
                    text = text + " {}={}".format( key, overlap[key] )
            text = text + ' / \n'
        return text + '&OVERLAP / \n'
#----------------------------------------------------------------------------
class Fresco_Couplings():
    """
    One Coupling part of Fresco input
    """
    str_variables=[]
    int_variables=['ICTO','ICFROM','KIND','IP1','IP2','IP3','IP4','IP5',
                   'KFRAG','KCORE','IN','IB','IA','KN','A']
    float_variables=['P1','P2','JMAX','RMX']

    CFP_variables = [ 'IN','IB','IA','KN','A','KEEP','NO','K','QSCALE']
    def __init__(self):
        self.data = []

    def add_a_coupling(self,**attributes):
        temp_dictionary={}
        for item,value in attributes.items():
            if (item in (self.str_variables
                         +self.int_variables
                         +self.float_variables) ):
                temp_dictionary[item]=value
        self.data.append(temp_dictionary)

    def change_a_coupling(self,coupling_index,**attributes):
        for item,value in attributes.items():
            if (item in (self.str_variables
                         +self.int_variables
                         +self.float_variables) ):
                self.data[coupling_index][item] = value
    def remove_coupling(self,coupling_index):
        del self.data[coupling_index]

    def add_cfp(self,coupling_index,**attributes):
        # If CFP does not exist yet
        if not('CFP' in self.data[coupling_index]):
            self.data[coupling_index]['CFP']=[]
        temp_dictionary={}
        for item,value in attributes.items():
            if (item in self.CFP_variables ):
                temp_dictionary[item] = value
        self.data[coupling_index]['CFP'].append(temp_dictionary)

    def change_cfp(self,coupling_index, cfp_index,**attributes):
        for item,value in attributes.items():
            if (item in self.CFP_variables ):
                self.data[coupling_index]['CFP'][cfp_index][item]=value

    def remove_cfp(self,coupling_index, cfp_index):
        del self.data[coupling_index]['CFP'][cfp_index]

    def write(self):
        text=""
        for coup in self.data :
            text = text+'&COUPLING'
            keys =  list(coup.keys())
            try:
                keys.remove('CFP') # CFP will be treated later
            except:
                pass

            for item in keys:
                if item in self.str_variables:
                    text = text+" {}=\'{}\'".format(item,coup[item])
                elif item in self.float_variables:
                    text = text +' {}={:.5f}'.format(item,coup[item])    
                else:
                    text = text+" {}={}".format(item,coup[item])
            text = text +' /\n'
            try: # CFP case
                for cfp in coup['CFP'] :
                    text = text +' &CFP '
                    for item in cfp.keys():
                        text = text +' {}={}'.format(item,cfp[item])
                    text = text +' /\n'
                text = text +' &CFP /\n'
            except: # no CFP
                pass
        return text+'&COUPLING / \n'
#-------------------------------------------------------------------------
class Fresco_input():
    """
    Class to store informations to run Fresco code

    they are stored in

    self.head
    self.partitions
    self.potentials
    self.overlaps
    self.couplings

    each have its methods of add/remove/change

    """
    def __init__(self,title='',head=None,partitions=None
                 ,potentials=None,overlaps=None,couplings=None):
        """
        Initialize the class.
        If objects are not given in initialization  create it.
        Parameters
        ----------
        title : string
            short description of the reaction
        head : Fresco_Head object, optional
        partitions : Fresco_partiiton object, optional
        potentials : Fresco_Potential object, optional
        overlaps : Fresco_overlaps object , optional
        couplings : Fresco_coupling object , optional

        """
        if head :
            self.head = head
        else :
            self.head = Fresco_Head(title)
        if partitions :
            self.partitions = partitions
        else :
            self.partitions = Fresco_Partitions()
        if potentials :
            self.potentials = potentials
        else:
            self.potentials = Fresco_Potentials()
        if overlaps :
            self.overlaps = overlaps
        else:
            self.overlaps = Fresco_Overlaps()
        if couplings :
            self.couplings = couplings
        else :
            self.couplings = Fresco_Couplings()

    def write(self,):
        """
        write the input text of FRESCO

        But, it is important to delete or not fort.4 files
        before writing !!

        """
        #delete previous fort.4 file
        try:
            os.remove('fort.4')
            #delete_command='del' # OS_dependence ?
            #print('To do: generalize delete command for os dependent')
            #proc = Popen(delete_command+' '+' fort.4',shell=True )
            #proc.wait()
            #proc.terminate()
        except:
            pass
        text = ''
        text = text+self.head.write()
        text = text+self.partitions.write()
        text = text+self.potentials.write() # write fort.4 for external potential
        text = text+self.overlaps.write()
        text = text+self.couplings.write()
        return text

    def example_CCBA(self,):
        """
        setting an example of CRC calculation
        'CCBA 28Si(19F,16O)31P , cluster form factor'
        """
        self.head = Fresco_Head('CCBA 28Si(19F,16O)31P , cluster form factor')
        self.head.set_items(HCM=0.10, RMATCH=25.0, RINTP=0.50, HNL=0.100, RNL=3.0,
                           JTMIN=0., JTMAX=80.0, ABSEND=0.01,
                           KQMAX=1,THMIN=0.0,THMAX=60.0, THINC=2.5,
                           ITER=1, IBLOCK=2, NNU=30,
                           CHANS=1, LISTCC=2, SMATS=1,
                           ELAB= 60.0)
        self.partitions = Fresco_Partitions() #first partition
        self.partitions.change_partition(0,
                        NAMEP='19F',MASSP=19.0,ZP=9,
                        NAMET='28Si',MASST=28.0,ZT=14)
        self.partitions.add_state(0,
                        JP=0.5,BANDP=1,EP=0.0,
                        CPOT=1,
                        JT=0.0,BANDT=1,ET=0.0)
        self.partitions.add_state(0,
                        JP=2.5,BANDP=1,EP=0.2)
        self.partitions.add_partition(NAMEP='16O',MASSP=16.0,ZP=8,
                                      NAMET='31-P',MASST=31.0,ZT=15,
                                      QVAL=6.199)
        self.partitions.add_state(1,
                        JP=0.0,BANDP=1,EP=0.0,
                        CPOT=2,
                        JT=0.5,BANDT=1,ET=0.0)
        self.potentials = Fresco_Potentials()
        self.potentials.add_potential() #KP=1
        self.potentials.change_term(0,TYPE=0,AT=28.0,AP=19.0,RC=1.2)
        self.potentials.add_potential() #KP=2
        self.potentials.add_term(1,TYPE=0,AT=31.0,AP=16.0,RC=1.35)
        self.potentials.add_term(1,TYPE=1,P1=31.2,P2=1.45,P3=0.47,
                                          P4=15.10,P5=1.27,P6=0.31)
        self.potentials.add_term(1,TYPE=3,P1=0.75,P2=1.24,P3=0.37)
        self.potentials.add_potential() #KP=3
        self.potentials.add_term(2,TYPE=0,AT=19.0,RC=1.25,AC=0.65)
        self.potentials.add_term(2,TYPE=1,P1=115.0,P2=1.25,P3=0.65)
        self.potentials.add_term(2,TYPE=2,P1=6.30,P2=1.25,P3=0.65)



    def example_elastic(self,):
        self.head = Fresco_Head('p + 112Cd elastic')
        self.head.set_items(HCM= 0.100,RMATCH=  20.000,JTMIN=  0.0,JTMAX=200.0,
                            ABSEND=0.001,
                     THMIN=  0.00,THMAX=180.00,THINC=2.00,XSTABL=1,
                     CHANS=1, SMATS=2,
                     ELAB=27.9)
        self.partitions = Fresco_Partitions()
        self.partitions.change_partition(0,
                    NAMEP='Proton  ',MASSP=  1.0, ZP=  1,
                    NAMET='112Cd   ',MASST=112.0, ZT= 48, QVAL=  0)
        self.partitions.add_state(0,
                    JP= 0.5, BANDP= 1, EP=  0.0,
                    CPOT=  1,
                    JT= 0.0, BANDT= 1, ET=  0.0)
        self.potentials = Fresco_Potentials()
        self.potentials.add_potential()
        self.potentials.change_term(0,0,TYPE=0,SHAPE=0,P1=112,P3=1.2)
        self.potentials.add_term(0,TYPE=1,SHAPE=0,
                       P1=52.500,P2=1.17,P3=0.75,
                       P4=3.5000,P5=1.32,P6=0.61 )
        self.potentials.add_term(0,TYPE=2,SHAPE=0,
                       P4= 8.5000,P5=1.3200,P6=0.6100 )
        self.potentials.add_term(0,TYPE=3,SHAPE=0,
                       P1=6.200,P2=1.0100,P3=0.7500 )
        self.overlaps = Fresco_Overlaps()
        self.couplings = Fresco_Couplings()

    def example_inelastic_rotor(self,):
        self.head = Fresco_Head('alpha+12C->alpha+12C* @ 100 MeV; nuc deform')
        self.head.set_items(HCM=0.050, RMATCH=20.0,
                            JTMIN=0.0, JTMAX=40.0,
                            THMIN=0.0, THMAX=180.0,THINC=1.0,
                            CHANS=1, SMATS=2,XSTABL=1,
                            ELAB = 100.0,
                            ITER=1, IPS=0.0, IBLOCK=0, # new for inelastic
                            RINTP=0.20, ABSEND=0.01
                            )
        self.partitions = Fresco_Partitions()
        self.partitions.change_partition(0,
                       NAMEP='alpha',MASSP=4.0, ZP=2,
                       NAMET='12C', MASST=12.0, ZT=6, QVAL=0)
        self.partitions.add_state(0,
                      JP=0.0, BANDP=1,EP=0.0,
                      JT=0.0, BANDT=1,ET=0.0,
                      CPOT=1 )
        self.partitions.add_state(0,
                       COPYP=1,
                       JT=2.0,BANDT=1,ET=4.4390, CPOT=1)
        self.potentials = Fresco_Potentials()
        self.potentials.add_potential()
        self.potentials.change_term(0,0,TYPE=0,P1=4.0,P2=12.0,P3=1.2 )
        self.potentials.add_term(0,TYPE=1,P1=40.0,P2=1.2,P3=0.65,P4=10.0,P5=1.2,P6=0.5)
        self.potentials.add_term(0,TYPE=11,SHAPE=10,P2=1.3)
        self.potentials.add_term(0,TYPE=2, P1=0.0,P2=1.2,P3=0.65,P4=6.0,P5=1.2,P6=0.5 )
        self.potentials.add_term(0,TYPE=11,SHAPE=10,P2=1.3)
        self.overlaps = Fresco_Overlaps()
        self.couplings = Fresco_Couplings()

    def example_dfold_elastic(self,):
        """
        Example input for the case of double folding potential
        This requires external potential in 'fort.4' file.

        Parameters in this example is not realistic.
        Just for testing format.
        """
        #----external potential------
        R=np.arange(0.,12.,0.02 ); R0=1.25*12**(1./3.);a0=0.65;
        test_folded_potential={'R': R,
                               'Isoscalar': 1./(1+np.exp(-(R-R0)/a0)),
                               'Isovector': None }
        #------------------------------
        self.head = Fresco_Head('alpha+12C->alpha+12C* @ 100 MeV; double folding M3Y')
        self.head.set_items(HCM=0.050, RMATCH=20.0,
                            JTMIN=0.0, JTMAX=40.0,
                            THMIN=0.0, THMAX=180.0,THINC=1.0,
                            CHANS=1, SMATS=2,XSTABL=1,
                            ELAB = 100.0,
                            ITER=1, IPS=0.0, IBLOCK=0, # new for inelastic
                            RINTP=0.20, ABSEND=0.01
                            )
        self.partitions = Fresco_Partitions()
        self.partitions.change_partition(0,
                       NAMEP='alpha',MASSP=4.0, ZP=2,
                       NAMET='12C', MASST=12.0, ZT=6, QVAL=0)
        self.partitions.add_state(0,
                      JP=0.0, BANDP=1,EP=0.0,
                      JT=0.0, BANDT=1,ET=0.0,
                      CPOT=1 )
        self.potentials = Fresco_Potentials()
        self.potentials.add_potential()
        self.potentials.change_term(0,0,TYPE=0,P1=4.0,P2=12.0,P3=1.2 )
        self.potentials.add_term(0,TYPE=1,SHAPE=9,P1=1.0,P2=1.0)
        #--- prepare external input file-------------------------------
        temp1=np.arange(0.01,11.0,0.05) #for test external potential
        temp2=np.exp(-temp1)/temp1
        self.potentials.set_external_potential(0,1,temp1,temp2)
        #-------------------------------------------------------------
        self.overlaps = Fresco_Overlaps()
        self.couplings = Fresco_Couplings()

#----------------------------------------------------------------------------
class SFresco_Input():
    """
    class to store experimental data information
    it can be used as an input to SFRESCO
    """
    str_variables=['NAME','NOPT','DATA_FILE']
    int_variables=['KIND','KP','PLINE','COL','NAFRAC',
                   'TERM','PAR','CHANNEL',
                   'DATASET',
                   'TYPE','POINTS' ]
    float_variables=['STEP','VALMIN','VALMAX','NULL',
                     'POTENTIAL','DATANORM',
                     'AFRAC',
                     'ENERGY','JTOT','WIDTH',
                     'DELTA']

    Var_keywords=['NAME','KIND','STEP','VALMIN','VALMAX',
                   'NULL','KP','PLINE','COL',
                   'POTENTIAL','DATASET','DATANORM',
                   'NAFRAC','AFRAC',
                   'ENERGY','JTOT','PAR','CHANNEL','WIDTH']
    Data_keywords=['TYPE','DATA_FILE','POINTS',
                    'DELTA','XMIN','LAB','ENERGY',
                    'ANGLE','IDIR','ISCALE',
                    'ABSERR','IC','IA','K',
                    'Q','JTOT','PAR','CHANNEL',
                    'VALUE','ERROR','PEL','EXL',
                    'LABE','LIN','LEX' ]

    def __init__(self,fresco_input_object=None ):
        self.fresco_input = fresco_input_object
        self.search_var=[ ]
        self.data = "" # do not use self.data directly!  Use set_exp_data
        self.data_input = {'TYPE': 1,'ISCALE': 0,'IDIR': 0,
                           'LAB': 'F','ABSERR': 'T' }

    def set_data_info(self,**attributes):
        """
        set &DATA line of search file
        """
        for item,value in attributes.items():
            if item in self.Data_keywords:
                self.data_input[item]=value

    def remove_data_info(self,item):
        if item in self.Data_keywords:
            del self.data_input[item]

    def set_exp_data(self,exp_data):
        # exp data text
        # sfreco does not take comments in exp_data 
        # remove comments  
        new = ''
        exp_data = exp_data.strip()
        ll =exp_data.split('\n') 
        for i in ll :
            if i: # not  empty line 
                if not(i.strip()[0] in ['#','!','\n']):
                    new = new + i+'\n'
        self.data = new         

    def write_data(self,):
        """
         write SFresco input data part
        """
        txt ='&DATA'
        txt = txt +' TYPE={}'.format(self.data_input['TYPE'])
        txt = txt +' ISCALE={}'.format(self.data_input['ISCALE'])
        txt = txt +' IDIR={}'.format(self.data_input['IDIR'] )
        txt = txt +" LAB={}".format(self.data_input['LAB'])
        txt = txt +" ABSERR={}".format(self.data_input['ABSERR'])
        txt = txt +' / \n'
        # experimental data
        txt = txt + self.data
        txt = txt+'&\n'
        return txt

    def write_minuit(self,fname='_test.min',search_file='_test.search'):
        """
        write the input file for minuit
        to fname

        plot files are generated with
        fname_head+'_init.plot'
        fname_head+'_fit.plot'
        """
        fname_head = fname.split('.')[0]
        ff=open(fname,'w')
        txt ='\'{}\' \n'.format( search_file ) #card1
        txt = (txt
               +'q\n'
               +'plot {}_init.plot\n'.format(fname_head)
               +'min\n'                            #card2
               +'migrad \n'
               +'end\n'
               +'plot {}_fit.plot\n'.format(fname_head)
               +'ex\n')
        ff.write(txt)
        ff.close()

    def write_search(self,search_fname='_test.search',
                     frin_fname='_test.in',
                     frout_fname = '_test.out'):
        # use self.search_var information....
        num_var= len(self.search_var)
        num_exp=1
        txt=""
        txt=txt +"\'{}\' \'{}\' {} {} \n".format(frin_fname, frout_fname, num_var, num_exp)
        for var in self.search_var:
            txt= txt+'&VARIABLE '
            for item in var:
                if item in self.str_variables:
                    txt= txt +' {}=\'{}\' '.format(item, var[item] )
                else:
                    txt= txt +' {}={} '.format(item, var[item] )
            txt = txt +'/\n'
        txt=txt+self.write_data()
        ff=open(search_fname,'w')
        ff.write(txt)
        ff.close()
        return txt

    def add_search_variable(self,**attributes):
        """
        add search variables to SFRESCO with keyword arguments
        """
        var_input ={}
        for item,value in attributes.items():
            if item in self.Var_keywords:
                var_input[item]=value
        self.search_var.append(var_input)
        
    def reset_search_variables(self,):
        self.search_var =[] 

#=============================================================================
def get_fresco_result():
    """ read fort.16 file results and return dictionary
    with arrays
    """
    clean_comm('fort.16') #remove comments
    out = read_fresco_res('fort.16x')
    for i in out:
        out[i] = np.array(out[i])
    return out

def chck_fresco_out(fname=''):
  """
    test whether the FRESCO ended normally
    return 0 if okay
             1 if error

    Originally checked
      'Total CPU '
      and
      'Recommended RNL: non-local width' as 'OK'

    However this seems to be not correct in Windows...
  """
  ff=open(fname,'r')
  ll=ff.readlines()
  ff.close()
  chck=0
  for i in ll[:-20:-1]: #last 9 lines
    if 'Total CPU ' in i:
       chck=2
    if 'ACCURACY ANALYSIS at ' in i:
       chck=2

  if chck==0 :
    print('ERROR in %s file'%fname)
    return 1
  elif chck==1 :
    print('ERROR in %s file'%fname)
    return 1
  elif chck==2 :
    print('ok')
    return 0

def get_elastic_result_from_fresco_out(fname='_output.out'):
    """
    Read Fresco output file and read elastic channel X-S.
    for both mb/sr and ratio values.

    return array = [angle, mb/sr, ratio]
 
    In case of chargeless neutron, there is no ratio! 
    needs a special treatment 
    """
    ff=open(fname,'r')
    ll=ff.readlines()
    ff.close()
    elastic=[]
    for i in range(len(ll)):
        if ' CROSS SECTIONS FOR OUTGOING ' in ll[i]: # found elastic x-s
            i = i+1 # one line skip
            while True:
                i = i+1
                if (' Integrated' in ll[i]
                    or 'CROSS SECTIONS FOR OUTGOING' in ll[i]
                    or 'Finished all xsecs' in ll[i]): # end of elastic x-s
                    #print('end of elastic')
                    break
                if 'deg.: X-S' in ll[i]:
                    ww = ll[i].split()
                    theta, mb = float(ww[0]),float(ww[4])
                    try: # usual case , ratio is printed 
                        ww = ll[i+1].split()
                        ratio = float(ww[3])
                        elastic.append([theta,mb,ratio])
                        i=i+1
                    except: # neutral particle case 
                        elastic.append([theta,mb,0])
                    #print(theta,mb,ratio)
            break #end of search 
    return np.array(elastic)

#------------------------------------------------------------------------------
def run_fresco_from_input_txt(fresco_input_txt,
                              fresco_path='fresco.exe',
                              fresco_input_path='_test.in',
                              fresco_output_path='_test.out'):
    """
    run fresco by using the fresco input text. 
    """
    #----need to delete fort before calculation
    import glob 
    fort_files = glob.glob('fort.*')    
    for i in fort_files:
        if not(i == 'fort.4'):
            os.remove(i) 
    # remove previous input/ouputs
    try:
        os.remove(fresco_input_path)    
        os.remove(fresco_output_path)
    except:
        print('Error removing tempoary files')
    
    #----prepare input file and fort.4 file 
    ff= open(fresco_input_path,'w')
    ff.write(fresco_input_txt)
    ff.close() 
    
    proc = Popen(fresco_path +" < "+fresco_input_path
                 +" > "+fresco_output_path,shell=True)
    proc.wait() 
    proc.terminate() 
             
    if chck_fresco_out(fname=fresco_output_path)==0:
        out = get_fresco_result()
        return out 
    else: 
        print('There is an Error!!. Check input and output')
        return 
    
            

    
    