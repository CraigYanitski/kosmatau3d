import Shape
import VoxelGrid
import Orientation
import Observations
class Model():
  '''
  This is the highest class in the hierarchy of the KOSMA-tau^3 simulation.
  It contains all of the information needed to properly model a PDR (I think).
  '''
  # PRIVATE
  def __init__(self, x, y, z, modelType=''):
    self.__type = modelType
    self.__shape = Shape(x, y, z, modelType=modelType)
    self.__grid = VoxelGrid(self.__shape.getDimensions())
    self.__orientation = Orientation()
    self.__observation = Observations
    self.__species = []
    self.__molecules = []   #list of molecules to include in model
    self.__moleculeNames = []
    self.__moleculeNumber = []
    self.__moleculeTransitions = {}
    self.__moleculeFrequencies = {}
    self.__dust = []        #list of dust to include in model
    self.__dustNames = []
    self.__dustNumber = []
    self.__dustTransitions = {}
    self.__dustFrequencies = {}
    return

  # PUBLIC
  def initialiseModel(self):
    self.__species.append(self.__molecules)
    self.__species.append(self.__dust)
    self.__grid.initialiseVoxels(self.__species, self.__observations)
    return
  def addDust(self, transition):
    (numbers,species,transitions,frequencies) = self.__observations.__speciesData
    if len(self.__dust)==0:
      self.__dust.append(Dust(transition))
      self.__dustNames.append('dust')
      self.__dustNumber.append(numbers[species=='dust' and transitions==transition])
      self.__dustTransitions['dust'].append(transition)
      self.__dustFrequencies['dust'].append(frequencies[species=='dust' and transitions==transition])
    else:
      self.__dust[0].addTransition(transition)
      self.__dustNumber.append(numbers[species=='dust' and transitions==transition])
      self.__dustTransitions['dust'].append(transition)
      self.__dustFrequencies['dust'].append(frequencies[species=='dust' and transitions==transition])
    return
  def addMolecule(self, molecule, transition):
    (numbers,species,transitions,frequencies) = self.__observations.__speciesData
    if len(self.__dust)==0:
      self.__molecules.append(Dust(transition))
      self.__moleculeNames.append(molecule)
      self.__moleculeNumber.append(numbers[species==molecule and transitions==transition])
      self.__moleculeTransitions[molecule].append(transition)
      self.__moleculeFrequencies[molecule].append(frequencies[species=='dust' and transitions==transition])
    else:
      self.__molecules[self.__molecules==molecule].addTransition(transition)
      self.__moleculeNumber.append(numbers[species==molecule and transitions==transition])
      self.__moleculeTransitions[molecule].append(transition)
      self.__moleculeFrequencies[molecule].append(frequencies[species==molecule and transitions==transition])
    return
  def addSpecies(self, species, transition):
    # find transition number as defined in
    # SetUpOrionBarModelEnvironment.v.1.1.nb
    #from PDR import _globals  
    (number,species,transition,frequency) = self.__observations.__speciesData
    #gbl._globals['compound']['number'] = []
    #gbl._globals['compound']['frequency'] = []
    for element in species:
      if element=='Dust':
        self.addDust(transition)
        return
      else:
        self.addMolecule(species, transition)
        return







    for i in range(gbl._globals['compound']['nspe']):
      # for each species 
      number = -1
      species = gbl._globals['compound']['species'][i]
      transition = gbl._globals['compound']['transition'][i]
      if species == 'C+': number = 1
      elif species == 'C':  
        if transition > 3 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 1 + transition
      elif species == 'O':
        if transition > 3 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 4 + transition
      elif species == 'CO': 
        if transition > 49 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 7 + transition
      elif species == '13CO': 
        if transition > 49 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 56 + transition
      elif species == '13C+': number = 106
      elif species == '13C':
        if transition > 3 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 106 + transition
      elif species == 'HCO+':
        if transition > 15 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 109 + transition
      elif species == 'H13CO+':
        if transition > 30 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 124 + transition
      elif species == 'H3O+':
        if transition > 17 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 154 + transition
      elif species == 'C18O':
        if transition > 49 or transition < 1:
          import sys    
          sys.exit('transition does not exist...exiting...')
        else: number = 171 + transition
      elif species == 'Dust':
        if transition > 51 or transition < 1:
          import sys
          sys.exit('PDR_functions:this Dust transition is not implemented...exiting...')
        else: number = 191 + transition
  #    elif species == 'CH+': 
  #        if transition > 10 or transition < 1:
  #            import sys    
  #            sys.exit('transition does not exist...exiting...')
  #        else: number = 171 + transition
      if number == -1:
        print('species: ', species)
        print('transition: ', transition)
        import sys    
        sys.exit('species/transition combination does not exist...exiting...')            
      gbl._globals['compound']['number'].append(number) 
      with open(gbl.INPUTPATH+'frequencies.dat','r') as f:
        for i, line in enumerate(f):
          if i == number:
            try: freq = float(line.split()[3]) 
            # read frequency from frequency.dat    
            except IndexError:
              import sys
              sys.exit('Problem :( frequency information needs to be added to file frequencies.dat ...exiting...')          
            gbl._globals['compound']['frequency'].append(freq)
            #print 'hier'
            
          #elif i > number:
            #print 'break'
            #break
      #gbl._globals['compound']['frequency'].append(freq)
    if len(gbl._globals['compound']['transition']) != len( gbl._globals['compound']['frequency']):
      print('frequencies found ' , len( gbl._globals['compound']['frequency']))
      print('frequencies seeked ', len(gbl._globals['compound']['transition']))
      import sys
      sys.exit('Problem :( Frequency information does not match. Maybe needs to be added to file frequencies.dat ...exiting...')
    
      # print 'number', gbl._globals['compound']['number']
      # print 'frequency',  gbl._globals['compound']['frequency'] 
      # pause = input('..ok?..remove again..in PDRfunctions..')
      return gbl._globals 
        #gbl._globals['compound']['frequency'].append(freq)
    # print 'number', gbl._globals['compound']['number']
    # print 'frequency',  gbl._globals['compound']['frequency'] 
    # pause = input('..ok?..remove again..in PDRfunctions..')
    return gbl._globals