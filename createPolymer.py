import argparse as a
from time import sleep
import decimal
from itertools import combinations
import sys

def cutString (inputString, subString1, subString2):
	startIndex = inputString.find (subString1) + len (subString1)
	endIndex = inputString.find (subString2)

	return inputString [startIndex:endIndex]

def readConfig (configFile):
	# Reads the config file to get user defined input parameters
	with open (configFile, "r") as file:
		fileString = file.read ()

	fileLines = fileString.split ("\n")

	for line in fileLines:
		if ("translateAxis" in line):
			translateAxis = line.split (":")
		if ("endGroup1" in line):
			endGroup1 = line.split (":")
		if ("endGroup2" in line):
			endGroup2 = line.split (":")
		if ("repeatMonomers" in line):
			repeatMonomers = line.split (":")
		if ("translateDistance" in line):
			translateDistance = line.split (":")
		if ("mainChain" in line):
			mainChain = line.split (":")
		if ("pendantGroup" in line):
			pendantGroup = line.split (":")
		if ("improper" in line):
			improper = line.split (":")

	return translateAxis [1], int (endGroup1 [1].replace ("a", "")), int (endGroup2 [1].replace ("a", "")), int (repeatMonomers [1]), float (translateDistance [1]), mainChain [1], pendantGroup [1], int (improper [1])

def readCharges (chargesFile):
	chargesDict = {}
	with open (chargesFile, "r") as file:
		fileString = file.read ()

	fileString2 = fileString.split ("\n")

	for line in fileString2:
		fileString3 = line.replace (" ", "").replace ("a", "")
		fileString4 = fileString3.split (":")
		if (len (fileString4) == 2):
			chargesDict [int (fileString4[0])] = fileString4[1]

	return chargesDict

def createPolymer (cmlFile, atomTypesFile, configFile, chargesFile):
	# Inputs monomer atomic coordinates from 'cmlFile'
	# Replicates the monomer based on information given in 'configFile'
	# Bond information are also read from 'cmlFile'
	translateAxis, endGroup1, endGroup2, repeatMonomers, translateDistance, mainChain, pendantGroup, improper = readConfig (configFile)

	chargesDict = readCharges (chargesFile)

	isXTranslate = 0
	isYTranslate = 0
	isZTranslate = 0
	if ('x' in translateAxis):
		isXTranslate = 1
	if ('y' in translateAxis):
		isYTranslate = 1
	if ('z' in translateAxis):
		isZTranslate = 1

	atomTypes = {}
	with open (atomTypesFile, "r") as file:
		for line in file:
			lineArray = line.split (" ")
			atomTypes [lineArray [0]] = int (lineArray [1].replace ("\n", ""))
			
	with open (cmlFile, "r") as file:
		cmlInput = file.read ()

	atomEntries = []
	bondEntries = []
	polymerEntries = []

	# String processing is done to extract the necessary information from input 'cmlFile'
	atomInput = cutString (cmlInput, "<atomArray>", "</atomArray>")
	atomArray = atomInput.split ("\n")
	nAtomsPerMonomer = 0

	# Extracting atomic coordinates
	for atoms in atomArray:
		if (len (atoms) > 10):
			currentID = int (cutString (atoms, "<atom id=\"", "\" elementType=\"").replace ("a", ""))
			currentElementType = cutString (atoms, "\" elementType=\"", "\" x3=\"")
			currentX = float (cutString (atoms, "\" x3=\"", "\" y3=\""))
			currentY = float (cutString (atoms, "\" y3=\"", "\" z3=\""))
			currentZ = float (cutString (atoms, "\" z3=\"", "\"/>"))
			try:
				atomEntries.append ({'id': currentID, 'elementType': currentElementType, 'x': currentX, 'y': currentY, 'z': currentZ, 'charge': chargesDict [currentID], 'bondAtom1': 0, 'bondAtom2': 0, 'bondAtom3': 0, 'bondAtom4': 0})
			except:
				print ("Check all the input files. Cannot add atom entries. Code stopped at line: 99")
				sys.exit (1)
			nAtomsPerMonomer = int (cutString (atoms, "<atom id=\"", "\" elementType=\"").replace ("a", ""))

	bondInput = cutString (cmlInput, "<bondArray>", "</bondArray>")
	bondArray = bondInput.split ("\n")
	
	# Extracting bond information
	for bonds in bondArray:
		if (len (bonds) > 10):
			bondAtoms = cutString (bonds, "atomRefs2=\"", "\" order=\"").replace ("a", "")
			atom1, atom2 = bondAtoms.split (" ")
			bondEntries.append ({'atom1': atom1, 'atom2': atom2, 'bondType': cutString (bonds, "\" order=\"", "\"/>")})
			if (atomEntries [int (atom1) - 1]['bondAtom1'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom1'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom2'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom2'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom3'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom3'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom4'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom4'] = int (atom2)

			if (atomEntries [int (atom2) - 1]['bondAtom1'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom1'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom2'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom2'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom3'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom3'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom4'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom4'] = int (atom1)

	# Replicating the initial monomeric structure
	endGroupArr = []
	for x in range (0, repeatMonomers, 1):
		for atoms in atomEntries:
			# atoms ['bondAtom1'] are edited if they contain non-zero value
			if (atoms ['bondAtom1'] > 0):
				bondAtom1_polymerEntries = atoms ['bondAtom1'] + (x * nAtomsPerMonomer)
			else:
				bondAtom1_polymerEntries = 0
			if (atoms ['bondAtom2'] > 0):
				bondAtom2_polymerEntries = atoms ['bondAtom2'] + (x * nAtomsPerMonomer)
			else:
				bondAtom2_polymerEntries = 0
			if (atoms ['bondAtom3'] > 0):
				bondAtom3_polymerEntries = atoms ['bondAtom3'] + (x * nAtomsPerMonomer)
			else:
				bondAtom3_polymerEntries = 0
			if (atoms ['bondAtom4'] > 0):
				bondAtom4_polymerEntries = atoms ['bondAtom4'] + (x * nAtomsPerMonomer)
			else:
				bondAtom4_polymerEntries = 0

			# Each monomer has endGroup1 and endGroup2
			# Bonds are defined between these endGroups while replicating the monomer
			if (x == 0 and atoms ['id'] == endGroup2):
				bondAtom4_polymerEntries = endGroup1 + ((x + 1) * nAtomsPerMonomer)
			elif (x == (repeatMonomers - 1) and atoms ['id'] == endGroup1):
				bondAtom4_polymerEntries = endGroup2 + ((x - 1) * nAtomsPerMonomer)
			elif (x > 0 and x < (repeatMonomers - 1)):
				if (atoms ['id'] == endGroup1):
					bondAtom4_polymerEntries = endGroup2 + ((x - 1) * nAtomsPerMonomer)
				if (atoms ['id'] == endGroup2):
					bondAtom4_polymerEntries = endGroup1 + ((x + 1) * nAtomsPerMonomer)

			# Above info is added to polymerEntries
			# polymerEntries array is the return type
			polymerEntries.append ({'sino': atoms ['id'] + (x * nAtomsPerMonomer), 'monomerAtomID': atoms ['id'], 'elementType': atoms ['elementType'], 'atomType': atomTypes [atoms ['elementType']], 'molType': 1, 'charge': atoms ['charge'], 'x': round (atoms ['x'] + (isXTranslate * x * translateDistance), 4), 'y': round (atoms ['y'] + (isYTranslate * x * translateDistance), 4), 'z': round (atoms ['z'] + (isZTranslate * x * translateDistance), 4), 'bondAtom1': bondAtom1_polymerEntries, 'bondAtom2': bondAtom2_polymerEntries, 'bondAtom3': bondAtom3_polymerEntries, 'bondAtom4': bondAtom4_polymerEntries})

	return polymerEntries, atomEntries, atomTypes

def createPolymer2 (cmlFile, atomTypesFile, configFile):
	# 
	# ~~~~~~~~~~~~~~~~~
	# OUTDATED FUNTION
	# ~~~~~~~~~~~~~~~~~
	# This function is no longer used currently
	# 
	translateAxis, endGroup1, endGroup2, repeatMonomers, translateDistance = readConfig (configFile)

	isXTranslate = 0
	isYTranslate = 0
	isZTranslate = 0
	if ('x' in translateAxis):
		isXTranslate = 1
	if ('y' in translateAxis):
		isYTranslate = 1
	if ('z' in translateAxis):
		isZTranslate = 1

	atomTypes = {}
	with open (atomTypesFile, "r") as file:
		for line in file:
			lineArray = line.split (" ")
			atomTypes [lineArray [0]] = int (lineArray [1].replace ("\n", ""))
			
	with open (cmlFile, "r") as file:
		cmlInput = file.read ()

	atomEntries = []
	bondEntries = []
	polymerEntries = []
	atomInput = cutString (cmlInput, "<atomArray>", "</atomArray>")
	atomArray = atomInput.split ("\n")
	nAtomsPerMonomer = 0
	for atoms in atomArray:
		if (len (atoms) > 10):
			atomEntries.append ({'id': int (cutString (atoms, "<atom id=\"", "\" elementType=\"").replace ("a", "")), 'elementType': cutString (atoms, "\" elementType=\"", "\" x3=\""), 'x': float (cutString (atoms, "\" x3=\"", "\" y3=\"")), 'y': float (cutString (atoms, "\" y3=\"", "\" z3=\"")), 'z': float (cutString (atoms, "\" z3=\"", "\"/>")), 'bondAtom1': 0, 'bondAtom2': 0, 'bondAtom3': 0, 'bondAtom4': 0})
			nAtomsPerMonomer = int (cutString (atoms, "<atom id=\"", "\" elementType=\"").replace ("a", ""))

	bondInput = cutString (cmlInput, "<bondArray>", "</bondArray>")
	bondArray = bondInput.split ("\n")
	for bonds in bondArray:
		if (len (bonds) > 10):
			bondAtoms = cutString (bonds, "atomRefs2=\"", "\" order=\"").replace ("a", "")
			atom1, atom2 = bondAtoms.split (" ")
			bondEntries.append ({'atom1': atom1, 'atom2': atom2, 'bondType': cutString (bonds, "\" order=\"", "\"/>")})
			if (atomEntries [int (atom1) - 1]['bondAtom1'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom1'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom2'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom2'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom3'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom3'] = int (atom2)
			elif (atomEntries [int (atom1) - 1]['bondAtom4'] == 0):
				atomEntries [int (atom1) - 1]['bondAtom4'] = int (atom2)

			if (atomEntries [int (atom2) - 1]['bondAtom1'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom1'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom2'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom2'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom3'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom3'] = int (atom1)
			elif (atomEntries [int (atom2) - 1]['bondAtom4'] == 0):
				atomEntries [int (atom2) - 1]['bondAtom4'] = int (atom1)

	for x in range (0, repeatMonomers, 1):
		for atoms in atomEntries:
			polymerEntries.append ({'sino': atoms ['id'] + (x * nAtomsPerMonomer), 'elementType': atoms ['elementType'], 'atomType': atomTypes [atoms ['elementType']], 'molType': 1, 'charge': 0, 'x': round (atoms ['x'] + (isXTranslate * x * translateDistance), 4), 'y': round (atoms ['y'] + (isYTranslate * x * translateDistance), 4), 'z': round (atoms ['z'] + (isZTranslate * x * translateDistance), 4), 'bondAtom1': atoms ['bondAtom1'] + (x * nAtomsPerMonomer), 'bondAtom2': atoms ['bondAtom2'] + (x * nAtomsPerMonomer), 'bondAtom3': atoms ['bondAtom3'] + (x * nAtomsPerMonomer), 'bondAtom4': atoms ['bondAtom4'] + (x * nAtomsPerMonomer)})

	return polymerEntries, atomTypes

def extract_numbers (line):
	# Extracts numbers from a line of type 'string'
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def readAtomInfo (inputFileName, atomInfo):
	# 
	# ~~~~~~~~~~~~~~~~~
	# OUTDATED FUNTION
	# ~~~~~~~~~~~~~~~~~
	# This function is no longer used currently
	# 
	atomTypeArr = []
	with open ("atomEntries.testing", "r") as inputFile:
		for line in inputFile:
			lineArray = list (extract_numbers (line))
			atomInfo.append ({'sino': int (lineArray[0]), 'molType': int (lineArray[1]), 'atomType': int (lineArray[2]), 'charge': int (lineArray[3]), 'x': float (lineArray[4]), 'y': float (lineArray[5]), 'z': float (lineArray[6]), 'bondAtom1': int (lineArray[7]), 'bondAtom2': int (lineArray[8]), 'bondAtom3': int (lineArray[9]), 'bondAtom4': int (lineArray[10])})

	for atomLine in atomInfo:
		if (atomLine ['atomType'] not in atomTypeArr):
			atomTypeArr.append (atomLine ['atomType'])

	return atomInfo, atomTypeArr

def createBonds (polymerEntries, bondInfo):
	# Based on the list of connected atoms, this function
	# creates bond information in LAMMPS format
	lineArray = []
	sino_bond = 1

	# bondTypeArr array is used to determine the bondType
	bondTypeArr = []
	bondTypeArr.append({'atom1': 0, 'atom2': 0, 'bondType': 0})

	# Check if the bond is already present
	# This function returns 0 is the bond is already present
	# and returns 1 if otherwise
	def bondCheck (atom1, atom2, bondInfo):
		for bond in bondInfo:
			if ((atom1 == bond['bondAtom1'] or atom1 == bond['bondAtom2']) and (atom2 == bond['bondAtom1'] or atom2 == bond['bondAtom2'])):
				return 0

		return 1

	# Determine the bondType
	# This function takes two atoms of a bond and checks for old records
	# If record is not present, then a new bondType is assigned
	def findBondType (atom1, atom2, bondTypeArr):
		returnBondType = 0

		# Atoms are arranged in ascending order
		# To make sure atom1---atom2 is the same as atom2---atom1
		if (atom1 < atom2):
			ascAtom1 = atom1
			ascAtom2 = atom2
		else:
			ascAtom1 = atom2
			ascAtom2 = atom1

		for items in bondTypeArr:
			if (ascAtom1 == items ['atom1'] and ascAtom2 == items ['atom2']):
				returnBondType = items ['bondType']

		if (returnBondType == 0):
			existingLength = len (bondTypeArr)
			bondTypeArr.append ({'atom1': ascAtom1, 'atom2': ascAtom2, 'bondType': existingLength})
			returnBondType = existingLength

		return returnBondType, bondTypeArr

	# Add bonds to bondInfo array
	def addBond (atomLine, polymerEntries, dictString, bondInfo, sino_bond, bondTypeArr):
		if (atomLine [dictString] and bondCheck (atomLine ['sino'], atomLine [dictString], bondInfo)):
			try:
				secondAtomType = polymerEntries [int (atomLine [dictString]) - 1]['atomType']
				bondType, bondTypeArr = findBondType (atomLine ['atomType'], secondAtomType, bondTypeArr)
				if (atomLine ['sino'] < atomLine [dictString]):
					bondInfo.append ({'sino': sino_bond, 'bondType': bondType, 'bondAtom1': atomLine ['sino'], 'bondAtom2': atomLine [dictString], 'bondAtom1Type': atomLine ['atomType'], 'bondAtom2Type': secondAtomType})
					sino_bond += 1
				else:
					bondInfo.append ({'sino': sino_bond, 'bondType': bondType, 'bondAtom1': atomLine [dictString], 'bondAtom2': atomLine ['sino'], 'bondAtom1Type': secondAtomType, 'bondAtom2Type': atomLine ['atomType']})
					sino_bond += 1
			except:
				pass

		return bondInfo, bondTypeArr, sino_bond

	# Iterate through the polymer to add bonds
	# Bond info are stored in 'bondInfo' array
	# Bond types are stored in 'bondTypeArr'
	for atomLine in polymerEntries:
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom1', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom2', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom3', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom4', bondInfo, sino_bond, bondTypeArr)

	return bondInfo, bondTypeArr

def createAngles (polymerEntries, angleInfo, bondInfo):
	# angleInfo contains sino, angleType, angleAtom1, angleAtom2, angleAtom3, angleAtom1Type, angleAtom2Type, angleAtom3Type
	# This function uses bondInfo to define angleInfo
	# angleInfo array is used to create LAMMPS angle entries
	angleTypeArr = []
	sino_angle = 1

	# Initializing connectedAtoms dict of lists
	connectedAtoms = {}
	for atomLine in polymerEntries:
		connectedAtoms [atomLine ['sino']] = []

	angleTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'angleType': 0})

	def findConnectedAtoms (concernedBondAtom, primaryConnect, bondInfo, connectedAtoms):
		# concernedBondAtom = bondAtom1 from bondInfo
		# primaryConnect = bondAtom2 from bondInfo

		# Adding concernedBondAtom and primaryConnect to their respective 'connectedAtoms[]' dict of lists
		if (primaryConnect not in connectedAtoms [concernedBondAtom]):
			connectedAtoms [concernedBondAtom].append (primaryConnect)
		if (concernedBondAtom not in connectedAtoms [primaryConnect]):
			connectedAtoms [primaryConnect].append (concernedBondAtom)

		# Iterating through bondInfo list
		# Checking for all the atoms connected to concernedBondAtom and primaryConnect
		# All the connected atoms are stored in 'connectedAtoms[]' dict of lists
		for bondLine in bondInfo:
			if ((bondLine ['bondAtom1'] == concernedBondAtom) and (bondLine ['bondAtom2'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom2'])
				connectedAtoms [concernedBondAtom].sort ()
			if ((bondLine ['bondAtom2'] == concernedBondAtom) and (bondLine ['bondAtom1'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom1'])
				connectedAtoms [concernedBondAtom].sort ()

		return connectedAtoms

	# Iterates through bondInfo array
	# All connected atoms are used to create 'angle' list
	for bondLine in bondInfo:
		connectedAtoms = findConnectedAtoms (bondLine ['bondAtom1'], bondLine ['bondAtom2'], bondInfo, connectedAtoms)
	
	def findAngleType (firstAtomType, secondAtomType, thirdAtomType, angleTypeArr):
		if (firstAtomType <= thirdAtomType):
			ascFirstAtomType = firstAtomType
			ascThirdAtomType = thirdAtomType
		else:
			ascFirstAtomType = thirdAtomType
			ascThirdAtomType = firstAtomType

		assignedType = 0

		for types in angleTypeArr:
			if (ascFirstAtomType == types ['atom1'] and secondAtomType == types ['atom2'] and ascThirdAtomType == types ['atom3']):
				assignedType = types ['angleType']

		if (assignedType == 0):
			assignedType = len (angleTypeArr)
			angleTypeArr.append ({'atom1': ascFirstAtomType, 'atom2': secondAtomType, 'atom3': ascThirdAtomType, 'angleType': len (angleTypeArr)})

		return assignedType, angleTypeArr

	# Creating angleInfo
	for atom in connectedAtoms:
		if (len (connectedAtoms [atom]) > 1):
			comb = combinations (connectedAtoms [atom], 2)
			for i in list (comb):
				firstAtomType = polymerEntries [int (i [0]) - 1]['atomType']
				secondAtomType = polymerEntries [int (atom) - 1]['atomType']
				thirdAtomType = polymerEntries [int (i [1]) - 1]['atomType']
				angleType, angleTypeArr = findAngleType (firstAtomType, secondAtomType, thirdAtomType, angleTypeArr)
				angleInfo.append ({'sino': sino_angle, 'angleType': angleType, 'angleAtom1': i [0], 'angleAtom2': atom, 'angleAtom3': i [1], 'angleAtom1Type': firstAtomType, 'angleAtom2Type': secondAtomType, 'angleAtom3Type': thirdAtomType})
				sino_angle += 1

	return angleInfo, angleTypeArr

def createDihedrals (polymerEntries, dihedralInfo, angleInfo, bondInfo):
	# dihInfo contains sino, dihType, dihAtom1, dihAtom2, dihAtom3, dihAtom4, dihAtom1Type, dihAtom2Type, dihAtom3Type, dihAtom4Type
	dihTypeArr = []
	sino_dih = 1

	dihTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'atom4': 0, 'dihType': 0})

	def findDihType (firstAtomType, secondAtomType, thirdAtomType, fourthAtomType, dihTypeArr):
		if (firstAtomType < fourthAtomType):
			ascFirstAtomType = firstAtomType
			ascSecondAtomType = secondAtomType
			ascThirdAtomType = thirdAtomType
			ascFourthAtomType = fourthAtomType
		else:
			ascFirstAtomType = fourthAtomType
			ascSecondAtomType = thirdAtomType
			ascThirdAtomType = secondAtomType
			ascFourthAtomType = firstAtomType

		assignedType = 0

		for types in dihTypeArr:
			if (ascFirstAtomType == types ['atom1'] and ascSecondAtomType == types ['atom2'] and ascThirdAtomType == types ['atom3'] and ascFourthAtomType == types ['atom4']):
				assignedType = types ['dihType']

		if (assignedType == 0):
			assignedType = len (dihTypeArr)
			dihTypeArr.append ({'atom1': ascFirstAtomType, 'atom2': ascSecondAtomType, 'atom3': ascThirdAtomType, 'atom4': ascFourthAtomType, 'dihType': len (dihTypeArr)})

		return assignedType, dihTypeArr

	for x in range (0, len (angleInfo), 1):
		for y in range (x + 1, len (angleInfo), 1):
			if (angleInfo [x]['angleAtom2'] == angleInfo [y]['angleAtom1'] and angleInfo [x]['angleAtom3'] == angleInfo [y]['angleAtom2']):
				firstAtomType = polymerEntries [int (angleInfo [x]['angleAtom1']) - 1]['atomType']
				secondAtomType = polymerEntries [int (angleInfo [x]['angleAtom2']) - 1]['atomType']
				thirdAtomType = polymerEntries [int (angleInfo [x]['angleAtom3']) - 1]['atomType']
				fourthAtomType = polymerEntries [int (angleInfo [y]['angleAtom3']) - 1]['atomType']
				dihType, dihTypeArr = findDihType (firstAtomType, secondAtomType, thirdAtomType, fourthAtomType, dihTypeArr)
				dihedralInfo.append ({'sino': sino_dih, 'dihType': dihType, 'dihAtom1': angleInfo [x]['angleAtom1'], 'dihAtom2': angleInfo [x]['angleAtom2'], 'dihAtom3': angleInfo [x]['angleAtom3'], 'dihAtom4': angleInfo [y]['angleAtom3'], 'dihAtom1Type': firstAtomType, 'dihAtom2Type': secondAtomType, 'dihAtom3Type': thirdAtomType, 'dihAtom4Type': fourthAtomType})
				sino_dih += 1

	return dihedralInfo, dihTypeArr

def createImpropers (polymerEntries, atomEntries, improperInfo, configFile):
	# sino, impType, impAtom1, impAtom2, impAtom3, impAtom4 (impAtom2 is the central atom from main chain)
	# The following code assumes there is only one type of improper angle
	# Because the most common use for an improper angle is to maintain tacticity of the molecule

	# Temp dict to store atom IDs (sino)
	improperEntries = {}

	translateAxis, endGroup1, endGroup2, repeatMonomers, translateDistance, mainChain, pendantGroup, improper = readConfig (configFile)

	mainChain = mainChain.replace ("a", "").replace (" ", "")
	mainChainAtoms = mainChain.split (",")
	mainChainAtoms = [int(i) for i in mainChainAtoms]
	pendantGroupAtoms = pendantGroup.replace ("a", "").replace (" ", "")

	# Find the atom (from main chain) connected to the pendant atom
	mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
	mainChainAtomsArray = []
	pendantAtomsArray = []

	# Classifying the atoms into two groups, (i) mainChainAtoms, and (ii) pendantGroupAtoms
	for atoms in polymerEntries:
		if (str (atoms ['monomerAtomID']) in mainChainAtoms):
			mainChainAtomsArray.append (atoms)
		if (str (atoms ['monomerAtomID']) in pendantGroupAtoms):
			pendantAtomsArray.append (atoms)

	# Iterate through the pendant atoms
	# Find the main chain atom directly connected to the pendant atoms
	# Then find the two main chain atoms connected to the original main chain atom
	for atoms in pendantAtomsArray:
		bondAtom1 = atoms ['bondAtom1']
		bondAtom2 = atoms ['bondAtom2']
		bondAtom3 = atoms ['bondAtom3']
		bondAtom4 = atoms ['bondAtom4']
		currentSino = atoms['sino']

		# improperEntries contain the final information to be returned
		improperEntries [currentSino] = []

		# Checking if the bondAtom1 is not zero and if the bondAtom1 is present in main chain (as specified in the *.config file)
		# This 'if' statement gives the main chain atom that is connected to pendant group (impAtom2)
		if (bondAtom1 and polymerEntries [bondAtom1 - 1]['monomerAtomID'] in mainChainAtoms):
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []

			# All atoms connected to that specific main chain atom (impAtom2)
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [bondAtom1 - 1]['bondAtom1'], polymerEntries [bondAtom1 - 1]['bondAtom2'], polymerEntries [bondAtom1 - 1]['bondAtom3'], polymerEntries [bondAtom1 - 1]['bondAtom4']]

			# Iterating through the atoms connected to impAtom2
			# sino of all connected atoms are added to improperEntries dict
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					improperEntries [currentSino].append (polymerEntries [bondAtom1 - 1]['sino'])
					improperEntries [currentSino].append (polymerEntries [atoms1 - 1]['sino'])

		# Checking if the bondAtom1 is not zero and if the bondAtom1 is present in main chain (as specified in the *.config file)
		# This 'if' statement gives the main chain atom that is connected to pendant group (impAtom2)
		if (bondAtom2 and polymerEntries [bondAtom2 - 1]['monomerAtomID'] in mainChainAtoms):
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []

			# All atoms connected to that specific main chain atom (impAtom2)
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [bondAtom2 - 1]['bondAtom1'], polymerEntries [bondAtom2 - 1]['bondAtom2'], polymerEntries [bondAtom2 - 1]['bondAtom3'], polymerEntries [bondAtom2 - 1]['bondAtom4']]

			# Iterating through the atoms connected to impAtom2
			# sino of all connected atoms are added to improperEntries dict
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					improperEntries [currentSino].append (polymerEntries [bondAtom2 - 1]['sino'])
					improperEntries [currentSino].append (polymerEntries [atoms1 - 1]['sino'])

		# Checking if the bondAtom1 is not zero and if the bondAtom1 is present in main chain (as specified in the *.config file)
		# This 'if' statement gives the main chain atom that is connected to pendant group (impAtom2)
		if (bondAtom3 and polymerEntries [bondAtom3 - 1]['monomerAtomID'] in mainChainAtoms):
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []

			# All atoms connected to that specific main chain atom (impAtom2)
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [bondAtom3 - 1]['bondAtom1'], polymerEntries [bondAtom3 - 1]['bondAtom2'], polymerEntries [bondAtom3 - 1]['bondAtom3'], polymerEntries [bondAtom3 - 1]['bondAtom4']]

			# Iterating through the atoms connected to impAtom2
			# sino of all connected atoms are added to improperEntries dict
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					improperEntries [currentSino].append (polymerEntries [bondAtom3 - 1]['sino'])
					improperEntries [currentSino].append (polymerEntries [atoms1 - 1]['sino'])

		# Checking if the bondAtom1 is not zero and if the bondAtom1 is present in main chain (as specified in the *.config file)
		# This 'if' statement gives the main chain atom that is connected to pendant group (impAtom2)
		if (bondAtom4 and polymerEntries [bondAtom4 - 1]['monomerAtomID'] in mainChainAtoms):
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []

			# All atoms connected to that specific main chain atom (impAtom2)
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [bondAtom4 - 1]['bondAtom1'], polymerEntries [bondAtom4 - 1]['bondAtom2'], polymerEntries [bondAtom4 - 1]['bondAtom3'], polymerEntries [bondAtom4 - 1]['bondAtom4']]

			# Iterating through the atoms connected to impAtom2
			# sino of all connected atoms are added to improperEntries dict
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					improperEntries [currentSino].append (polymerEntries [bondAtom4 - 1]['sino'])
					improperEntries [currentSino].append (polymerEntries [atoms1 - 1]['sino'])

	# Internal 'sino' count for improperInfo array
	impCount = 0

	# Iterating through the improperEntries dict (temporary variable)
	for entries in improperEntries:
		impAtom1 = 0
		impAtom2 = 0
		impAtom3 = 0
		impAtom4 = 0
		impCount += 1

		# impAtom1 is set as the pendant atom by default
		impAtom1 = entries
		for entries1 in improperEntries [entries]:
			# impAtom2 is the center atom
			# It is the main chain atom directly connected to the pendant atom
			# In the dict of array, it is repeated twice (extracted from bond connection information)
			if (improperEntries [entries].count (entries1) == 2):
				impAtom2 = entries1
			# Remaining atoms are added to impAtom3 and impAtom4
			elif (impAtom3 == 0):
				impAtom3 = entries1
			elif (impAtom4 == 0):
				impAtom4 = entries1

		# The length is less than 4 for end groups
		# No improper is set for end group monomers
		# So they are atactic in the setup.
		# Tacticity of remaining segments are set
		if (len (improperEntries [entries]) == 4):
			# improperInfo array is returned
			# Only one type of improper is considered
			# It is assumed that improper is set only to control tacticity
			# and the polymer contains only one pendant group (simple polymer structures)
			# For complex polymers with multiple pendant groups (or groups within groups)
			# this code must be modified and several impTypes must be considered
			improperInfo.append ({'sino': impCount, 'impType': 1, 'impAtom1': impAtom1, 'impAtom2': impAtom2, 'impAtom3': impAtom3, 'impAtom4': impAtom4})

	return improperInfo

def printDataFile (polymerEntries, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, outputdir):

	nAtoms = len (polymerEntries)
	nBonds = len (bondInfo)
	nAngles = len (angleInfo)
	nDihedrals = len (dihedralInfo)
	nImpropers = len (improperInfo)

	nAtomTypes = len (atomTypeArr)
	nBondTypes = len (bondTypeArr) - 1
	nAngleTypes = len (angleTypeArr) - 1
	nDihedralTypes = len (dihTypeArr) - 1
	nImproperTypes = 1

	coords_x = []
	coords_y = []
	coords_z = []

	with open (outputdir + "atomEntries.output", "w") as file:
		for atomLine in polymerEntries:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (atomLine ['sino'], atomLine ['molType'], atomLine ['atomType'], atomLine ['charge'], atomLine ['x'], atomLine ['y'], atomLine ['z']))
			coords_x.append (atomLine ['x'])
			coords_y.append (atomLine ['y'])
			coords_z.append (atomLine ['z'])

	coords_xlo = min (coords_x) - 0.5 * min (coords_x)
	coords_xhi = max (coords_x) + 0.5 * max (coords_x)
	coords_ylo = min (coords_y) - 0.5 * min (coords_y)
	coords_yhi = max (coords_y) + 0.5 * max (coords_y)
	coords_zlo = min (coords_z) - 0.5 * min (coords_z)
	coords_zhi = max (coords_z) + 0.5 * max (coords_z)

	with open (outputdir + "bondEntries.output", "w") as file:
		for bondLine in bondInfo:
			file.write ("\t{}\t{}\t{}\t{}\n".format (bondLine ['sino'], bondLine ['bondType'], bondLine ['bondAtom1'], bondLine ['bondAtom2']))

	with open (outputdir + "angleEntries.output", "w") as file:
		for angleLine in angleInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\n".format (angleLine ['sino'], angleLine ['angleType'], angleLine ['angleAtom1'], angleLine ['angleAtom2'], angleLine ['angleAtom3']))

	with open (outputdir + "dihedralEntries.output", "w") as file:
		for dihedralLine in dihedralInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (dihedralLine ['sino'], dihedralLine ['dihType'], dihedralLine ['dihAtom1'], dihedralLine ['dihAtom2'], dihedralLine ['dihAtom3'], dihedralLine ['dihAtom4']))

	with open (outputdir + "output.data", "w") as file:
		file.write ("Created by you v1.8.1 on today, this month, this year, current time.\n\n\t{}\tatoms\n\t{}\tbonds\n\t{}\tangles\n\t{}\tdihedrals\n\t{}\timpropers\n\n\t{} atom types\n\t{} bond types\n\t{} angle types\n\t{} dihedral types\n\t{} improper types\n\n\t{}\t{}\txlo xhi\n\t{}\t{}\tylo yhi\n\t{}\t{}\tzlo zhi\n\nMasses\n\n\t1\t13.0907\t#CG311 CH\n\t2\t14.1707\t#CG321 CH2\n\t3\t15.2507\t#CG331 CH3\n\nAtoms\n\n".format (nAtoms, nBonds, nAngles, nDihedrals, nImpropers, nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes, round (coords_xlo, 4), round (coords_xhi, 4), round (coords_ylo, 4), round (coords_yhi, 4), round (coords_zlo, 4), round (coords_zhi, 4)))

		for atomLine in polymerEntries:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (atomLine ['sino'], atomLine ['molType'], atomLine ['atomType'], atomLine ['charge'], atomLine ['x'], atomLine ['y'], atomLine ['z']))

		file.write ("\nBonds\n\n")
		for bondLine in bondInfo:
			file.write ("\t{}\t{}\t{}\t{}\n".format (bondLine ['sino'], bondLine ['bondType'], bondLine ['bondAtom1'], bondLine ['bondAtom2']))

		file.write ("\nAngles\n\n")
		for angleLine in angleInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\n".format (angleLine ['sino'], angleLine ['angleType'], angleLine ['angleAtom1'], angleLine ['angleAtom2'], angleLine ['angleAtom3']))

		file.write ("\nDihedrals\n\n")
		for dihedralLine in dihedralInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (dihedralLine ['sino'], dihedralLine ['dihType'], dihedralLine ['dihAtom1'], dihedralLine ['dihAtom2'], dihedralLine ['dihAtom3'], dihedralLine ['dihAtom4']))

	with open (outputdir + "bondTypes.output", "w") as file:
		file.write ("atomType1 atomType2 bondType\n")
		for bondType in bondTypeArr [1:]:
			file.write ("{} {} {}\n".format (bondType ['atom1'], bondType ['atom2'], bondType ['bondType']))

	with open (outputdir + "angleTypes.output", "w") as file:
		file.write ("atomType1 atomType2 atomType3 angleType\n")
		for angleType in angleTypeArr [1:]:
			file.write ("{} {} {} {}\n".format (angleType ['atom1'], angleType ['atom2'], angleType ['atom3'], angleType ['angleType']))

	with open (outputdir + "dihTypes.output", "w") as file:
		file.write ("atomType1 atomType2 atomType3 atomType4 dihType\n")
		for dihType in dihTypeArr [1:]:
			file.write ("{} {} {} {} {}\n".format (dihType ['atom1'], dihType ['atom2'], dihType ['atom3'], dihType ['atom4'], dihType ['dihType']))

if __name__ == '__main__':
	parser = a.ArgumentParser (formatter_class = a.RawDescriptionHelpFormatter, description = "~~~~~~~~~~~~~~~~~\nABOUT THE PROGRAM\n~~~~~~~~~~~~~~~~~\n\ncreatePolymer.py can be used to generate LAMMPS data file for polymeric materials. Follow the steps below,\n\n   1. Create the monomer structure using Avogadro software\n\n   2. Carry out energy minimization of the structure\n\n   3. Orient the monomer along a prefered axis\n\n   4. Save the structure in CML format\n\nUse the following arguments to run the program\n\n   1. --cml Specify the input CML file containing energy minimized monomer structure\n\n   2. --atomTypes Mention the atom type number for the corresponding 'element type' in CML input\n\n   3. --config Config files can be used to input additional parameters such as\n\n\t(i) translateAxis: Same as the axis of orientation used in Avogadro software while creating the monomer\n\n\t(ii) translateDistance: Equal to the distance between identical atoms in two adjacent monomeric units\n\n\t(iii) endGroup1 and endGroup2: Defines the connection point for the monomers\n\n\t(iv) repeatMonomers: Number of monomers required in the polymer\n\n\t(v) mainChain: If the polymer contains side groups, then specify the main chain atoms\n\n\t(vi) pendantGroup: Define the first atom of the pendant group. This information is used to define improper dihedral angles to maintain the necessary tacticity\n\timproper: Number of impropers required. Currently set as 1 by default.\n\n\t~~~~~~~~~~~~~~~~~~~~\n\tExample config file:\n\t~~~~~~~~~~~~~~~~~~~~\n\n\ttranslateAxis: x\n\ttranslateDistance: 2.6\n\tendGroup1: a1\n\tendGroup2: a2\n\trepeatMonomers: 10\n\tmainChain: a1, a2\n\tpendantGroup: a3\n\timproper: 1\n\n\t~~~~~~~~~~~~~~~~~~~~\n\n")
	parser.add_argument ("--cml", "-c", type = str, required = True)
	parser.add_argument ("--atomtypes", "-at", type = str, required = True)
	parser.add_argument ("--config", "-co", type = str, required = True)
	parser.add_argument ("--charges", "-ch", type = str, required = True)
	parser.add_argument ("--outputdir", "-od", type = str, default = "")
	args = parser.parse_args ()
	polymerEntries, atomEntries, atomTypes = createPolymer (args.cml, args.atomtypes, args.config, args.charges)

	atomTypeArr = []
	# polymerEntries = []
	bondInfo = []
	angleInfo = []
	dihedralInfo = []
	improperInfo = []

	for key in atomTypes:
		atomTypeArr.append (atomTypes [key])

	# polymerEntries, atomTypeArr = readAtomInfo ("atomEntries.testing", polymerEntries)
	bondInfo, bondTypeArr = createBonds (polymerEntries, bondInfo)
	angleInfo, angleTypeArr = createAngles (polymerEntries, angleInfo, bondInfo)
	dihedralInfo, dihTypeArr = createDihedrals (polymerEntries, dihedralInfo, angleInfo, bondInfo)
	improperInfo = createImpropers (polymerEntries, atomEntries, improperInfo, args.config)

	printDataFile (polymerEntries, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, args.outputdir)
