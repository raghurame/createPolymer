import argparse as a
from time import sleep
import decimal
from itertools import combinations

def cutString (inputString, subString1, subString2):
	startIndex = inputString.find (subString1) + len (subString1)
	endIndex = inputString.find (subString2)

	return inputString [startIndex:endIndex]

def readConfig (configFile):
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

def createPolymer (cmlFile, atomTypesFile, configFile):
	translateAxis, endGroup1, endGroup2, repeatMonomers, translateDistance, mainChain, pendantGroup, improper = readConfig (configFile)

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

	endGroupArr = []
	for x in range (0, repeatMonomers, 1):
		for atoms in atomEntries:
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

			if (x == 0 and atoms ['id'] == endGroup2):
				bondAtom4_polymerEntries = endGroup1 + ((x + 1) * nAtomsPerMonomer)
			elif (x == (repeatMonomers - 1) and atoms ['id'] == endGroup1):
				bondAtom4_polymerEntries = endGroup2 + ((x - 1) * nAtomsPerMonomer)
			elif (x > 0 and x < (repeatMonomers - 1)):
				if (atoms ['id'] == endGroup1):
					bondAtom4_polymerEntries = endGroup2 + ((x - 1) * nAtomsPerMonomer)
				if (atoms ['id'] == endGroup2):
					bondAtom4_polymerEntries = endGroup1 + ((x + 1) * nAtomsPerMonomer)

			polymerEntries.append ({'sino': atoms ['id'] + (x * nAtomsPerMonomer), 'monomerAtomID': atoms ['id'], 'elementType': atoms ['elementType'], 'atomType': atomTypes [atoms ['elementType']], 'molType': 1, 'charge': 0, 'x': round (atoms ['x'] + (isXTranslate * x * translateDistance), 4), 'y': round (atoms ['y'] + (isYTranslate * x * translateDistance), 4), 'z': round (atoms ['z'] + (isZTranslate * x * translateDistance), 4), 'bondAtom1': bondAtom1_polymerEntries, 'bondAtom2': bondAtom2_polymerEntries, 'bondAtom3': bondAtom3_polymerEntries, 'bondAtom4': bondAtom4_polymerEntries})

	return polymerEntries, atomEntries, atomTypes

def createPolymer2 (cmlFile, atomTypesFile, configFile):
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
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def readAtomInfo (inputFileName, atomInfo):
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
	lineArray = []
	sino_bond = 1
	bondTypeArr = []
	bondTypeArr.append({'atom1': 0, 'atom2': 0, 'bondType': 0})

	def bondCheck (atom1, atom2, bondInfo):
		for bond in bondInfo:
			if ((atom1 == bond['bondAtom1'] or atom1 == bond['bondAtom2']) and (atom2 == bond['bondAtom1'] or atom2 == bond['bondAtom2'])):
				return 0

		return 1

	def findBondType (atom1, atom2, bondTypeArr):
		returnBondType = 0
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

	for atomLine in polymerEntries:
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom1', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom2', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom3', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, polymerEntries, 'bondAtom4', bondInfo, sino_bond, bondTypeArr)

	return bondInfo, bondTypeArr

def createAngles (polymerEntries, angleInfo, bondInfo):
	# angleInfo contains sino, angleType, angleAtom1, angleAtom2, angleAtom3, angleAtom1Type, angleAtom2Type, angleAtom3Type
	angleTypeArr = []
	sino_angle = 1

	# Initializing connectedAtoms dict of lists
	connectedAtoms = {}
	for atomLine in polymerEntries:
		connectedAtoms [atomLine ['sino']] = []

	angleTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'angleType': 0})

	def findConnectedAtoms (concernedBondAtom, primaryConnect, bondInfo, connectedAtoms):
		if (primaryConnect not in connectedAtoms [concernedBondAtom]):
			connectedAtoms [concernedBondAtom].append (primaryConnect)
		if (concernedBondAtom not in connectedAtoms [primaryConnect]):
			connectedAtoms [primaryConnect].append (concernedBondAtom)

		for bondLine in bondInfo:
			if ((bondLine ['bondAtom1'] == concernedBondAtom) and (bondLine ['bondAtom2'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom2'])
				connectedAtoms [concernedBondAtom].sort ()
			if ((bondLine ['bondAtom2'] == concernedBondAtom) and (bondLine ['bondAtom1'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom1'])
				connectedAtoms [concernedBondAtom].sort ()

		return connectedAtoms

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
	# sino impType impAtom1 impAtom2 impAtom3 impAtom4 (impAtom2 is the central atom from main chain)
	# The following code assumes there is only one type of improper angle
	# Because the most common use for an improper angle is to maintain tacticity of the molecule
	improperEntries = {}
	improperTypeArr = []
	improperTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'atom4': 0, 'improperType': 0})

	translateAxis, endGroup1, endGroup2, repeatMonomers, translateDistance, mainChain, pendantGroup, improper = readConfig (configFile)

	mainChain = mainChain.replace ("a", "").replace (" ", "")
	mainChainAtoms = mainChain.split (",")
	mainChainAtoms = [int(i) for i in mainChainAtoms]
	pendantGroupAtoms = pendantGroup.replace ("a", "").replace (" ", "")
	# print (mainChainAtoms)
	# print (pendantGroupAtoms)
	# print (improper)

	# Find the atom (from main chain) connected to the pendant atom
	mainChainAtomConnectedToPendant = 0
	mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
	mainChainAtomsArray = []
	pendantAtomsArray = []

	for atoms in polymerEntries:
		if (str (atoms ['monomerAtomID']) in mainChainAtoms):
			mainChainAtomsArray.append (atoms)
		if (str (atoms ['monomerAtomID']) in pendantGroupAtoms):
			pendantAtomsArray.append (atoms)

	# print ("mainChainAtoms", mainChainAtoms)
	# print ("pendantGroupAtoms", pendantGroupAtoms)
	# print ("mainChainAtomsArray", mainChainAtomsArray)
	# print ("pendantAtomsArray", pendantAtomsArray)

	for atoms in pendantAtomsArray:
		improperEntries [atoms['sino']] = []
		# print ("\n\n")
		if (atoms ['bondAtom1'] and polymerEntries [atoms ['bondAtom1'] - 1]['monomerAtomID'] in mainChainAtoms):
			# polymerEntries [atoms ['bondAtom1'] - 1]
			mainChainAtomConnectedToPendant = 0
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
			mainChainAtomConnectedToPendant = polymerEntries [atoms ['bondAtom1'] - 1]['monomerAtomID']
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [atoms ['bondAtom1'] - 1]['bondAtom1'], polymerEntries [atoms ['bondAtom1'] - 1]['bondAtom2'], polymerEntries [atoms ['bondAtom1'] - 1]['bondAtom3'], polymerEntries [atoms ['bondAtom1'] - 1]['bondAtom4']]
			# print (mainChainAtomConnectedToPendant_connectedMainChainAtoms)
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					# print ("Pendant atoms: ", atoms)
					# print ("first main chain atom: ", polymerEntries [atoms ['bondAtom1'] - 1])
					# print ("second main chain atom: ", polymerEntries [atoms1 - 1])
					improperEntries [atoms['sino']].append (polymerEntries [atoms ['bondAtom1'] - 1]['sino'])
					improperEntries [atoms['sino']].append (polymerEntries [atoms1 - 1]['sino'])

		if (atoms ['bondAtom2'] and polymerEntries [atoms ['bondAtom2'] - 1]['monomerAtomID'] in mainChainAtoms):
			# polymerEntries [atoms ['bondAtom2'] - 1]
			mainChainAtomConnectedToPendant = 0
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
			mainChainAtomConnectedToPendant = polymerEntries [atoms ['bondAtom2'] - 1]['monomerAtomID']
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [atoms ['bondAtom2'] - 1]['bondAtom1'], polymerEntries [atoms ['bondAtom2'] - 1]['bondAtom2'], polymerEntries [atoms ['bondAtom2'] - 1]['bondAtom3'], polymerEntries [atoms ['bondAtom2'] - 1]['bondAtom4']]
			# print (mainChainAtomConnectedToPendant_connectedMainChainAtoms)
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					# print ("Pendant atoms: ", atoms)
					# print ("first main chain atom: ", polymerEntries [atoms ['bondAtom2'] - 1])
					# print ("second main chain atom: ", polymerEntries [atoms1 - 1])
					improperEntries [atoms['sino']].append (polymerEntries [atoms ['bondAtom2'] - 1]['sino'])
					improperEntries [atoms['sino']].append (polymerEntries [atoms1 - 1]['sino'])

		if (atoms ['bondAtom3'] and polymerEntries [atoms ['bondAtom3'] - 1]['monomerAtomID'] in mainChainAtoms):
			# polymerEntries [atoms ['bondAtom3'] - 1]
			mainChainAtomConnectedToPendant = 0
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
			mainChainAtomConnectedToPendant = polymerEntries [atoms ['bondAtom3'] - 1]['monomerAtomID']
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [atoms ['bondAtom3'] - 1]['bondAtom1'], polymerEntries [atoms ['bondAtom3'] - 1]['bondAtom2'], polymerEntries [atoms ['bondAtom3'] - 1]['bondAtom3'], polymerEntries [atoms ['bondAtom3'] - 1]['bondAtom4']]
			# print (mainChainAtomConnectedToPendant_connectedMainChainAtoms)
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					# print ("Pendant atoms: ", atoms)
					# print ("first main chain atom: ", polymerEntries [atoms ['bondAtom3'] - 1])
					# print ("second main chain atom: ", polymerEntries [atoms1 - 1])
					improperEntries [atoms['sino']].append (polymerEntries [atoms ['bondAtom3'] - 1]['sino'])
					improperEntries [atoms['sino']].append (polymerEntries [atoms1 - 1]['sino'])

		if (atoms ['bondAtom4'] and polymerEntries [atoms ['bondAtom4'] - 1]['monomerAtomID'] in mainChainAtoms):
			# polymerEntries [atoms ['bondAtom4'] - 1]
			mainChainAtomConnectedToPendant = 0
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = []
			mainChainAtomConnectedToPendant = polymerEntries [atoms ['bondAtom4'] - 1]['monomerAtomID']
			mainChainAtomConnectedToPendant_connectedMainChainAtoms = [polymerEntries [atoms ['bondAtom4'] - 1]['bondAtom1'], polymerEntries [atoms ['bondAtom4'] - 1]['bondAtom2'], polymerEntries [atoms ['bondAtom4'] - 1]['bondAtom3'], polymerEntries [atoms ['bondAtom4'] - 1]['bondAtom4']]
			# print (mainChainAtomConnectedToPendant_connectedMainChainAtoms)
			for atoms1 in mainChainAtomConnectedToPendant_connectedMainChainAtoms:
				if (polymerEntries [atoms1 - 1]['monomerAtomID'] in mainChainAtoms):
					# print ("Pendant atoms: ", atoms)
					# print ("first main chain atom: ", polymerEntries [atoms ['bondAtom4'] - 1])
					# print ("second main chain atom: ", polymerEntries [atoms1 - 1])
					improperEntries [atoms['sino']].append (polymerEntries [atoms ['bondAtom4'] - 1]['sino'])
					improperEntries [atoms['sino']].append (polymerEntries [atoms1 - 1]['sino'])

	# for atoms in polymerEntries:
	# 	if (atoms ['monomerAtomID'] == mainChainAtomConnectedToPendant):
	# 		polymerEntries [atoms ['bondAtom1'] - 1]['monomerAtomID'] polymerEntries [atoms ['bondAtom2'] - 1]['monomerAtomID'] polymerEntries [atoms ['bondAtom3'] - 1]['monomerAtomID'] polymerEntries [atoms ['bondAtom4'] - 1]['monomerAtomID']

	# for atoms in polymerEntries:
	# 	print (atoms)
	# 	sleep (1)

	print (improperEntries)

	for entries in improperEntries:
		print (entries, improperEntries [entries])
		sleep (1)

	return improperInfo, improperTypeArr

def printDataFile (polymerEntries, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, improperTypeArr):

	nAtoms = len (polymerEntries)
	nBonds = len (bondInfo)
	nAngles = len (angleInfo)
	nDihedrals = len (dihedralInfo)
	nImpropers = len (improperInfo)

	nAtomTypes = len (atomTypeArr)
	nBondTypes = len (bondTypeArr) - 1
	nAngleTypes = len (angleTypeArr) - 1
	nDihedralTypes = len (dihTypeArr) - 1
	nImproperTypes = 0

	coords_x = []
	coords_y = []
	coords_z = []

	with open ("atomEntries.output", "w") as file:
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

	with open ("bondEntries.output", "w") as file:
		for bondLine in bondInfo:
			file.write ("\t{}\t{}\t{}\t{}\n".format (bondLine ['sino'], bondLine ['bondType'], bondLine ['bondAtom1'], bondLine ['bondAtom2']))

	with open ("angleEntries.output", "w") as file:
		for angleLine in angleInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\n".format (angleLine ['sino'], angleLine ['angleType'], angleLine ['angleAtom1'], angleLine ['angleAtom2'], angleLine ['angleAtom3']))

	with open ("dihedralEntries.output", "w") as file:
		for dihedralLine in dihedralInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (dihedralLine ['sino'], dihedralLine ['dihType'], dihedralLine ['dihAtom1'], dihedralLine ['dihAtom2'], dihedralLine ['dihAtom3'], dihedralLine ['dihAtom4']))

	with open ("output.data", "w") as file:
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

if __name__ == '__main__':
	parser = a.ArgumentParser (description = "Read CML files")
	parser.add_argument ("--cml", "-c", type = str, required = True)
	parser.add_argument ("--atomtypes", "-at", type = str, required = True)
	parser.add_argument ("--config", "-co", type = str, required = True)
	args = parser.parse_args ()
	polymerEntries, atomEntries, atomTypes = createPolymer (args.cml, args.atomtypes, args.config)

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
	improperInfo, improperTypeArr = createImpropers (polymerEntries, atomEntries, improperInfo, args.config)

	printDataFile (polymerEntries, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, improperTypeArr)
