--Exercise 1: Consider the distribution of distances between CA's 
--  separated in sequence by k=2 in a given PDB file. Explain the 
--  bimodal character in PDB file 2P9R. Correlate the distance with 
--  phi and psi angles.

-- Get the distances between k=2 alpha Carbons
main = do
    p <- pr "/Users/mortal/Desktop/2P9R.pdb"
    let cas = getAtomName p "CA"
    let chaina = getChain cas 'A'
    let chainb = getChain cas 'B'
    let adist = lindist chaina 2
    print adist
    print $ length adist
    let bdist = lindist chainb 2
    print bdist
    print $ length bdist

-- Get the phi and psi angles
main = do
    p <- pr "/Users/mortal/Desktop/2P9R.pdb"
    let chaina = getChain p 'A'
    let chainb = getChain p 'B'
    let resa = groupRes chaina
    let resb = groupRes chainb
    let psisa = psiangles resa
    let psisb = psiangles resb
    print psisa
    print $ length psisa
    print psisb
    print $ length psisb
    let phisa = phiangles resa
    let phisb = phiangles resb
    print phisa
    print $ length phisa
    print phisb
    print $ length phisb

-- Correlate the distances between CA's separated in sequence by k
--  with the secondary structure call in a given PDB file. That is,
--  give the (three, separate) distributions of distances as a 
--  function of k for helices and sheets and loops 
--  (that is, for residues not in either structure).
    p <- pr "/Users/mortal/Desktop/3O4P.pdb"
    -- Get generic data about protein
    -- Get rid of the B altlocs
    let theas = getAlt p ['A', ' ']
    --print theas
    --print $ length theas
    
    -- Get the sheet and helix data structures
    helices <- hr "/Users/mortal/Desktop/3O4P.pdb"
    sheets <- sr "/Users/mortal/Desktop/3O4P.pdb"
    --print $ length helices
    --print $ length sheets
    --print helices
    --print sheets
    
    let helixAtomsList = map (getHelixAtoms theas) helices
    let sheetAtomsList = map (getSheetAtoms theas) sheets
    --print helixAtomsList
    --print sheetAtomsList
    let helixAtoms = concat helixAtomsList
    let sheetAtoms = concat sheetAtomsList
    let loopAtoms = getLoopAtoms theas helices sheets
--    print $ length helixAtoms
--    print $ length sheetAtoms
--    print $ length loopAtoms
    --let resHelix = groupRes helixAtoms
    --let resSheet = groupRes sheetAtoms
    --let resLoop = groupRes loopAtoms
    --print $ length resHelix
    --print $ length resSheet
    --print $ length resLoop

    let hcas = getAtomName helixAtoms "CA"
    let scas = getAtomName sheetAtoms "CA"
    let lcas = getAtomName loopAtoms "CA"
    --let hdist = concat $ map (lindist hcas) [1,2,3,4]
    --print hdist
    --let sdist = concat $ map (lindist scas) [1..152]
    --print sdist
    let ldist = concat $ map (lindist lcas) [1..155]
    print ldist

-- Determine whether the choice of atom chi_delta for Asp or chi_epsilon
--  for Glu in equation 5.6 matters by examining high-resolution 
--  structures. For each Asp (and Glu), compute the difference of dihedral
--  angles for the two choices of terminal atoms, and plot the 
--  distribution of angle distances, in degrees.
    p <- pr "/Users/mortal/Desktop/3O4P.pdb"
    -- Get generic data about protein
    let theas = getAlt p ['A', ' ']
    print theas
    print $ length theas
    
    -- Get specifically interesting data
    let asps = getAmino p "ASP"
    let glus = getAmino p "GLU"
    let aspsA = getAlt asps ['A', ' ']
    let glusA = getAlt glus ['A', ' ']
    let resAsps = groupRes aspsA
    let resGlus = groupRes glusA
    print $ resAsps
    print $ resGlus
    let asp1 = map (dihedral "OD1") resAsps
    print asp1
    let asp2 = map (dihedral "OD2") resAsps
    print asp2
    let glu1 = map (dihedral "OE1") resGlus
    print glu1
    let glu2 = map (dihedral "OE2") resGlus
    print glu2
    print $ map length [asp1, asp2, glu1, glu2]