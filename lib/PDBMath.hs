-- A library of general mathematics functions
module PDBMath (cd, dihedral, phiangles, psiangles) where

import Numeric.Container

-- Custom libraries
import Meta
import PDBDataStructures
import Math

-- calculate 3D Cartesian distance
-- Calculate the 3D distance between two atoms
cd :: Atom -> Atom -> Double
cd a1 a2 = sqrt ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    where
        [x1, y1, z1] = toList $ atomPos a1
        [x2, y2, z2] = toList $ atomPos a2

-- Take a residue and calculate the sidechain rotamer dihedral 
--      angle relative to a given fourth atom name
--      Output is in radians
-- Returns a value of 0 if something screwy happens
dihedral :: String -> [Atom] -> Double
dihedral _ [] = 0.0
dihedral str p
    | xAtom /= [] = 180*dirad/pi
    | otherwise = 0.0 where
    nAtom = getAtomName p "N"
    caAtom = getAtomName p "CA"
    cbAtom = getAtomName p "CB"
    xAtom = getAtomName p str
    a = atomPos $ head nAtom
    b = atomPos $ head caAtom
    c = atomPos $ head cbAtom
    d = atomPos $ head xAtom
    n1 = (a `sub` b) `cross` (c `sub` b)
    n2 = (b `sub` c) `cross` (d `sub` c)
    dirad = acos $ (n1 <.> n2)/( norm2 n1 * norm2 n2)
    --residue = atomRes $ head xAtom

-- Take a protein and output the phi angle as defined in Scott, page 72
--  equation (5.5)
--  Input should be a protein run through groupRes
--  Output angle values in degrees
phiangles :: [[Atom]] -> [Double]
phiangles (p:ps)
    -- Check that there are at least two residues in the list
    | length p >= 2 = ((acos (n1 <.> n2)/(norm2 n1 * norm2 n2)) * 180 / pi) : phiangles ps
    -- If there are not enough residues left in the list, then done
    | otherwise = [] where
    -- Get the interesting atoms from each residue
    c0Atom = getAtomName p "C"
    n1Atom = getAtomName (ps!!1) "N"
    ca1Atom = getAtomName (ps!!1) "CA"
    c1Atom = getAtomName (ps!!1) "C"
    -- Get the coordinates of each intersting atom
    a = atomPos $ head c0Atom
    b = atomPos $ head n1Atom
    c = atomPos $ head ca1Atom
    d = atomPos $ head c1Atom
    -- Define the math quantities for the final equation
    n1 = (a `sub` b) `cross` (c `sub` b)
    n2 = (b `sub` c) `cross` (d `sub` c)

-- Take a protein and output the phi angle as defined in Scott, page 72
--  equation (5.5)
--  Input should be a protein run through groupRes
--  Output angle values in degrees
psiangles :: [[Atom]] -> [Double]
psiangles p
    -- Check that there are at least two residues in the list
    | length p >= 2 = ((acos (n1 `dot` n2)/(norm2 n1 * norm2 n2)) * 180 / pi) : psiangles (tail p)
    -- If there are not enough residues left in the list, then done
    | otherwise = [] where
    -- Get the first and second residues
    a0 = head p
    a1 = head $ tail p
    -- Get the interesting atoms from each residue
    n0Atom = getAtomName a0 "N"
    ca0Atom = getAtomName a0 "CA"
    c0Atom = getAtomName a0 "C"
    n1Atom = getAtomName a1 "N"
    -- Get the coordinates of each intersting atom
    a = atomPos $ head n0Atom
    b = atomPos $ head ca0Atom
    c = atomPos $ head c0Atom
    d = atomPos $ head n1Atom
    -- Define the math quantities for the final equation
    n1 = (a `sub` b) `cross` (c `sub` b)
    n2 = (b `sub` c) `cross` (d `sub` c)
