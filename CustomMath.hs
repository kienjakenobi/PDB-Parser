-- Description: A library of general mathematics functions

module CustomMath (cd, l2norm, scalar, dihedral, phiangles, psiangles) where

-- Custom libraries
import DataStructures
import Meta

-- calculate 3D Cartesian distance
-- Calculate the 3D distance between two atoms
cd :: Atom -> Atom -> Double
cd a1 a2 = sqrt ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    where
        Point3D x1 y1 z1 = atomPos a1
        Point3D x2 y2 z2 = atomPos a2

-- Calculate the l2-Norm == Euclidean Norm of a Point3D
l2norm :: Point3D -> Double
l2norm p  = sqrt(x**2 + y**2 + z**2) where
    Point3D x y z = p

scalar :: Point3D -> Point3D -> Double
scalar p1 p2 = (x1*x2) + (y1*y2) + (z1*z2) where
    Point3D x1 y1 z1 = p1
    Point3D x2 y2 z2 = p2
    
-- Take a residue and calculate the sidechain rotamer dihedral 
--      angle relative to a given fourth atom name
--      Output is in radians
-- Returns a value of 0 if something screwy happens
dihedral :: String -> Protein -> Double
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
    n1 = (a - b)*(c - b)
    n2 = (b - c)*(d - c)
    dirad = acos $ (scalar n1 n2)/((l2norm n1)*(l2norm n2))
    --residue = atomRes $ head xAtom

-- Take a protein and output the phi angle as defined in Scott, page 72
--  equation (5.5)
--  Input should be a protein run through groupRes
--  Output angle values in degrees
phiangles :: [Protein] -> [Double]
phiangles p
    -- Check that there are at least two residues in the list
    | length p >= 2 = [(acos ((scalar n1 n2)/((l2norm n1)*(l2norm n2))))*180/pi] ++ phiangles (tail p)
    -- If there are not enough residues left in the list, then done
    | otherwise = [] where
    -- Get the first and second residues
    r0 = p !! 0
    r1 = p !! 1
    -- Get the interesting atoms from each residue
    c0Atom = getAtomName r0 "C"
    n1Atom = getAtomName r1 "N"
    ca1Atom = getAtomName r1 "CA"
    c1Atom = getAtomName r1 "C"
    -- Get the coordinates of each intersting atom
    a = atomPos $ head c0Atom
    b = atomPos $ head n1Atom
    c = atomPos $ head ca1Atom
    d = atomPos $ head c1Atom
    -- Define the math quantities for the final equation
    n1 = (a - b)*(c - b)
    n2 = (b - c)*(d - c)

-- Take a protein and output the phi angle as defined in Scott, page 72
--  equation (5.5)
--  Input should be a protein run through groupRes
--  Output angle values in degrees
psiangles :: [Protein] -> [Double]
psiangles p
    -- Check that there are at least two residues in the list
    | length p >= 2 = [(acos ((scalar n1 n2)/((l2norm n1)*(l2norm n2))))*180/pi] ++ psiangles (tail p)
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
    n1 = (a - b)*(c - b)
    n2 = (b - c)*(d - c)