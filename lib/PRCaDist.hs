-- Description: Functions for calculating the distances between alpha Carbons

-- Usage:
-- p <- pr "/path/to/file.pdb"
-- To get atom distances from left to right of atom: lindist p
-- To get only backbone atoms: backbn p
-- To calculate distance between all 

-- Top level functions (Those not referenced by any other function)
--  lindist (I=Protein and sequential distnace number k, O=distances between k
--        separated atoms)
--  backbn (I=Protein, O=List of backbone atoms N, CA, C, O)
--  alldist (I=Protein, sequential distance number k, and atom a; O=
--      list of distances between atom a and all atoms at least k away
--      in the given protein)
--  alphasdist (I=Protein out of pr, sequential number k; O= list of distnaces
--      between all pairs of alpha carbons at least k residues away)

module PRCaDist (dist, lindist, backbn, alldist, alphasdist) where

import qualified Data.ByteString.Char8 as BC 

-- Custom libraries
import PDBDataStructures
import Meta
import PDBMath

--Always start with ind = 0 to start with the first element in the 
--  given list of atoms p
--  Input=Protein from pr, Output=List of distances each atom has from its
--      neighbors that are k away in sequence tupled with atom id
distM :: [Atom] -> Int -> Int -> [(Int, Double)]
distM p k ind
-- Check to see if there exist entries in the list p k places away from
    | ind == ed = [(aid, cd ca la)]
    | (ind - k) >= 0 && (ind + k) <= ed = [(aid, cd ca la), (aid, cd ca ra)] ++ distM p k (ind + 1)
    | (ind - k) < 0 && (ind + k) <= ed = [(aid, cd ca ra)] ++ distM p k (ind + 1)
    | (ind - k) >= 0 && (ind + k) > ed = [(aid, cd ca la)] ++ distM p k (ind + 1)
    | otherwise = []
    where
        -- Calculate the end of the list
        ed = (length p) - 1
        -- left attom
        la = p !! (ind - k)
        -- right atom
        ra = p !! (ind + k)
        -- current atom
        ca = p !! ind
        aid = atomRef ca

-- Call distM with ind = 0, because the index should always start at 0
dist :: [Atom] -> Int -> [(Int, Double)]
dist p k = distM p k 0

-- A linear version of distM.  That is, only a single point for a given x value
--  Useful for calculating averages
--  Input=Protein from pr, Output=List of distances each atom has from 
--      its neighbor to the right that is k away in sequence
lindistM :: [Atom] -> Int -> Int -> [(Int, Double)]
lindistM p k ind
    | ind + (k - 1) == ed = []
    | otherwise = [(k, cd ca ra)] ++ lindistM p k (ind + 1)
    where
        -- Calculate the end of the list
        ed = (length p) - 1
        --right atom
        ra = p !! (ind + k)
        --current atom
        ca = p !! ind

-- Call lindist with ind = 0 because indicator should always start 0
lindist :: [Atom] -> Int -> [(Int, Double)]
lindist p k = lindistM p k 0

-- Take a Protein and output a list of atoms in the protein which 
--      produces the N, CA, C, O sequence, which is signature of 
--      the protein's backbone
backbn :: [Atom] -> [Atom]
backbn [] = []
backbn p
    | ana == n && anb == ca && anc == cstring && ane == o = [a, b, c, e] ++ backbn (tail p)
    | otherwise = backbn (tail p)
    where
        a = head p
        ana = atomName a
        b = p !! 1
        anb = atomName b
        c = p !! 2
        anc = atomName c
        e = p !! 3
        ane = atomName e
        -- Conver these atom names from strings to ByteStrings for comparison
        n = BC.pack "N"
        ca = BC.pack "CA"
        cstring = BC.pack "C"
        o = BC.pack "O"

-- TODO: Write a function which calculates the distance between a given atom and all
--  the atoms in a given list of atoms
-- TODO: Then modify that function so that it can exclude distance measurements if 
--  the other atom is within k residues of the given atom
alldist :: [Atom] -> Int -> Atom -> [Double]
alldist [] _ _ = []
alldist p k a
    | h /= a && abs((atomRes h) - atomRes a) > k = [cd a h] ++ alldist (tail p) k a
    | otherwise = alldist (tail p) k a
    where
        h = head p

-- Take a protein straight out of pr and calculate the distances between all alpha carbons
--  that are k residues away
alphasdist :: [Atom] -> Int -> [Double]
alphasdist p k = concat $ map (alldist alphas k) alphas
    where
        --Get the backbone and then remove the N, C, and O's
        alphas = delAtomName (backbn p) ["N", "C", "O"]

----single point distance function
----calculate the distance between pairs of atoms in a list
--spdf :: Protein -> [(Int, Double)]
--spdf [] = []
--spdf p = [(aid, cd a1 a2)] ++ spdf (tail (tail p))
--    where
--        a1 = p !! 0
--        a2 = p !! 1
--        (aid, _, _, _, _, _) = a1