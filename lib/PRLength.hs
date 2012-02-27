-- Description: Functions to determine the longest "side" of residues

-- Usage:
-- p <- pr [PDB file]
-- let q = furthestCalc p
-- let a = avgDist q

-- Top Level Functions (Those not referenced by any other function)
--  furthestCalc (I=Protein, O=List of furthest atom from alpha
--         carbon grouped by amino acid and residue)
--  avgDist (I=furthestCalc output, O=List of average largest distances grouped
--         by amino acid)

module PRLength (furthestCalc, avgDist) where

import qualified Data.ByteString.Char8 as BC
import qualified Data.ByteString as B
import Data.List -- \\

-- Custom
import PDBDataStructures
import Meta
import PDBMath

-- furthest atom locator
-- Locate the atom furthest from the alpha Carbon in a given residue
fal :: Residue -> (Double, Atom)
fal res = maximum dAtoms
    where
        -- Atoms in the residue
        ats = resAtoms res
        -- I safely assume that a given residue group has only a single alpha Carbon
        caAtom = head $ [x | x <- ats, atomName x == BC.pack "CA"]
        -- Atoms that are not the alpha Carbon
        otherAtoms = [x | x <- ats, atomName x /= BC.pack "CA"]
        -- List of distances from the alpha Carbon to the other atoms
        distances = map (cd caAtom) otherAtoms
        --List of all tuples atoms and their distances from caAtom
        dAtoms = zip distances otherAtoms

-- TODO: This function fails on some PDB files, such as 1KDF
-- Take a list of a residues grouped by amino acid type and output a list of 
-- furthest lengths from the alpha Carbon within each residue
-- Feed this something out of resGrouper
furthestCalcM :: [[Residue]] -> [(B.ByteString, [(B.ByteString, Double)])]
furthestCalcM [] = []
furthestCalcM (p:ps) =  [(atomAmino ref, together)] ++ furthestCalcM ps
    where
        -- List of the distances between the 
        distances = fst $ unzip $ furthest
        -- List of the atoms
        theAtoms = snd $ unzip $ furthest
        -- List of the atomNames for each atom in the above list
        names = map atomName theAtoms
        -- 2-Tuples of the atnomNames and distances
        together = zip names distances
        -- List of the atoms from each residue group within a given amino acid
        --  which are furthest from the alpha of that residue group
        furthest = map fal p
        -- Reference atom
        ref = snd $ head furthest

-- Interface to the furthestCalcM function
-- Feed this something out of pr
furthestCalc :: Protein -> [(B.ByteString, [(B.ByteString, Double)])]
furthestCalc pr = furthestCalcM $ groupResAmino $ head $ proModels pr

-- Calculate the average distance for each atomName in a given amino acid group from furthestCalc
avgDistM :: [(B.ByteString, Double)] -> [(B.ByteString, Int, Double)]
avgDistM [] = []
avgDistM p = [(interAtomName, numberOfThem, average)] ++ avgDistM (p \\ interestingTuples)
    where
        h = head p
        interAtomName = fst h
        interestingTuples = [x | x <- p, fst x == interAtomName]
        quantity = fromIntegral (length interestingTuples) :: Double
        numbers = snd $ unzip p
        average = (sum numbers) / quantity
        numberOfThem = length interestingTuples

-- Take output from furthestCalcs and calculate average distances
avgDist :: [(B.ByteString, [(B.ByteString, Double)])] -> [(B.ByteString, [(B.ByteString, Int, Double)])]
avgDist [] = []
avgDist p = [(amino, (averages))] ++ avgDist (tail p)
    where
        h = head p
        amino = fst h
        datas = snd h
        averages = avgDistM datas
        