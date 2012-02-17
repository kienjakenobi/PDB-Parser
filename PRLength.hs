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
import DataStructures
import Meta
import CustomMath

-- furthest atom locator
-- Locate the atom furthest from the alpha Carbon in a given residue group
fal :: Protein -> (Double, Atom)
fal p = maximum dAtoms
    where
        -- I safely assume that a given residue group has only a single alpha Carbon
        caAtom = head $ [x | x <- p, atomName x == BC.pack "CA"]
        -- Atoms that are not the alpha Carbon
        otherAtoms = [x | x <- p, atomName x /= BC.pack "CA"]
        -- List of distances from the alpha Carbon to the other atoms
        distances = map (cd caAtom) otherAtoms
        --List of all tuples atoms and their distances from caAtom
        dAtoms = zip distances otherAtoms

-- TODO: This function fails on some PDB files, such as 1KDF
--Take a list sorted by amino acid and residue group and output a list of
-- Feed this something out of resGrouper
furthestCalcM :: [[[Atom]]] -> [(B.ByteString, [(B.ByteString, Double)])]
furthestCalcM [] = []
furthestCalcM p =  [(atomAmino ref, together)] ++ furthestCalcM (tail p)
    where
        -- The group of atoms of the same amino acid type grouped by residue
        --  that we are currently working on
        h = head p
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
        furthest = map fal h
        -- Reference atom
        ref = snd $ head furthest

-- Interface to the furthestCalcM function
furthestCalc :: Protein -> [(B.ByteString, [(B.ByteString, Double)])]
furthestCalc p = furthestCalcM $ groupResAmino p

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
        