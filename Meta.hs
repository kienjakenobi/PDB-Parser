-- Description: Functions dealing with meta attributes, like finding atoms with 
-- certain properties

module Meta (getAmino, aminosList, groupRes, fBE, delAtomName, 
    getAlt, getAtomName, getAtomID, groupResAmino, getChain,
    getHelixAtoms, getSheetAtoms, getLoopAtoms) where

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.ByteString.Internal --c2w
import Data.List -- \\

-- Custom libraries
import DataStructures

-- This file contains three types of functions, as indicated by the name of each function:
-- get = Take a list of atoms and return that list with only those elements matching some query
-- del = Take a list of atoms and return that list with only those that do NOT match some query
-- group = Take a list of atoms and create further lists out of the atoms based on some query

--find aminos
-- Return a list of all atoms matching an amino acid
getAmino :: Protein -> String -> Protein
getAmino [] _ = []
getAmino p bs
    | abb == BC.pack bs = [h] ++ getAmino (tail p) bs
    | otherwise = [] ++ getAmino (tail p) bs
    where
        h = head p
        abb = atomAmino h

--Take a list of atoms and turn it into a list of residue groups, grouped by amino acid
-- Feed this something out of pr
groupResAmino :: Protein -> [[[Atom]]]
groupResAmino p = resGroups
    where
        -- List of groups of atoms of the same amino acid type
        aminoGroups = map (getAmino p) aminosList -- :t [Protein], length = 20
        -- List of groups of atoms of the same residue chain, within the same amino acid type
        resGroups = map groupRes aminoGroups -- :t [[[Atom]]], length = 20

--List of all amino acid abbreviations
aminosList :: [String]
aminosList = list
    where
        list = ["ALA","ARG","ASN"
               ,"ASP","CYS","GLN"
               ,"GLU","GLY","HIS"
               ,"ILE","LEU","LYS"
               ,"MET","PHE","PRO"
               ,"SER","THR","TRP"
               ,"TYR","VAL"]

-- residue seperator
-- Separate a list of atoms into a list of residue groups
groupRes :: Protein -> [[Atom]]
groupRes [] = []
groupRes p = [taken] ++ groupRes (p \\ taken)
    where
        h = head p
        nm = atomRes h
        cm = atomChain h
        taken = [x | x <- p, atomRes x == nm, atomChain x == cm]

--TODO: Create a single function which performs all four function: fBE, foBE, fBN, foBN

--filter by elements
--filter a Protein by the atoms' elements
-- Take a protein from pr and an element by character and return a list of
--  atoms which are of that element
-- Get only those atoms who match the given atom element symbol
fBE :: Protein -> String -> Protein
fBE [] _ = []
fBE p str
    | symbol == query = [x] ++ fBE (tail p) str
    | otherwise = fBE (tail p) str
    where
        x = head p
        symbol = atomSymbol x
        query = BC.pack str

--get all atoms in protein with given atom name, such as "CA" or "OG1"
getAtomName :: Protein -> String -> Protein
getAtomName [] _ = []
getAtomName p str
    | atomName h == query = [h] ++ getAtomName (tail p) str
    | otherwise = getAtomName (tail p) str where
    h = head p
    query = BC.pack str

-- Remove those atoms which match the given atom name, such as "CA" or "OG1"
delAtomName :: Protein -> [String] -> Protein
delAtomName [] _ = []
delAtomName p s
    | e `elem` converted = delAtomName (tail p) s
    | otherwise = [x] ++ delAtomName (tail p) s where
    x = head p
    e = atomName x
    converted = map BC.pack s

-- Take a protein of atoms and a helix and return only those atoms
--  which are in that helix
getHelixAtoms :: Protein -> Helix -> Protein
getHelixAtoms [] _ = []
getHelixAtoms p h
    | startLoc <= atomRes a && atomRes a <= endLoc = [a] ++ getHelixAtoms (tail p) h
    | otherwise = getHelixAtoms (tail p) h where
    a = head p
    (startLoc, endLoc) = helixResList h

-- Take a protein of atoms and a sheet and return only those atoms
--  which are in that sheet
getSheetAtoms :: Protein -> Sheet -> Protein
getSheetAtoms [] _ = []
getSheetAtoms p s
    | startLoc <= atomRes a && atomRes a <= endLoc = [a] ++ getSheetAtoms (tail p) s
    | otherwise = getSheetAtoms (tail p) s where
    a = head p
    (startLoc, endLoc) = sheetResList s

-- Take a protein, a list of all its helices, a list of all its helices, 
--  a list of all its sheets, and return a list of atoms which are not
--  in those helices or sheets
getLoopAtoms :: Protein -> [Helix] -> [Sheet] -> Protein
getLoopAtoms p h s = [i | i <- p, i `notElem` helixAtoms, i `notElem` sheetAtoms] where
    helixAtomsList = map (getHelixAtoms p) h
    sheetAtomsList = map (getSheetAtoms p) s
    helixAtoms = concat helixAtomsList
    sheetAtoms = concat sheetAtomsList

-- filter by altloc
-- Return a list of all atoms of a certain altloc type
-- Example: getAlt p ['A', ' ']
getAlt :: Protein -> [Char] -> Protein
getAlt [] _ = []
getAlt p l
    -- Take not only those atoms with an identical altloc indicator, but
    --      also those with no altloc indocator, that is, a space
    | atomAlt a `elem` values = [a] ++ getAlt (tail p) l
    | otherwise = getAlt (tail p) l where
    a = head p
    values = map c2w l

-- Take a protein and return a list of atoms with the given atom reference ID #
getAtomID :: Protein -> Int -> Protein
getAtomID [] _ = []
getAtomID p i
    | atomRef h == i = [h] ++ getAtomID (tail p) i
    | otherwise = getAtomID (tail p) i where
    h = head p

getChain :: Protein -> Char -> Protein
getChain [] _ = []
getChain p c
    | atomChain a == wc = [a] ++ getChain (tail p) c
    | otherwise = getChain (tail p) c where
    a = head p
    wc = c2w c