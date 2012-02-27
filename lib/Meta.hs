-- Meta.hs
-- Functions dealing with meta attributes, like finding atoms with 
-- certain properties

module Meta (getAmino, aminosList, groupRes, fBE, delAtomName, 
    getAlt, getAtomName, getAtomID, groupResAmino, getChain,
    getHelixAtoms, getSheetAtoms, getLoopAtoms, getResID) where

import qualified Data.ByteString.Char8 as BC
import Data.ByteString.Internal --c2w
import Data.List -- \\

-- Custom libraries
import PDBDataStructures

-- This file contains three types of functions, as indicated by the name of each function:
-- get = Take a list of atoms and return that list with only those elements matching some query
-- del = Take a list of atoms and return that list with only those that do NOT match some query
-- group = Take a list of atoms and create further lists out of the atoms based on some query

-- TODO: Should these functions take and produce Protein or [Atom]?

--find aminos
-- Return a list of all atoms matching an amino acid
getAmino :: [Atom] -> String -> [Atom]
getAmino [] _ = []
getAmino (p:ps) bs
    | atomAmino p == BC.pack bs = p : getAmino ps bs
    | otherwise = [] ++ getAmino ps bs

--Take a list of atoms and turn it into a list of residue groups, grouped by amino acid
-- Feed this something out of pr
groupResAmino :: [Atom] -> [[Residue]]
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
groupRes :: [Atom] -> [Residue]
groupRes [] = []
groupRes pr = res : groupRes (pr \\ takenAtoms)
    where
        p = head pr
        nm = atomRes p
        cm = atomChain p
        takenAtoms = [x | x <- pr, atomRes x == nm, atomChain x == cm]
        res = Residue
            { resName = atomAmino p
            , resID = atomRes p
            , resAtoms = takenAtoms
            }

-- Take a list of residues and a desired ID
--  Return the 
getResID :: [Residue] -> Int -> Maybe Residue
getResID [] _ = Nothing
getResID (r:rs) i
    --
    | resID r == i = Just r
    | otherwise = getResID rs i where
    

--TODO: Create a single function which performs all four function: fBE, foBE, fBN, foBN

--filter by elements
--filter a Protein by the atoms' elements
-- Take a protein from pr and an element by character and return a list of
--  atoms which are of that element
-- Get only those atoms who match the given atom element symbol
fBE :: [Atom] -> String -> [Atom]
fBE [] _ = []
fBE (p:ps) str
    -- If atomSymbol same as given string, then take it
    | atomSymbol p == BC.pack str = p : fBE ps str
    | otherwise = fBE ps str

-- Get all atoms in protein with given atom name, such as "CA" or "OG1"
-- Example: getAtomName ps 
getAtomName :: [Atom] -> String -> [Atom]
getAtomName [] _ = []
getAtomName (p:ps) str
    | atomName p == query = p : getAtomName ps str
    | otherwise = getAtomName ps str where
    query = BC.pack str

-- Remove those atoms which match the given atom name, such as "CA" or "OG1"
delAtomName :: [Atom] -> [String] -> [Atom]
delAtomName [] _ = []
delAtomName (p:ps) s
    | e `elem` converted = delAtomName ps s
    | otherwise = p : delAtomName ps s where
    e = atomName p
    converted = map BC.pack s

-- Take a protein of atoms and a helix and return only those atoms
--  which are in that helix
getHelixAtoms :: [Atom] -> Helix -> [Atom]
getHelixAtoms [] _ = []
getHelixAtoms (p:ps) h
    | startLoc <= atomRes p && atomRes p <= endLoc = p : getHelixAtoms ps h
    | otherwise = getHelixAtoms ps h where
    (startLoc, endLoc) = helixResList h

-- Take a protein of atoms and a sheet and return only those atoms
--  which are in that sheet
getSheetAtoms :: [Atom] -> Sheet -> [Atom]
getSheetAtoms [] _ = []
getSheetAtoms (p:ps) s
    -- If the atom's location in the chain is within the sheet's range,
    --  then take it
    | startLoc <= atomRes p && atomRes p <= endLoc = p : getSheetAtoms ps s
    | otherwise = getSheetAtoms ps s where
    (startLoc, endLoc) = sheetResList s

-- Take a protein, a list of all its helices, a list of all its helices, 
--  a list of all its sheets, and return a list of atoms which are not
--  in those helices or sheets
getLoopAtoms :: [Atom] -> [Helix] -> [Sheet] -> [Atom]
getLoopAtoms p h s = [i | i <- p, i `notElem` helixAtoms, i `notElem` sheetAtoms] where
    helixAtomsList = map (getHelixAtoms p) h
    sheetAtomsList = map (getSheetAtoms p) s
    helixAtoms = concat helixAtomsList
    sheetAtoms = concat sheetAtomsList

-- filter by altloc
-- Return a list of all atoms of a certain altloc type
-- Example: getAlt p ['A', ' ']
getAlt :: [Atom] -> String -> [Atom]
getAlt [] _ = []
getAlt p l
    -- Take not only those atoms with an identical altloc indicator, but
    --      also those with no altloc indocator, that is, a space
    | atomAlt a `elem` values = a : getAlt (tail p) l
    | otherwise = getAlt (tail p) l where
    a = head p
    values = map c2w l

-- Take a protein and return a list of atoms with the given atom reference ID #
getAtomID :: [Atom] -> Int -> [Atom]
getAtomID [] _ = []
getAtomID p i
    | atomRef h == i = h : getAtomID (tail p) i
    | otherwise = getAtomID (tail p) i where
    h = head p

-- Get only atoms of a certain chain
-- Example: getChain p 'A'
getChain :: [Atom] -> Char -> [Atom]
getChain [] _ = []
getChain (p:ps) c
    -- If the atom's chain is the same as the given char, then take it
    | atomChain p == c2w c = p : getChain ps c
    | otherwise = getChain ps c where
