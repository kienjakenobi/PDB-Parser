-- Author: Alex Dunn
-- Last Modified: 25 January 2012
-- Description: Data structures and functions to import and store
--  information from a Protein Database file

-- Usage:
--  Open a PDB file: p <- 'pr [file]'

-- Top Level Functions (Those functions not referenced by any other function)
--  pr (I=PDB file location as a string, O=List of Atoms)

module Main where

import qualified Data.ByteString as B -- ByteString class for efficient string operations
import qualified Data.ByteString.Char8 as BC -- unpack, readInt
import Data.ByteString.Lex.Double -- readDouble, from the bytestring-lexing module
import Data.Maybe -- fromJust
import Data.List -- \\

-- Custom libraries:
import DataStructures
import PRLength
import PRCaDist
import Meta
import CustomMath
-- import System.Environment --getArgs

--Read a ByteString into an Int
ri :: B.ByteString -> Int
ri bs = fst $ fromJust $ BC.readInt bs

--Read a ByteString into a Double
rd :: B.ByteString -> Double
rd bs = fst $ fromJust $ readDouble bs

-- ByteString -> [Word8]
-- Data.ByteString.unpack ByteString = B.unpack String

-- Word8 -> Char
-- Data.ByteString.Internal.w2c Word8

-- Char -> Word8
-- Data.ByteString.Internal.c2w Char

-- String -> ByteString
--Data.ByteString.Char8.pack String = BC.pack String

-- ByteString -> String
-- Data.ByteString.Char8.unpack string = BC.unpack String

-- Take substring of a larger ByteString
-- Take start index a and end index b, indices starting at 1
-- Also remove spaces
get :: Int -> Int -> B.ByteString -> B.ByteString
get a b str = B.filter (/=32) $ B.drop (a - 1) $ B.take b str

-- atomize
-- This is the actual parser
-- Take a single line of 
at :: B.ByteString -> Maybe Atom
at str
    -- Test that it is an ATOM entry
    | recName == B.pack a8 = Just Atom 
        {atomRef = serial
        , atomName = name
        , atomAlt = altloc
        , atomAmino = amino
        , atomChain = chain
        , atomRes = resNum
        , atomPos = Point3D x y z
        , atomSymbol = symbol} 
    | otherwise = Nothing where
    -- Get the title of the record
    recName = get 1 6 str
    --Attributes of the atom:
    serial = ri $ get 7 11 str
    name = get 13 16 str
    altloc = B.index str 16
    amino = get 18 20 str
    chain = B.index str 21
    resNum = ri $ get 23 26 str
    x = rd $ get 31 38 str
    y = rd $ get 39 46 str
    z = rd $ get 47 54 str
    symbol = get 77 78 str
    -- "ATOM" as a [Word8]
    a8 = [65,84,79,77]

-- Get all the helices out of a PDB file
ht :: B.ByteString -> Maybe Helix
ht str
    -- Test that it is a HELIX entry
    | recName == B.pack h8 = Just Helix 
        {helixResList = (startLoc, endLoc)} 
    | otherwise = Nothing where
    -- Get the title of the record
    recName = get 1 6 str
    --Attributes of the sheet:
    startLoc = ri $ get 22 25 str
    endLoc = ri $ get 34 37 str
    -- "HELIX" as a [Word8]
    h8 = [72,69,76,73,88]

-- Get all the strips of sheets out of a PDB file
st :: B.ByteString -> Maybe Sheet
st str
    -- Check if it is a SHEET entry
    | recName == B.pack s8 = Just Sheet
        {sheetResList = (startLoc, endLoc)}
    | otherwise = Nothing where
    -- Get the title of the record
    recName = get 1 6 str
    --Attributes of the sheet:
    startLoc = ri $ get 23 26 str
    endLoc = ri $ get 34 37 str
    -- "SHEET" as a [Word8]
    s8 = [83,72,69,69,84]

-- proteinize
-- Perform the listize, furtherListize, and atoms functions on a file
pr :: FilePath -> IO Protein
pr s = do
    pdb <- B.readFile s
    -- Split the file on line breaks
    let split = B.splitWith (==10) pdb
    return $ mapMaybe at split

-- Get all the helices which are in the PDB
hr :: FilePath -> IO [Helix]
hr s = do
    -- Open the file
    pdb <- B.readFile s
    -- Split the file on line breaks
    let split = B.splitWith (==10) pdb
    return $ mapMaybe ht split
    
    
-- Get all sheet strips which are in the PDB
sr :: FilePath -> IO [Sheet]
sr s = do
    -- Open the file
    pdb <- B.readFile s
    -- Split the file on line breaks
    let split = B.splitWith (==10) pdb
    return $ mapMaybe st split

-- Get all atoms which are inside a secondary structure sheet
getSheets :: Protein -> [[Atom]]
getSheets = undefined

-- Get all atoms which are inside a secondary structure loop
getLoops :: Protein -> [[Atom]]
getLoops = undefined

-- TODO: Define this such that it takes a PDB file as its argument and
--  flags to perform operations on the PDB file
main :: IO ()
main = do
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
    