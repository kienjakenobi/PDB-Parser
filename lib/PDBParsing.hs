-- PDBParsing.hs
-- Functions to take PDB files and output usable data structures

module PDBParsing (pr, folderFiles) where

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC -- unpack, readInt
import Data.ByteString.Lex.Double -- readDouble, from the bytestring-lexing module
import Data.Maybe -- fromJust
import Numeric.Container
import System.Directory

-- Custom libraries
import PDBDataStructures


--Read a ByteString into an Int
ri :: B.ByteString -> Int
ri bs = fst $ fromJust $ BC.readInt bs

--Read a ByteString into a Double
rd :: B.ByteString -> Double
rd bs = fst $ fromJust $ readDouble bs

-- Take substring of a larger ByteString
-- Take start index a and end index b, indices starting at 1
-- Also remove spaces
get :: Int -> Int -> B.ByteString -> B.ByteString
get a b str = B.filter (/=32) $ B.drop (a - 1) $ B.take b str

-- TODO: Get HETATM's
-- atomize
-- This is the actual parser
-- Take a single line of 
at :: B.ByteString -> Atom
at str = Atom 
        {atomRef = serial
        , atomName = name
        , atomAlt = altloc
        , atomAmino = amino
        , atomChain = chain
        , atomRes = resNum
        , atomPos = fromList [x,y, z]
        , atomSymbol = symbol} where
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

-- TODO: Clean this up.  It is a mess.  What are n-tuples doing here?
-- Take a list of PDB file lines in ByteString format and iterate through
--  each one.  For each, determine 
parse :: [B.ByteString] -> [Atom] -> ([Model], [(Int, Int)], [(Int, Int)], B.ByteString) -> ([Model], [(Int, Int)], [(Int, Int)], B.ByteString)
parse [] curr_atoms (models, sheets, helices, experiment) = (models ++ [curr_atoms], sheets, helices, experiment)
parse (b:bs) curr_atoms (models, sheets, helices, experiment)
    -- Turn an ATOM record into an Atom data type
    | recName == atom = parse bs (curr_atoms ++ [at b]) (models, sheets, helices, experiment)
    | recName == model && curr_atoms /= [] = parse bs [] (models ++ [curr_atoms], sheets, helices, experiment)
    | recName == helix = parse bs curr_atoms (models, sheets, helices ++ [ht b], experiment)
    | recName == sheet = parse bs curr_atoms (models, sheets ++ [st b], helices, experiment)
    | otherwise = parse bs curr_atoms (models, sheets, helices, experiment) where
    recName = get 1 6 b
    atom = B.pack [65,84,79,77]
    model = B.pack [77, 79, 68, 69, 76]
    helix = B.pack [72,69,76,73,88]
    sheet = B.pack [83,72,69,69,84]

-- Turn a HELIX entry line into a Helix data structure
ht :: B.ByteString -> (Int, Int)
ht str = (startLoc, endLoc) where
    --Attributes of the sheet:
    startLoc = ri $ get 22 25 str
    endLoc = ri $ get 34 37 str

-- Get all the strips of sheets out of a PDB file
st :: B.ByteString -> (Int, Int)
st str = (startLoc, endLoc) where
    --Attributes of the sheet:
    startLoc = ri $ get 23 26 str
    endLoc = ri $ get 34 37 str

-- proteinize
-- Perform the listize, furtherListize, and atoms functions on a file
pr :: FilePath -> IO Protein
pr s = do
    pdb <- B.readFile s
    -- Split the file on line breaks
    let split = B.splitWith (==10) pdb
    -- Parse the lines
    let parsed = parse split [] ([], [], [], B.pack [])
    let (models, sheets, helices, experiment) = parsed
    let prot = Protein { proMethod = experiment
        , proModels = models
        , proSheets = sheets
        , proHelices = helices
        , fullPDB = True
        }
    return prot
    
-- Take a string that is a path to a folder and return a list
--  of complete paths to all files in that folder
-- Make sure that the given file path ends in '/'
folderFiles :: FilePath -> IO [FilePath]
folderFiles s = do
    -- Get all the file names in the directory
    filenames <- getDirectoryContents s
    -- Remove hidden files
    let noHiddens = [x | x <- filenames, head x /= '.']
    -- Append the path to each filename
    let fullPaths = map (s ++ ) noHiddens
    return fullPaths
