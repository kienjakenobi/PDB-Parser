-- PDBDataStructures.hs
-- Definition of the data structures used for storing PDB data

module PDBDataStructures where

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.ByteString.Internal
import Data.Word
import Data.Packed.Vector

-- Custom imports
import MathDataStructures

-- Note: I think the correct method of using the new data structures is to
--  propogate the highest-level data structure for as long as possible.
--  That is, a function which takes a Protein should also output Protein.
--  This is optimal because it tells me where the data are coming from, 
--  and it maintains information such as multiple MODEL's from NMR files.

-- Record data structure for storing data on each atom
data Atom = Atom
    { atomRef       :: !Int             -- atom reference number: 1, 2, 3, 4, ...
    , atomName      :: B.ByteString     -- atom name: CA, CB, CG, ...
    , atomAlt       :: Word8            -- atom alternate location: A, B, C, ...
    , atomAmino     :: B.ByteString     -- amino acid that the atom belongs to: ALA, ARG, ...
    , atomChain     :: Word8            -- chain of the protein that the atom belongs to: A, B, C, ...
    , atomRes       :: !Int             -- Residue that the atom belongs to: 1, 1, 1, 2, ...
    , atomPos       :: !(Vector Double) -- x, y, z coordinates of the atom's position
    , atomSymbol    :: B.ByteString     -- Symbol of the atom's element: C, N, O, ...
    } deriving (Eq, Ord)

-- A MODEL is a single data set of atoms in a PDB
--  X-Ray Crystallogrophy produces 1 MODEL and NMR produces many models
type Model = [Atom]

-- A protein contians meta information from the PDB
--  and a list of all models in the file
data Protein = Protein
    { proMethod     :: B.ByteString     -- Describe the method used to collect the data (X-Ray, NMR, ...)
    , proModels     :: [Model]          -- A list of all the atoms in the protein
    , proSheets     :: ![(Int, Int)]    -- List of start and stop sequences numbers of the sheets
    , proHelices    :: ![(Int, Int)]    -- List of start and stop sequence numbers of the helices
    , fullPDB       :: Bool             -- 1 if proModels contains all atoms in the originating PDB, 0 otherwise
    } deriving (Show, Eq, Ord)

data Residue = Residue
    { resName   :: B.ByteString -- Name of the residue's amino acid type: ALA, ARG, ...
     , resID    :: !Int        -- ID of the residue: 1, 2, 3, 4, ...
     , resAtoms :: [Atom]      -- List of atoms contained in the residue
    } deriving (Eq, Ord)

-- Print Residues with just enough information to idenitfy them
instance Show Residue where
    show a = "Residue " ++
             show (resID a) ++ " " ++
             show (resName a)

-- Must convert all values to print-able strings
-- Format Atoms into PDB-like format
instance Show Atom where
    show a = show (atomRef a) ++ "\t" ++ 
             BC.unpack (atomName a) ++ "\t" ++ 
             [w2c (atomAlt a)] ++ "\t" ++
             BC.unpack (atomAmino a) ++ "\t" ++
             [w2c (atomChain a)] ++ "\t" ++
             show (atomRes a) ++ "\t" ++
             show (x) ++ "\t" ++
             show (y) ++ "\t" ++
             show (z) ++ "\t" ++
             BC.unpack (atomSymbol a) ++ "\n" where
        [x, y, z] = toList $ atomPos a

-- Data structure to store the start and stop sequence numbers of a helix
data Helix = Helix
    { helixResList :: !(Int, Int)
    } deriving (Show, Eq, Ord)
    
-- Data structure to store the start and stop sequence numbers of a sheet
data Sheet = Sheet
    { sheetResList :: !(Int, Int)
    } deriving (Show, Eq, Ord)

-- Store the information for an aromatic ring from Phe, Tyr, or Trp
data Aromatic = Aromatic
    { aroResName    :: B.ByteString     -- amino acid that the atom belongs to: ALA, ARG, ...
    , aroResID      :: !Int             -- Residue that the aromatic ring is from: 1, 2, 3, 4, ...
    , aroCenter     :: Vector Double    -- Position of the center point of the aromatic ring
    , aroPlane      :: !Plane3P         -- The ring's plane
    }

-- Print just enough to identify where the Aromtic is in the PDB
instance Show Aromatic where
    show a = "Aromatic " ++ 
             show (aroResID a) ++ " " ++
             show (aroResName a) ++ " " ++
             show (aroPlane a) ++ " " ++
             show (aroCenter a)
