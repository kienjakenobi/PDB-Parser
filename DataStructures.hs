-- Description: Definition of the data structures used for storing PDB atom data

module DataStructures where

import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.ByteString.Internal
import Data.Word

-- TODO: Create a data structure to store information by amino acid

-- NOTE: Any time you have more than a 2-tuple, rethink how you are storing
--  data and instead consider using a record

-- TODO: Instead of defining teh Point3D data type, use a more general vector/matrix
--      data type so that you can use a general linear algebra library to do math with
--      the atom points

--  data structure to contain a 3D point in space
data Point3D = Point3D !Double !Double !Double
    deriving (Eq, Ord)

instance Show Point3D where
    show p = "(" ++ show(x) ++ ", " ++ show (y) ++ ", " ++ show (z) ++ ")" where
        Point3D x y z = p

-- TODO: May need to take into consideration Models in the future.  For example, 
--      1CFC has various models within its PDB file

instance Num Point3D where
    -- Addition
    Point3D x1 y1 z1 + Point3D x2 y2 z2 = Point3D (x1 + x2) (y1 + y2) (z1 + z2)
    
    -- Subtraction
    Point3D x1 y1 z1 - Point3D x2 y2 z2 = Point3D (x1 - x2) (y1 - y2) (z1 - z2)
    
    -- Multiplication, scalar == Dot product
    --Point3D x1 y1 z1 * Point3D x2 y2 z2 = x1*x2 + y1*y2 + z1*z2
    
    -- Multiplication, cross
    Point3D x1 y1 z1 * Point3D x2 y2 z2 = Point3D (y1*z2 - z1*y2) (z1*x2 - x1*z2) (x1*y2 - y1*x2)
    
    -- Absolute Value == l2 Norm == Euclidean Norm
    --abs (Point3D x y z) = 1
    
    -- Cross Product
    --Point3D x1 y1 z1 x Point3D x2 y2 z2 = Point3D (y1*z2 - z1*y2) (y1 - y2) (z1 - z2)

-- Record data structure for storing data on each atom
data Atom = Atom 
    { atomRef :: !Int           -- atom reference number: 1, 2, 3, 4, ...
    , atomName :: B.ByteString  -- atom name: CA, CB, CG, ...
    , atomAlt :: Word8          -- atom alternate location: A, B, C, ...
    , atomAmino :: B.ByteString -- amino acid that the atom belongs to: ALA, ARG, ...
    , atomChain :: Word8        -- chain of the protein that the atom belongs to: A, B, C, ...
    , atomRes :: !Int           -- Residue that the atom belongs to: 1, 1, 1, 2, ...
    , atomPos :: !Point3D       -- x, y, z coordinates of the atom's position
    , atomSymbol :: B.ByteString-- Symbol of the atom's element: C, N, O, ...
    } deriving (Eq, Ord)

-- Must convert all values to print-able strings
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
        Point3D x y z = atomPos a

data Helix = Helix
    { helixResList :: (Int, Int)
    } deriving (Show, Eq, Ord)
    
data Sheet = Sheet
    { sheetResList :: (Int, Int)
    } deriving (Show, Eq, Ord)

-- TODO: Consider altering this type declaration into a data type that
--  can store helix, sheet, and loop information, as well as the list of
--  atoms

-- A protein is a list of atoms
type Protein = [Atom]