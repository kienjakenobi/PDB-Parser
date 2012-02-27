module MathDataStructures where

import Data.Packed.Vector

-- A line defined by two points
data Line2P = Line2P
    { lp1 :: !(Vector Double) -- 3-Vector (x1, y1, z2)
    , lp2 :: !(Vector Double) -- 3-Vector (x2, y2, z2)
    }
    deriving (Show, Eq, Ord)

-- A plane defined by three points
data Plane3P = Plane3P
    { pp1 :: !(Vector Double)   -- 3-Vector (x1, y1, z1)
    , pp2 :: !(Vector Double)   -- 3-Vector (x2, y2, z2)
    , pp3 :: !(Vector Double)   -- 3-Vector (z3, y3, z3)
    }
    deriving (Show)