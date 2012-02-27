-- Math.hs
-- A library of general mathematics functions

module Math (dist, cross, multC, appLineLineInter, linePlaneInter) where

-- Custom libraries
import MathDataStructures
import Numeric.LinearAlgebra.Algorithms

import Numeric.Container
-- scalar product: v1 <.> v2
-- l2norm: norm2 v2
-- Custom made cross product: v1 `cross` v2
-- Subtract: v1 `sub` v2
-- Add: v1 `add` v2
-- Multiply element by element: v1 `mul` v2
-- Add constant to every element: addConstant c v1
-- Custom made multiply constant to every element: c `multC` v1
-- Get the second element of a vector: @>1
-- Determinant of Matrix: det m

-- Other reminders:
-- Exponent of a double: 4.0**2 == 16.0

-- mapM for Matrix types
--mapMatrixM :: Monad m => (a -> m b) -> Matrix a -> m (Matrix b)
--mapMatrixM f m = do
--    v <- mapM f m
--    return $ liftMatrix v

-- calculate 3D Cartesian distance
-- Calculate the 3D distance between two atoms
dist :: Vector Double -> Vector Double -> Double
dist v1 v2 = sqrt ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    where
        [x1, y1, z1] = toList v1
        [x2, y2, z2] = toList v2

-- 3-Vector cross product
cross :: Vector Double -> Vector Double -> Vector Double
cross v1 v2 = fromList [y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2] where
      [x1, y1, z1] = toList v1
      [x2, y2, z2] = toList v2

-- Multiply all elements of a vector v by the constant c
multC :: Double -> Vector Double -> Vector Double
multC c = mapVector (c *)

-- Take two lines in R^3 in the two-point form and produce
--  the point which minimizes the distance between the two
--  lines.  This is used to find the approximate intersection
--  of two lines
-- Idea source: http://math.stackexchange.com/questions/61719/
--  finding-the-intersection-point-of-many-lines-in-3d-point-closest-to-all-lines
-- Note: This does not check for the chase of equal or parallel lines
appLineLineInter :: Line2P -> Line2P -> Vector Double
appLineLineInter l1 l2 = c where
    -- Get the points from the line structures
    a1 = lp1 l1
    b1 = lp2 l1
    a2 = lp1 l2
    b2 = lp2 l2
    -- Difference vectors
    d1 = b1 `sub` a1
    d2 = b2 `sub` a2
    -- Square of the l2-norm of the difference vectors
    p = norm2 d1 ** 2
    q = norm2 d2 ** 2
    -- Do the algebra
    x = appLineLineInterh (a1@>0) (a2@>0) (d1@>0) (d2@>0) p q
    y = appLineLineInterh (a1@>1) (a2@>1) (d1@>1) (d2@>1) p q
    z = appLineLineInterh (a1@>2) (a2@>2) (d1@>2) (d2@>2) p q
    -- Solution vector
    c = fromList[x, y, z]

-- Helper function which actually does the alebgra for appIntersection
appLineLineInterh :: Double -> Double -> Double -> Double -> Double -> Double -> Double
appLineLineInterh a1 a2 d1 d2 p q = num/denom where
    num = p*q*a1 + p*q*a2 - q*a1*(d1**2) - p*a2*(d2**2)
    denom = 2*p*q - q*(d1**2) - p*(d2**2)

-- Take a line and plane and give the point of their intersection
--  Note: This does not check for the case that the line is embedded in
--  the plane or parallel to the plane
-- Method from http://mathworld.wolfram.com/Line-PlaneIntersection.html
linePlaneInter :: Plane3P -> Line2P -> Vector Double
linePlaneInter pl l = c where
    -- Get the points defining the line and plane
    a1 = pp1 pl
    a2 = pp2 pl
    a3 = pp3 pl
    a4 = lp1 l
    a5 = lp2 l
    -- Calculate the parameter t
    num = (4><4) [ 1, 1, 1, 1,
                   a1@>0, a2@>0, a3@>0, a4@>0,
                   a1@>1, a2@>1, a3@>1, a4@>1, 
                   a1@>2, a2@>2, a3@>2, a4@>2 ]
    denom = (4><4) [ 1, 1, 1, 0,
                     a1@>0, a2@>0, a3@>0, a5@>0 - a4@>0,
                     a1@>1, a2@>1, a3@>1, a5@>1 - a4@>1,
                     a1@>2, a2@>2, a3@>2, a5@>2 - a4@>2 ]
    t = (-(det num)) / det denom
    -- Find solutions
    x = (a4@>0) + (a5@>0 - a4@>0)*t
    y = (a4@>1) + (a5@>1 - a4@>1)*t
    z = (a4@>2) + (a5@>2 - a4@>2)*t
    -- Solution point stored in Vector
    c = fromList [x,y,z]
