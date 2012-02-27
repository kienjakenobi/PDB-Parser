-- PRAromatics.hs
-- Functions to deal with the aromatic rings in Phe, Tyr, Trp

module PRAromatics (mainH, sideH, getAromatics, hBondMFinder, 
    hBondSFinder) where

import qualified Data.ByteString as B
import Data.Maybe

import Numeric.Container

-- Custom libraries
import Math
import PDBDataStructures
import MathDataStructures
import Meta

-- TODO: Think of a clever way to 
-- Take any residue and return a list of the Hydrogens which should be
--  attached to the mainchain.  Predict location of Hydrogens.
mainH :: Residue -> Maybe Atom
mainH res
    | nits /= [] && oxys /= [] && cars /= [] = Just newH
    | otherwise = Nothing where
    ats = resAtoms res
    a = head ats
    newH = Atom 
        { atomRef = 0
        , atomName = B.pack [72, 65]
        , atomAlt = atomAlt a
        , atomAmino = resName res
        , atomChain = atomChain a
        , atomRes = atomRes a
        , atomPos = hpos
        , atomSymbol = B.pack [72]
        }
    -- Get the atom positions required by the equation
    nits = getAtomName ats "N"
    cars = getAtomName ats "C"
    oxys = getAtomName ats "O"
    nit = atomPos $ head nits
    car = atomPos $ head cars
    oxy = atomPos $ head oxys
    covec = car `sub` oxy
    inorm = 1 / norm2 covec
    hpos = nit `add` (inorm `multC` covec)

-- Take any residue and return a list of the Hydrogens which should be
--  attached to the sidechains.  Predict location of Hydrogens.
--  This is intended to produce the Hydrogens Asn and Gln
sideH :: Residue -> Maybe (Atom, Atom)
sideH res
    -- Make sure that there are the necessary atoms
    | nits /= [] && cars /= [] && oxys /= [] = Just (newH1, newH2)
    | otherwise = Nothing where
    (a:as) = resAtoms res
    newH1 = Atom 
        { atomRef = 0
        , atomName = B.pack [72, 65]
        , atomAlt = atomAlt a
        , atomAmino = resName res
        , atomChain = atomChain a
        , atomRes = atomRes a
        , atomPos = h1pos
        , atomSymbol = B.pack [72]
        }
    newH2 = Atom 
        { atomRef = 0
        , atomName = B.pack [72, 65]
        , atomAlt = atomAlt a
        , atomAmino = resName res
        , atomChain = atomChain a
        , atomRes = atomRes a
        , atomPos = h2pos
        , atomSymbol = B.pack [72]
        }
    -- Get the atom positions required by the equation
    nits = getAtomName (a:as) "N"
    cars = getAtomName (a:as) "C"
    oxys = getAtomName (a:as) "O"
    nit = atomPos $ head nits
    car = atomPos $ head cars
    oxy = atomPos $ head oxys
    covec = car `sub` oxy
    ocvec = oxy `sub` car
    ncvec = nit `sub` car
    inorco = 1 / norm2 covec
    inornc = 1 / norm2 ncvec
    inoroc = 1 / norm2 ocvec
    -- Calculate the hydrogen positions
    h1pos = nit `add` (inorco `multC` covec)
    h2pos = nit `add` ( 0.5 `multC` ((inoroc`multC` ocvec) `add` (inornc `multC` ncvec)))

-- Take a residue and a list of atoms out of which to compose the aromatic ring
--  and return the corresponding aromatic ring
--  Test if the expected atoms actually exist, return Nothing if they don't
makeAromatic :: Residue -> (String, String, String, String) -> Maybe Aromatic
makeAromatic r (s1, s2, s3, s4) = case (c1, c2, c3, c4) of
    -- Test that all the atoms exist in the residue
    (_:_, _:_, _:_, _:_) -> Just aro 
    _ -> Nothing
    where
    -- Get the atoms
    c1 = getAtomName ats s1
    c2 = getAtomName ats s2
    c3 = getAtomName ats s3
    c4 = getAtomName ats s4
    -- Get atom positions to define lines and plane
    c1p = atomPos $ head c1
    c2p = atomPos $ head c2
    c3p = atomPos $ head c3
    c4p = atomPos $ head c4
    -- Get the atoms from the residue
    ats = resAtoms r
    -- Create the output aromatic
    aro = Aromatic
        { aroResName = resName r
        , aroResID = resID r
        , aroCenter = appLineLineInter Line2P {lp1 = c1p, lp2 = c2p} Line2P {lp1 = c3p, lp2 = c4p}
        , aroPlane = Plane3P {pp1 = c1p, pp2 = c2p, pp3 = c3p}
        }

-- Take a residue and return the aromatic ring contains in it, provided it contains
--  an aromatic ring at all (PHE, TYR, and TRP)
getAromatic :: Residue -> Maybe Aromatic
getAromatic r
    -- Get aromatics for PHE and TYR
    | resName r `elem` [phe, tyr] = case makeAromatic r ("CG", "CZ", "CE1", "CD2") of
        -- Add to list if something is returned, otherwise skip
        Just aro -> Just aro
        _ -> Nothing
    -- Get aromatics for TRP
    | resName r == trp = case makeAromatic r ("CD2", "CH2", "CZ2", "CE3") of
        -- Add to list if something is returned, otherwise skip
        Just aro -> Just aro
        _ -> Nothing
    -- Pass everything which does not have any aromatics to get
    | otherwise = Nothing where
    -- Define the residue names as ByteStrings
    tyr = B.pack [84, 89, 82]
    phe = B.pack [80, 72, 69]
    trp = B.pack [84, 82, 80]

-- Map getAromatics over a list of Residues
getAromatics :: [Residue] -> Maybe [Aromatic]
getAromatics resl = case x of 
    y:ys -> Just (y:ys)
    _ -> Nothing
    where
    x = mapMaybe getAromatic resl

-- Given a list of residues and an aromatic, find instances where there
--  is a Hydrogen bond between the aromatic and the Hydrogen on the
--  mainchain of any residue
hBondMFinder :: [Residue] -> Aromatic -> [(Double, Double)]
hBondMFinder [] _ = []
hBondMFinder (r:rs) aro
    -- Requirements: Hydrogen placement succeeds, 
    --  H in the N-H vector is closer to the aromatic than the N, 
    --  R is less than 7.0 and r is less than 6.0
    | isJust newH && backboneNPos `dist` inter > newHPos `dist` inter &&
        r1 <= 6.0 && bR1 <= 7.0  = 
        (r1, bR1) : hBondMFinder rs aro
    | otherwise = hBondMFinder rs aro where
    newH = mainH r
    -- Get the location of the Hydrogen on the mainchain
    newHPos = atomPos $ fromJust newH
    -- Get the location of the Nitrogen on the mainchain of the Residue
    ats = resAtoms r
    backboneNPos = atomPos $ head $ getAtomName ats "N"
    hNVect = Line2P {lp1 = newHPos, lp2 = backboneNPos}
    inter = linePlaneInter (aroPlane aro) hNVect
    -- Center position of the aromatic
    cp = aroCenter aro
    bR1 = cp `dist` newHPos
    r1 = inter `dist` cp

-- Given a list of residues and an aromatic, find instances where there
--  is a Hydrogen bond between the aromatic and the Hydrogens on Gln or Asn
hBondSFinder :: [Residue] -> Aromatic -> [(Double, Double)]
hBondSFinder [] _ = []
hBondSFinder (r:rs) aro
    -- Requiments: Hydrogen placement prediction must have succeeded, 
    -- H in the N-H vector must be closer to the aromatic then the N, 
    --  the residue being considered must be GLN or ASN,
    --  r must be less than 6.0, and R must be less than 7.0
    -- Case where both sidechain Hydrogens appears to be bonding
    | commonBool && h1Bool && h2Bool = 
        [h1sol, h2sol] ++ hBondSFinder rs aro
    -- Case where only h1 appears to be bonding
    | commonBool && h1Bool = 
         h1sol : hBondSFinder rs aro
    -- Case where only h2 appears to be bonding
    | commonBool && h2Bool = 
         h2sol : hBondSFinder rs aro
    | otherwise = hBondSFinder rs aro where
    -- Tests which should pass for any situation
    commonBool = (resName r == gln || resName r == asn) && 
        isJust newHs
    h1Bool = backboneNPos `dist` inter1 > newH1Pos `dist` inter1 &&
        r1 <= 6.0 && bR1 <= 7.0
    h2Bool = backboneNPos `dist` inter2 > newH2Pos `dist` inter2 &&
        r2 <= 6.0 && bR2 <= 7.0
    -- Sidechain identifying strings
    gln = B.pack [71, 76, 78]
    asn = B.pack [65, 83, 78]
    -- Get the location of the Hydrogens on the residue's sidechain
    newHs = sideH r
    newH1Pos = atomPos $ fst $ fromJust newHs
    newH2Pos = atomPos $ snd $ fromJust newHs
    -- Get the location of the Nitrogen on the mainchain of the Residue
    ats = resAtoms r
    backboneNPos = atomPos $ head $ getAtomName ats "N"
    -- Lines formed by the N-H vectors
    nHVect1 = Line2P {lp1 = backboneNPos, lp2 = newH1Pos}
    nHVect2 = Line2P {lp1 = backboneNPos, lp2 = newH2Pos}
    -- Points where the N-H vector intersects with the plane of the aromatic
    inter1 = linePlaneInter (aroPlane aro) nHVect1
    inter2 = linePlaneInter (aroPlane aro) nHVect2
    -- Calculate R
    bR1 = cp `dist` newH1Pos
    bR2 = cp `dist` newH2Pos
    -- Calculate r
    cp = aroCenter aro
    r1 = inter1 `dist` cp
    r2 = inter2 `dist` cp
    -- Data points
    h1sol = (r1, bR1)
    h2sol = (r2, bR2)
