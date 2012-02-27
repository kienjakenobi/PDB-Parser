--Exercise 12.3: Catalog hydrogen bonds having aromatics as 
--  an acceptor in the PDB. Determine the distribution of distances
--  from the donor to the face of the other aromatic, as well as the
--  distance from the intersection of the donor's N-H vector with the plane
--  and the center of the aromatic ring.
--  Do this for all mainchain Hydrogens and for two sidechains Hydrogens
--  of GLN and ASN.

-- Used to test the output values at every step
main :: IO ()
main  = do
    -- Get the file
    p <- pr "/Users/mortal/Desktop/2ACX.pdb"
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    let resAs = groupRes theas
    
    -- Get the interesting aromatic
    let aroRes = fromJust $ getResID resAs 173
    let aro = head $ getAromatics [aroRes]
    print aro
    -- Get the interesting donor residue
    let donRes = fromJust $ getResID resAs 525
    
    let newHs = sideH donRes
    
    print $ fromJust $ newHs
    
    print $ hBondSFinder [donRes] aro
    
    
    -- Print the points on the plane where the NH vectors intersect it
    let ats = resAtoms donRes
    let backboneNPos = atomPos $ head $ getAtomName ats "N"
    let newH1Pos = atomPos $ fst $ fromJust newHs
    let newH2Pos = atomPos $ snd $ fromJust newHs
    let nHVect1 = Line2P {lp1 = backboneNPos, lp2 = newH1Pos}
    let nHVect2 = Line2P {lp1 = backboneNPos, lp2 = newH2Pos}
    let inter1 = linePlaneInter (aroPlane aro) nHVect1
    let inter2 = linePlaneInter (aroPlane aro) nHVect2
    print $ inter1
    print $ inter2
    
    print $ show (backboneNPos `dist` inter1) ++ " vs. " ++ show (newH1Pos `dist` inter1)
    print $ show (backboneNPos `dist` inter2) ++ " vs. " ++ show (newH2Pos `dist` inter2)

-- These are all of the helper functions which are called by main()

pheMain :: FilePath -> FilePath -> IO ()
pheMain dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    -- Sort the entire list of proteins into residues
    let resAs = groupRes theas

    -- Get the Phe residues
    let phe = getAmino theas "PHE"
    -- Group the Phe atoms by residue
    let resPhe = groupRes phe
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resPhe
    -- Get all HBond situations
    let hbonds = map (hBondMFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)

pheSide :: FilePath -> FilePath -> IO ()
pheSide dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    -- Sort the entire list of proteins into residues
    let resAs = groupRes theas

    -- Get the Phe residues
    let phe = getAmino theas "PHE"
    -- Group the Phe atoms by residue
    let resPhe = groupRes phe
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resPhe
    -- Get all HBond situations
    let hbonds = map (hBondSFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)
    
tyrMain :: FilePath -> FilePath -> IO ()
tyrMain dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    -- Sort the entire list of proteins into residues
    let resAs = groupRes theas

    -- Get the Phe residues
    let tyr = getAmino theas "TYR"
    -- Group the Phe atoms by residue
    let resTyr = groupRes tyr
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resTyr
    -- Get all HBond situations
    let hbonds = map (hBondMFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)

tyrSide :: FilePath -> FilePath -> IO ()
tyrSide dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    -- Sort the entire list of proteins into residues
    let resAs = groupRes theas

    -- Get the Phe residues
    let tyr = getAmino theas "TYR"
    -- Group the Phe atoms by residue
    let resTyr = groupRes tyr
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resTyr
    -- Get all HBond situations
    let hbonds = map (hBondSFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)

trpMain :: FilePath -> FilePath -> IO ()
trpMain dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    -- Sort the entire list of proteins into residues
    let resAs = groupRes theas

    -- Get the Phe residues
    let trp = getAmino theas "TRP"
    -- Group the Phe atoms by residue
    let resTrp = groupRes trp
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resTrp
    -- Get all HBond situations
    let hbonds = map (hBondMFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)

trpSide :: FilePath -> FilePath -> IO ()
trpSide dest s = do
    -- Status report
    print $ "Working on " ++ s
    -- Get the file
    p <- pr s
    -- Get only the A altlocs
    let theas = getAlt (head $ proModels p) ['A', ' ']
    let resAs = groupRes theas

    -- Get the TRP residues
    let trp = getAmino theas "TRP"
    -- Group the TRP atoms by residue
    let resTrp = groupRes trp
    -- Get a list of all the corresponding aromatics
    let aromatics = getAromatics resTrp
    -- Get all HBond situations
    let hbonds = map (hBondSFinder resAs) aromatics
    
    -- Only save those which have content
    case (concat hbonds) of
        [] -> print $ "None found in " ++ s
        _ -> appendFile dest (show $ concat hbonds)

main :: IO ()
main = do
    -- Get a list of all PDB files to iterate over
    pdbs <- folderFiles "/Users/mortal/Downloads/nrpdb/"
    -- Set the location to save data
    let dest = "/Users/mortal/Desktop/pheMain_R=ToCenter.txt"
    -- Make sure the destination file is empty
    writeFile dest ""
    -- Map the function across all the PDB files, throwing out the result,
    --  which is printed
    mapM_ (pheMain dest) pdbs