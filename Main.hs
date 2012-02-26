-- Main.hs
-- main function to execute for PDB Parser

module Main where

-- Custom libraries:
import PDBDataStructures
--import PRLength
--import PRCaDist
import Meta
import PRAromatics
import PDBParsing

-- TODO: Remove all instances of head and fromJust via pattern matching

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

-- TODO: Define this such that it takes a PDB file as its argument and
--  flags to perform operations on the PDB file, then the program will
--  be ready to be used in compiled form
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
