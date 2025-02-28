import sys

def rename_pdb_file(filename, resNameMap={'A':'ADE','G':'GUA','C':'CYT','U':'URA','T':'THY'}):
    # PDB format version 2.3
    last_serial = None
    last_chainId = None
    chainId_map = {}
    conect = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                serial = int(line[6:12])
                name = line[12:16].strip()
                altLoc = line[16:17]
                if altLoc=="":
                    altLoc=" "
                resName = line[17:21].strip()
                # rename!!
                if resName in resNameMap:
                    resName = resNameMap[resName]
                chainId = line[21:22]
                if chainId=="":
                    chainId="A" # default chainId if chainId is not defined
                resSeq = int(line[22:26])
                iCode = line[26:27]
                if iCode=="":
                    iCode=" "
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occupancy = float(line[54:60])
                tempFactor = float(line[60:66])
                segId = line[72:76]
                element = line[76:78]
                charge = line[78:80]

                if line.startswith('ATOM'):
                    print('ATOM  ', end='')
                else:
                    print('HETATM', end='')
                print('{:5d} {:4s}{}{:3s} {}{:4d}{}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:4s}{:2s}{:2s}'.format(
                    serial, name, altLoc, resName, 
                    chainId, resSeq, iCode, x, y, z, occupancy, tempFactor, segId, element, charge
                ))
                
            if line.startswith('TER'):
                serial = int(line[6:12])
                name = line[12:16].strip()
                altLoc = line[16:17]
                if altLoc=="":
                    altLoc=" "
                resName = line[17:21].strip()
                chainId = line[21:22]
                if chainId=="":
                    chainId=" "
                resSeq = int(line[22:26])
                iCode = line[26:27]
                if iCode=="":
                    iCode=" "
                
                print('TER   ', end='')
                print('{:5d} {:4s}{}{:3s} {}{:4d}{}'.format(serial, name, altLoc, resName, chainId, resSeq, iCode))

    
for infile in sys.argv[1:]:
    rename_pdb_file(infile)
print('END')
