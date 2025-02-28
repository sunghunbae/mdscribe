""" merge multiple PDBs - unique serial & chain """
import sys

def read_pdb_file(filename, serial_ended=0, chainId_ended=None):
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

                if chainId_ended and (chainId not in chainId_map):
                    chainId_ended = chr(ord(chainId_ended)+1)
                    chainId_map[chainId] = chainId_ended
                
                # override serial and chainId
                serial = serial + serial_ended
                if chainId_map:
                    chainId = chainId_map[chainId]

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
                
                # override serial and chainId
                serial = serial + serial_ended
                if chainId_map:
                    chainId = chainId_map[chainId]

                print('TER   ', end='')
                print('{:5d} {:4s}{}{:3s} {}{:4d}{}'.format(serial, name, altLoc, resName, chainId, resSeq, iCode))

            if line.startswith('CONECT'):
                # keep the conect lines
                numbers = [ int(n) + serial_ended for n in line[6:].strip().split() ]
                conect_line = 'CONECT '
                for n in numbers:
                    conect_line += '{:4d}'.format(n)
                conect.append(conect_line)
            
            if not last_serial or serial > last_serial:
                last_serial = serial

            if not last_chainId or chainId > last_chainId:
                last_chainId = chainId
    
    # print("filename", filename)
    # print("chainId_map", chainId_map)
    # print("last_serial", last_serial)
    # print("last_chainId", last_chainId)
    # print()
    return (last_serial, last_chainId, conect) # last serial and chainId
    
last_serial = None
conect_lines = []
for infile in sys.argv[1:]:
    if last_serial:
        (last_serial, last_chainId, conect_) = read_pdb_file(infile, last_serial, last_chainId)
        if conect_:
            conect_lines.extend(conect_)
    else:
        (last_serial, last_chainId, conect_) = read_pdb_file(infile)
        if conect_:
            conect_lines.extend(conect_)
for line in conect_lines:
    print(line)
print('END')
