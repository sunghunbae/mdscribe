#!/usr/bin/bash

awk '
# Main Processing
{
    if ($1=="ATOM") 
    {
    # PDB format version 2.3
        serial =  substr($0,7,5);
        name = substr($0,13,4);
        altLoc = substr($0,17,1);if(altLoc=="") altLoc=" ";
        resName = substr($0,18,3);
        chainId = substr($0,22,1);if(chainId=="") chainId=" ";
        resSeq = substr($0,23,4);
        iCode = substr($0,27,1);if(iCode=="") iCode=" ";
        x = substr($0,31,8);
        y = substr($0,39,8);
        z = substr($0,47,8);
        occupancy = substr($0,55,6);
        tempFactor = substr($0,61,6);
        segId = substr($0,73,4);
        element = substr($0,77,2);
        charge = substr($0,79,2);

        segId1 = substr(segId,0,1);
	if(segId1 == "A") {
            occupancy = 1
        }
        else if (segId1 == "B") {
            occupancy = 2
        }
        else if (segId1 == "C") {
            occupancy = 3
        }

        # output (according to PDB version 2.3)
        printf("ATOM  %5d %-4s%c%3s ",serial,name,altLoc,resName);
        printf("%c%4d%c   ",chainId,resSeq,iCode);
        printf("%8.3f%8.3f%8.3f",x,y,z);
        printf("%6.2f%6.2f      ",occupancy,tempFactor);
        printf("%4s%2s%2s\n",segId,element,charge);

    } 
    else print $0; # do not touch other lines
}' $*
