#!/usr/bin/perl 
#
# written by Sung-Hun Bae Oct. 19, 2010
#

$debug = 0;

if ($#ARGV != 2) {
  printf("\n\tusage: copyxyz SOURCE_PDB SOURCE_PDB_OFFSET TARGET_PDB\n");
  printf("\tCopy XYZ coordinates from source to target\n");
  printf("\n\n");
  exit;
  }

# conversion table
%conversion = (
	" OP1" => " O1P" ,
	" OP2" => " O2P" , 
	"1H2 " => " H21" , "2H2 " => " H22" ,
	"1H4 " => " H41" , "2H4 " => " H42" ,
	"1H6 " => " H61" , "2H6 " => " H62" ,
	"1H5'" => " H5'" , "2H5'" => "H5''" , 
	"HO2'" => " H2'" , " H2'" => "H2''");


@PDB_SOURCE = ();

# store SOURCE PDB
$filename = $ARGV[0];
$offset = $ARGV[1];
open(pdbfile, $filename) or die "cannot open $filename";
	foreach $line (<pdbfile>) {
	if ($line =~ m/^ATOM/) {
    	chomp($line);
    	$name 							= substr($line,12,4);
    	$resSeq 						= substr($line,22,4);
		# apply residue number offset
		$resSeq = 0 + $resSeq + $offset;
		# convert atom name 
		if (defined $conversion{$name}) {
			$name = $conversion{$name};
			}
	    $PDB_SOURCE{$resSeq}{$name}{x}	= substr($line,30,8);
	    $PDB_SOURCE{$resSeq}{$name}{y}	= substr($line,38,8);
	    $PDB_SOURCE{$resSeq}{$name}{z}	= substr($line,46,8);
		}
	}
close(pdbfile);

$filename = $ARGV[2];
open(pdbfile, $filename) or die "cannot open $filename";
foreach $line (<pdbfile>) {
	if ($line =~ m/^ATOM/) {
    	chomp($line);
		# PDB format version 2.3
		$serial =  substr($line,6,5);
		$name = substr($line,12,4);
		$altLoc = substr($line,16,1);
		$resName = substr($line,17,3);
		$chainId = substr($line,21,1);
		$resSeq = substr($line,22,4) + 0;
		$iCode = substr($line,26,1);
		if (defined $PDB_SOURCE{$resSeq}{$name}{x} and
			defined $PDB_SOURCE{$resSeq}{$name}{y} and
			defined $PDB_SOURCE{$resSeq}{$name}{z}) {
			$x = $PDB_SOURCE{$resSeq}{$name}{x};
			$y = $PDB_SOURCE{$resSeq}{$name}{y};
			$z = $PDB_SOURCE{$resSeq}{$name}{z};
			}
		else {
			if ($debug) {
				printf ("REMARK WARNING : $resSeq $name not found in source\n",
					$resSeq, $name);
				}
	    	$x = substr($line,30,8);
	    	$y = substr($line,38,8);
	    	$z = substr($line,46,8) + 100.0;
			}
		$occupancy = substr($line,54,6);
		$tempFactor = substr($line,60,6);
		$segId = substr($line,72,4);
		$element = substr($line,76,2);
		$charge = substr($line,78,2);
		# output (according to PDB version 2.3)
		printf("ATOM  %5d %-4s%1s%3s ",$serial,$name,$altLoc,$resName);
		printf("%1s%4d%1s   ",$chainId,$resSeq,$iCode);
		printf("%8.3f",$x);
		printf("%8.3f",$y);
		printf("%8.3f",$z);
		printf("%6.2f%6.2f      ",$occupancy,$tempFactor);
		printf("%4s%2s%2s\n",$segId,$element,$charge);
		}
	else {
    	# keep other information
		print $line;
		}
  	}
