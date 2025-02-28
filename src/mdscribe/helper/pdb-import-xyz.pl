#!/usr/bin/perl 
#
# written by Sung-Hun Bae Aug. 15, 2008
#
# Import XYZ coordinates from RCSB or CYANA format PDB files
# Atom names and other PDB contents of the XPLOR_PDB remains intact
# Some atoms may not be found in the imported PDB files
# In that case, the coordinates of those atoms will be assigned to (0,0,0)
# These missing atoms can be treated by separate XPLOR script
# opt_cyana.inp or opt_rcsb.inp 
#
# verbose option provides details of importing process

if ($#ARGV != 2) {
  printf("\n\tusage: importxyz [[r|c]v] XPLOR_PDB RCSB_or_CYANA_PDB\n");
  printf("\tImport XYZ coordinates from RCSB or CYANA format for XPLOR\n");
  printf("\n\t\tr or c : RCSB or CYANA format");
  printf("\n\t\tv : Verbose [option]");
  printf("\n\n");
  exit;
  }

$import_rcsb = 0;
$import_cyana = 0;
$verbose = 0;

if ($ARGV[0] =~ m/c/) {
  $import_cyana = 1;
  }
if ($ARGV[0] =~ m/r/) {
  $import_rcsb = 1;
  }
if ($ARGV[0] =~ m/v/) {
  $verbose = 1;
  }
if ($import_rcsb and $import_cyana) {
  printf("error: choose RCSB or CYANA\n");
  exit;
  }

# converting RCSB to XPLOR
# string length should be 4
# second character is atom code [H/N/C/O/S] 
# Note: XPLOR has additional atoms
#   HT1,HT2,HT3 @ N terminal
#   OT1,OT2 @ C terminal

%rcsb = (
  " H  " => " HN " , 
  "1HA " => " HA1" , "2HA " => " HA2" ,
  "1HB " => " HB1" , "2HB " => " HB2" , "3HB " => " HB3" , 
  "1HG " => " HG1" , "2HG " => " HG2" , 
  "1HG1" => "HG11" , "2HG1" => "HG12" , "3HG1" => "HG13" ,
  "1HG2" => "HG21" , "2HG2" => "HG22" , "3HG2" => "HG23" ,
  "1HD " => " HD1" , "2HD " => " HD2" , 
  "1HD1" => "HD11" , "2HD1" => "HD12" , "3HD1" => "HD13" ,
  "1HD2" => "HD21" , "2HD2" => "HD22" , "3HD2" => "HD23" ,
  "1HE " => " HE1" , "2HE " => " HE2" , "3HE " => " HE3" , 
  "1HE2" => "HE21" , "2HE2" => "HE22" ,
  "1HH1" => "HH11" , "2HH1" => "HH12" , "1HH2" => "HH21" , "2HH2" => "HH22" 
  );

# converting CYANA to XPLOR
# string length should be 4
# second character is atom code [H/N/C/O/S] 
# Note: XPLOR has additional atoms
#   HT1,HT2,HT3 @ N terminal
#   OT1,OT2 @ C terminal
#   HE2 @ histidine

%cyana = (
  THR => {
    " H  " => " HN " , 
    "1HG " => "HG1 " ,
    "1HG2" => "HG21" , "2HG2" => "HG22" , "3HG2" => "HG23" ,
    },
  LYS => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    "2HD " => " HD1" , "3HD " => " HD2" , 
    " HD2" => " HD1" , " HD3" => " HD2" , 
    "2HE " => " HE1" , "3HE " => " HE2" , 
    " HE2" => " HE1" , " HE3" => " HE2" , 
    "1HZ " => " HZ1" , "2HZ " => " HZ2" , "3HZ " => " HZ3",
    },
  GLY => {
    " H  " => " HN " , 
    "2HA " => " HA1" , "3HA " => " HA2" ,
    " HA2" => " HA1" , " HA3" => " HA2" ,
    },
  ALA => {
    " H  " => " HN " , 
    "1HB " => " HB1" , "2HB " => " HB2" , "3HB " => " HB3" ,
    },
  VAL => {
    " H  " => " HN " , 
    "1HB " => " HB1" , "2HB " => " HB2" , "3HB " => " HB3" ,
    "1HG1" => "HG11" , "2HG1" => "HG12" , "3HG1" => "HG13" ,
    "1HG2" => "HG21" , "2HG2" => "HG22" , "3HG2" => "HG23" ,
    },
  LEU => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    "2HG1" => "HG11" , "3HG1" => "HG12" ,
    "1HG2" => "HG21" , "2HG2" => "HG22" , "3HG2" => "HG23" ,
    "1HD1" => "HD11" , "2HD1" => "HD12" , "3HD1" => "HD13" ,
    "1HD2" => "HD21" , "2HD2" => "HD22" , "3HD2" => "HD23" ,
    },
  ILE => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    "2HG1" => "HG11" , "3HG1" => "HG12" ,
    "1HG2" => "HG21" , "2HG2" => "HG22" , "3HG2" => "HG23" ,
    "1HD1" => "HD11" , "2HD1" => "HD12" , "3HD1" => "HD13" ,
    },
  SER => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  CYS => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  MET => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    "1HE " => " HE1" , "2HE " => " HE2" , "3HE " => " HE3" , 
    },
  ASP => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" ,
    " HB2" => " HB1" , " HB3" => " HB2" ,
    },
  ASN => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" ,
    " HB2" => " HB1" , " HB3" => " HB2" ,
    "1HD2" => "HD21" , "2HD2" => "HD22" ,
    },
  GLN => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" ,
    " HB2" => " HB1" , " HB3" => " HB2" ,
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    "1HE2" => "HE21" , "2HE2" => "HE22" , 
    },
  HIS => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  GLU => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" ,
    " HB2" => " HB1" , " HB3" => " HB2" ,
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    },
  ARG => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" ,
    " HB2" => " HB1" , " HB3" => " HB2" ,
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    "2HD " => " HD1" , "3HD " => " HD2" , 
    " HD2" => " HD1" , " HD3" => " HD2" , 
    "1HH1" => "HH11" , "2HH1" => "HH12" , 
    "1HH2" => "HH21" , "2HH2" => "HH22" , 
    },
  PHE => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  TYR => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  TRP => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    },
  PRO => {
    " H  " => " HN " , 
    "2HB " => " HB1" , "3HB " => " HB2" , 
    " HB2" => " HB1" , " HB3" => " HB2" , 
    "2HG " => " HG1" , "3HG " => " HG2" , 
    " HG2" => " HG1" , " HG3" => " HG2" , 
    "2HD " => " HD1" , "3HD " => " HD2" , 
    " HD2" => " HD1" , " HD3" => " HD2" , 
    },
  );

open(INPUT_A,$ARGV[1]) or die "cannot open $ARGV[1]";
@PDB_A = <INPUT_A>;
close(INPUT_A);

open(INPUT_B,$ARGV[2]) or die "cannot open $ARGV[2]";
@PDB_B = <INPUT_B>;
close(INPUT_B);

@PDB = ();
@NOTE = ();

foreach $line_a (@PDB_A) {
  if ($line_a =~ m/^ATOM/) {
    chomp($line_a);
    $serial_a =  substr($line_a,6,5);
    $name_a = substr($line_a,12,4);
    $resName_a = substr($line_a,17,3);
    $resSeq_a = substr($line_a,22,4);

    if ($verbose) {
      printf("REMARK %s: %s %s %s <- ",$ARGV[1],$resName_a,$resSeq_a,$name_a);
      }

    # search a matching atom in PDB_B -----------------------------------v
    foreach $line_b (@PDB_B) {
      if ($line_b =~ m/^ATOM/) {
	chomp($line_b);
	$serial_b =  substr($line_b,6,5);
	$name_b = substr($line_b,12,4);
	$resName_b = substr($line_b,17,3);
	$resSeq_b = substr($line_b,22,4);
	if ($resName_a =~ $resName_b and $resSeq_a =~ $resSeq_b) {
	  # convert
	  if ($import_rcsb) {
	    $converted = $rcsb{$name_b};
	    }
	  if ($import_cyana) {
	    $converted = $cyana{$resName_b}{$name_b};
	    }
	  if ($name_a =~ $name_b or $name_a =~ $converted) {
	    if ($verbose) {
	      printf("%s: %s %s %s",$ARGV[2],$resName_b,$resSeq_b,$name_b);
	      }
	    $NOTE{$serial_b} = $name_a;
	    $PDB{$serial_a}{x} = substr($line_b,30,8);
	    $PDB{$serial_a}{y} = substr($line_b,38,8);
	    $PDB{$serial_a}{z} = substr($line_b,46,8);
	    }
	  }
	}
      }
    # search a matching atom in PDB_B -----------------------------------^
    if ($verbose) {
      printf("\n");
      }
    }
  }

foreach $line_a (@PDB_A) {
  # modify ATOM lines
  if ($line_a =~ m/^ATOM/) {
    chomp($line_a);
    # PDB format version 2.3
    $serial =  substr($line_a,6,5);
    $name = substr($line_a,12,4);
    $altLoc = substr($line_a,16,1);
    $resName = substr($line_a,17,3);
    $chainId = substr($line_a,21,1);
    $resSeq = substr($line_a,22,4);
    $iCode = substr($line_a,26,1);
    # XYZ coordinate from PDB_B
    # XYZ coordinate from PDB_B
    # XYZ coordinate from PDB_B
    $occupancy = substr($line_a,54,6);
    $tempFactor = substr($line_a,60,6);
    $segId = substr($line_a,72,4);
    $element = substr($line_a,76,2);
    $charge = substr($line_a,78,2);
    # output (according to PDB version 2.3)
    printf("ATOM  %5d %-4s%1s%3s ",$serial,$name,$altLoc,$resName);
    printf("%1s%4d%1s   ",$chainId,$resSeq,$iCode);
    printf("%8.3f",$PDB{$serial}{x});
    printf("%8.3f",$PDB{$serial}{y});
    printf("%8.3f",$PDB{$serial}{z});
    printf("%6.2f%6.2f      ",$occupancy,$tempFactor);
    printf("%4s%2s%2s\n",$segId,$element,$charge);
    }
  else {
    # keep other information
    print $line_a;
    }
  }
