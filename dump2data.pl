#!/usr/bin/perl
use Getopt::Long;
#use Data::Dumper qw(Dumper);
#Getopt::Long::Configure ('bundling');
#use Math::MatrixDecomposition::Eigen;
#use Cwd qw(abs_path);
use Cwd;

$logfile = "perl.log";
$gpu_num = 0;
#$step_step = 400000;
$step_step = 100000;
#$step_total = 26;
$step_total = 51;

$outlog = getcwd.">".$0;
for($i=0;$i<scalar(@ARGV);$i++)
{
 $outlog .= " ".$ARGV[$i];
}
$outlog .= "\n";


GetOptions (
'f|file=s' => \$in_file, 
'd|data=s' => \$data_file, 
'o|out=s' => \$outfile, 
'b|begin=i' => \$start_time, 
'e|end=i' => \$end_time,
'v|every=i' => \$every_time,
#'n|name=s' => \$col_name,
'h|help' => \$help
);

if($help) {
  print "Options for this program:
-f --file for input file default= dump.all
-d --data for data file default= data
-o --output for output file default= 'input file'_r2.dat
-b --begin start time
-e --end end time
-n --name column to mol
";
  exit;
  }

if(!col_name)
 {
 die "no column name defined";
 }

open(LOG, ">>${logfile}") or die "Cannot open ${logfile} : $!";
print LOG $outlog;
close(LOG);


if(!$in_file) {
    $in_file = "dump.all";
}

if(!$data_file) {
    $data_file = "data";
}

if(!$outfile)
 {
 $outfile = $in_file."_r2.dat";
 }

for($i=0;$i<$step_total;$i++)
 {
 $step_start[$i] = sprintf "%d", ($i+0.5)*$step_step;
 $step_end[$i]   = ($i+1)*$step_step;
 $step_show[$i] = "eps".($i+1);
 }

sub in_step
 {
 my ($st) = @_;
 my $i;
 for($i=0;$i<scalar(@step_end);$i++)
  {
#  if(($st>=$step_start[$i])&&($st<=$step_end[$i]))
   if(($st==$step_start[$i])||($st==$step_end[$i])) # not "in step" now but only on borders
   {
   return $i;
   }
  }
 return -1;
 }

sub in_step
{
return 1;
}

sub in_step
{
my $step = shift;
if((!defined($start_time)||($step>=$start_time))&&(!defined($end_time)||($step<=$end_time)))
 {
 return 1;
 }
 else
 {
 return -1;
 }
}

open DATA, "<$data_file" or die "Cannot open $data_file: $!";

$Natoms = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\d+\s+atoms\s+$/)
  {
  ($Natoms) = ($str =~ /^\s*(\d+)\s+atoms\s+$/);
  last;
  }
 }
if($Natoms==0){die "cannot find number of atoms in $data_file";}
$Nbonds = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\d+\s+bonds\s+$/)
  {
  ($Nbonds) = ($str =~ /^\s*(\d+)\s+bonds\s+$/);
  last;
  }
 }
if($Nbonds==0){die "cannot find number of bonds in $data_file";}

$Xhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+xlo\s+xhi\s+$/)
  {
  ($Xlo,$Xhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+xlo\s+xhi\s+$/);
  last;
  }
 }
if($Xhi==0){die "cannot find xhi in $data_file";}
$Yhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+ylo\s+yhi\s+$/)
  {
  ($Ylo,$Yhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+ylo\s+yhi\s+$/);
  last;
  }
 }
if($Yhi==0){die "cannot find yhi in $data_file";}
$Zhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+zlo\s+zhi\s+$/)
  {
  ($Zlo,$Zhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+zlo\s+zhi\s+$/);
  last;
  }
 }
if($Zhi==0){die "cannot find zhi in $data_file";}

$atoms_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Atoms\s+/)
  {
  $atoms_flag = 1;
  last;
  }
 }
if($atoms_flag==0){die "cannot find atoms section in $data_file";}

for($i=0;$i<$Natoms;$i++)
 {
 $str = <DATA>;
 while($str =~ /^\s+$/)
   {
   $str = <DATA>;
   }
 $str =~ s/^\s+//;
 my @arr = split /\s+/, $str;
# (($i+1)==$arr[0]) or die "atom numeration error";
# $atom_mol[$i]  = $arr[1];
# $atom_type[$i] = $arr[2];
 $atom_mol[$arr[0]-1]  = $arr[1];
 $atom_type[$arr[0]-1] = $arr[2];
 }

$bonds_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Bonds\s+/)
  {
  $bonds_flag = 1;
  last;
  }
 }
if($bonds_flag==0){die "cannot find bonds section in $data_file";}

for($i=0;$i<$Nbonds;$i++)
 {
 $str = <DATA>;
 while($str =~ /^\s+$/)
   {
   $str = <DATA>;
   }
 $str =~ s/^\s+//;
 my @arr = split /\s+/, $str;
 (($i+1)==$arr[0]) or die "bond numeration error";
 $bond_type[$i] = $arr[1];
 $bond_i[$i]    = $arr[2];
 $bond_j[$i]    = $arr[3];
 }

close(DATA);


#open INPUT, "<$in_file" or die "Cannot open $in_file: $!";
if($in_file =~ /\.gz$/)
 {
 open INPUT, "zcat $in_file |" or die "Cannot open $in_file: $!";
 }
 else
 {
open INPUT, "<$in_file" or die "Cannot open $in_file: $!";
 }
#open(OUTPUT,">".$outfile) or die "cant open ".$outfile;


$flag = 0;

OUTCYCLE: while(1)
 {
 while($str=<INPUT>)
  {
  if($str =~ /^ITEM: TIMESTEP/)
   {
   last;
   }
  }

  if(eof(INPUT))
   {
   close(INPUT);
#   close(OUTPUT);

   die "OK";
   }

 $str=<INPUT>;
 $step = sprintf "%d", $str;

print $step."\n";

 my $index = in_step($step);
  if($index == -1)
   {
   next OUTCYCLE;
   }  

 while($str=<INPUT>)
  {
  if($str =~ /^ITEM: NUMBER OF ATOMS/)
   {
   last;
   }
  }
 $str=<INPUT>;
 $Natoms2 = sprintf "%d", $str;
 ($Natoms == $Natoms2) or die "number of atoms differ in data file ( $Natoms ) and trajectory file ( $Natoms2 ) (timestep $step )";

 while($str=<INPUT>)
  {
  if($str =~ /^ITEM: BOX BOUNDS/)
   {
   last;
   }
  }
 ($rest) = $str =~ /^ITEM: BOX BOUNDS(.+)$/;
 $rest =~ s/^\s+//;
 @periodic = split /\s+/, $rest;

 $str=<INPUT>;
 @arr = split /\s+/, $str;
 $lox = $arr[0];
 $hix = $arr[1];
 $str=<INPUT>;
 @arr = split /\s+/, $str;
 $loy = $arr[0];
 $hiy = $arr[1];
 $str=<INPUT>;
 @arr = split /\s+/, $str;
 $loz = $arr[0];
 $hiz = $arr[1];
# end ITEM: BOX BOUNDS


 while($str=<INPUT>)
  {
  if($str =~ /^ITEM: ATOMS/)
   {
   ($str1) = $str =~ /^ITEM: ATOMS\s*(.*)$/;
   @column_names = split /\s+/, $str1;
   for($itmp=0;$itmp<scalar(@column_names);$itmp++)
    {
    $column_index{$column_names[$itmp]}=$itmp;
    }
   last;
   }
  }
 $molecules_number = 0;

  defined($column_index{'id'}) || die "no id ";
  defined($column_index{'type'}) || die "no type ";
  defined($column_index{'mol'}) || die "no mol ";
#  defined($column_index{$col_name}) || die "no col found ";
if(defined($column_index{'xu'}))
{
  defined($column_index{'xu'}) || die "no xu ";
  defined($column_index{'yu'}) || die "no yu ";
  defined($column_index{'zu'}) || die "no zu ";
  $coord_flag = "u";
}
else
{
 if(defined($column_index{'xs'}))
  {
  defined($column_index{'xs'}) || die "no xs ";
  defined($column_index{'ys'}) || die "no ys ";
  defined($column_index{'zs'}) || die "no zs ";
  $coord_flag = "s";
  }
 else
 {
 die "no coords ";
 }
}
#  defined($column_index{'ix'}) || die "no ix ";
#  defined($column_index{'iy'}) || die "no iy ";
#  defined($column_index{'iz'}) || die "no iz ";

# start calculating statistics

@atom2_mol = ();
@atom2_type = ();
@x = ();
@y = ();
@z = ();

@x_c = ();
@y_c = ();
@z_c = ();
@num = ();
@rx = ();
@ry = ();
@rz = ();
@rg2x = ();
@rg2y = ();
@rg2z = ();
@num2 = ();
@rg = ();
@rg2 = ();
$rh2 = 0.0;
$rh = 0.0;

 for($iatom=0;$iatom<$Natoms;$iatom++)
  {
  $str=<INPUT>;
  $str =~ s/^\s+//;
  @arr = split /\s+/, $str;

  $atom2_mol[$arr[$column_index{'id'}]]  = $arr[$column_index{'mol'}];
  $atom2_type[$arr[$column_index{'id'}]] = $arr[$column_index{'type'}];
  $atom2_newmol[$arr[$column_index{'id'}]] = $arr[$column_index{$col_name}];

if($coord_flag eq "u")
 {
  $x[$arr[$column_index{'id'}]] = $arr[$column_index{'xu'}];
  $y[$arr[$column_index{'id'}]] = $arr[$column_index{'yu'}];
  $z[$arr[$column_index{'id'}]] = $arr[$column_index{'zu'}];
 }
 else
 {
  $x[$arr[$column_index{'id'}]] = $arr[$column_index{'xs'}]*($hix-$lox);
  $y[$arr[$column_index{'id'}]] = $arr[$column_index{'ys'}]*($hiy-$loy);
  $z[$arr[$column_index{'id'}]] = $arr[$column_index{'zs'}]*($hiz-$loz);
 }
  $ix[$arr[$column_index{'id'}]] = $arr[$column_index{'ix'}];
  $iy[$arr[$column_index{'id'}]] = $arr[$column_index{'iy'}];
  $iz[$arr[$column_index{'id'}]] = $arr[$column_index{'iz'}];
#  print "i x y z $arr[$column_index{'id'}] $x[$arr[$column_index{'id'}]] $y[$arr[$column_index{'id'}]] $z[$arr[$column_index{'id'}]]\n";
  }
  $lx = ($Xhi-$Xlo);
  $ly = ($Yhi-$Ylo);
  $lz = ($Zhi-$Zlo);


if((!defined($start_time)||($step>=$start_time))&&(!defined($end_time)||($step<=$end_time))&&((!defined($every_time))||(($step-$start_time)%$every_time)==0))
#if(in_step($step)>0)
{
#start calculate rg via r_ij

print "calculating step $step \n";


$new_out_file = $outfile."_s".$step.".data";

print "writing to $new_out_file \n";

open(OUT,">".$new_out_file) or die "cant create file ".$new_out_file;

print OUT $desc."\n\n";
printf OUT "%11d atoms\n", $Natoms;
printf OUT "%11d bonds\n", $Nbonds;

print OUT <<END;

         4 atom types
         1 bond types

END

printf OUT "      %f %f xlo xhi\n", $lox, $hix;
printf OUT "      %f %f ylo yhi\n", $loy, $hiy;
printf OUT "      %f %f zlo zhi\n", $loz, $hiz;

print OUT "\nMasses\n\n";
print OUT "1 1.0\n";
print OUT "2 1.0\n";
print OUT "3 1.0\n";
print OUT "4 1.0\n";

#print OUT "\nBond Coeffs\n\n";
#print OUT "1 100 1.0\n";

print OUT "\nAtoms\n\n";

for($iatom=1;$iatom<=$Natoms;$iatom++)
 {
 printf OUT "%6d %6d %6d %8.2f %8.2f %8.2f\n",$iatom,$atom2_mol[$iatom],$atom2_type[$iatom],$x[$iatom],$y[$iatom],$z[$iatom];
 }

print OUT "\nVelocities\n\n";

for($iatom=1;$iatom<=$Natoms;$iatom++)
 {
 printf OUT "%6d %8.2f %8.2f %8.2f\n",$iatom,0,0,0;
 }

print OUT "\nBonds\n\n";

for($ibond=1;$ibond<=$Nbonds;$ibond++)
 {
 printf OUT "%6d %6d %6d %6d\n",$ibond,$bond_type[$ibond],$bond_i[$ibond],$bond_j[$ibond];
 }

close(OUT);


} #start time end time

#finish calculate rg via r_ij

 } # OUTCYCLE

