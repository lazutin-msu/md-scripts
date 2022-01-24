#!/usr/bin/perl
use Getopt::Long;
#Getopt::Long::Configure ('bundling');

@cutoffs = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100);
$Ncutoffs = scalar(@cutoffs);

GetOptions ('f|file=s' => \$in_file, 'd|data=s' => \$data_file, 'o|out=s' => \$outfile, 'b|begin=i' => \$start_time, 'e|end=i' => \$end_time, 'h|help' => \$help);

if($help) {
  print "Options for this program:\n-f --file for input file default= dump.all\n-d --data for data file default= data\n-t --dt time to output\n-s --slice slice height to use\n-o --output for output file default=dump.vtf\n-l --log log file with aggregates data @t\n";
  exit;
  }


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


open DATA, "<$data_file" or die "Cannot open $data_file: $!";

$Natoms = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s+\d+\s+atoms\s+$/)
  {
  ($Natoms) = ($str =~ /^\s+(\d+)\s+atoms\s+$/);
  last;
  }
 }
if($Natoms==0){die "cannot find number of atoms in $data_file";}
$Nbonds = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s+\d+\s+bonds\s+$/)
  {
  ($Nbonds) = ($str =~ /^\s+(\d+)\s+bonds\s+$/);
  last;
  }
 }
if($Nbonds==0){die "cannot find number of bonds in $data_file";}

$Xhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s+\d+\.?\d*\s+\d+\.?\d*\s+xlo\s+xhi\s+$/)
  {
  ($Xlo,$Xhi) = ($str =~ /^\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+xlo\s+xhi\s+$/);
  last;
  }
 }
if($Xhi==0){die "cannot find xhi in $data_file";}
$Yhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s+\d+\.?\d*\s+\d+\.?\d*\s+ylo\s+yhi\s+$/)
  {
  ($Ylo,$Yhi) = ($str =~ /^\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+ylo\s+yhi\s+$/);
  last;
  }
 }
if($Yhi==0){die "cannot find yhi in $data_file";}
$Zhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s+\d+\.?\d*\s+\d+\.?\d*\s+zlo\s+zhi\s+$/)
  {
  ($Zlo,$Zhi) = ($str =~ /^\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+zlo\s+zhi\s+$/);
  last;
  }
 }
if($Zhi==0){die "cannot find zhi in $data_file";}

$atoms_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Atoms\s+$/)
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
 (($i+1)==$arr[0]) or die "atom numeration error";
 $atom_mol[$i]  = $arr[1];
 $atom_type[$i] = $arr[2];
 }

$bonds_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Bonds\s+$/)
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

# output
open(OUT,">".$outfile) or die "cant open ".$outfile;
print OUT "r p2\n";
for($icut=0;$icut<$Ncutoffs;$icut++)
 {
 if($p2_norm[$icut]>0)
  { 
  printf OUT "%g %g %g\n",$cutoffs[$icut], $p2_hist[$icut]/$p2_norm[$icut], $p2_norm[$icut];
  }
  else
  {
  printf OUT "%g %g %g\n",$cutoffs[$icut], $p2_hist[$icut], $p2_norm[$icut];
  }
 }
 printf OUT "inf %g %g\n",$p2_hist_all/$p2_norm_all, $p2_norm_all;
close(OUT);
# end output

   die "OK";
   }

 $str=<INPUT>;
 $step = sprintf "%d", $str;

print $step."\n";
 
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
#  defined($column_index{'c_clust_z'.$slice}) || die "no c_clust_z${slice} ";
#  defined($column_index{'c_cluster1'}) || die "no c_cluster1 ";
  defined($column_index{'xu'}) || die "no xu ";
  defined($column_index{'yu'}) || die "no yu ";
  defined($column_index{'zu'}) || die "no zu ";
  defined($column_index{'ix'}) || die "no ix ";
  defined($column_index{'iy'}) || die "no iy ";
  defined($column_index{'iz'}) || die "no iz ";

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
#  (($iatom+1) == $arr[0]) or die "atom order is wrong in $in_file timestep $step";
#  $atom_id[$iatom]
#  $x[$arr[0]] = $arr[2];
#  $y[$arr[0]] = $arr[3];
#  $z[$arr[0]] = $arr[4];

#  $x[$arr[0]] = $arr[3];
#  $y[$arr[0]] = $arr[4];
#  $z[$arr[0]] = $arr[5];
  $atom2_mol[$arr[$column_index{'id'}]]  = $arr[$column_index{'mol'}];
  $atom2_type[$arr[$column_index{'id'}]] = $arr[$column_index{'type'}];
  $x[$arr[$column_index{'id'}]] = $arr[$column_index{'xu'}];
  $y[$arr[$column_index{'id'}]] = $arr[$column_index{'yu'}];
  $z[$arr[$column_index{'id'}]] = $arr[$column_index{'zu'}];
#  print "i x y z $arr[$column_index{'id'}] $x[$arr[$column_index{'id'}]] $y[$arr[$column_index{'id'}]] $z[$arr[$column_index{'id'}]]\n";
  }
  $lx = ($Xhi-$Xlo);
  $ly = ($Yhi-$Ylo);
  $lz = ($Zhi-$Zlo);


# orientation vectors calculation
  for($iatom=1;$iatom<=$Natoms;$iatom++)
  {
  if($atom2_type[$iatom]==2)
   {
   $jatom = $iatom - 1;
        my $dx = ($x[$iatom] - $x[$jatom]);
        my $dy = ($y[$iatom] - $y[$jatom]);
        my $dz = ($z[$iatom] - $z[$jatom]);
        while($dx>$lx/2.0)
         {
         $dx = $dx - $lx;
         }
        while($dx<-$lx/2.0)
         {
         $dx = $dx + $lx;
         }
        while($dy>$ly/2.0)
         {
         $dy = $dy - $ly;
         }
        while($dy<-$ly/2.0)
         {
         $dy = $dy + $ly;
         }
   $atom_vector_x[$iatom] = $dx;
   $atom_vector_x[$jatom] = $dx;
   $atom_vector_y[$iatom] = $dy;
   $atom_vector_y[$jatom] = $dy;
   $atom_vector_z[$iatom] = $dz;
   $atom_vector_z[$jatom] = $dz;
   $atom_vector_norm[$iatom] = $dx*$dx + $dy*$dy + $dz*$dz;
   $atom_vector_norm[$jatom] = $dx*$dx + $dy*$dy + $dz*$dz;
#   if($atom_vector_norm[$iatom]==0){die "step $step atom $iatom $jatom vectornorm 0";}
   }
  }

if((!defined($start_time)||($step>=$start_time))&&(!defined($end_time)||($step<=$end_time)))
{
#start calculate rg via r_ij
print "calculating step $step \n";
  for($iatom=1;$iatom<=$Natoms;$iatom++)
  {
#     for($jatom=1;$jatom<=$Natoms;$jatom++)
     for($jatom=1;$jatom<=$iatom;$jatom++)
     {
      if($iatom!=$jatom)
       {
        my $dx = ($x[$iatom] - $x[$jatom]);
        my $dy = ($y[$iatom] - $y[$jatom]);
        my $dz = ($z[$iatom] - $z[$jatom]);
        while($dx>$lx/2.0)
         {
         $dx = $dx - $lx;
         }
        while($dx<-$lx/2.0)
         {
         $dx = $dx + $lx;
         }
        while($dy>$ly/2.0)
         {
         $dy = $dy - $ly;
         }
        while($dy<-$ly/2.0)
         {
         $dy = $dy + $ly;
         }
#        $rg2x[$atom_aggr[$iatom]] += $dx*$dx;
#        $rg2y[$atom_aggr[$iatom]] += $dy*$dy;
#        $rg2z[$atom_aggr[$iatom]] += $dz*$dz;
#        $rh2 += 1.0/($dx*$dx + $dy*$dy + $dz*$dz);
#        $rh += 1.0/sqrt($dx*$dx + $dy*$dy + $dz*$dz);
#        $num2[$atom_aggr[$iatom]] += 1;
         my $r2ij = $dx*$dx + $dy*$dy + $dz*$dz;
         my $multij = ( $atom_vector_x[$iatom]*$atom_vector_x[$jatom] + $atom_vector_y[$iatom]*$atom_vector_y[$jatom] + $atom_vector_z[$iatom]*$atom_vector_z[$jatom] ) ;
         my $cos2ij = $multij * $multij / ($atom_vector_norm[$iatom] * $atom_vector_norm[$jatom]);
         my $p2 = 0.5 * (3.0*$cos2ij - 1.0);
         if($p2>1.0){print "time $step i $iatom j $jatom p $p2 cos $cos2ij mult $multij norm ".$atom_vector_norm[$iatom]." ". $atom_vector_norm[$jatom]." \n";}
         for($icut=0;$icut<$Ncutoffs;$icut++)
           {
           if($r2ij<=$cutoffs[$icut]*$cutoffs[$icut])
            {
            $p2_hist[$icut] = $p2_hist[$icut] + $p2;
            $p2_norm[$icut] = $p2_norm[$icut] + 1.0;
            }
           }
          $p2_hist_all = $p2_hist_all + $p2;
          $p2_norm_all = $p2_norm_all + 1.0;

#         $num2 += 1;

       }
     }
  }

} #start time end time

#finish calculate rg via r_ij

 } # OUTCYCLE

