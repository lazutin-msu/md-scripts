#!/usr/bin/perl
use Getopt::Long;
#use Data::Dumper qw(Dumper);
#Getopt::Long::Configure ('bundling');
#use Math::MatrixDecomposition::Eigen;

#@cutoffs = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100);
#$Ncutoffs = scalar(@cutoffs);
#$Npoly=100;
#$delta_r = 0.1;

GetOptions (
'f|file=s' => \$in_file, 
'd|data=s' => \$data_file, 
'o|out=s' => \$outfile, 
'b|begin=i' => \$start_time, 
'e|end=i' => \$end_time, 
't|type=i' => \$out_type, 
'r|cutoff=f' => \$cutoff,
'g|grid=f' => \$grid_size,
'c|chains=i' => \$chains_per_particle,
'p|particles=i' => \$particles_number,
'h|help' => \$help
);

if($help) {
  print "Options for this program:\n-f --file for input file default= dump.all\n-d --data for data file default= data\n-t --dt time to output\n-s --slice slice height to use\n-o --output for output file default=dump.vtf\n-t --type atom type to output\n-l --log log file with aggregates data @t\n";
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

#if(!$out_type)
# {
# $out_type = 2;
# }

#if(!$cutoff)
# {
# $cutoff = 1.3;
# }
#$cutoff2 = $cutoff * $cutoff;

sub sphere_bins
{
  # number of bins should be even
  my $phi    = shift;
  my $theta  = shift;
  my $Nphi   = shift;
  my $Ntheta = shift;
  my $iphi   = myround( $Nphi   * $phi   / pi / 2.0);
  my $itheta = myround( $Ntheta * cos( $theta ) / 2.0 );
  my @return = ($iphi, $itheta);
  return @return;
}

sub mypolymer
{
my $iatom = shift;
my $jatom = shift;
my $chains_per_particle = shift;
my $particles_number = shift;
my $ci = int(($iatom-1) / ($chains_per_particle * 100 + 1) );
my $cj = int(($jatom-1) / ($chains_per_particle * 100 + 1) );
return ($ci==$cj);
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
 if($str =~ /^\s*\-?\d+\.?\d*e?[\+\-]?\d+\s+\d+\.?\d*e?[\+\-]?\d+\s+xlo\s+xhi\s*$/)
  {
  ($Xlo,$Xhi) = ($str =~ /^\s*(\-?\d+\.?\d*e?[\+\-]?\d+)\s+(\d+\.?\d*e?[\+\-]?\d+)\s+xlo\s+xhi\s*$/);
  last;
  }
 }
if($Xhi==0){die "cannot find xhi in $data_file";}
$Yhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\-?\d+\.?\d*e?[\+\-]?\d+\s+\d+\.?\d*e?[\+\-]?\d+\s+ylo\s+yhi\s*$/)
  {
  ($Ylo,$Yhi) = ($str =~ /^\s*(\-?\d+\.?\d*e?[\+\-]?\d+)\s+(\d+\.?\d*e?[\+\-]?\d+)\s+ylo\s+yhi\s*$/);
  last;
  }
 }
if($Yhi==0){die "cannot find yhi in $data_file";}
$Zhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\-?\d+\.?\d*e?[\+\-]?\d+\s+\d+\.?\d*e?[\+\-]?\d+\s+zlo\s+zhi\s*$/)
  {
  ($Zlo,$Zhi) = ($str =~ /^\s*(\-?\d+\.?\d*e?[\+\-]?\d+)\s+(\d+\.?\d*e?[\+\-]?\d+)\s+zlo\s+zhi\s*$/);
  last;
  }
 }
if($Zhi==0){die "cannot find zhi in $data_file";}

$atoms_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Atoms\s*/)
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
open(OUT,">".$outfile) or die "cant open ".$outfile;


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

################################################################################################################################

print OUT "i r hist_my norm_my hist_oth norm_oth\n";

for($i=0;$i<=$max_ir;$i++)
{
printf OUT "%d %g %g %g %g %g\n", $i, $i*$grid_size, $densitymy[$i], $densitymy[$i]/$denmy_norm, $densityoth[$i], $densityoth[$i]/$denoth_norm;
}

close(OUT);
# end output

   die "OK";
   }

 $str=<INPUT>;
 $step = sprintf "%d", $str;

print $step."\n";

if(!((!defined($start_time)||($step>=$start_time))&&(!defined($end_time)||($step<=$end_time))))
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
#  defined($column_index{'mol'}) || die "no mol ";
#  defined($column_index{'c_clust_z'.$slice}) || die "no c_clust_z${slice} ";
#  defined($column_index{'c_cluster1'}) || die "no c_cluster1 ";
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


# orientation vectors calculation

#   if($atom_vector_norm[$iatom]==0){die "step $step atom $iatom $jatom vectornorm 0";}
#   }
#  }

if((!defined($start_time)||($step>=$start_time))&&(!defined($end_time)||($step<=$end_time)))
{
#start calculate rg via r_ij

print "calculating step $step \n";

# generate golden spiral

sub distance
{
my $x1 = shift;
my $y1 = shift;
my $z1 = shift;
my $x2 = shift;
my $y2 = shift;
my $z2 = shift;

my $distance = sqrt(($x1-$x2)*($x1-$x2) + ($y1-$y2)*($y1-$y2) + ($z1-$z2)*($z1-$z2));

return $distance;
}

sub myint
{
    my ($x, $lx) = @_;
    
    while($x<0)
    {$x= $x + $lx;
    }
    while($x>=$lx)
    {$x= $x - $lx;
    }
    return $x;
}

$lx2 = $hix - $lox;
$ly2 = $hiy - $loy;
$lz2 = $hiz - $loz;

$nx = int( $lx2 / $grid_size) ;
$ny = int( $ly2 / $grid_size) ;
$nz = int( $lz2 / $grid_size) ;

print "check grid size ".$lx2." ".($nx*$grid_size)."\n";
print "check grid size ".$ly2." ".($ny*$grid_size)."\n";
print "check grid size ".$lz2." ".($nz*$grid_size)."\n";

print "check beads number ".$Natoms." ".($chains_per_particle*$particles_number*100)." = ".$chains_per_particle." * ".$particles_number." * 100\n";

 for($iatom=1;$iatom<=$Natoms;$iatom++)
 {
 if($atom2_type[$iatom]==4)
  {
  my $xi = $x[$iatom];
  my $yi = $y[$iatom];
  my $zi = $z[$iatom];
  
  while($xi<0) {$xi= $xi + $lx2;}
  while($xi>=$lx2) {$xi= $xi - $lx2;}
  while($yi<0) {$yi= $yi + $ly2;}
  while($yi>=$ly2) {$yi= $yi - $ly2;}
  while($zi<0) {$zi= $zi + $lz2;}
  while($zi>=$lz2) {$zi= $zi - $lz2;}

  for($jatom=1;$jatom<$Natoms;$jatom++)
   {
   
   if($atom2_type[$jatom]!=4)
    {

    my $xj = $x[$jatom];
    my $yj = $y[$jatom];
    my $zj = $z[$jatom];
  
    while($xj<0) {$xj= $xj + $lx2;}
    while($xj>=$lx2) {$xj= $xj - $lx2;}
    while($yj<0) {$yj= $yj + $ly2;}
    while($yj>=$ly2) {$yj= $yj - $ly2;}
    while($zj<0) {$zj= $zj + $lz2;}
    while($zj>=$lz2) {$zj= $zj - $lz2;}
    
    my $dx = $xi - $xj;
    my $dy = $yi - $yj;
    my $dz = $zi - $zj;
    
    while($dx<-0.5*$lx2) {$dx= $dx + $lx2;}
    while($dx>=0.5*$lx2) {$dx= $dx - $lx2;}
    while($dy<-0.5*$ly2) {$dy= $dy + $ly2;}
    while($dy>=0.5*$ly2) {$dy= $dy - $ly2;}
    while($dz<-0.5*$lz2) {$dz= $dz + $lz2;}
    while($dz>=0.5*$lz2) {$dz= $dz - $lz2;}
    
    my $dr2 = $dx * $dx + $dy * $dy + $dz * $dz;
    my $dr = sqrt($dr2);

    my $ir = int ($dr / $grid_size);
    
    if(mypolymer($iatom,$jatom,$chains_per_particle,$particles_number))
     {
     $densitymy[$ir] += 1;
     $denmy_norm += 1;
     }
     else
     {
     $densityoth[$ir] += 1;
     $denoth_norm += 1;
     }
    if($ir>$max_ir){$max_ir = $ir;}
    
    }
   } # jatom

  
  }
 } # iatom


} #start time end time

#finish calculate rg via r_ij

 } # OUTCYCLE

