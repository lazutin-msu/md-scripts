#!/usr/bin/perl
use Getopt::Long;
use File::Copy;
use Math::Trig;
use POSIX;
    
$dens = 8;
$N = 50;
$chains = 10; #10x10
$time = 2000000;
$repeat = 1;
$R_sphere = 30;

$Ngpu = 1;
#$Ncpu = 6;
#$path_to_dump = "/home/lazutin/NAS-Personal/lazutin/l1_BB_change_pot_re";
$path_to_dump = "";
$r_aggr = 1.3;

open(LOG, ">>perl.log");
print LOG  $0;
for($i=0;$i<scalar(@ARGV);$i++)
{
 print LOG " ".$ARGV[$i];
}
print LOG "\n";
close(LOG);

GetOptions (
'datafile|f=s' => \$datafile,
'density|d=f' => \$dens,
'length|n=i' => \$N,
'chains|c=i' => \$chains,
'rsphere|a=f' => \$R_sphere,
'repeat|r=i' => \$repeat,
'time|t=i' => \$time,
'step|s=f' => \$step,
'g|gpuid=i' => \$gpu_num,
'p|spheres=i' => \$spheres,
'x|dx=f' => \$spheres_dx,
'pressure=f' => \$pressure,
'structure=s' => \$structure,
'stiff=f' => \$k_stiff,
'fix_particle' => \$fix_particle,
'lomo' => \$lomo_flag,
'fix_eps=f{3}' => \@fix_eps,
'step_dump=i' => \$step_dump,
'comm_mod=f' => \$comm_mod,
'rbond=f' => \$rbond,
'kbond=f' => \$kbond,
'np=i' => \$np,
'Ncpu=i' => \$Ncpu,
'help|h' => \$help);
if($help)
 {
 die("Usage make_data.pl options
  -f datafile initial configuration
  -d density = 8
  -c chains = 10
  -n length = 50
  -a rsphere = 30
  -r repeat = 1
  -t time = 2000000
  -s step 
  -p spheres
  -x dx spheres dx
  --rbond (0)
  --kbond (40)
  -g gpuid
  --pressure
  --structure HP or PH
  --lomo to calc on Lomonosov-2
  --fix_particle to fix particle mass center
  --fix_eps epsAA epsAB epsBB run at this energy parameters 
  --comm_mod communication cutoff distance
  --np number of cores to run on
  --Ncpu number of tasks which run in parallel
  -h help\n");
 }

if($structure ne "PH") {$structure = "HP";}

if(!$Ncpu){$Ncpu = 6;}

if(!$rbond){$rbond = 0.0;}
if(!$kbond){$kbond = 40.0;}
if(!$k_stiff){$k_stiff = 20;}

print STDERR "eps = ".$epsilon[0]." ".$epsilon[1]." ".$epsilon[2]."\n";

$time_lastdump = $time / 2;
$time_dump = $time / 20;

if(scalar(@fix_eps))
 {
 if($step_dump)
   {
   $time_dump = $step_dump;
   }
 }

sub myround
{ # for positive x
my $x = shift;
return int($x+0.5);
}



for($idirs=0;$idirs<$repeat;$idirs++)
{

$gpuid = 0;
$cpuid = 0;

if(-f ".make_data.conf")
 {
 open(CONFIG, "<.make_data.conf");
 my $str = <CONFIG>;
 print STDERR "'".$str."' '".chomp($str)."'\n";
# $gpuid = int(chomp($str));
 $str =~ s/\s$//;
 $str =~ s/^\s//;
 $gpuid = int($str);
 my $str = <CONFIG>;
 print STDERR "'".$str."' '".chomp($str)."'\n";
# $gpuid = int(chomp($str));
 $str =~ s/\s$//;
 $str =~ s/^\s//;
 $cpuid = int($str);
 close(CONFIG);
 }

 print STDERR $gpuid."\n";
 $gpuid++;
 print STDERR $gpuid."\n";
 if( $gpuid >= $Ngpu )
  {
  $gpuid=0;
  }
 print STDERR $gpuid."\n";
 print STDERR $cpuid."\n";
 $cpuid++;
 print STDERR $cpuid."\n";
 if( $cpuid >= $Ncpu )
  {
  $cpuid=0;
  }
 print STDERR $cpuid."\n";
 open(CONFIG, ">.make_data.conf");
 printf CONFIG "%d\n%d\n",$gpuid,$cpuid;
 close(CONFIG);

if($gpu_num)
{
$gpuid = $gpu_num;
}


#if($gpuid==1){$gpuid=2;}
printf "R_sphere = %f\nchains = %d\nlength = %d \ngpuid= %d\n",$R_sphere,$chains,$N,$gpuid;

$Ndir = 1;

$desc = sprintf "brush_p%d_c%d_R%f_N%d_str%s_press%.2f",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];

if($datafile)
{
$is_datafile = 1;
}
else
{
$is_datafile = 0;
$datafile = sprintf "data_p%d_c%d_R%f_N%d_str%s_press%.2f",$spheres,$chains,$R_sphere,$N,$structure,$pressure;
}

if(scalar(@fix_eps))
 {
$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  $dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  $output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
  }
 }
 }
 else
 {
$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  $dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  $output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
  }
 }
 
 } #fix_eps

`mkdir $dirname`;

$dirname = $dirname."/";

#$dx = sqrt(1.0 / $dens);
$dx = $dens; # dens now size of grafting points' cell
$radius = 0.5;
#$bond = 2.0*$radius;
$bond = 0.5;
#if($dx<(2.0*$bond)){$bond=$dx/2.0;}

$delta_rho = $bond;

$natoms = 2 * $chains * $N + 1;
if(!$k_stiff)
 {
$nbonds = (2*$N - 1) * $chains;
 }
 else
 {
$nbonds = (2*$N - 1) * $chains + $chains ;
 }

$totalatoms = $natoms * $spheres;
$totalbonds = $nbonds * $spheres;

#$lx =  ( 2 * $N + $R_sphere ) +2;
#$ly =  ( 2 * $N + $R_sphere ) +2;
#$lz =  ( 2 * $N + $R_sphere ) +2;
$lx =  ( $N * $bond + $R_sphere ) +2;
$ly =  ( $N * $bond + $R_sphere ) +2;
$lz =  ( $N * $bond + $R_sphere ) +2;


print "$desc \n";
$ibond = 1;

$Nspheres = ceil( $spheres ** (1/3) );
#$totallx = $Nspheres * $spheres_dx + 2* $lx;
#$totally = $Nspheres * $spheres_dx + 2* $ly;
#$totallz = $Nspheres * $spheres_dx + 2* $lz;
$totallx = ($Nspheres - 1) * $spheres_dx + 2* $lx;
$totally = ($Nspheres - 1) * $spheres_dx + 2* $ly;
$totallz = ($Nspheres - 1) * $spheres_dx + 2* $lz;

$isphere = 0;

SPHERES: for($iz=0;$iz<$Nspheres;$iz++)
 {
 for($iy=0;$iy<$Nspheres;$iy++)
  {
  for($ix=0;$ix<$Nspheres;$ix++)
   {
   

$sp_dlong = pi*(3.0-sqrt(5.0)); #  /* ~2.39996323 */
$sp_dz    = 2.0/$chains;
$sp_long = 0;
$z    = 1.0 - $sp_dz/2.0;

$delta_phi = 1/ $R_sphere; 

  $x[$natoms * $isphere + 1] = $ix * $spheres_dx + $lx +    0;
  $y[$natoms * $isphere + 1] = $iy * $spheres_dx + $ly +    0;
  $z[$natoms * $isphere + 1] = $iz * $spheres_dx + $lz +    0;
  $type[$natoms * $isphere + 1] = 4;
  $mol[$natoms * $isphere + 1] = ($chains + 1) * $isphere + 1;
  $coreid[$natoms * $isphere + 1] = $isphere + 1;
  $partnum[$natoms * $isphere + 1] = $isphere + 1;

for($k=0;$k<$chains;$k++)
   {
    $r    = sqrt(1.0-$z*$z);
    ($xchain, $ychain, $zchain) = (cos($long)*$r, sin($long)*$r, $z);
   my $rho   = sqrt($r*$r+$z*$z)*$R_sphere;
   my $phi   = $sp_long;
   my $theta = atan2($r,$z);
    $z    = $z - $sp_dz;
    $sp_long = $sp_long + $sp_dlong;
   $rho = $rho + $delta_rho * 0.5;

  for($imono=0;$imono<$N;$imono++)
   {
   $iatom = 2*$imono + 2*$k * $N + 2;

   $x[$natoms * $isphere + $iatom] = $ix * $spheres_dx + $lx + $rho * sin($theta) * cos($phi);
   $y[$natoms * $isphere + $iatom] = $iy * $spheres_dx + $ly + $rho * sin($theta) * sin($phi);
   $z[$natoms * $isphere + $iatom] = $iz * $spheres_dx + $lz + $rho * cos($theta) ;
   if($z[$iatom]>0)
    {
    $my_delta_theta = $bond / ($R_sphere + $imono * $delta_rho);
    }
    else
    {
    $my_delta_theta = -$bond / ($R_sphere  + $imono * $delta_rho);
    }

   $type[$natoms * $isphere + $iatom] = 1;
   if($imono==0){$type[$natoms * $isphere + $iatom] = 3; $coreid[$natoms * $isphere + $iatom] = $isphere + 1;}
   $mol[$natoms * $isphere + $iatom] = ($chains + 1) * $isphere + $k + 2;
   $partnum[$natoms * $isphere + $iatom] = $isphere + 1;


   $x[$natoms * $isphere + $iatom+1] = $ix * $spheres_dx + $lx + $rho * sin($theta + $my_delta_theta) * cos($phi);
   $y[$natoms * $isphere + $iatom+1] = $iy * $spheres_dx + $ly + $rho * sin($theta + $my_delta_theta) * sin($phi);
   $z[$natoms * $isphere + $iatom+1] = $iz * $spheres_dx + $lz + $rho * cos($theta + $my_delta_theta);
   $rho = $rho + $delta_rho; 
   $type[$natoms * $isphere + $iatom+1] = 2;
   $mol[$natoms * $isphere + $iatom+1] = ($chains + 1) * $isphere + $k + 2;
   $partnum[$natoms * $isphere + $iatom+1] = $isphere + 1;
if(!$k_stiff)
 {
   if($imono!=0)
    {
    $at1_bond[$ibond] = $natoms * $isphere + $iatom - 2;
    $at2_bond[$ibond] = $natoms * $isphere + $iatom;
    $bond_type[$ibond] = 1;
    $ibond++;
    }
 }
 else
 {
   if($imono!=0)
    {
    $at1_bond[$ibond] = $natoms * $isphere + $iatom - 2;
    $at2_bond[$ibond] = $natoms * $isphere + $iatom;
    $bond_type[$ibond] = 1;
    $ibond++;
    }
    else
    {
    $at1_bond[$ibond] = $natoms * $isphere + 1;
    $at2_bond[$ibond] = $natoms * $isphere + $iatom;
    $bond_type[$ibond] = 2;
    $ibond++;
    }
 }
   $at1_bond[$ibond] = $natoms * $isphere + $iatom;
   $at2_bond[$ibond] = $natoms * $isphere + $iatom+1;
   $bond_type[$ibond] = 1;
   $ibond++;
   }
   }

  if($isphere<$spheres-1)
   {
   $isphere++;
   }
   else
   {
   last SPHERES;
   }

  } # ix
 } # iy
} # iz

(($ibond-1)==$totalbonds) or die "ibond $ibond != totalbonds $totalbonds";

if(!$is_datafile)
{
print "writing to $datafile \n";
open(OUT,">".$dirname.$datafile) or die "cant create file ".$dirname.$datafile;

print OUT $desc."\n\n";
printf OUT "%11d atoms\n", $totalatoms;
printf OUT "%11d bonds\n", $totalbonds;

print OUT <<END;

         5 atom types
END
if(!$k_stiff)
{
print OUT <<END;
         1 bond types

END
}
else
{
print OUT <<END;
         2 bond types

END
}
#printf OUT "      %f %f xlo xhi\n", -$lx, $lx;
#printf OUT "      %f %f ylo yhi\n", -$ly, $ly;
#printf OUT "      %f %f zlo zhi\n", -$lz, $lz;
printf OUT "      %f %f xlo xhi\n", 0, $totallx;
printf OUT "      %f %f ylo yhi\n", 0, $totally;
printf OUT "      %f %f zlo zhi\n", 0, $totallz;

print OUT "\nMasses\n\n";
print OUT "1 1.0\n";
print OUT "2 1.0\n";
print OUT "3 1.0\n";
print OUT "4 1.0\n";
print OUT "5 1.0\n";

#print OUT "\nBond Coeffs\n\n";
#print OUT "1 100 1.0\n";

print OUT "\nAtoms\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %6d %6d 0.0 %8.2f %8.2f %8.2f\n",$iatom,$mol[$iatom],$type[$iatom],$x[$iatom],$y[$iatom],$z[$iatom];
 }

print OUT "\nVelocities\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %8.2f %8.2f %8.2f\n",$iatom,$vx[$iatom],$vy[$iatom],$vz[$iatom];
 }

print OUT "\nBonds\n\n";

for($ibond=1;$ibond<=$totalbonds;$ibond++)
 {
if(!$k_stiff)
{
 printf OUT "%6d %6d %6d %6d\n",$ibond,1,$at1_bond[$ibond],$at2_bond[$ibond];
}
else
{
 printf OUT "%6d %6d %6d %6d\n",$ibond,$bond_type[$ibond],$at1_bond[$ibond],$at2_bond[$ibond];
}
 }

close(OUT);
}
else
{
copy($datafile,$dirname.$datafile) or die "Copy failed: $!";
}

#@epsAA = (0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0, -6.5, -7.0);
#@epsAB = (0.0, 0.25,  0.5, 0.75,  1.0, 1.25,  1.5, 1.75,  2.0, 2.25,  2.5, 2.75,  3.0, 3.25,  3.5);
#@epsBB = (0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0);
# A - main chain, B - addons
#@epsAB = (0.0, 0.25,  0.5, 0.75,  0.875, 1.0, 1.125, 1.25, 1.375,  1.5, 1.625,  1.75);
#@epsBB = (0.0,  0.5,  1.0,  1.5,  1.75,  2.0,  2.25,  2.5,  2.75,  3.0,  3.25,  3.5);

#@epsAA = (0.0, -1.0, -2.0, -3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0, -6.5, -7.0);
#@epsAB = (0.0,  0.5,  1.0,  1.5,  1.75,  2.0,  2.25,  2.5,  2.75,  3.0,  3.25,  3.5);
#@epsBB = (0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0);

#@eps = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0);

#@eps = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10, 11, 12, 13, 14, 15);
#@time_mult = (  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    2,   2,    2,   2,    2,   2,    2,   2,    2,   2,    2,   2,    2,   2,    2,   2,    2,   2,    2,  2);
#@eps =       (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5,   6.7 );

#@eps =       (6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0);
#@time_mult = (   1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,    1);

#@eps = 	     (   0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4.0, 4.125, 4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0);
#@time_mult = (   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,     1,    1,     1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1);

#@time_mult = (  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    2,   2,    2,   2,    2,   2,     2 );
#@time_mult = (  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,   1,    1,   1,    1,   1,     1 );

#@eps = (0.0);

#@eps =       (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10);
#@time_mult = (  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,   1,    1,  1);
#
#for($ip=0;$ip<scalar(@eps);$ip++)
# {
# $epsAA[$ip] = $eps[$ip]*$epsilon[0] + $epsilon0[0];
# $epsAB[$ip] = $eps[$ip]*$epsilon[1] + $epsilon0[1];
# $epsBB[$ip] = $eps[$ip]*$epsilon[2] + $epsilon0[2];
# }

if($structure eq "PH")
{
@epsAA = (  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0);
@epsAB = (  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 , 5.0, 5.0 , 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0);
@epsBB = ( -0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0,-4.5,-5.0,-5.25,-5.5,-5.75,-6.0,-6.25,-6.5,-6.75,-7.0,-7.25,-7.5,-7.75,-8.0,-8.25,-8.5,-8.75,-9.0,-9.25,-9.5,-9.75, -10);
}
elsif($structure eq "HP")
{
@epsBB = (  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0);
@epsAB = (  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 , 5.0, 5.0 , 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0);
@epsAA = ( -0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0,-4.5,-5.0,-5.25,-5.5,-5.75,-6.0,-6.25,-6.5,-6.75,-7.0,-7.25,-7.5,-7.75,-8.0,-8.25,-8.5,-8.75,-9.0,-9.25,-9.5,-9.75, -10);
}
else
{
die "--structure is not HP nor PH '$structure' \n";
}

#for($ip=0;$ip<$points;$ip++)
# {
# $epsAA[$ip] = $epsilon[0];
# $epsAB[$ip] = $epsilon[1];
## $epsBB[$ip] = $epsilon[2] + $ip*$step;
# $epsBB[$ip] = $epsilon[2] + $eps[$ip]*$step;
# }

#$gpuid = int(rand(3));
$seed = int(rand(9999999));
$seed2 = int(rand(9999999));
$seed_create = int(rand(9999999));
$seed_mol = int(rand(9999999));

$num_mols = myround(3*$totallx*$totally*$totallz - $totalatoms);
#$num_mols = myround(3*($totallx*$totally*$totallz - 4.0/3.0*3.14159*$R_sphere*$R_sphere*$R_sphere ) - $totalatoms);


$xlo_bound = 0;
$xhi_bound = $totallx;
$ylo_bound = 0;
$yhi_bound = $totally;
$zlo_bound = 0;
$zhi_bound = $totallz;


$time_energy = $time /4;
$step_energy = 100;
$num_energy = $time_energy / $step_energy;




$set_temp = 1.0;

$R_sphere_cutoff = 1.1224620 * $R_sphere;
$R_sphere2 = $R_sphere - 0.5;
$hi_bound = $R_sphere + $N;

open(SCRIPT,">".$dirname.$scriptname);

if($lomo_flag)
{
print SCRIPT <<END;
#package gpu 1 neigh no gpuID $gpuid $gpuid
#suffix gpu
END
}
else
{
print SCRIPT <<END;
package gpu 1 neigh yes gpuID $gpuid $gpuid
suffix gpu
END
}

print SCRIPT <<END;
units lj
END

if($k_stiff)
{
 if($comm_mod)
 {
print SCRIPT <<END;
comm_modify cutoff ${comm_mod}
newton off
END
 }
}

print SCRIPT <<END;
atom_style full
#neighbor        1.5 bin
#neigh_modify    every 5 delay 5 check yes one 10000
neighbor        1.0 bin
neigh_modify    delay 0 every 1 check yes cluster no 
comm_modify vel yes


#pair_style hybrid/overlay yukawa 1.2 4.0 lj/cut 1.1224620 lj/expand 1.1224620
END

print SCRIPT <<END;
boundary p p p
#special_bonds fene
special_bonds lj/coul 1.0 1.0 1.0

read_data $datafile

region mysphere sphere 0 0 0 ${R_sphere} side out
#fix sphwall all wall/region mysphere lj93 2.0 1.0 0.8583742

region          bxx block ${xlo_bound} ${xhi_bound} ${ylo_bound} ${yhi_bound} ${zlo_bound} ${zhi_bound}

region byy intersect 2 bxx mysphere

create_atoms 5 random ${num_mols} ${seed_create} bxx 
#create_atoms 5 random ${num_mols} ${seed_create} byy 

bond_style harmonic
bond_coeff 1 ${kbond} ${rbond}
END

if(!$k_stiff)
{
}
else
{
$R_bond = $R_sphere + 0.5;
print SCRIPT <<END;
bond_coeff 2 ${k_stiff} ${R_bond}
END
}



print SCRIPT <<END;

pair_style      dpd 1.0 1.0 ${seed2}

pair_coeff      * * 25  4.5 1

mass * 1.0

write_data      ${output_filename}.orig.data
minimize 0.000000001 0.000000001 10000 10000

group roots type 4
group rest type 1 2 3 5
group type2 type 2
group type1 type 1 3

END

if($fix_particle)
{
print SCRIPT <<END;
velocity roots zero linear
fix root roots move linear 0 0 0 
velocity rest create ${set_temp} $seed2
fix             1 rest nve
END
}
else
{
print SCRIPT <<END;
velocity all create ${set_temp} $seed2
fix             1 all nve
END
}

print SCRIPT <<END;
thermo_style    custom step temp pe epair etotal 
thermo_modify flush yes
thermo    1000

write_data      ${output_filename}.min.data

reset_timestep 0

timestep        0.03

dump            1 all atom 100000 ${output_filename}.lammpstrj

run 100000

write_data ${output_filename}_mix.data

pair_coeff      * * 25  4.5 1

pair_coeff      1 1 25  4.5 1
pair_coeff      2 2 25  4.5 1
pair_coeff      3 3 25  4.5 1
pair_coeff      4 4 25  4.5 1
pair_coeff      5 5 25  4.5 1

pair_coeff      2 5 200 4.5 1

pair_coeff      1 2 75  4.5 1
pair_coeff      2 3 75  4.5 1

pair_coeff      1 5 20  4.5 1
pair_coeff      2 5 20  4.5 1

pair_coeff      1 4 1000  4.5 1
pair_coeff      2 4 1000  4.5 1
pair_coeff      3 4 1000  4.5 1

#pair_coeff      1 2 25.23  4.5 1
#pair_coeff      1 3 31.7  4.5 1
#pair_coeff      1 4 27.5  4.5 1
#pair_coeff      1 5 32.5  4.5 1

#pair_coeff      2 3 29.4  4.5 1
#pair_coeff      2 4 26.2  4.5 1
#pair_coeff      2 5 30.0  4.5 1


#pair_coeff      3 4 26  4.5 1
#pair_coeff      3 5 25.1  4.5 1

#pair_coeff      4 5 26.3  4.5 1


#pair_coeff 1*3 1*3 lj/cut 1.0 1.0 1.1224620
#pair_coeff * 4 lj/expand 1.0 1.0 ${R_sphere2} 1.1224620
#pair_coeff 3 4 lj/expand 0.0 1.0 ${R_sphere2} 1.1224620
#pair_coeff 1 1 yukawa $epsAA[0]
#pair_coeff 1 3 yukawa $epsAA[0]
#pair_coeff 3 3 yukawa $epsAA[0]
#pair_coeff 1 2 yukawa $epsAB[0]
#pair_coeff 2 3 yukawa $epsAB[0]
#pair_coeff 2 2 yukawa $epsBB[0]


END


print SCRIPT <<END;

run ${time} 

write_data ${output_filename}_eq.data

END

close(SCRIPT);

open(SH,">>".$shellname);
#print SH "nohup /usr/mount-opt/lammps/lmp_ubuntu_1Feb14 <$scriptname  &";
if($lomo_flag)
{
print SH <<END;
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-10Feb15/src/lmp_cuda -sf gpu  -in $scriptname 
#/home/lazutin/src/lammps-11Aug17/src/lmp_cuda -sf gpu  -in $scriptname 
#/home/lazutin/src/lammps-16Mar18/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
#mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
sbatch -n 256 --time=1-23:30:00 -p compute_prio ompi /home/lazutin_2123/_scratch/src/lammps-29Oct20/build-nogpu/lmp  -in script
cd ..
END
}
else
{
if($np)
 {
 if($np>1)
  {
print SH <<END;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-10Feb15/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-11Aug17/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-16Mar18/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
#mpirun -np ${np} /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
mpirun -np ${np} /opt/lammps/lmp_ubuntu_24Dec20_kokkos_cuda_cmake -sf gpu  -in ${scriptname}
cd ..
END
  }
 else
  {
print SH <<END;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-10Feb15/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-11Aug17/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-16Mar18/src/lmp_cuda -sf gpu  -in $scriptname
/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
#mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
cd ..
END
  }
 }
 else
 {
print SH <<END;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-10Feb15/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-11Aug17/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-16Mar18/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
cd ..
END
 }
}

#print SH "sbatch -n 1 -p gputest ompi lmp_linux -in $scriptname \n";
#print SH "sbatch -n 1014 -t 4300 -p regular4 ompi lmp_linux  -in $scriptname \n";
close(SH);

#if($Ndir == 1)
#{
open(SH,">".$shellname_all);
for($icpu=0;$icpu<$Ncpu;$icpu++)
{
$shellname1 = sprintf "run_cpu%d.sh",$icpu;  #look at init block too
print SH <<END;
nohup sh $shellname1 &
END
}
close(SH);
#}

}
