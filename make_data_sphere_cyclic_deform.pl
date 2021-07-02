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
$Ncpu = 6;
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
#'step|s=f' => \$step,
'g|gpuid=i' => \$gpu_num,
'p|spheres=i' => \$spheres,
'x|dx=f' => \$spheres_dx,
'pressure=f' => \$pressure,
'time_deform=i' => \$time_deform,
'period_deform=i' => \$period_deform,
'amp_deform=f' => \$amp_deform,
'dir_deform=s' => \$dir_deform,
'eps|e=f{3}' => \@epsilon,
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
#  -s step 
  -p spheres
  -x dx spheres dx
  -g gpuid
  --pressure
  --time_deform=100
  --period_deform=5000 (tau)
  --amp_deform=0.05 (x l_box)
  --dir_deform=x (x, y or z)
  -e --eps epsAA epsAB epsBB
  -h help\n");
 }

if(!$time_deform) {$time_deform = 100;}
if(!$period_deform) {$period_deform = 5000;}
if(!$amp_deform) {$amp_deform = 0.05;}
if(($dir_deform ne "x")&&($dir_deform ne "y")&&($dir_deform ne "z")) {$dir_deform = "x";}

print STDERR "eps = ".$epsilon[0]." ".$epsilon[1]." ".$epsilon[2]."\n";

$time_lastdump = $time / 2;
$time_dump = $time / 20;

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


if($gpuid==1){$gpuid=2;}
printf "R_sphere = %f\nchains = %d\nlength = %d \ngpuid= %d\n",$R_sphere,$chains,$N,$gpuid;

$Ndir = 1;

$desc = sprintf "brush_p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];

if($datafile)
{
$is_datafile = 1;
}
else
{
$is_datafile = 0;
$datafile = sprintf "data_p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];
}

$dirname = sprintf "p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_period%d_amp_%.3f_n%d",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$period_deform,$amp_deform,$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = sprintf "p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_period%d_amp_%.3f_n%d",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$period_deform,$amp_deform,$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  $dirname = sprintf "p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_period%d_amp_%.3f_n%d",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$period_deform,$amp_deform,$Ndir;
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  $output_filename = sprintf "p%d_c%d_R%f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_period%d_amp_%.3f_n%d",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$period_deform,$amp_deform,$Ndir;
  }
 }

`mkdir $dirname`;

$dirname = $dirname."/";

#$dx = sqrt(1.0 / $dens);
$dx = $dens; # dens now size of grafting points' cell
$radius = 0.5;
$bond = 2.0*$radius;
if($dx<(2.0*$bond)){$bond=$dx/2.0;}

$delta_rho = $bond;

$natoms = 2 * $chains * $N + 1;
$nbonds = (2*$N - 1) * $chains;

$totalatoms = $natoms * $spheres;
$totalbonds = $nbonds * $spheres;

#$lx =  ( 2 * $N + $R_sphere ) +2;
#$ly =  ( 2 * $N + $R_sphere ) +2;
#$lz =  ( 2 * $N + $R_sphere ) +2;
$lx =  ( $N + $R_sphere ) +2;
$ly =  ( $N + $R_sphere ) +2;
$lz =  ( $N + $R_sphere ) +2;


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
   if($imono!=0)
    {
    $at1_bond[$ibond] = $natoms * $isphere + $iatom - 2;
    $at2_bond[$ibond] = $natoms * $isphere + $iatom;
    $ibond++;
    }
   $at1_bond[$ibond] = $natoms * $isphere + $iatom;
   $at2_bond[$ibond] = $natoms * $isphere + $iatom+1;
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

if(!$is_datafile)
{
print "writing to $datafile \n";
open(OUT,">".$dirname.$datafile) or die "cant create file ".$dirname.$datafile;

print OUT $desc."\n\n";
printf OUT "%11d atoms\n", $totalatoms;
printf OUT "%11d bonds\n", $totalbonds;

print OUT <<END;

         4 atom types
         1 bond types

END

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

#print OUT "\nBond Coeffs\n\n";
#print OUT "1 100 1.0\n";

print OUT "\nAtoms\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %6d %6d %8.2f %8.2f %8.2f\n",$iatom,$mol[$iatom],$type[$iatom],$x[$iatom],$y[$iatom],$z[$iatom];
 }

print OUT "\nVelocities\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %8.2f %8.2f %8.2f\n",$iatom,$vx[$iatom],$vy[$iatom],$vz[$iatom];
 }

print OUT "\nBonds\n\n";

for($ibond=1;$ibond<=$totalbonds;$ibond++)
 {
 printf OUT "%6d %6d %6d %6d\n",$ibond,1,$at1_bond[$ibond],$at2_bond[$ibond];
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

#@epsAA = (  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0);
#@epsAB = (  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 , 5.0, 5.0 , 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0,  5.0, 5.0);
#@epsBB = ( -0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0,-4.5,-5.0,-5.25,-5.5,-5.75,-6.0,-6.25,-6.5,-6.75,-7.0,-7.25,-7.5,-7.75,-8.0,-8.25,-8.5,-8.75,-9.0,-9.25,-9.5,-9.75, -10);

@epsAA = ( -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,   0.0 );
@epsAB = (   5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,   5.0 );
@epsBB = (   0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0, -6.5, -7.0, -7.5, -8.0, -8.5, -9.0, -9.5, -10.0 );

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

$xlo_bound = 0 + 3;
$xhi_bound = $lx - 3;
$ylo_bound = 0 + 3;
$yhi_bound = $ly - 3;

$time_energy = $time /4;
$step_energy = 100;
$num_energy = $time_energy / $step_energy;

$time_msd = $time;
$step_msd = 1000;
$num_msd = $time_msd / $step_msd;

$time_rgyr = 1000000;
$step_rgyr = 1000;
$num_rgyr = $time_rgyr / $step_rgyr;

$time_clust = 1000000;
$step_clust = 1000;
$num_clust = $time_clust / $step_clust;

$timestep = 0.005; 

$time_modul = $time_deform;
$step_modul = 10;
$num_modul = $time_modul / $step_modul;


$set_temp = 1.0;

$R_sphere_cutoff = 1.1224620 * $R_sphere;
$R_sphere2 = $R_sphere - 0.5;
$hi_bound = $R_sphere + $N;

open(SCRIPT,">".$dirname.$scriptname);
print SCRIPT <<END;
package gpu 1 neigh no gpuID $gpuid $gpuid
suffix gpu
units lj
atom_style bond
neighbor        1.5 bin
neigh_modify    every 5 delay 5 check yes

pair_style hybrid/overlay yukawa 1.2 4.0 lj/cut 1.1224620 lj/expand 1.1224620
bond_style fene 
boundary p p p
#special_bonds fene
special_bonds lj/coul 1.0 1.0 1.0

read_data $datafile

END

if(!datafile)
{
print SCRIPT <<END;
velocity all create ${set_temp} $seed2
END
}
else
{
print SCRIPT <<END;
#velocity all create ${set_temp} $seed2
END
}
print SCRIPT <<END;

pair_coeff 1*3 1*3 lj/cut 1.0 1.0 1.1224620
pair_coeff * 4 lj/expand 1.0 1.0 ${R_sphere2} 1.1224620
pair_coeff 3 4 lj/expand 0.0 1.0 ${R_sphere2} 1.1224620
pair_coeff 1 1 yukawa $epsilon[0]
pair_coeff 1 3 yukawa $epsilon[0]
pair_coeff 3 3 yukawa $epsilon[0]
pair_coeff 1 2 yukawa $epsilon[1]
pair_coeff 2 3 yukawa $epsilon[1]
pair_coeff 2 2 yukawa $epsilon[2]


bond_coeff 1 30.0 1.5 0.0 1.0

group roots type 3 4
group rest type 1 2
group type2 type 2
group type1 type 1 3



#velocity roots zero linear

#fix root roots move linear 0 0 0 

fix lan rest langevin ${set_temp} ${set_temp} 100.0 $seed
fix finve rest nve 

END

#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT "group sphere".$igroup." id ".($natoms * $igroup + 1).":".($natoms * ($igroup + 1))."\n";
#}
#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT "group core".$igroup." intersect roots sphere".$igroup."\n";
#}
#print SCRIPT "fix firig roots rigid/nve group ".$spheres;

open(CORE,">".$dirname."coreid.txt");

print CORE $totalatoms."\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
# if($coreid[$iatom]>0)
#  {
  printf CORE "%d %d\n", $iatom,$coreid[$iatom];
#  }
 }

close(CORE);

open(PART,">".$dirname."partnum.txt");

print PART $totalatoms."\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
# if($coreid[$iatom]>0)
#  {
  printf PART "%d %d\n", $iatom,$partnum[$iatom];
#  }
 }

close(PART);


print SCRIPT "variable coreid atomfile coreid.txt\n";

print SCRIPT "fix firig roots rigid/nve custom v_coreid langevin ${set_temp} ${set_temp} 100.0 $seed \n";

#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT " core".$igroup;
#}
#print SCRIPT "\n";

print SCRIPT <<END;
#fix lanri roots langevin ${set_temp} ${set_temp} 100.0 $seed

run 0


#min_style quickmin 
#timestep 0.0001

#dump dumpmin all atom 1 ${output_filename}.minimize.lammpstrj.gz

#minimize 0.0 1.0 1000 1000

END

if(!$is_datafile)
{
print SCRIPT <<END;
timestep 0.002

dump dump_cluster3 all custom 10000 ${output_filename}.rel2.dump.gz id type mol xu yu zu ix iy iz

run 100000
timestep 0.003
run  50000

#fix pres all press/berendsen iso 10.0 10.0 100.0

#fix pres all press/berendsen iso 100.0 100.0 100.0
fix pres all press/berendsen iso ${pressure} ${pressure} 100.0

thermo 1000
thermo_style    custom step time temp ke pe evdwl ebond etotal press density lx vol

run  300000

unfix pres

undump dump_cluster3

reset_timestep 0
END
}

print SCRIPT <<END;


variable partnum atomfile partnum.txt
compute part all chunk/atom v_partnum
compute msd all msd/chunk part 
fix fixmsd all ave/time ${step_msd} 1 ${step_msd} c_msd[1] file msd1.txt mode vector 
fix fixmsd2 all ave/time ${step_msd} 1 ${step_msd} c_msd[2] file msd2.txt mode vector 
fix fixmsd3 all ave/time ${step_msd} 1 ${step_msd} c_msd[3] file msd3.txt mode vector 
fix fixmsd4 all ave/time ${step_msd} 1 ${step_msd} c_msd[4] file msd4.txt mode vector 

compute gyr all gyration/chunk part
fix fixgyr all ave/time ${step_rgyr} 1 ${step_rgyr} c_gyr file gyr1.txt mode vector
fix fixgyr2 all ave/time ${step_rgyr} ${num_rgyr} ${time_rgyr} c_gyr file gyr2.txt mode vector
#fix gyrhist all ave/histo 1000 1 1000 0 50 500 c_gyr  mode vector ave one  beyond end file gyrhist.txt
#fix gyrhist2 all ave/histo 1000 500 1000000 0 50 500 c_gyr  mode vector ave one  beyond end file gyrhist2.txt

compute cluster3 all aggregate/atom 1.0
compute cc3 all chunk/atom c_cluster3 compress yes
compute size3 all property/chunk cc3 count
#fix cluhist3 all ave/histo 1000 1 1000 0 100000 100000 c_size3 mode vector ave one beyond ignore file newcluster3.txt
fix cluhist3b all ave/histo ${step_clust} ${num_clust} ${time_clust} 0 100000 100000 c_size3 mode vector ave one beyond ignore file newcluster3b.txt

compute yukpair all pair yukawa
compute ljpair all pair lj/cut
thermo 100
thermo_style    custom step time temp ke pe evdwl c_yukpair c_ljpair ebond  etotal press lx ly lz pxx pyy pzz vol density

compute binsphere all chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
compute binspherea type1 chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
compute binsphereb type2 chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
fix spatall all ave/chunk ${step_energy} ${num_energy} ${time_energy} binsphere density/number file ${output_filename}_spat_all.txt ave one

fix spat all ave/chunk ${step_energy} ${num_energy} ${time} binsphere density/number file ${output_filename}_spat.txt ave one
fix spata type1 ave/chunk ${step_energy} ${num_energy} ${time} binspherea density/number file ${output_filename}_spatA.txt ave one
fix spatb type2 ave/chunk ${step_energy} ${num_energy} ${time} binsphereb density/number file ${output_filename}_spatB.txt ave one

dump dump_cluster all custom ${time_dump} ${path_to_dump}${output_filename}.dump.gz id type mol xu yu zu ix iy iz

# heat capacity
compute peall all pe 
compute keall all ke 
variable natoms equal atoms
variable settemp equal ${set_temp}
variable enall equal c_peall+c_keall
variable pe2 equal c_peall*c_peall
variable en2 equal v_enall*v_enall
#fix ener all ave/time ${step_energy} ${num_energy} ${time_energy} c_peall v_pe2 v_enall v_en2 file  ${output_filename}_energy.txt ave one
fix enout all ave/time ${step_energy} 1 ${step_energy} c_peall c_keall v_enall file  ${output_filename}_raw_en.txt ave one


END

print SCRIPT <<END;
fix ener1 all ave/time ${step_energy} ${num_energy} ${time_energy} c_peall v_pe2 v_enall v_en2 ave one
variable peavg1 equal f_ener1[1]/v_natoms
variable cvp1 equal sqrt(f_ener1[2]-f_ener1[1]*f_ener1[1])/v_natoms/(v_settemp*v_settemp)
variable enavg1 equal f_ener1[3]/v_natoms
variable cve1 equal sqrt(f_ener1[4]-f_ener1[3]*f_ener1[3])/v_natoms/(v_settemp*v_settemp)
END

print SCRIPT "fix eneravg all ave/time ${time_energy} 1 ${time_energy} v_peavg1 v_cvp1 v_enavg1 v_cve1 file ${output_filename}_energy.txt ave one\n";

print SCRIPT <<END;

timestep 0.005

fix bal all balance 50000 1.0 shift xyz 10 1.1

dump dumplast all atom ${time_lastdump} ${output_filename}.last.lammpstrj.gz

variable tmp equal "lx"
variable L0 equal \${tmp}
variable Amp equal ${amp_deform}
variable AmpL equal \${tmp}*\${Amp}
print "Initial Length, L0: \${L0}"
print "Amplitude stain, Amp: \${Amp}"
print "Amplitude Length : \${AmpL}"

#fix def all deform 100000 x erate 0.00002 y volume z volume
#fix def all deform 50000 x wiggle 0.05 5000 y volume z volume units lattice
#fix def all deform 50000 x wiggle ${AmpL} 5000 units box
END

$str_deform = "fix def all deform ${time_deform}";

foreach $dir ("x", "y", "z")
{
if($dir eq $dir_deform)
 {
 $str_deform = $str_deform . " ".$dir." wiggle \${AmpL} ${period_deform}";
 }
 else
 {
 $str_deform = $str_deform . " ".$dir." volume";
 }
}
$str_deform = $str_deform . " units box\n";

print SCRIPT $str_deform;

print SCRIPT <<END;
#fix press all press/berendsen y 0 0 5 z 0 0 5

variable strain equal "(l${dir_deform} - v_L0)/v_L0"
variable length equal "l${dir_deform}"

thermo_style    custom step time temp ke pe evdwl c_yukpair c_ljpair ebond  etotal press v_strain pxx pyy pzz

variable px equal "pxx"
variable py equal "pyy"
variable pz equal "pzz"

#fix modul all ave/time 100000 500 100 v_strain pxx pyy pzz file module.txt
#fix modul all ave/time 100000 100000 100 v_strain v_px v_py v_pz file module.txt
#fix modul all ave/time 100000 100 1000 pxx pyy pzz  file module.txt
#fix modul all ave/time 100 1 100 v_length v_strain v_px v_py v_pz file module.txt
fix modul all ave/time ${step_modul} 1 ${step_modul} v_length v_strain v_px v_py v_pz file ${output_filename}_module.txt
fix module all ave/time ${step_modul} ${num_modul} ${time_modul} v_length v_strain v_px v_py v_pz file ${output_filename}_module2.txt
#fix modulee all ave/time 1 100000 100000 v_length v_strain v_px v_py v_pz file module3.txt


run $time
write_data ${output_filename}_deform.data

END

close(SCRIPT);

open(SH,">>".$shellname);
#print SH "nohup /usr/mount-opt/lammps/lmp_ubuntu_1Feb14 <$scriptname  &";
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
