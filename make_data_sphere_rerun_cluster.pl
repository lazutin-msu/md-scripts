#!/usr/bin/perl
use Getopt::Long;
use File::Copy;
use Math::Trig;
use POSIX;
    
$dens = 8;
$N = 50;
$chains = 10; #10x10
$time = 10000000;
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
#'density|d=f' => \$dens,
#'length|n=i' => \$N,
#'chains|c=i' => \$chains,
'rsphere|a=f' => \$R_sphere,
'dump|d=s' => \$dumpfile,
'repeat|r=i' => \$repeat,
#'time|t=i' => \$time,
#'step|s=f' => \$step,
'g|gpuid=i' => \$gpu_num,
'm|max_cluster=i' => \$max_cluster,
#'p|spheres=i' => \$spheres,
#'x|dx=f' => \$spheres_dx,
#'pressure=f' => \$pressure,
#'time_deform=i' => \$time_deform,
#'rate_deform=f' => \$rate_deform,
#'dir_deform=s' => \$dir_deform,
'eps|e=f{3}' => \@epsilon,
'ncpu=i' => \$Ncpu,
'np=i' => \$num_np,
'prefix=s' => \$prefix,
'cut=f' => \$cut,
'help|h' => \$help);
if($help)
 {
 die("Usage make_data.pl options
  -f datafile initial configuration
  -d dump trajectory
#  -d density = 8
#  -c chains = 10
#  -n length = 50
  -a rsphere = 30
  -r repeat = 1
#  -t time = 10000000
#  -s step 
#  -p spheres
#  -x dx spheres dx
  -g gpuid
  -m max_cluster
#  --pressure
#  --time_deform=100000
#  --rate_deform=0.00002 
#  --dir_deform=x (x, y or z)
  -e --eps epsAA epsAB epsBB
  --prefix prefix to output files
  --cut default 1.05
  -h help\n");
 }

if(!$time_deform) {$time_deform = 100000;}
#if(!$period_deform) {$period_deform = 5000;}
if(($dir_deform ne "x")&&($dir_deform ne "y")&&($dir_deform ne "z")) {$dir_deform = "x";}
if(($num_np!=1)&&($num_np!=2)&&($num_np!=4)&&($num_np!=8)) {die "--np $num_np : should be 1, 2, 4, or 8\n";}
if(!$rate_deform){$rate_deform = 0.00002;}
if(!$cut){$cut = 1.05}

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


#if($gpuid==1){$gpuid=2;}
printf "R_sphere = %f\nchains = %d\nlength = %d \ngpuid= %d\n",$R_sphere,$chains,$N,$gpuid;

$Ndir = 1;

$desc = sprintf "brush_p%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];

if($datafile)
{
$is_datafile = 1;
}
else
{
$is_datafile = 0;
$datafile = sprintf "data_%sp%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f",$prefix,$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];
}

$dirname = sprintf "%seps_%.2f_%.2f_%.2f_cluster_n%d",$prefix,$epsilon[0],$epsilon[1],$epsilon[2],$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = sprintf "%seps_%.2f_%.2f_%.2f_cluster_n%d",$prefix,$epsilon[0],$epsilon[1],$epsilon[2],$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
#  $dirname = sprintf "%sp%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_rate_%.3e_n%d",$prefix,$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$rate_deform,$Ndir;
  $dirname = sprintf "%seps_%.2f_%.2f_%.2f_cluster_n%d",$prefix,$epsilon[0],$epsilon[1],$epsilon[2],$Ndir;

#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
#  $output_filename = sprintf "%sp%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_rate_%.3e_n%d",$prefix,$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$rate_deform,$Ndir;
  $output_filename = sprintf "%seps_%.2f_%.2f_%.2f_cluster_n%d",$prefix,$epsilon[0],$epsilon[1],$epsilon[2],$Ndir;
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
#one 100000 page 1000000

pair_style hybrid/overlay yukawa 1.2 4.0 lj/cut 1.1224620 lj/expand 1.1224620
bond_style fene 
boundary p p p
#special_bonds fene
special_bonds lj/coul 1.0 1.0 1.0

read_data $datafile

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

compute clusterA type1 aggregate/atom ${cut}
compute ccA type1 chunk/atom c_clusterA compress yes
compute sizeA type1 property/chunk ccA count
fix cluhistA type1 ave/histo 50000 10 1000000 0 ${max_cluster} ${max_cluster} c_sizeA mode vector ave one beyond ignore file ${output_filename}_cluA.txt
#fix cluhistA type1 ave/histo 50000 1 50000 0 ${max_cluster} ${max_cluster} c_sizeA mode vector ave one beyond ignore file ${output_filename}_cluA.txt

compute clusterB type2 aggregate/atom ${cut}
compute ccB type2 chunk/atom c_clusterB compress yes
compute sizeB type2 property/chunk ccB count
fix cluhistB type2 ave/histo 50000 10 1000000 0 ${max_cluster} ${max_cluster} c_sizeB mode vector ave one beyond ignore file ${output_filename}_cluB.txt
#fix cluhistB type2 ave/histo 50000 1 50000 0 ${max_cluster} ${max_cluster} c_sizeB mode vector ave one beyond ignore file ${output_filename}_cluB.txt

compute clusterall all aggregate/atom ${cut}
compute ccall all chunk/atom c_clusterall compress yes
compute sizeall all property/chunk ccall count
fix cluhistall all ave/histo 50000 10 1000000 0 ${max_cluster} ${max_cluster} c_sizeall mode vector ave one beyond ignore file ${output_filename}_cluall.txt
#fix cluhistall all ave/histo 50000 1 50000 0 ${max_cluster} ${max_cluster} c_sizeall mode vector ave one beyond ignore file ${output_filename}_cluall.txt

rerun ../${dumpfile}  dump x y z wrapped no scaled no

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
END
if($num_np==1)
{
print SH <<END;
/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
END
}
else
{
print SH <<END;
mpirun -np $num_np /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
END
}
print SH <<END;
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
