#!/usr/bin/perl
use Getopt::Long;
use File::Copy;
use Math::Trig;
use POSIX;
    
$dens = 3.0; # hex lattice size
$N = 50;
#$chains = 10; #10x10
$time = 2000000;
$repeat = 1;
$R_sphere = 5;
$aspect_ratio = 1;

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
'aspect=i' => \$aspect_ratio,
'repeat|r=i' => \$repeat,
'time|t=i' => \$time,
'step|s=f' => \$step,
'g|gpuid=i' => \$gpu_num,
#'p|spheres=i' => \$spheres,
'x|dx=f' => \$spheres_dx,
'pressure=f' => \$pressure,
'structure=s' => \$structure,
'stiff=f' => \$k_stiff,
'fix_particle' => \$fix_particle,
'lomo' => \$lomo_flag,
'fix_eps=f{3}' => \@fix_eps,
'step_dump=i' => \$step_dump,
'comm_mod=f' => \$comm_mod,
'np=i' => \$np,
'structure_step=f' => \$structure_step,
'short' => \$short_flag,
'Ncpu=i' => \$Ncpu,
'help|h' => \$help);
if($help)
 {
 die("Usage make_data.pl options
  -f datafile initial configuration
  -d density - hex lattice size = 3
  -c chains = 10
  -n length = 50
  -a rsphere = 5
  -r repeat = 1
  -t time = 2000000
  -s step 
#  -p spheres
  -x dx spheres dx
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

#if($structure ne "PH") {$structure = "HP";}

if(!$Ncpu){$Ncpu = 6;}

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

$R_sphere_wall = $R_sphere;

#$R_sphere = $R_sphere + 0.5;

if($k_stiff)
 {
 $name_prefix = "mov"
 }
 else
 {
 $name_prefix = "fix"
 }

$m_real = pi / asin( 0.5* $dens/$R_sphere );
$m_int = int($m_real + 0.5);

$drad =  2.0 * pi / $m_int;
$dx = 2.0 * $R_sphere * sin( $drad * 0.5 );

$dz = $dx * sqrt(3.0) * 0.5;
$m_z_real = $R_sphere / $dz;
$m_z_int = $aspect_ratio * 2 * int($m_z_real + 0.5);


#if($gpuid==1){$gpuid=2;}
#printf "R_sphere = %f\nchains = %d\nlength = %d \ngpuid= %d\n",$R_sphere,$chains,$N,$gpuid;
printf "R_sphere = %f\nchain_z = %d\nchain_rad = %d\nlength = %d \ngpuid= %d\n",$R_sphere_wall,$m_z_int,$m_int,$N,$gpuid;

$Ndir = 1;

#$desc = sprintf "brush_p%d_c%d_R%f_N%d_str%s_press%.2f",$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2];
$desc = sprintf "brush_cz%d_crad%d_R%f_N%d_str%s_press%.2f",$m_z_int,$m_int,$R_sphere_wall,$N,$epsilon[0],$epsilon[1],$epsilon[2];

if($datafile)
{
$is_datafile = 1;
}
else
{
$is_datafile = 0;
#$datafile = sprintf "data_p%d_c%d_R%f_N%d_str%s_press%.2f",$spheres,$chains,$R_sphere,$N,$structure,$pressure;
$datafile = sprintf "data_cz%d_crad%d_R%f_N%d_str%s_press%.2f",$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure;
}

if(scalar(@fix_eps))
 {
#$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
$dirname = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
#$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
$output_filename = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  #$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
  $dirname = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  #$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
  $output_filename = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_fixeps_%.2f_%.2f_%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$fix_eps[0],$fix_eps[1],$fix_eps[2],$Ndir;
  }
 }
 }
 else
 {
#$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
$dirname = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
#$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
$output_filename = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  #$dirname = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
  $dirname = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$Ndir;
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  #$output_filename = sprintf "p%d_c%d_R%f_N%d_str%s_press%.2f_n%d",$spheres,$chains,$R_sphere,$N,$structure,$pressure,$Ndir;
  $output_filename = sprintf "%s_cz%d_crad%d_R%f_N%d_str%s_press%.2f_n%d",$name_prefix,$m_z_int,$m_int,$R_sphere_wall,$N,$structure,$pressure,$Ndir;
  }
 }
 
 } #fix_eps

`mkdir $dirname`;

$dirname = $dirname."/";

#$dx = sqrt(1.0 / $dens);
#$dx = $dens; # dens now size of grafting points' cell



$radius = 0.5;
$bond = 2.0*$radius;
if($dx<(2.0*$bond)){$bond=$dx/2.0;}

$delta_rho = $bond;


print "m_int $m_int m_z_int $m_z_int dx $dx dz $dz \n";

$chains = $m_z_int * $m_int;

#$natoms = 2 * $chains * $N + 1;
$natoms = 2 * $chains * $N ;

if(!$k_stiff)
 {
$nbonds = (2*$N - 1) * $chains;
 }
 else
 {
#$nbonds = (2*$N - 1) * $chains + $chains ;
$nbonds = (2*$N - 1) * $chains ;
 }

#$totalatoms = $natoms * $spheres;
#$totalbonds = $nbonds * $spheres;
$totalatoms = $natoms ;
$totalbonds = $nbonds ;

#$lx =  ( 2 * $N + $R_sphere ) +2;
#$ly =  ( 2 * $N + $R_sphere ) +2;
#$lz =  ( 2 * $N + $R_sphere ) +2;
$lx =  ( $N + $R_sphere ) +2;
$ly =  ( $N + $R_sphere ) +2;


#$lz =  ( $N + $R_sphere ) +2;
$lz = $m_z_int * $dz;

print "$desc \n";
$ibond = 1;

#$Nspheres = ceil( $spheres ** (1/3) );
$Nspheres = 1;
#$totallx = $Nspheres * $spheres_dx + 2* $lx;
#$totally = $Nspheres * $spheres_dx + 2* $ly;
#$totallz = $Nspheres * $spheres_dx + 2* $lz;
#$totallx = ($Nspheres - 1) * $spheres_dx + 2* $lx;
#$totally = ($Nspheres - 1) * $spheres_dx + 2* $ly;
#$totallz = ($Nspheres - 1) * $spheres_dx + 2* $lz;

$totallx = 2* $lx;
$totally = 2* $ly;
$totallz = $lz;


$isphere = 0;

SPHERES: for($iz=0;$iz<$Nspheres;$iz++)
 {
 for($iy=0;$iy<$Nspheres;$iy++)
  {
  for($ix=0;$ix<$Nspheres;$ix++)
   {
   

#$sp_dlong = pi*(3.0-sqrt(5.0)); #  /* ~2.39996323 */
#$sp_dz    = 2.0/$chains;
#$sp_long = 0;
#$z    = 1.0 - $sp_dz/2.0;

#$delta_phi = 1/ $R_sphere; 

#  $x[$natoms * $isphere + 1] = $ix * $spheres_dx + $lx +    0;
#  $y[$natoms * $isphere + 1] = $iy * $spheres_dx + $ly +    0;
#  $z[$natoms * $isphere + 1] = $iz * $spheres_dx + $lz +    0;
#  $type[$natoms * $isphere + 1] = 4;
#  $mol[$natoms * $isphere + 1] = ($chains + 1) * $isphere + 1;
#  $coreid[$natoms * $isphere + 1] = $isphere + 1;
#  $partnum[$natoms * $isphere + 1] = $isphere + 1;

#for($k=0;$k<$chains;$k++)
for($kz=0;$kz<$m_z_int;$kz++)
 {
 if($kz % 2 == 0)
  {
  $phi0 = 0.;
  }
  else
  {
  $phi0 = $drad * 0.5;
  }
  
 for($krad=0;$krad<$m_int;$krad++)
   {
#    $r    = sqrt(1.0-$z*$z);
#    ($xchain, $ychain, $zchain) = (cos($long)*$r, sin($long)*$r, $z);
    
    $phi = $phi0 + $drad * $krad;
    $xchain = $R_sphere * cos( $phi );
    $ychain = $R_sphere * sin( $phi );
    $zchain = $kz * $dz;
#   my $rho   = sqrt($r*$r+$z*$z)*$R_sphere;
#   my $phi   = $sp_long;
#   my $theta = atan2($r,$z);
#    $z    = $z - $sp_dz;
#    $sp_long = $sp_long + $sp_dlong;
#   $rho = $rho + $delta_rho * 0.5;

  for($imono=0;$imono<$N;$imono++)
   {
#   $iatom = 2*$imono + 2*$k * $N + 1;
   $iatom = 2*$imono + ($kz*$m_int + $krad) * 2 * $N + 1;

#   $x[$natoms * $isphere + $iatom] = $ix * $spheres_dx + $lx + $rho * sin($theta) * cos($phi);
#   $y[$natoms * $isphere + $iatom] = $iy * $spheres_dx + $ly + $rho * sin($theta) * sin($phi);
#   $z[$natoms * $isphere + $iatom] = $iz * $spheres_dx + $lz + $rho * cos($theta) ;
   $rho_mono = $R_sphere + $bond * $imono;
   $x[$natoms * $isphere + $iatom] = $lx + $rho_mono * cos($phi);
   $y[$natoms * $isphere + $iatom] = $ly + $rho_mono * sin($phi);
#   $z[$natoms * $isphere + $iatom] = $lz + $zchain ;
   $z[$natoms * $isphere + $iatom] = $zchain ;

#   if($z[$iatom]>0)
#    {
#    $my_delta_theta = $bond / ($R_sphere + $imono * $delta_rho);
#    }
#    else
#    {
#    $my_delta_theta = -$bond / ($R_sphere  + $imono * $delta_rho);
#    }
   $my_delta_phi = atan( $bond / $rho_mono );

   $type[$natoms * $isphere + $iatom] = 1;
   if($imono==0){$type[$natoms * $isphere + $iatom] = 3; $coreid[$natoms * $isphere + $iatom] = $isphere + 1;}
   $mol[$natoms * $isphere + $iatom] = ($chains + 1) * $isphere + ($kz*$m_int + $krad) + 1;
   $partnum[$natoms * $isphere + $iatom] = $isphere + 1;


#   $x[$natoms * $isphere + $iatom+1] = $ix * $spheres_dx + $lx + $rho * sin($theta + $my_delta_theta) * cos($phi);
#   $y[$natoms * $isphere + $iatom+1] = $iy * $spheres_dx + $ly + $rho * sin($theta + $my_delta_theta) * sin($phi);
#   $z[$natoms * $isphere + $iatom+1] = $iz * $spheres_dx + $lz + $rho * cos($theta + $my_delta_theta);

   $x[$natoms * $isphere + $iatom+1] = $lx + $rho_mono * cos($phi + $my_delta_phi);
   $y[$natoms * $isphere + $iatom+1] = $ly + $rho_mono * sin($phi + $my_delta_phi);
#   $z[$natoms * $isphere + $iatom+1] = $lz + $zchain;
   $z[$natoms * $isphere + $iatom+1] = $zchain;

#   $rho = $rho + $delta_rho; 
   $type[$natoms * $isphere + $iatom+1] = 2;
   $mol[$natoms * $isphere + $iatom+1] = ($chains + 1) * $isphere + ($kz*$m_int + $krad) + 2;
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
#    else
#    {
#    $at1_bond[$ibond] = $natoms * $isphere + 1;
#    $at2_bond[$ibond] = $natoms * $isphere + $iatom;
#    $bond_type[$ibond] = 2;
#    $ibond++;
#    }
 }
   $at1_bond[$ibond] = $natoms * $isphere + $iatom;
   $at2_bond[$ibond] = $natoms * $isphere + $iatom+1;
   $bond_type[$ibond] = 1;
   $ibond++;
   } # imono
   } # krad
   } # kz

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

         3 atom types
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
         1 bond types

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
#print OUT "4 1.0\n";

#print OUT "\nBond Coeffs\n\n";
#print OUT "1 100 1.0\n";

print OUT "\nAtoms # bond\n\n";

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
elsif($structure eq "PH2")
{
@epsAA = ();
@epsAB = ();
@epsBB = ();
for($istep=0;$istep<=int(10.0/$structure_step+0.5);$istep++)
 {
 push @epsAA, 0.0;
 push @epsAB, 5.0;
 push @epsBB, -1*$istep*$structure_step;
 }
}
elsif($structure eq "PH3")
{
@epsAA = ();
@epsAB = ();
@epsBB = ();
for($istep=0;$istep<=int(10.0/$structure_step+0.5);$istep++)
 {
 push @epsAA, 0.0;
 push @epsAB, 5.0;
 push @epsBB, -1*$istep*$structure_step;
 }
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

$xlo_bound = 0 + 3;
$xhi_bound = $lx - 3;
$ylo_bound = 0 + 3;
$yhi_bound = $ly - 3;

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
package gpu 1 neigh no gpuID $gpuid $gpuid
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
#comm_modify cutoff ${comm_mod}
#newton off
END
 }
}

print SCRIPT <<END;
atom_style bond
neighbor        1.5 bin
neigh_modify    every 5 delay 5 check yes one 10000

#pair_style hybrid/overlay yukawa 1.2 4.0 lj/cut 1.1224620 lj/expand 1.1224620
pair_style hybrid/overlay yukawa 1.2 4.0 lj/cut 1.1224620 
END

if(!$k_stiff)
{
print SCRIPT <<END;
bond_style fene 
END
}
else
{
#print SCRIPT <<END;
#bond_style hybrid fene harmonic 
#END

print SCRIPT <<END;
bond_style fene  
END

}

print SCRIPT <<END;
boundary p p p
#special_bonds fene
special_bonds lj/coul 1.0 1.0 1.0


read_data $datafile

velocity all create ${set_temp} $seed2

pair_coeff 1*3 1*3 lj/cut 1.0 1.0 1.1224620
#pair_coeff * 4 lj/expand 1.0 1.0 ${R_sphere2} 1.1224620
#pair_coeff 3 4 lj/expand 0.0 1.0 ${R_sphere2} 1.1224620
pair_coeff 1 1 yukawa $epsAA[0]
pair_coeff 1 3 yukawa $epsAA[0]
pair_coeff 3 3 yukawa $epsAA[0]
pair_coeff 1 2 yukawa $epsAB[0]
pair_coeff 2 3 yukawa $epsAB[0]
pair_coeff 2 2 yukawa $epsBB[0]


END
if(!$k_stiff)
{
print SCRIPT <<END;
bond_coeff 1 30.0 1.5 0.0 1.0
END
}
else
{
$R_bond = $R_sphere + 0.5;
#print SCRIPT <<END;
#bond_coeff 1 fene 30.0 1.5 0.0 1.0
#bond_coeff 2 harmonic ${k_stiff} ${R_bond}
#END
print SCRIPT <<END;
bond_coeff 1 30.0 1.5 0.0 1.0
END
}

if($k_stiff)
 {
print SCRIPT <<END;

#group roots type 4
group rest type 1 2 3
group type2 type 2
group type1 type 1 3
group walled type 3

region mywall cylinder z ${lx} ${ly} ${R_sphere2} INF INF side out
region mywall2 cylinder z ${lx} ${ly} ${R_sphere2} INF INF side out
fix repwall all wall/region mywall lj126 1.0 1.0 1.1224620
fix bondwall walled  wall/region mywall2 harmonic ${k_stiff} 0.5 5.5
fix_modify energy yes

END
 }
 else
 {
print SCRIPT <<END;

#group roots type 3 4
group roots type 3 
group rest type 1 2
group type2 type 2
group type1 type 1 3

region mywall cylinder z ${lx} ${ly} ${R_sphere2} 0.0 ${totallz} side out
fix repwall all wall/region mywall lj126 1.0 1.0 1.1224620

END
 }

#if($k_stiff)
#{
#print SCRIPT <<END;
#compute bty all property/local batom1 batom2 btype
#compute dis all bond/local dist
#compute dis_min all reduce min c_dis[*]
#compute dis_max all reduce max c_dis[*]
#
#thermo 1000
#thermo_style    custom step time temp ke pe evdwl ebond  etotal press c_dis_min c_dis_max
#END
#}

print SCRIPT <<END;
#velocity roots zero linear

#fix root roots move linear 0 0 0 
END

if(!$k_stiff)
{
print SCRIPT <<END;
fix lan rest langevin ${set_temp} ${set_temp} 100.0 $seed
fix finve rest nve 
velocity roots zero linear

fix root roots move linear 0 0 0

END
}
else
{
print SCRIPT <<END;
fix lan rest langevin ${set_temp} ${set_temp} 100.0 $seed
fix finve rest nve 

END

}

#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT "group sphere".$igroup." id ".($natoms * $igroup + 1).":".($natoms * ($igroup + 1))."\n";
#}
#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT "group core".$igroup." intersect roots sphere".$igroup."\n";
#}
#print SCRIPT "fix firig roots rigid/nve group ".$spheres;

if(!$k_stiff)
{
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
}

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

#if(!$k_stiff)
#{
#print SCRIPT "variable coreid atomfile coreid.txt\n";
#
#print SCRIPT "fix firig roots rigid/nve custom v_coreid langevin ${set_temp} ${set_temp} 100.0 $seed \n";
#}

#for($igroup=0;$igroup<$spheres;$igroup++)
#{
#print SCRIPT " core".$igroup;
#}
#print SCRIPT "\n";

print SCRIPT <<END;
#fix lanri roots langevin ${set_temp} ${set_temp} 100.0 $seed

run 0


#fix zlowall all wall/lj93 zlo EDGE 2.0 1.0 0.8583742
#fix zhiall all wall/lj93 zhi EDGE 2.0 1.0 0.8583742

#region mysphere sphere 0 0 0 ${R_sphere} side out
#fix sphwall all wall/region mysphere lj93 2.0 1.0 0.8583742

END

if(!$is_datafile)
{
print SCRIPT <<END;

min_style quickmin 
#timestep 0.01
timestep 0.0001

#dump dumpmin all atom 1 ${output_filename}.minimize.lammpstrj.gz

minimize 0.0 1.0 1000 1000

#undump dumpmin

#write_data ${output_filename}_minimize.data

#bond_coeff 1 30.0 1.5 0.0 1.0

timestep 0.002

dump dump_cluster3 all custom 10000 ${output_filename}.rel2.dump.gz id type mol xu yu zu ix iy iz

run 100000
timestep 0.003
run  50000

END

if($spheres>1)
{
print SCRIPT <<END;

#fix pres all press/berendsen iso 10.0 10.0 100.0

#fix pres all press/berendsen iso 100.0 100.0 100.0
fix pres all press/berendsen iso ${pressure} ${pressure} 100.0

thermo 1000
thermo_style    custom step time temp ke pe evdwl ebond etotal press density lx vol 

run  300000

unfix pres
END
}

print SCRIPT <<END;

undump dump_cluster3

reset_timestep 0
END
}

print SCRIPT <<END;

#bond_coeff 1 30.0 1.5 0.0 1.0

#region cut block ${xlo_bound} ${xhi_bound} ${ylo_bound} ${yhi_bound} INF INF side in units lattice
#group gcut dynamic type1 region cut every ${time_dump}

variable partnum atomfile partnum.txt
compute part all chunk/atom v_partnum
compute msd all msd/chunk part 
fix fixmsd all ave/time 1000 1 1000 c_msd[1] file msd1.txt mode vector 
fix fixmsd2 all ave/time 1000 1 1000 c_msd[2] file msd2.txt mode vector 
fix fixmsd3 all ave/time 1000 1 1000 c_msd[3] file msd3.txt mode vector 
fix fixmsd4 all ave/time 1000 1 1000 c_msd[4] file msd4.txt mode vector 

compute gyr all gyration/chunk part
fix fixgyr all ave/time 1000 1 1000 c_gyr file gyr1.txt mode vector
fix fixgyr2 all ave/time 1000 500 1000000 c_gyr file gyr2.txt mode vector
#fix gyrhist all ave/histo 1000 1 1000 0 50 500 c_gyr  mode vector ave one  beyond end file gyrhist.txt
#fix gyrhist2 all ave/histo 1000 500 1000000 0 50 500 c_gyr  mode vector ave one  beyond end file gyrhist2.txt

compute cluster all aggregate/atom 1.1
compute cc1 all chunk/atom c_cluster compress yes
compute size all property/chunk cc1 count
fix cluhist all ave/histo  5000  50  250000 0 100000 100000 c_size mode vector ave one beyond ignore file newcluster.txt
fix cluhistb all ave/histo 5000 100 1000000 0 100000 100000 c_size mode vector ave one beyond ignore file newclusterb.txt

compute cluster2 all aggregate/atom 1.05
compute cc2 all chunk/atom c_cluster2 compress yes
compute size2 all property/chunk cc2 count
fix cluhist2 all ave/histo  5000  50  250000 0 100000 100000 c_size2 mode vector ave one beyond ignore file newcluster2.txt
fix cluhist2b all ave/histo 5000 100 1000000 0 100000 100000 c_size2 mode vector ave one beyond ignore file newcluster2b.txt

compute cluster3 all aggregate/atom 1.0
compute cc3 all chunk/atom c_cluster3 compress yes
compute size3 all property/chunk cc3 count
fix cluhist3  all ave/histo 5000  50  250000 0 100000 100000 c_size3 mode vector ave one beyond ignore file newcluster3.txt
fix cluhist3b all ave/histo 5000 100 1000000 0 100000 100000 c_size3 mode vector ave one beyond ignore file newcluster3b.txt


compute yukpair all pair yukawa
compute ljpair all pair lj/cut
END
if(!$k_stiff)
{
print SCRIPT <<END;
thermo 1000
thermo_style    custom step time temp ke pe evdwl c_yukpair c_ljpair ebond  etotal press 
END
}
else
{
#print SCRIPT <<END;
#thermo 1000
#thermo_style    custom step time temp ke pe evdwl c_yukpair c_ljpair ebond  etotal press  c_dis_min c_dis_max
#END
print SCRIPT <<END;
thermo 1000
thermo_style    custom step time temp ke pe evdwl c_yukpair c_ljpair ebond  etotal press 
END
}
print SCRIPT <<END;
#compute gyrtensor all gyration/molecule tensor 
#compute rgyr all gyration/molecule

#compute chunkmol all chunk/atom molecule nchunk once ids once compress yes
#compute rgyr all gyration/chunk chunkmol

#fix outgyr all ave/time 1000 1 1000 c_gyrtensor mode vector file ${output_filename}_tensor.txt
#fix outrgyr all ave/time 1000 1 1000 c_rgyr mode vector file ${output_filename}_rgyr.txt

#variable rgyravg equal ave(c_rgyr)
#fix outrgyravg all ave/time 1000 1 1000 v_rgyravg  file ${output_filename}_rgyravg.txt

#compute crdf all rdf ${step_energy} * *
#fix outcrdf all ave/time ${step_energy} ${num_energy} ${time_energy} c_crdf[*] file ${output_filename}_rdf.txt mode vector
#compute crdf11 all rdf ${step_energy} 1 1
#fix outcrdf11 all ave/time ${step_energy} ${num_energy} ${time_energy} c_crdf11[*] file ${output_filename}_rdf11.txt mode vector
#compute crdf12 all rdf ${step_energy} 1 2
#fix outcrdf12 all ave/time ${step_energy} ${num_energy} ${time_energy} c_crdf12[*] file ${output_filename}_rdf12.txt mode vector
#compute crdf22 all rdf ${step_energy} 2 2
#fix outcrdf22 all ave/time ${step_energy} ${num_energy} ${time_energy} c_crdf22[*] file ${output_filename}_rdf22.txt mode vector

#fix outrgyrhist all ave/histo 1000 250 250000 0 $N $N c_rgyr file ${output_filename}_rgyr_hist.txt ave one

#fix spat all ave/spatial ${step_energy} ${num_energy} ${time_energy} z lower 0.1 density/number file ${output_filename}_spat.txt ave one
#compute binz all chunk/atom bin/1d z lower 0.1 units box
compute binsphere all chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
compute binspherea type1 chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
compute binsphereb type2 chunk/atom bin/sphere 0 0 0 ${R_sphere} ${hi_bound} 1000 units box
#fix spat all ave/chunk ${step_energy} ${num_energy} ${time_energy} binz density/number file ${output_filename}_spat.txt ave one
fix spatall all ave/chunk ${step_energy} ${num_energy} ${time_energy} binsphere density/number file ${output_filename}_spat_all.txt ave one

fix spat all ave/chunk ${step_energy} ${num_energy} ${time} binsphere density/number file ${output_filename}_spat.txt ave one
fix spata type1 ave/chunk ${step_energy} ${num_energy} ${time} binspherea density/number file ${output_filename}_spatA.txt ave one
fix spatb type2 ave/chunk ${step_energy} ${num_energy} ${time} binsphereb density/number file ${output_filename}_spatB.txt ave one

#compute cluster1 type1 cluster/atom ${r_aggr}
#compute cluster2 type2 cluster/atom ${r_aggr}

#compute clust_cut gcut cluster/atom ${r_aggr}

#compute coord1 all coord/atom ${r_aggr} 1 2 3 *
#compute coord1 all coord/atom cutoff ${r_aggr} 1 2 3 *

#compute clustnum all cluster/atom ${r_aggr}
compute poten all pe/atom pair
compute bonen all pe/atom bond

#dump dump_cluster all custom ${time_dump} ${output_filename}.dump.gz id type mol xu yu zu ix iy iz c_cluster1 c_cluster2 c_coord1[1] c_coord1[2] c_coord1[3] c_coord1[4] c_clustnum c_poten c_bonen c_cluster1b c_cluster1bb c_cluster1b0 c_cluster1b1 c_cluster1b2 c_cluster1b3 c_cluster1b4 c_clust_z5 c_clust_z10 c_clust_z15 c_clust_z20 c_clust_z25 c_clust_z30 c_clust_z35 c_clust_z40 c_clust_z45 c_clust_z50 c_clust_z55 c_clust_z60 c_clust_z65 c_clust_z70 c_clust_z75 c_clust_z80
#dump dump_cluster all custom ${time_dump} ${path_to_dump}${output_filename}.dump.gz id type mol xu yu zu ix iy iz c_cluster1 c_cluster2 c_coord1[1] c_coord1[2] c_coord1[3] c_coord1[4] c_clustnum c_poten c_bonen c_clust_cut c_clust_z5 c_clust_z10 c_clust_z15 c_clust_z20 c_clust_z25 c_clust_z30 c_clust_z35 c_clust_z40 c_clust_z45 c_clust_z50 c_clust_z55 c_clust_z60 c_clust_z65 c_clust_z70 c_clust_z75 c_clust_z80
#dump dump_cluster all custom ${time_dump} ${path_to_dump}${output_filename}.dump.gz id type mol xu yu zu ix iy iz c_cluster1 c_cluster2 c_coord1[1] c_coord1[2] c_coord1[3] c_coord1[4] c_clustnum c_poten c_bonen c_clust_cut
END

if($lomo_flag)
{
print SCRIPT <<END;
dump dump_cluster all custom ${time_dump} ${path_to_dump}${output_filename}.lammpstrj id type mol xu yu zu ix iy iz c_poten c_bonen
END
}
else
{
print SCRIPT <<END;
dump dump_cluster all custom ${time_dump} ${path_to_dump}${output_filename}.lammpstrj.gz id type mol xu yu zu ix iy iz c_poten c_bonen
END
}

print SCRIPT <<END;
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
compute pebond all pe/atom bond
compute pepair all pe/atom pair
compute peatall all pe/atom 

compute pebondA type1 pe/atom bond
compute pepairA type1 pe/atom pair
compute peatallA type1 pe/atom 

compute pebondB type2 pe/atom bond
compute pepairB type2 pe/atom pair
compute peatallB type2 pe/atom 

fix peb all ave/histo ${step_energy} ${num_energy} ${time_energy}    0 100 1000 c_pebond   file  ${output_filename}_hist_bond.txt ave one mode vector
fix pep all ave/histo ${step_energy} ${num_energy} ${time_energy} -100 100 2000 c_pepair   file  ${output_filename}_hist_pair.txt ave one mode vector
fix pea all ave/histo ${step_energy} ${num_energy} ${time_energy} -100 200 3000 c_peatall  file  ${output_filename}_hist_all.txt ave one mode vector

fix peba type1 ave/histo ${step_energy} ${num_energy} ${time_energy}    0 100 1000 c_pebondA  file  ${output_filename}_hist_bondA.txt ave one mode vector
fix pepa type1 ave/histo ${step_energy} ${num_energy} ${time_energy} -100 100 2000 c_pepairA  file  ${output_filename}_hist_pairA.txt ave one mode vector
fix peaa type1 ave/histo ${step_energy} ${num_energy} ${time_energy} -100 200 3000 c_peatallA file  ${output_filename}_hist_allA.txt ave one mode vector

fix pebb type2 ave/histo ${step_energy} ${num_energy} ${time_energy}    0 100 1000 c_pebondB  file  ${output_filename}_hist_bondB.txt ave one mode vector
fix pepb type2 ave/histo ${step_energy} ${num_energy} ${time_energy} -100 100 2000 c_pepairB  file  ${output_filename}_hist_pairB.txt ave one mode vector
fix peab type2 ave/histo ${step_energy} ${num_energy} ${time_energy} -100 200 3000 c_peatallB file  ${output_filename}_hist_allB.txt ave one mode vector

END

print SCRIPT <<END;

timestep 0.005

fix bal all balance 50000 1.0 shift xyz 10 1.1

#dump dump all atom ${time_dump} ${path_to_dump}${output_filename}.lammpstrj.gz
END
if($lomo_flag)
{
print SCRIPT <<END;
dump dumplast all atom ${time_lastdump} ${output_filename}.last.lammpstrj
END
}
else
{
print SCRIPT <<END;
dump dumplast all atom ${time_lastdump} ${output_filename}.last.lammpstrj.gz
END
}

if(scalar @fix_eps)
 {
 print SCRIPT <<END;
pair_coeff 1 1 yukawa $fix_eps[0]
pair_coeff 1 3 yukawa $fix_eps[0]
pair_coeff 3 3 yukawa $fix_eps[0]
pair_coeff 1 2 yukawa $fix_eps[1]
pair_coeff 2 3 yukawa $fix_eps[1]
pair_coeff 2 2 yukawa $fix_eps[2]

run ${time}
write_data ${output_filename}.data
END
 }
 else
 {

print SCRIPT <<END;
#dump_modify dump scale no

run $time
write_data ${output_filename}_e0.data

END

for($ip=1;$ip<scalar(@epsAA);$ip++)
 {
#my $time_out = $time_mult[$ip]*$time;
if($structure eq "PH3")
{
if($epsBB[$ip]<=-5)
 {
$time_out = $time;
 }
 else
 {
$time_out = $time/10;
 }
}
else
{
$time_out = $time;
}

 print SCRIPT <<END;
pair_coeff 1 1 yukawa $epsAA[$ip]
pair_coeff 1 3 yukawa $epsAA[$ip]
pair_coeff 3 3 yukawa $epsAA[$ip]
pair_coeff 1 2 yukawa $epsAB[$ip]
pair_coeff 2 3 yukawa $epsAB[$ip]
pair_coeff 2 2 yukawa $epsBB[$ip]

run ${time_out}
write_data ${output_filename}_e$ip.data
END
 }

 }
#for($ip=scalar(@epsAA)-2;$ip>=0;$ip--)
# {
# print SCRIPT <<END;
#pair_coeff 1 1 yukawa $epsAA[$ip]
#pair_coeff 1 3 yukawa $epsAA[$ip]
#pair_coeff 3 3 yukawa $epsAA[$ip]
#pair_coeff 1 2 yukawa $epsAB[$ip]
#pair_coeff 2 3 yukawa $epsAB[$ip]
#pair_coeff 2 2 yukawa $epsBB[$ip]
#
#run $time
#END
# }

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
mpirun -np ${np} /opt/lammps/lmp_ubuntu_24Dec20_kokkos_cuda_cmake -in $scriptname
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
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
/opt/lammps/lmp_ubuntu_24Dec20_kokkos_cuda_cmake  -in $scriptname
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
