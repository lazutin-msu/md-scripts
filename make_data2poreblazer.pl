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
$Ncpu = 12;
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
'rsphere|a=f' => \$R_sphere,
'grid|z=f' => \$grid_size,
#'dump|d=s' => \$dumpfile,
'repeat|r=i' => \$repeat,
'g|gpuid=i' => \$gpu_num,
#'m|max_cluster=i' => \$max_cluster,
#'p|spheres=i' => \$spheres,
#'x|dx=f' => \$spheres_dx,
#'pressure=f' => \$pressure,
#'time_deform=i' => \$time_deform,
#'rate_deform=f' => \$rate_deform,
#'dir_deform=s' => \$dir_deform,
#'eps|e=f{3}' => \@epsilon,
'ncpu=i' => \$Ncpu,
#'np=i' => \$num_np,
'prefix=s' => \$prefix,
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
  -h help\n");
 }

if(!$time_deform) {$time_deform = 100000;}
#if(!$period_deform) {$period_deform = 5000;}
if(($dir_deform ne "x")&&($dir_deform ne "y")&&($dir_deform ne "z")) {$dir_deform = "x";}
#if(($num_np!=1)&&($num_np!=2)&&($num_np!=4)&&($num_np!=8)) {die "--np $num_np : should be 1, 2, 4, or 8\n";}
if(!$rate_deform){$rate_deform = 0.00002;}

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

$dirname = sprintf "%scluster_n%d",$prefix,$Ndir;
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = sprintf "%scluster_n%d",$prefix,$Ndir;

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
#  $dirname = sprintf "%sp%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_rate_%.3e_n%d",$prefix,$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$rate_deform,$Ndir;
  $dirname = sprintf "%scluster_n%d",$prefix,$Ndir;

#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
#  $output_filename = sprintf "%sp%d_c%d_R%.1f_N%d_eps_%.2f_%.2f_%.2f_deform_r%s_time%d_rate_%.3e_n%d",$prefix,$spheres,$chains,$R_sphere,$N,$epsilon[0],$epsilon[1],$epsilon[2],$dir_deform,$time_deform,$rate_deform,$Ndir;
  $output_filename = sprintf "%scluster_n%d",$prefix,$Ndir;
  }
 }

`mkdir $dirname`;

$dirname = $dirname."/";

#copy($datafile,$dirname.$datafile) or die "Copy failed: $!";


#########################


open DATA, "<$datafile" or die "Cannot open $datafile: $!";

$Natoms = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\d+\s+atoms\s+$/)
  {
  ($Natoms) = ($str =~ /^\s*(\d+)\s+atoms\s+$/);
  last;
  }
 }
if($Natoms==0){die "cannot find number of atoms in $datafile";}
$Nbonds = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*\d+\s+bonds\s+$/)
  {
  ($Nbonds) = ($str =~ /^\s*(\d+)\s+bonds\s+$/);
  last;
  }
 }
if($Nbonds==0){die "cannot find number of bonds in $datafile";}

$Xhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+xlo\s+xhi\s+$/)
  {
  ($Xlo,$Xhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+xlo\s+xhi\s+$/);
  last;
  }
 }
if($Xhi==0){die "cannot find xhi in $datafile";}
$Yhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+ylo\s+yhi\s+$/)
  {
  ($Ylo,$Yhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+ylo\s+yhi\s+$/);
  last;
  }
 }
if($Yhi==0){die "cannot find yhi in $datafile";}
$Zhi = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*-?\d+\.?\d*e?[+-]?\d*\s+\d+\.?\d*e?[+-]?\d*\s+zlo\s+zhi\s+$/)
  {
  ($Zlo,$Zhi) = ($str =~ /^\s*(-?\d+\.?\d*e?[+-]?\d*)\s+(\d+\.?\d*e?[+-]?\d*)\s+zlo\s+zhi\s+$/);
  last;
  }
 }
if($Zhi==0){die "cannot find zhi in $datafile";}

$atoms_flag = 0;
while($str=<DATA>)
 {
 if($str =~ /^\s*Atoms\s+/)
  {
  $atoms_flag = 1;
  last;
  }
 }
if($atoms_flag==0){die "cannot find atoms section in $datafile";}

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
 $x[$arr[0]-1] = $arr[3];
 $y[$arr[0]-1] = $arr[4];
 $z[$arr[0]-1] = $arr[5];
 $ix[$arr[0]-1] = $arr[6];
 $iy[$arr[0]-1] = $arr[7];
 $iz[$arr[0]-1] = $arr[8];
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
if($bonds_flag==0){die "cannot find bonds section in $datafile";}

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

open(XYZ,">".$dirname.$output_filename.".xyz");

print XYZ $Natoms."\n";
print XYZ "XYZ\n";
for($i=0;$i<$Natoms;$i++)
{
my $str_type;
if($atom_type[$i]==1)
 {
 $str_type="A";
 }
 elsif($atom_type[$i]==2)
 {
 $str_type="B";
 }
 elsif($atom_type[$i]==3)
 {
 $str_type="C";
 }
 elsif($atom_type[$i]==4)
 {
 $str_type="D";
 }
 else
 {
 die "atom type out of range - atom $i has type $atom_type[$i]";
 }
my $myx = $x[$i] - $Xlo;
my $myy = $y[$i] - $Ylo;
my $myz = $z[$i] - $Zlo;
printf XYZ "%-9s%16.9f%16.9f%16.9f\n", $str_type, $myx, $myy, $myz;
}

close(XYZ);

open(DAT,">".$dirname.$output_filename.".dat");

print DAT $output_filename.".xyz\n";

my $lx = $Xhi - $Xlo;
my $ly = $Yhi - $Ylo;
my $lz = $Zhi - $Zlo;

printf DAT "%16.9f%16.9f%16.9f\n", $lx, $ly, $lz;
printf DAT "%16d%16d%16d\n", 90, 90, 90;

close(DAT);

$seed = int(rand(9999999));

open(DEF,">".$dirname."defaults.dat");
print DEF <<END ;
my.atoms
1.00, 1.00, 100.0, 20.0
1.00
500
$grid_size
100.0, 3.00
$seed
0
END

close(DEF);

open(FF,">".$dirname."my.atoms");

my $diam = 2.0 * $R_sphere;

print FF <<END ;
4
A       1.000    1.0     1.0
B       1.000    1.0     1.0
C       1.000    1.0     1.0
D       ${diam}    1.0     1.0
END
close(FF);

#########################

#$gpuid = int(rand(3));
$seed2 = int(rand(9999999));


open(SH,">>".$shellname);
#print SH "nohup /usr/mount-opt/lammps/lmp_ubuntu_1Feb14 <$scriptname  &";
print SH <<END;
cd $dirname
/home/lazutin/scripts/PoreBlazer/src/poreblazer.exe < ${output_filename}.dat > ${output_filename}.txt
cd ..
END
close(SH);

open(SH,">".$shellname_all);
for($icpu=0;$icpu<$Ncpu;$icpu++)
{
$shellname1 = sprintf "run_cpu%d.sh",$icpu;  #look at init block too
print SH <<END;
nohup sh $shellname1 &
END
}
close(SH);

}
