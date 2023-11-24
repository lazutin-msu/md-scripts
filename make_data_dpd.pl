#!/usr/bin/perl
use Getopt::Long;
use File::Copy;
use Math::Trig;
use POSIX;
    
$time = 10000000;
$repeat = 1;

$Ngpu = 1;
$Ntypes = 4;

$path_to_dump = "";

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
#'sequence|s=s' => \$sequence,
'NA=i' => \$NA,
'NB=i' => \$NB,
'NB2=i' => \$NB2,
'phi=f' => \$phi,
'repeat|r=i' => \$repeat,
'time|t=i' => \$time,
'g|gpuid=i' => \$gpu_num,
'lomo' => \$lomo_flag,
'c|cell=f' => \$cell_size,
#'e|energy=s' => \$energyfile,
'chi' => \$chi,
'np=i' => \$np,
'Ncpu=i' => \$Ncpu,
'relax=i' => \$trelax,
'noanneal' => \$flag_noanneal,
'kbond=f' => \$kbond,
'rbond=f' => \$rbond,
'stiff=f' => \$stiff,
'help|h' => \$help);
if($help)
 {
 die("Usage make_data.pl options
  -f datafile initial configuration
  -s sequence O   C   B
  --NA length of A block
  --NB length of B block
  --NB2 length of B homopolymer
  --phi volume ratio of homopolymer
  -r repeat = 1
  -t time = 2000000
  -g gpuid
  --lomo to calc on Lomonosov-2
  -c --cell cell size = 40
  --chi chi_AB
  --np number of cores to run on
  --Ncpu number of tasks which run in parallel
  --noanneal no annealing
  --relax (100000) relaxation time
  --kbond (40) bond koefficient
  --rbond (1) bond length
  -h help\n");
 }

if(!$Ncpu){$Ncpu = 6;}

if(!$cell_size){$cell_size = 40;}

if(!$kbond){$kbond = 40.0;}
if(!$rbond){$rbond = 1.0;}

#if($sequence !~ /^[OCB]+$/) {die "sequence $sequence should contain only O C B \n";}

#$countO = $sequence =~ tr/O//;
#$countC = $sequence =~ tr/C//;
#$countB = $sequence =~ tr/B//;
#$countO = length( $str =~ s/[^\Q$char\E]//rg )
#$countC = $str =~ tr/C//;
#$countB = $str =~ tr/B//;

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

sub description
{
my $NA = shift
my $NB = shift
my $NB2 = shift
my $time = shift;
my $cell_size = shift;

my $desc = sprintf "#dpd_A%dB%d_B%d_c%d_t%d",$NA,$NB,$NB2,$cell_size,$time;
return $desc;
}

sub datafilename1
{
my $NA = shift
my $NB = shift
#my $NB2 = shift
my $time = shift;
my $N = shift;
my $cell_size = shift;

my $desc = sprintf "A%dB%d_t%d_c%d_n%d.ml",$NA,$NB,$time,$cell_size,$N;
return $desc;
}

sub datafilename2
{
#my $NA = shift
#my $NB = shift
my $NB2 = shift
my $time = shift;
my $N = shift;
my $cell_size = shift;

my $desc = sprintf "B%d_t%d_c%d_n%d.ml",$NB2,$time,$cell_size,$N;
return $desc;
}

sub dirname
{
my $NA = shift
my $NB = shift
my $NB2 = shift
my $time = shift;
my $N = shift;
my $cell_size = shift;

my $desc = sprintf "A%dB%d_B%d_t%d_c%d_n%d",$NA,$NB,$NB2,$time,$cell_size,$N;
return $desc;
}

sub myround
{ # for positive x
my $x = shift;
return int($x+0.5);
}

#if($gpuid==1){$gpuid=2;}
printf "R_sphere = %f\nchains = %d\nlength = %d \ngpuid= %d\n",$R_sphere,$chains,$N,$gpuid;

$Ndir = 1;

$desc = description($NA, $NB, $NB2, $time, $cell_size);

if($datafile)
{
$is_datafile = 1;
}
else
{
$is_datafile = 0;
$datafile1 = datafilename1($NA, $NB, $time, $Ndir, $cell_size);
$datafile2 = datafilename2($NB2, $time, $Ndir, $cell_size);
}

$dirname = dirname($NA, $NB, $NB2, $time, $Ndir, $cell_size);
$scriptname = "script";
$shellname_all = "run.sh";
#$shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
$shellname = sprintf "run_cpu%d.sh",$cpuid;  #look at run.sh too
$output_filename = dirname($NA, $NB, $NB2, $time, $Ndir, $cell_size);

if( -d $dirname)
 {
# die "$dirname exists";
 while(-d $dirname)
  {
  $Ndir++;
  $dirname = dirname($NA, $NB, $NB2, $time, $Ndir, $cell_size);
#  $shellname = sprintf "run_c%d_d%f_N%d_n%d.sh",$chains,$dens,$N,$Ndir;
  $output_filename = dirname($NA, $NB, $NB2, $time, $Ndir, $cell_size);
  }
 }
 

`mkdir $dirname`;

$dirname = $dirname."/";

$bond = 1.0;

#$natoms = 2 * $chains * $N + 1;

#$totalatoms = $natoms * $spheres;
#$totalbonds = $nbonds * $spheres;

$lx = $cell_size;
$ly = $cell_size;
$lz = $cell_size;

print "$desc \n";
$iatom = 1;
$ibond = 1;
$iangle = 1;

#$N = length($sequence);

  for($imono=0;$imono<$NA;$imono++)
   {
   $x[$iatom] = $imono * $bond;
   $y[$iatom] = 0;
   $z[$iatom] = 0;

   $type[$iatom] = 1; # A
   
   if($imono!=0)
    {
    $at1_bond[$ibond] = $iatom - 1;
    $at2_bond[$ibond] = $iatom;
    $bond_type[$ibond] = 1;
    $ibond++;
    }
   $iatom++;
   }

  for($imono=$NA;$imono<$NA+$NB;$imono++)
   {
   $x[$iatom] = $imono * $bond;
   $y[$iatom] = 0;
   $z[$iatom] = 0;

   $type[$iatom] = 2; # B
   
#   if($imono!=0)
#    {
    $at1_bond[$ibond] = $iatom - 1;
    $at2_bond[$ibond] = $iatom;
    $bond_type[$ibond] = 1;
    $ibond++;
#    }
   $iatom++;
   }



# beads were made

#(($ibond-1)==$totalbonds) or die "ibond $ibond != totalbonds $totalbonds";
$totalatoms = $iatom-1;
$totalbonds = $ibond-1;

$totalatoms1 = $totalatoms;

if(!$is_datafile)
{
print "writing to $datafile \n";
open(OUT,">".$dirname.$datafile1) or die "cant create file ".$dirname.$datafile1;

print OUT $desc."\n\n";
printf OUT "%11d atoms\n", $totalatoms;
printf OUT "%11d bonds\n", $totalbonds;

print OUT "\nCoords\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %8.2f %8.2f %8.2f\n",$iatom,$x[$iatom],$y[$iatom],$z[$iatom];
 }

print OUT "\nTypes\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %6d\n",$iatom,$type[$iatom];
 }

print OUT "\nBonds\n\n";

for($ibond=1;$ibond<=$totalbonds;$ibond++)
 {
 printf OUT "%6d %6d %6d %6d\n",$ibond,$bond_type[$ibond],$at1_bond[$ibond],$at2_bond[$ibond];
 }

close(OUT);

$iatom = 1;
$ibond = 1;
$iangle = 1;

#$N = length($sequence);

  for($imono=0;$imono<$NB2;$imono++)
   {
   $x[$iatom] = $imono * $bond;
   $y[$iatom] = 0;
   $z[$iatom] = 0;

   $type[$iatom] = 2; # B
   
   if($imono!=0)
    {
    $at1_bond[$ibond] = $iatom - 1;
    $at2_bond[$ibond] = $iatom;
    $bond_type[$ibond] = 1;
    $ibond++;
    }
   $iatom++;
   }

# beads were made

#(($ibond-1)==$totalbonds) or die "ibond $ibond != totalbonds $totalbonds";
$totalatoms = $iatom-1;
$totalbonds = $ibond-1;
$totalangles = $iangle-1;

$totalatoms2 = $totalatoms;

if(!$is_datafile)
{
print "writing to $datafile \n";
open(OUT,">".$dirname.$datafile2) or die "cant create file ".$dirname.$datafile2;

print OUT $desc."\n\n";
printf OUT "%11d atoms\n", $totalatoms;
printf OUT "%11d bonds\n", $totalbonds;
if($totalangles>0)
{
printf OUT "%11d angles\n", $totalangles;
}

print OUT "\nCoords\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %8.2f %8.2f %8.2f\n",$iatom,$x[$iatom],$y[$iatom],$z[$iatom];
 }

print OUT "\nTypes\n\n";

for($iatom=1;$iatom<=$totalatoms;$iatom++)
 {
 printf OUT "%6d %6d\n",$iatom,$type[$iatom];
 }

print OUT "\nBonds\n\n";

for($ibond=1;$ibond<=$totalbonds;$ibond++)
 {
 printf OUT "%6d %6d %6d %6d\n",$ibond,$bond_type[$ibond],$at1_bond[$ibond],$at2_bond[$ibond];
 }

if($totalangles>0)
{
print OUT "\nAngles\n\n";

for($iangle=1;$iangle<=$totalangles;$iangle++)
 {
 printf OUT "%6d %6d %6d %6d %6d\n",$iangle,$ang_type[$iangle],$at1_angle[$iangle],$at2_angle[$iangle],$at3_angle[$iangle];
 }
}

close(OUT);

}
else
{
copy($datafile,$dirname.$datafile) or die "Copy failed: $!";
}

# write script

#$gpuid = int(rand(3));
$seed = int(rand(9999999));
$seed2 = int(rand(9999999));
$seed_create1 = int(rand(9999999));
$seed_mol1 = int(rand(9999999));
$seed_create2 = int(rand(9999999));
$seed_mol2 = int(rand(9999999));

#$num_mols = myround(3*$lx*$ly*$lz/$totalatoms);

$num_mols1 = myround(3*$lx*$ly*$lz*(1.0-$phi)/$totalatoms1);
$num_mols2 = myround(3*$lx*$ly*$lz*$phi/$totalatoms2);


$xlo_bound = 0;
$xhi_bound = $lx;
$ylo_bound = 0;
$yhi_bound = $ly;
$zlo_bound = 0;
$zhi_bound = $lz;

$time_energy = $time /4;
$step_energy = 100;
$num_energy = $time_energy / $step_energy;

$set_temp = 1.0;

open(SCRIPT,">".$dirname.$scriptname);

if($lomo_flag)
{
print SCRIPT <<END;
#package gpu 1 neigh yes gpuID $gpuid $gpuid
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
boundary        p p p
atom_style full

region          bxx block ${xlo_bound} ${xhi_bound} ${ylo_bound} ${yhi_bound} ${zlo_bound} ${zhi_bound}

END

print SCRIPT <<END;
create_box      6 bxx bond/types 1 extra/bond/per/atom 10 extra/special/per/atom 20
END

print SCRIPT <<END;
molecule 1 ${datafile1}
molecule 2 ${datafile2}

create_atoms 0 random ${num_mols1} ${seed_create1} bxx mol 1 ${seed_mol1}
create_atoms 0 random ${num_mols2} ${seed_create2} bxx mol 2 ${seed_mol2}

bond_style harmonic
bond_coeff 1 ${kbond} ${rbond} 

END

print SCRIPT <<END;

# type(A) = 1 
# type(B) = 2 

pair_style      dpd 1.0 1.0 ${seed2}

pair_coeff      * * 25  4.5 1

mass * 1.0

special_bonds lj   1.0 1.0 1.0

neighbor        1.0 bin
neigh_modify    delay 0 every 1 check yes cluster no 
comm_modify vel yes cutoff 3.5

write_data      ${output_filename}.data
minimize 0.000000001 0.000000001 10000 10000

velocity all create 1.0 ${seed}
fix             1 all nve

thermo_style    custom step temp pe epair etotal 
thermo_modify flush yes
thermo    1000

write_data      ${output_filename}.min.data

reset_timestep 0

timestep        0.03

restart  500000  restart.${output_filename}

dump            1 all atom 100000 ${output_filename}.lammpstrj

END

if($trelax)
{
print SCRIPT <<END;
run ${trelax}
END
}
else
{
print SCRIPT <<END;
run 100000
#run 10000
END
}

print SCRIPT <<END;
write_data ${output_filename}_mix.data

pair_coeff      1 1 25  4.5 1
pair_coeff      2 2 25  4.5 1
pair_coeff      3 3 25  4.5 1
pair_coeff      4 4 25  4.5 1
pair_coeff      5 5 25  4.5 1

END

$N = 10;
#$time1 = $time/$N;
$time1 = $time/$N/2;
$time2 = $time/2;

if($flag_noanneal)
{

  printf SCRIPT "pair_coeff      1  2 %.2f %.1f %d\n",$chi*3.27+25,4.5,1;
printf SCRIPT "run %d\n",$time;

}
else
{

for($k=1;$k<=$N;$k++)
{

  printf SCRIPT "pair_coeff      1  2 %.2f %.1f %d\n",$chi/$N*$k*3.27+25,4.5,1;

printf SCRIPT "run %d\n",$time1;
}

print SCRIPT <<END;

#run ${time} 
run ${time2} 
END
}
print SCRIPT <<END;

write_data ${output_filename}_eq.data

unfix 1


print "All done"


END


close(SCRIPT);

open(SH,">>".$shellname);

if($lomo_flag) # lomo
{
print SH <<END;
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
#mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
sbatch -n 256 --time=1-23:30:00 -p compute_prio ompi /home/lazutin_2123/_scratch/src/lammps-29Oct20/build-nogpu/lmp  -in ${scriptname}
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
/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
#mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp -sf gpu  -in ${scriptname}
cd ..
END
  }
 }
 else # no np
 {
print SH <<END;
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-6.5/targets/x86_64-linux/lib
cd $dirname
#/home/lazutin/src/lammps-10Feb15/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-11Aug17/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-16Mar18/src/lmp_cuda -sf gpu  -in $scriptname
#/home/lazutin/src/lammps-29Oct20/build/lmp  -in $scriptname
mpirun -np 4 /home/lazutin/src/lammps-29Oct20/build/lmp -sf gpu  -in ${scriptname}
cd ..
END
 }
}

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
