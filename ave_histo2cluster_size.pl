#!/usr/bin/perl
use Getopt::Long;
#print "step rgyr1 rgyr2\n";

#$N=2000;
#$N=500;
#$Ncolumn = 3;
$x_col = "Coord";
$y_col = "Count";
$Npart = 27;

GetOptions ('infile|f=s' => \$in_file,
'outfile|o=s' => \$out_file,
'xcol|x=s' => \$x_col,
'ycol|y=s' => \$y_col,
'every|e=i' => \$every,
'nmol|n=i' => \$nmol,
'npart=i' => \$Npart,
'npoly|p=i' => \$Npoly,
'help|h' => \$help);
if($help)
 {
 die("Usage split.pl options
  -f in_file 
  -o out_file = in_file.dat
  -x xcol = Coord
  -y ycol = density/number
  -e every = read only every x-th step
  -n nmol = number of monomers in particle
  --npart = 
  -p npoly = polymerization degree
  -h help\n");
 }

if(!defined($in_file))
 {
 die "in_file should be defined\n-h for help\n";
 }
if(!defined($out_file))
 {
 $out_file = $in_file . ".dat";
 }

open(IN,"<".$in_file) || die "cant open $in_file \n";

$str=<IN>;
$str=<IN>;
$str =~ /^# TimeStep Number-of-bins/ or die "wrong format: second line is not '# Timestep Number-of-bins'";
#$str =~ /^# Timestep Number-of-chunks/ or die "wrong format: second line is not '# Timestep Number-of-chunks'";
#$str =~ /^# TimeStep Number-of-rows/ or die "wrong format: second line is not '# TimeStep Number-of-rows' '${str}'  ";
$str=<IN>;
@arr = split /\s+/, $str;
$Ncol = scalar(@arr)-1;
for($i=0;$i<$Ncol;$i++)
 {
 $col_name[$i] = $arr[$i+1];
# print "$i ".$col_name[$i]."\n";
 }

%rev_col_names = {};
for($i=0;$i<scalar(@col_name);$i++)
 {
 $rev_col_names{$col_name[$i]}=$i;
 }

defined($rev_col_names{$x_col}) or die "can not find x column $x_col \n";
defined($rev_col_names{$y_col}) or die "can not find y column $y_col \n";
$x_num = $rev_col_names{$x_col} - 1;
$y_num = $rev_col_names{$y_col} - 1;



$n=0;
while(!eof(IN))
 {
 $str=<IN>;
 $str =~ s/^\s+//;
 @arr = split /\s+/, $str;
 $step[$n] = $arr[0];
 $num = $arr[1];

# print "step".$step[$n]." num".$num."\n";

if(!defined($every)||($step[$n]%$every==0))
{
 if(!defined($nmol))
 {
  for($i=0;$i<$num;$i++)
  {
  $str=<IN>;
  $str =~ s/^\s+//;
  @arr = split /\s+/, $str;
  if(($i+1)!=$arr[0]){die "step $step[$n] line ".($i+1)." contains line number $arr[0] \n";}
#  for($j=0;$j<scalar(@arr)-1;$j++)
  for($j=0;$j<$Ncol;$j++)
   {
   $data[$n][$i][$j] = $arr[$j+1];
#   print "step $n x $i col $j ".$data[$n][$i][$j]."\n";
   }  #for j
  } #for i
  $n++;
 } #nmol
 else
 {
  for($i=0;$i<$num;$i++)
  {
  $str=<IN>;
  $str =~ s/^\s+//;
  @arr = split /\s+/, $str;
  if(($i+1)!=$arr[0]){die "step $step[$n] line ".($i+1)." contains line number $arr[0] \n";}

   if($i%$nmol==0)
   {
   $imod = $i/$nmol;
   for($j=0;$j<$Ncol;$j++)
    {
#    $data[$n][$i][$j] = $arr[$j+1];
    $data[$n][$imod][$j] = $arr[$j+1];
    }  #for j
   }
   else
   {
   if($i==1)
    {
     if($arr[$y_num+1]!=$Npart)
     {
     die "step $step[$n] i = 1 is $arr[$y_num+1] not $Npart \n";
     }
    }
    else
    {
     if($arr[$y_num+1]!=0)
     {
     die "step $step[$n] i = $i is $arr[$y_num+1] not 0 \n";
     }
    }
   }# nmol
  } #for i
  $n++;
 
 }
}
else
{
 for($i=0;$i<$num;$i++)
  {
  $str=<IN>;
  }
} # $every
 } #while eof

close(IN);

#if(defined($nmol))
# {
# $num = $num / $nmol;
# }

#print "$x_num $y_num\n";

open(OUT,">".$out_file);

#print OUT $x_col;
for($a=0;$a<$n;$a++)
 {
 my $num = $step[$a]/1000;
# print OUT " ".$y_col.$num."k";
 }
#print OUT "\n";

if(defined($nmol))
{
for($i=0;$i<$num/$nmol;$i++)
 {
# print OUT $data[0][$i][$x_num]/$nmol;

# print OUT $i;
 for($a=0;$a<$n;$a++)
   {
#   print OUT " ".$data[$a][$i][$y_num];
   }
# print OUT "\n";
 }
}
else
{
if(!defined($Npoly))
 {
  for($i=0;$i<$num;$i++)
  {
#  print OUT $data[0][$i][$x_num];
  for($a=0;$a<$n;$a++)
   {
#   print OUT " ".$data[$a][$i][$y_num];
   }
#  print OUT "\n";
  }
 }
 else
 {
  for($i=0;$i<$num;$i++)
  {
  my $myx = int($data[0][$i][$x_num]);
  
  
  if(($myx!=0)&&($myx%$Npoly==0))
   {
#   print OUT  $myx/$Npoly;
   for($a=0;$a<$n;$a++)
    {
    $sum_N2[$a]    += $data[$a][$i][$y_num] * $myx/$Npoly * $myx/$Npoly;
    $sum_N[$a]    += $data[$a][$i][$y_num] * $myx/$Npoly;
    $sum_norm[$a] += $data[$a][$i][$y_num] ;
#    print OUT " ".$data[$a][$i][$y_num];
    }
#   print OUT "\n";
   }
   else
   {
   if($data[$a][$i][$y_num]!=0)
    {
    die "non zero y in timestep $step[$a] clustersize $data[0][$i][$x_num] y $data[$a][$i][$y_num] \n";
    }   
   }
  } # for i
 } # defined Npoly
} # no nmol

print OUT "Step Nnum Nweight K\n";
for($a=0;$a<$n;$a++)
 {
printf OUT "%d %g %g %g\n", $step[$a], $sum_N[$a]/$sum_norm[$a], $sum_N2[$a]/$sum_N[$a], ($sum_N2[$a]/$sum_N[$a])/($sum_N[$a]/$sum_norm[$a]);
 }

close(OUT);

