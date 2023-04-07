#!/usr/bin/perl
use Getopt::Long;

GetOptions (
'f|file=s' => \$in_file, 
'o|out=s' => \$outfile, 
'h|help' => \$help
);

if($help) {
  print <<END;
Options for this program:
-f --file for input file default= dump.all
-o --output for output file 
END
  exit;
  }

open IN, "<$in_file" or die "Cannot open $in_file: $!";
open OUT, ">$outfile" or die "Cannot open $outfile: $!" ;

while($str=<IN>)
 {
 if($str =~ /reset_timestep/) 
   {
   #print "run found\n"; 
   last;
   }
 }


while($str=<IN>)
 {
# if($str =~ /run 1000000/) 
 if($str =~ /run /) 
   {
   #print "run found\n"; 
   last;
   }
 }

while($str=<IN>)
 {
 if($str =~ /Step/) 
   {
   print OUT $str; 
   #print "title found\n"; 
   last;}
 }

while($str=<IN>)
 {
 if($str =~ /Loop/) 
   {
   last;
   }
   else
   {
   print OUT $str;
   }
 }

while(!eof(IN))
{
while($str=<IN>)
 {
 if($str =~ /Step/) { last;}
 }

while($str=<IN>)
 {
 if($str =~ /Loop/) 
   {
   last;
   }
   else
   {
   print OUT $str;
   }
 }
}


close(IN);
close(OUT);
