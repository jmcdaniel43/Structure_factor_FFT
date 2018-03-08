#!/usr/bin/perl
use POSIX ;

# this script sums up all contributions to S(k) at |kvec|=k over S(kvec)
$ifile=$ARGV[0];
open ifile, $ifile, or die ;
$maxk=300000;
my@Sk;my@countk;

for(my$i=0;$i<$maxk;$i++)
{ $Sk[$i] = 0  ; $countk[$k]=0 }
$dk=0.001;
$endk=0;
while(<ifile>)
{
    @array=split;
    my@kvec;
    $kvec[0]=$array[0];
    $kvec[1]=$array[1];
    $kvec[2]=$array[2];
    $Skmag = $array[3];

    $kmag = sqrt( $kvec[0]**2 + $kvec[1]**2 + $kvec[2]**2  );
    $index = floor ( $kmag / $dk );

    if ( $index > $endk ) { $endk = $index }

    $Sk[$index] = $Sk[$index] + $Skmag ; 
    $countk[$index] = $countk[$index] + 1 ; 
}

# print
for( my$k=0;$k<$endk;$k++)
{  $kmag = $k * $dk ;
if ( $countk[$k] > 0 ){
 $avgSk = $Sk[$k] / $countk[$k];
}
else { $avgSk = 0 }
print " $kmag   $avgSk \n"
}
