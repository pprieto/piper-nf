#!/usr/bin/env perl
use strict;
use warnings;

# overlap.pl bestHit query_gtf
my $best_result = shift;
my $query_gtf = shift;

#TAKE REAL QUERY COORDINATES
open (QGTF, "<$query_gtf") or die "Cannot open query gtf $query_gtf";
my @queryName_gtfLine_allExons = <QGTF>;
if (! @queryName_gtfLine_allExons){
  die "Error[overlap.pl]! the $query_gtf gtf file does not contain query gtf lines\n";
}

open (BH, "<$best_result") or die "Cannot open blast result";
my $blastOut = <BH>;
chomp $blastOut;

#COMPARE
my ($b_chr , $b_q_start, $b_q_end, $b_s_start , $b_s_end , $b_q_frame , $b_s_frame , $b_strand);
if ($blastOut=~/^\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/){
  $b_chr     = $1;
  $b_q_frame = $2;
  $b_s_frame = $5;
  $b_s_start = $6;
  $b_s_end   = $7;

  $b_q_start = $3;
  $b_q_end   = $4;

  if($b_s_frame =~ /^[+-]\d+$/ ) {
    $b_strand = '+' if (($b_q_frame=~/\+/) && ($b_s_frame=~/\+/));
    $b_strand = '-' if (($b_q_frame=~/\+/) && ($b_s_frame=~/-/));
    $b_strand = '-' if (($b_q_frame=~/-/) && ($b_s_frame=~/\+/));
    $b_strand = '+' if (($b_q_frame=~/-/) && ($b_s_frame=~/-/));
  } else {
    $b_strand = '+' if (($b_q_start < $b_q_end) && ($b_s_frame eq 'plus'));
    $b_strand = '-' if (($b_q_start < $b_q_end) && ($b_s_frame eq 'minus'));
    $b_strand = '-' if (($b_q_start > $b_q_end) && ($b_s_frame eq 'plus'));
    $b_strand = '+' if (($b_q_start > $b_q_end) && ($b_s_frame eq 'minus'));
  }
}
else {
  die "Error[overlap.pl]! impossible to parse $blastOut\n";
}
foreach my $queryName_gtfLine (@queryName_gtfLine_allExons){
  my ($g_chr , $g_start , $g_end , $g_strand);
  if ($queryName_gtfLine=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/){
    $g_chr    = $1;
    $g_start  = $2;
    $g_end    = $3;
    $g_strand = $4;
  }
  else {
    die "Error[overlap.pl]! $queryName_gtfLine is a wrong gtf format\n";
  }
  if (($g_chr eq $b_chr) && ($g_strand eq $b_strand) ){
    verifyOverlap($b_s_start , $b_s_end , $g_start , $g_end );
  }
}
print "0";
exit;

sub verifyOverlap {
  my ($b_start , $b_end , $g_start , $g_end ) = @_;
  if ((($b_start <= $g_end) && ($b_start >= $g_start)) || (($b_end >= $g_start) && ($b_end <= $g_end)) || (($b_start <= $g_start) && ($b_end >= $g_end)) ||   (($b_start >= $g_start) && ($b_end <= $g_end))){
    print "1";
    exit;
  }
}
