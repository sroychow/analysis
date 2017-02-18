#!/usr/bin/env perl

use strict;
use warnings;

#{run}:{lumi}:{event}:{mass4l:.2f}:{mZ1:.2f}:{mZ2:.2f}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
#{njets30:d}:{jet1pt:.2f}:{jet2pt:.2f}:{category} 
sub readFile
{
  my ($file, $h) = @_;
  print ">>> Reading $file\n";
  open INPUT, $file or die qq|Failed to open $file, died|;
  while (<INPUT>) {
    chop;
    my ($run, $lumi, $event, $mass4l, $massz1, $massz2, $dbkg_kin, $dbkg, $dgg, $dhjj, $d0, $njets, $jet1pt, $jet2pt, $cat) = (split /:/);
      #= (split /:/)[0,1,2,3,4,5,-4,-1]; 
#    my $key = join(":", $run, $lumi, $event, $mass4l, $massz1, $massz2, $dbkg_kin, $dbkg, $dgg, $dhjj, $d0, $njets, $jet1pt, $jet2pt, $cat);
    my $key = join(":", $run, $lumi, $event,$mass4l, $massz1, $massz2);
    $h->{$key} = $_;
  }
  close INPUT;
}
sub _main {
  my ($infile1, $infile2) = @_;
  my $dict1 = {};
  readFile($infile1, $dict1);

  my $dict2 = {};
  readFile($infile2, $dict2);

  print join("", ">>> Lines present in <", $infile1, "> not in <", $infile2, ">"), "\n";
  for my $key (sort keys %$dict1) {
    print $dict1->{$key}, "\n" unless exists $dict2->{$key};
  }

 # print join("", ">>> Lines present in <", $infile2, "> not in <", $infile1, ">"), "\n";
 # for my $key (sort keys %$dict2) {
 #   print $dict2->{$key}, "\n" unless exists $dict1->{$key};
 # }
}
_main($ARGV[0], $ARGV[1]);
__END__
