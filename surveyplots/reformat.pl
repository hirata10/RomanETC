open(IN, 'summary-bao.txt') or die;
while ($line = <IN>) {
  if ($line !~ m/^\#/) {
    if (length($line)>2) {
      @data = split ' ', $line;
      $N = scalar @data;
      for $i (0..$N-2) {print "$data[$i],";}
      print $data[$N-1];
      print "\n";
    }
  }
}
close IN;
