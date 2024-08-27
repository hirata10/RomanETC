system "gcc ../exptimecalc.c -lm -Wall -O3 -o baoetc-o3.exe -DBAO_MODE -DIN_DLON -DOIII_GAL -DFOWLER8";
system "gcc ../exptimecalc.c -lm -Wall -O3 -o baoetc-o2.exe -DBAO_MODE -DIN_DLON -DOII_GAL -DFOWLER8";
system "gcc ../exptimecalc.c -lm -Wall -O3 -o baoetc.exe -DBAO_MODE -DIN_DLON -DUSE_NII -DFOWLER8";

@exec = ('baoetc.exe', 'baoetc-o3.exe', 'baoetc-o2.exe');
@mod = (62, 1992, 2002);

$t1 = 4.03*10; $cost1 = 40.3/1.7;
$t2 = 4.03*100; $cost2 = 213.6/1.7;

open(SUMMARY, '>summary-bao.txt') or die;
for $lat (0,15,30,45,60,75) {
  print OUT (sprintf "# ecliptic latitude %2d degrees\n", $lat);
  for ($t=66.000001; $t<600; $t+=8.25) {

    for $L (0..2) {
      # loop over lines
      open(K, qq+| ./$exec[$L] > temp-$L.out+) or die;
      print K "1\n../config_jun2024/T24_Grism.dat\n.07\n2\n.05\n1.00\n1.93\n";
      print K ".6\n100\n";
      print K "$t\n5\n.015\n$lat\n120\n.025\n";
      print K "6\n";
      print K "5\n";
      print K "$mod[$L]\n";
      close K;

      # pull out the data
      $d = `cat temp-$L.out`;
      ($d1,$d2) = split "\n\n", $d;
      $e1 = (split /\|\n/m, $d1)[-1];
      @data = split "\n", $e1;
      $Nz[$L] = scalar @data;
      for $i (0..$Nz[$L]-1) {
        @data2 = split ' ', $data[$i];
        $ncol = scalar @data2;
        for $j (0..$ncol-1) {$val[$L][$i][$j] = $data2[$j];}
      }

      # get other info
      $density[$L] = (split ' ', `sed -n -e/Available/p temp-$L.out`)[-2];
    }

    $cost = $cost1 + ($cost2-$cost1)/($t2-$t1)*($t-$t1);
    print SUMMARY (sprintf "%2d %6.2f %8.3f   %8.2f %8.2f %8.2f  %9.5f\n", $lat, $t, $cost,
      $density[0], $density[1], $density[2], $val[0][15][4]/1e-19);

    $outfile = sprintf "datafiles/grs-$lat-6x%.2f.dat", $t;
    open(OUT, ">$outfile") or die;
    for $L (0..2) {
      for $i (0..$Nz[$L]-1) {
        for $j (0..$ncol-1) {print OUT (sprintf " %13.6E", $val[$L][$i][$j]);}
        print OUT "\n";
      }
      print OUT "\n";
    }
    close OUT;
  }
  print SUMMARY "\n";
}
close SUMMARY;
