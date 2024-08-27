for $N (3,5) {

$Np = $N+1;
$file = "summary-wl-$N.txt";

@bands = ('F213', 'F184', 'F158', 'F129', 'F106', 'F087', 'F062', 'F146');
@RGB = (
'ff0000', 'd08000', '80a000', '00d020', '0060a0', '0020ff', '6000d0', '404040'
);

# Generate WL yields versus time
for $lat (0,15,30,45,60,75) {
  open(G, q-|tee x.dat|gnuplot-) or die;
  print G qq|set term postscript enhanced 18 eps color\n|;
  print G qq|set output "plots/yields-$N-$lat.eps"\n|;
  print G qq|set xlabel "survey time [days per filter 1000 deg^2]"\n|;
  print G qq|set ylabel "WL sources n_{eff} [galaxies per arcmin^2]"\n|;
  print G qq|set xrange [0:70]; set yrange [0:75]\n|;
  print G qq|set grid; set key at 18,72\n|;
  print G qq+set title "$N exposures, |{/Symbol b}|=$lat deg"\n+;
  for $ib (0..7) {
    $ibp = $ib+1;
    $lw = 1.5;
    if ($ib==7) {$lw=2.5;}
    print G qq|set style line $ibp lt 1 lw $lw dt $ibp lc rgb "#$RGB[$ib]"\n|;
  }
  @ord = (7,2,3,4,1,5,6,0);
  print G 'plot ';
  for $i (0..7) {
    $ib = $ord[$i];
    $ii = $ib+1;
    $c = int(substr $bands[$ib], 1, 3);
    print G qq|"$file" using (\$1==$lat?\$8:1/0):(\$10==$c?\$4:1/0) with lines title "$bands[$ib]" ls $ii|;
    if ($i<7) {print G q:, :;}
  }
  close G;
}

# Generate WL yields versus exposure time
for $lat (0,15,30,45,60,75) {
  open(G, q-|tee x.dat|gnuplot-) or die;
  print G qq|set term postscript enhanced 18 eps color\n|;
  print G qq|set output "plots/galrate-$N-$lat.eps"\n|;
  print G qq|set xlabel "exposure time [s]"\n|;
  print G qq|set ylabel "effective shapes [10^6 galaxies/day]"\n|;
  print G qq|set xrange [0:300]; set yrange [0:9.2]\n|;
  print G qq|set grid; set key at 295,9.05\n|;
  print G qq+set title "$N exposures, |{/Symbol b}|=$lat deg"\n+;
  for $ib (0..7) {
    $ibp = $ib+1;
    $lw = 1.5;
    if ($ib==7) {$lw=2.5;}
    print G qq|set style line $ibp lt 1 lw $lw dt $ibp lc rgb "#$RGB[$ib]"\n|;
  }
  @ord = (7,2,3,4,1,5,6,0);
  print G 'plot ';
  for $i (0..7) {
    $ib = $ord[$i];
    $ii = $ib+1;
    $c = int(substr $bands[$ib], 1, 3);
    print G qq|"$file" using (\$1==$lat?\$2:1/0):(\$10==$c?\$4/\$8*3.6:1/0) with lines title "$bands[$ib]" ls $ii|;
    if ($i<7) {print G q:, :;}
  }
  close G;
}

# Generate depths versus time
for $lat (0,15,30,45,60,75) {
  open(G, q-|tee x.dat|gnuplot-) or die;
  print G qq|set term postscript enhanced 18 eps color\n|;
  print G qq|set output "plots/depth-$N-$lat.eps"\n|;
  print G qq|set xlabel "survey time [days per filter 1000 deg^2]"\n|;
  print G qq|set ylabel "depth (5{/Symbol s} pt src AB)"\n|;
  print G qq|set xrange [0:70]; set yrange [28.5:23]\n|;
  print G qq|set grid; set key at 67.5,23.1\n|;
  print G qq+set title "$N exposures, |{/Symbol b}|=$lat deg"\n+;
  for $ib (0..7) {
    $ibp = $ib+1;
    $lw = 1.5;
    if ($ib==7) {$lw=2.5;}
    print G qq|set style line $ibp lt 1 lw $lw dt $ibp lc rgb "#$RGB[$ib]"\n|;
  }
  @ord = (0..3,5,4,6,7);
  print G 'plot ';
  for $i (0..7) {
    $ib = $ord[$i];
    $ii = $ib+1;
    $c = int(substr $bands[$ib], 1, 3);
    print G qq|"$file" using (\$1==$lat?\$8:1/0):(\$10==$c?\$6:1/0) with lines title "$bands[$ib]" ls $ii|;
    if ($i<7) {print G q:, :;}
  }
  close G;
}

# redshift plot
open(G, q-|tee x.dat|gnuplot-) or die;
print G qq|set term postscript enhanced 18 eps color\n|;
print G qq|set output "plots/zmed-$N.eps"\n|;
print G qq|set xlabel "WL sources n_{eff} [galaxies per arcmin^2]"\n|;
print G qq|set ylabel "source median redshift"\n|;
print G qq|set xrange [0:75]; set yrange [0.2:1.3]\n|;
print G qq|set grid; set key at 67.5,.8\n|;
print G qq+set title "$N exposures, |{/Symbol b}|=$lat deg"\n+;
for $ib (0..7) {
  $ibp = $ib+1;
  $lw = 1.5;
  if ($ib==7) {$lw=2.5;}
  print G qq|set style line $ibp lt 1 lw $lw dt $ibp pt $ibp ps .7 lc rgb "#$RGB[$ib]"\n|;
}
@ord = (0..3,5,4,6,7);
print G 'plot ';
for $i (0..7) {
  $ib = $ord[$i];
  $ii = $ib+1;
  $c = int(substr $bands[$ib], 1, 3);
  print G qq|"$file" using (\$4):(\$10==$c?\$5:1/0) with points title "$bands[$ib]" ls $ii|;
  if ($i<7) {print G q:, :;}
}
close G;

system 'rm x.dat';

}
