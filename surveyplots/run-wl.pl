# number of exposures
($Ne) = @ARGV;

# Setup
system "rm wletc.exe";
# system "rm datafiles/WL_*";
system "gcc ../exptimecalc.c -lm -Wall -o wletc.exe -DWL_MODE -DWL_CUT_DEFAULT -DOUT_WL_CAT -DIN_DLON -DOUT_EXTRA_PSF_PROPERTIES";

# bands
@bands = ('F213', 'F184', 'F158', 'F129', 'F106', 'F087', 'F062', 'F146');
# these are the nominal wavelength limits
@lmin = (1.95, 1.68, 1.38, 1.13, 0.93, 0.76, 0.48, 0.93);
@lmax = (2.30, 2.00, 1.77, 1.45, 1.19, 0.98, 0.76, 2.00);

# survey times (days per 1000 deg**2 per filter)
if ($Ne==5) {
  $t1 =  33*3.08; $treq1 = 192.4/1.7/4;
  $t2 = 100*3.08; $treq2 = 484.8/1.7/4;
}

if ($Ne==3) {
  $t1 =  33*3.08; $treq1 = 252.8/1.7/8;
  $t2 = 100*3.08; $treq2 = 642.7/1.7/8;
}

open(OUT, qq:>summary-wl-$Ne.txt:) or die;

print OUT "# Columns:\n";
print OUT "#   ecliptic latitude (deg)\n";
print OUT "#   time per exposure (sec)\n";
print OUT "#   WL n (arcmin^-2)\n";
print OUT "#   WL n_eff (arcmin^-2)\n";
print OUT "#   median z for WL sources\n";
print OUT "#   5 sigma pt src mag\n";
print OUT "#   5 sigma ext src mag (r_eff=0.3 arcsec, exp profile)\n";
print OUT "#   survey time (days per 1000 deg^2)\n";
print OUT "#   number of exposures\n";
print OUT "#   filter\n";
print OUT "\n";

$ind = 0;
for $band (0..7) {
  print OUT "# ### $Ne exposures, band $bands[$band] ###\n\n";

  for $lat (0,15,30,45,60,75) {
    print OUT (sprintf "# ecliptic latitude %2d degrees\n", $lat);
    for ($t=33.000001; $t<303; $t+=8.25) {

      $file = sprintf "datafiles/WL_$bands[$band]_lat%02d_%1dx%06.2f.out", $lat, $Ne, $t;

      # Write configuration file
      open(IN, q->IN-) or die;
      print IN "1\n";
      print IN "../config_jun2024/T24_$bands[$band].dat\n";
      print IN ".07\n"; # RMS wavefront error
      print IN "2\n"; # detector type
      print IN ".014\n"; # jitter
      $la = 0.95*$lmin[$band];
      $lb = 1.05*$lmax[$band];
      print IN "$la\n$lb\n1\n";
      print IN "$t\n";
      print IN "5.\n"; # read noise floor
      print IN ".015\n"; # dark current
      print IN "$lat\n120.0\n"; # longitude relative to Sun = 120 deg
      print IN ".03\n"; # E(B-V)
      print IN "$Ne\n"; # number of exposures
      print IN ".4\n.2\n../data/2K.dat\n$file.temp\n"; # WL inputs
      close IN;

      # Run ETC
      print "($ind) > $file "; system q-date-;
      if ($ind>=0) {
        (system "./wletc.exe < IN > $file") and die;
      }
      $ind++;

      # Get information
      $str = `tail -n 3 $file | head -n 1`;
      $str =~ m/\ ([\d\.E\+\-]+)\ /;
      $n = $1;
      $str = `tail -n 1 $file`;
      $str =~ m/\ ([\d\.E\+\-]+)\ /;
      $neff = $1;

      $str = `sed -n -e/5\\\ sigma\\\ pt\\\ src\\\ /p $file`;
      $str =~ m/\ ([\d\.]+)\ mag\ AB/;
      $maglim = $1;
      $str = `sed -n -e/5\\\ sigma\\\ ext\\\ src\\\ /p $file`;
      $str =~ m/\ ([\d\.]+)\ mag\ AB/;
      $maglimExt = $1;

      # redshift information
      @z = ();
      $N=0;
      open(CAT, "sort -nk2 $file.temp |") or die;
      while ($line=<CAT>) {
        $z[$N] = (split ' ', $line)[1];
        $N++;
      }
      close CAT;
      $zmed = $z[int($N/2)];      

      # survey time requirement (days per 1000 deg^2 per filter)
      $treq = $treq1 + ($treq2-$treq1)/($t2-$t1)*($t-$t1);

      $c = substr $bands[$band], 1, 3;
      $line = sprintf "%2d %6.2f   %5.2f %5.2f %5.3f   %6.3f %6.3f   %6.2f   %1d $c\n",
        $lat, $t, $n, $neff, $zmed, $maglim, $maglimExt, $treq, $Ne;
      print $line;
      print OUT $line;
    }
    print OUT "\n";
  }

}

close OUT;

# Cleanup
