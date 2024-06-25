# RomanETC

This is the exposure time calculator used for Roman (previously WFIRST) 
computations. The current configuration files are included in the 
config_* directories.

## Compilation instructions

You should be able to compile with any standard C compiler, e.g.,:

gcc exptimecalc.c -lm -Wall -O3 -o my_executable.exe options

Options are described in the manual but may include the "big two" 
(-DBAO_MODE, -DWL_MODE) as well as the many options to configure the ETC 
differently.

## References

- Code: C. Hirata et al., "The WFIRST Galaxy Survey Exposure Time 
Calculator", arXiv:1204.5151

- Weak lensing source catalogs: data/2K.dat is based on the CANDELS 
GOODS-S data set. Y. Guo et al., ApJS 207:24 (2013) and L.-T. Hsu et 
al., ApJ 796:60 (2014)

- Luminosity function models for galaxy redshift survey mode are as 
cited in the Manual.

- Configuration files are based on the June 2024 update from the Roman 
Project Office.
