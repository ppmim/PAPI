#Some examples for doing some common tests about PAPI
dfits *.fits |fitsort OBJECT
./calDark.py -s /disk-a/caha/panic/DATA/SIMU_PANIC_3/darks50.txt -o /tmp/mdark.fits
./calDomeFlat.py -s /disk-a/caha/panic/DATA/SIMU_PANIC_3/dflats.txt -o /tmp/dflats.fits
./calBPM_2.py -f /disk-a/caha/panic/DATA/SIMU_PANIC_3/dflats.txt -d /tmp/mdark_50_10.fits  
./calBPM.py -f /disk-a/caha/panic/DATA/SIMU_PANIC_3/dflats_on.txt
swarp q1.fits q2.fits q3.fits q4.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE q1.weight.fits,q2.weight.fits,q3.weight.fits,q4.weight.fits -COMBINE_TYPE CHI2 
swarp double.fits -COPY_KEYWORDS RA,DEC,OBJECT,PIXSCALE,ROT-RTA,INSTRUME,TELESCOPE

