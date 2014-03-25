./papi_v1 /tmp/Q1/; cp double.fits split8/q1.fits ; cp double.weight.fits split8/q1.weight.fits ; cp double.head split8/q1.head
./papi_v1 /tmp/Q2/; cp double.fits split8/q2.fits ; cp double.weight.fits split8/q2.weight.fits ; cp double.head split8/q2.head
./papi_v1 /tmp/Q3/; cp double.fits split8/q3.fits ; cp double.weight.fits split8/q3.weight.fits ; cp double.head split8/q3.head
./papi_v1 /tmp/Q4/; cp double.fits split8/q4.fits ; cp double.weight.fits split8/q4.weight.fits ; cp double.head split8/q4.head
cd split8
swarp q1.fits q2.fits q3.fits q4.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE q1.weight.fits,q2.weight.fits,q3.weight.fits,q4.weight.fits

