# AWK: files for imcombine
BEGIN {
a = 0
}
{
	printf ("%s\n", $1) > "im.list"
	printf ("%f %f\n", -$2, -$3) > "off.list"
}
