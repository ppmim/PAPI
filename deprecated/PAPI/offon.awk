# AWK: select even rows
BEGIN {
a = 0
}
{
	a = a+1
	if (int(a/2)!=(a/2)) {
		printf ("%s\n", $1)
	}
}
