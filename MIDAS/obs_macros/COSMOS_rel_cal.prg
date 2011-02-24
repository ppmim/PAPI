! Relative calibration of COSMOS fields (see manual Sec. 13.9):


$cmd_o2000 filter H_OLD

! NE corner of 4c

$cmd_o2000 tele abs  10  0 58.91  2 19 55.06 2000
$cmd_o2000 sync

o2k/relcal COSMOS_NE = 2,15

! NW corner of 4c

$cmd_o2000 tele abs  9 59 58.31  2 19 55.21 2000
$cmd_o2000 sync

o2k/relcal COSMOS_NW = 2,15

! SE corner of 4c

$cmd_o2000 tele abs  10  0 58.89  2  4 46.76 2000
$cmd_o2000 sync

o2k/relcal COSMOS_SE = 2,15

! SW corner of 4c

$cmd_o2000 tele abs 9 59 58.3   2  4 46.89  2000
$cmd_o2000 sync

o2k/relcal COSMOS_SW = 2,15

return
