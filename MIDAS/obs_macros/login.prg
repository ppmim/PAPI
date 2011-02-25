!		LOGIN_quick.PRG        for OMEGA2000 quicklook


WRITE/KEY ERR_CTRL/I/1/3 0,2,1
write/key mid$sys/c/21/10 "$debu     "
write/key user/c/1/20 "OMEGA2000"
write/key INSTR_ID/C/1/8 OMEGA2k
write/key DETECTOR/C/1/8 OMEGA2k



set/midas prompt=QUICK
write/key mid_session/i/1/1 31

return
