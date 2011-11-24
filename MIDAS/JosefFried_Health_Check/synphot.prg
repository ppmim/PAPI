DEFI/PAR P1 gsatlas C "input table =?"
DEFI/PAR P2 test C "output table=?"
DEFI/PAR P3 10,5 N "#objects, # filters =?"
DEFI/PAR P4 photdata C "filtertable=?"
WRITE/KEYW IN_A/C/1/8 'P1'
WRITE/KEYW OUT_A/C/1/8 'P2'
WRITE/KEYW INPUTI/I/1/2 'P3'
WRITE/KEYW IN_B/C/1/8 'P4'
!
!
write/keyw filcol/i/1/10 45,46,47,48,49,50,51,52,53,54,55,56
write/keyw obcol/i/100/1 -1
write/keyw obcol/i/1/1 -1
!write/keyw inputi/i/3/4 'p4'
!WRITE/KEYW INPUTC/C/1/3 'P5'
RUN FMP:synphot
!
!correct for color indices of vega
! alte Korrektur
!COMP/TAB 'P2' :U_B = :JOHNSN_U-:JOHNSN_B-0.86
!COMP/TAB 'P2' :B_V = :JOHNSN_B-:JOHNSN_V+0.119
!COMP/TAB 'P2' :V_R = :JOHNSN_V-:JOHNSN_R+0.04
!COMP/TAB 'P2' :R_I = :JOHNSN_R-:JOHNSN_I+0.03
!
COMP/TAB 'P2' :U_B = :JOHNSN_U-:JOHNSN_B-1.16
COMP/TAB 'P2' :B_V = :JOHNSN_B-:JOHNSN_V+0.09
COMP/TAB 'P2' :V_R = :JOHNSN_V-:JOHNSN_R+0.213
COMP/TAB 'P2' :R_I = :JOHNSN_R-:JOHNSN_I+0.175
name/col 'p2' #7 "mag" f6.2
name/col 'p2' #8 "mag" f6.2
name/col 'p2' #9 "mag" f6.2
name/col 'p2' #10 "mag" f6.2
