DEFINE/PAR P1 ? ? "catalog.dat"
DEFINE/PAR P2 0 N "alpha0 in deg"
DEFINE/PAR P3 0 N "delta0 in deg"


CREATE/TAB cat 8 500 {P1} fmtad2xy.fmt

COMPUTE/TAB cat :x_rad = -COS(:DEC) * SIN(:R_A -{P2})/(COS({P3})*COS(:DEC)*COS(:R_A - {P2})+SIN({P3})*SIN(:DEC))
COMPUTE/TAB cat :y_rad = ( COS({P3}) * SIN(:DEC) - SIN({P3})*COS(:DEC)*COS(:R_A - {P2}))
COMPUTE/TAB cat :y_rad = :y_rad /(COS({P3})*COS(:DEC)*COS(:R_A -{P2})+SIN({P3})*SIN(:DEC))
