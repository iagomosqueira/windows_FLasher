# test.R - DESC
# /home/mosqu003/Sandbox/DG/test.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}



 f0 <- fwd(om, control=fwdControl(quant='fbar', value=0, year=seq(iy + 1, fy)))
f0_ss3 <- fwd(om_ss3, control=fwdControl(quant='fbar', value=0, year=seq(iy + 1, fy)))

units(srpars_ss3)
units(srpars)
units(srpars_ss3) <- '1000'
f0_ss3 <- fwd(om_ss3, control=fwdControl(quant='fbar', value=0, year=seq(iy + 1, fy)))
