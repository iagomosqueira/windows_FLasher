# test.R - DESC
# /home/mosqu003/Sandbox/ARS_18_19_20/test.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

library(mse)

# -

ctrl <- fwdControl(year=2025:2029, quant="fbar", value=20)

fut <- fwd(om, ctrl)

unitMeans(fbar(fut))[, ac(2024:2029)]

unitSums(ssb(fut))[, ac(2024:2029)]

# - 

ctrl.om <- readRDS('ctrl.om.rds')

out <- do.call('mpDispatch', ctrl.om)

fut <- fwd(stock(ctrl.om$om), sr=sr(ctrl.om$om), deviances=ctrl.om$deviances,
  control=ctrl.om$ctrl)
