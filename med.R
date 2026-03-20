# med.R - DESC
# /home/mosqu003/Sandbox/mse/med.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

install.packages('mse', repos="https://flr.r-universe.dev")

library(mse)

load('fwd.om_args.rda')

tes <- fwd(args$object, control=args$control, deviances=args$deviances[,,1]) 

unitSums(catch(tes)[,'2029'])

args$control

do.call(fwd.om, args)

mpDispatch(floval, args)
