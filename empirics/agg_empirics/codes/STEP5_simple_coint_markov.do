/* change this */
cd ~/GitHub/inventory_dynamics/empirics/agg_empirics/codes

/* Housekeeping */
clear all
set matsize 5000
set more off

/* load rawly processed dataset */
use ../data/quarterly_FRED, clear

/* gen some variables */
gen ln_rSalesGoods = log(rSalesGoods)

/* firstly let's check for the cointegration in I and S*/
/* doesn't seem to be cointegrated at all */
kpss ln_rInventory, maxlag(4)
kpss ln_rISratio if date > yq(1984,1), maxlag(4)
kpss ln_rISratio if date < yq(1984,1), maxlag(4)
kpss ln_rISratio, maxlag(4)

vec ln_rInventory ln_rSalesGoods, trend(constant) lags(4)
vecrank rInventory rSalesGoods, trend(constant) lags(4) levela

/* let's check out CIPI and GDP */
vecrank rGDP rCIPI, trend(constant) lags(4) levela
vec rGDP rCIPI, trend(constant) lags(4)
kpss rCIPI, maxlag(4)

/* now let's focus on inventory stock first */
mswitch ar D.ln_rInventory, states(3) ar(1) difficult level(95)
est sto m1
predict pr_state*, pr smethod(smooth)
nbercycles pr_inv_state1, file("crap1.do") replace
nbercycles pr_inv_state2, file("crap2.do") replace
nbercycles pr_inv_state3, file("crap3.do") replace

mswitch ar D.ln_rInventory, states(2) ar(1) vce(robust) difficult
est sto m2
predict pr_state*, pr smethod(smooth)
nbercycles pr_state1, file("crap1.do") replace
nbercycles pr_state2, file("crap2.do") replace
drop pr_state*

mswitch dr D.ln_rInventory, states(3) varswitch vce(robust) difficult
est sto m3
predict pr_state*, pr smethod(smooth)
nbercycles pr_state1, file("crap1.do") replace
nbercycles pr_state2, file("crap2.do") replace
nbercycles pr_state3, file("crap3.do") replace
drop pr_state*

mswitch dr D.ln_rInventory, states(3) vce(robust) difficult
est sto m4
predict pr_state*, pr smethod(smooth)
nbercycles pr_state1, file("crap1.do") replace
nbercycles pr_state2, file("crap2.do") replace
nbercycles pr_state3, file("crap3.do") replace

drop pr_state*

/* now take a look at sales, do the regimes align? */
mswitch ar D.ln_rSales, states(2) ar(1) difficult level(95)
predict pr_state*, pr smethod(smooth)
nbercycles pr_state1, file("crap1.do") replace
nbercycles pr_state2, file("crap2.do") replace
drop pr_state*

/* the joint outcome of inventory and sales */
mswitch dr D.ln_rISratio, states(2) varswitch difficult level(95)
predict pr_state*, pr smethod(smooth)
nbercycles pr_state1, file("crap1.do") replace
nbercycles pr_state2, file("crap2.do") replace
nbercycles pr_state3, file("crap3.do") replace
drop pr_state*

/* peak at share_rCIPI */
mswitch dr share_CIPI L.share_CIPI L2.share_CIPI, states(2)
capture {
  predict pr_share_state*, pr smethod(smooth)
}
if _rc != 0 {
  drop pr_state*
  predict pr_share_state*, pr smethod(smooth)
}
nbercycles pr_share_state1, file("crap1.do") replace
nbercycles pr_share_state2, file("crap2.do") replace

