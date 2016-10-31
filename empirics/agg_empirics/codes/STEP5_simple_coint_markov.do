/* change this */
cd ~/GitHub/inventory_dynamics/empirics/agg_empirics/codes

/* Housekeeping */
clear all
set matsize 5000
set more off

/* save rawly processed dataset */
use ../data/quarterly_FRED, clear

/* gen some variables */
gen ln_rSalesGoods = log(rSalesGoods)

/* firstly let's check for the cointegration */
vec ln_rInventory ln_rSalesGoods, trend(constant) lags(4)
vecrank ln_rInventory ln_rSalesGoods, trend(constant) lags(4)
