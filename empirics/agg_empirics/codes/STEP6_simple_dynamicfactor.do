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

/* simple factor model on quarterly series */
dfactor (D.ln_rInventory D.ln_rSales D.ln_rGDP D.ln_rInvestment =) (f =, ar(1))
predict factor, factor
tsline factor
