/* change this */
cd ~/GitHub/inventory_dynamics/empirics/agg_empirics/codes

/* Housekeeping */
clear all
set matsize 5000
set more off

/* load rawly processed dataset */
use ../data/quarterly_FRED, clear

/* gen some variables and stdze them*/
gen ln_rSalesGoods = log(rSalesGoods)
gen ln_rC = log(rC)
gen ln_rGrossInv = log(rGrossInv)
egen zSales = std(D.ln_rSalesGoods)
egen zStock = std(D.ln_rInventory)
egen zGDP = std(D.ln_rGDP)
egen zC = std(D.ln_rC)
egen zInv = std(D.ln_rGrossInv)

/* simple factor model on quarterly series */
dfactor (zStock zSales zGDP zC zInv=, noconstant) (f  =, ar(1/2)), difficult
predict factor, factorz smethod(smooth)
tsline factor

/* simple factor model on quarterly series number 2 */
dfactor (D.(ln_rInventory ln_rSalesGoods ln_rGDP ln_rC ln_rGrossInv) =) (f  =, constant ar(1)), difficult
predict factor, factor smethod(smooth)
tsline factor
estat ssdisplay /* secret command to display state space form */
