/* change this */
cd ~/GitHub/inventory_dynamics/empirics/agg_empirics/codes

/* Housekeeping */
clear all
set matsize 5000
set more off

/* Move and import */
local varlist  A812RC2Q027SBEA GDPC96 /* ISratio and real GDP*/
local varlist `varlist' A813RC2Q027SBEA /* inventory to goods struc. sale */
local varlist `varlist' B008RA3Q086SBEA /* Real private investment */
local varlist `varlist' A014RY2Q224SBEA /* Contribution to GDP growth all CIPI*/
local varlist `varlist' A191RL1Q225SBEA /* Real GDP growth */
local varlist `varlist' A015RY2Q224SBEA /* growth annualized nonfarm CIPI*/
local varlist `varlist' CBIC96 /* Real CIPI */
local varlist `varlist' A810RX1Q020SBEA /* Real final sale of goods/structure*/
local varlist `varlist' A810RC1Q027SBEA /* nominal sale goods struc*/
local varlist `varlist' A809RX1Q020SBEA /* Real final sales of domestic business*/
local varlist `varlist' A809RC1Q027SBEA /* nominal sales domestic business */
local varlist `varlist' A373RX1Q020SBEA /* Real nonfarm inventory*/
local varlist `varlist' A015RX1Q020SBEA /* real change in nonfarm inventories*/
local varlist `varlist' DHUTRX1Q020SBEA /* housing and utility services GDP */
local varlist `varlist' PRFIC1 /* residential fixed investment */
local varlist `varlist' A2009X1Q020SBEA /* real gross housing value added */
local varlist `varlist' A2009C1Q027SBEA /* nominal gross housing VA */
local varlist `varlist' GDP
local varlist `varlist' GDPDEF /* GDP deflator */
local varlist `varlist' GCEC96 /* real Government consumption */
local varlist `varlist' NETEXC96 /* real net export */
local varlist `varlist' FPIC96 /* real fixed private investment*/
local varlist `varlist' PCECC96 /* real consumption */
local varlist `varlist' GPDIC96 /* real gross private domestic investment */
local varlist `varlist' EXPGSC96 /* real export */
local varlist `varlist' IMPGSC96 /* real import */
local varlist `varlist' FINSLC96 /* real final sale of GDP */
local varlist `varlist' CQRMTSPL /* real sales of manufacturing and trade */
local varlist `varlist' A373RD3Q086SBEA /*deflator nonfarm inventory */
local varlist `varlist' LREM64TTUSQ156S /*employment rate 15-64*/
local varlist `varlist' GDPDEF /* GDP deflator */
local varlist `varlist' GDPPOT /* GDP deflator */
 
/* download and rename */
freduse `varlist'
rename A813RC2Q027SBEA ISratio
rename GDPC96 rGDP
rename FPIC96 rInvestment
rename GCEC96 rG
rename PCECC96 rC
rename GPDIC96 rGrossInv
rename EXPGSC96 rExport
rename IMPGSC96 rImport
rename CBIC96 rCIPI
rename A191RL1Q225SBEA GDPgrowth
rename A015RY2Q224SBEA CIPIgrowth
rename A014RY2Q224SBEA CIPIcontri
rename A809RX1Q020SBEA rSales
rename A810RX1Q020SBEA rSalesGoods
rename A810RC1Q027SBEA nSalesGoods
rename A809RC1Q027SBEA nSales
rename A373RX1Q020SBEA rInventory
rename A015RX1Q020SBEA rCIPI_nonfarm
rename DHUTRX1Q020SBEA rhousingservice
rename PRFIC1 resid_inv
rename A2009X1Q020SBEA rhousingVA
rename A2009C1Q027SBEA nhousingVA
rename FINSLC96 rSaleGDP
rename CQRMTSPL rMTSales
rename A373RD3Q086SBEA price_Inventory
rename GDPDEF deflator

/* get the dates right */
drop date
gen date = qofd(daten)
ds date, not
format date %tq
tsset date, quarterly
scalar ginf = yq(1970,1)
scalar gmod = yq(1983,4)
scalar grec = yq(2007,4)
scalar kt_begin = yq(1954,1)
scalar kt_end = yq(2002,1)

/* gen new vars */
gen ln_rGDP = log(rGDP)
gen ln_ISratio = log(ISratio)
gen ln_rInventory = log(rInventory)
gen ln_rSales = log(rSalesGoods)
gen ln_rInvestment = log(rInvestment)
gen rISratio = rInventory/rSalesGoods
gen ln_rISratio = log(rISratio)
gen f = 1/(1+rISratio)
gen ln_f = log(f)
gen g_CIPI = (rCIPI-L.rCIPI)/L.rCIPI
gen rX = rExport - rImport
gen rFI = rGrossInv - rCIPI
gen share_CIPI = rCIPI/rGDP
gen price_SalesGoods = 100*nSalesGoods/rSalesGoods
gen price_Sales = 100*nSales/rSales

/* reorder things */
order date, first

/* save rawly processed dataset */
save ../data/quarterly_FRED, replace 
export delimited using "../data/quarterly_FRED.txt", delimiter("|") datafmt replace
