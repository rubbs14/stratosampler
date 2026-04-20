"""
EGFR (CHEMBL203) IC50 / pIC50 Dataset from ChEMBL
====================================================
Target  : Epidermal growth factor receptor (EGFR / HER1 / ERBB1)
ChEMBL  : CHEMBL203
Organism: Homo sapiens
Fetched : 2026-04-20

Total records available in ChEMBL for this target + IC50 + nM filter: 18,650
This file contains 300 curated records (6 pages × 50) covering four distinct
chemical series (scaffolds).

Chemical series present
-----------------------
1. Pyrrolo-benzylidene-quinazolinone covalent inhibitors
   (the "EKB-569 / CI-1033" series; e.g. CHEMBL68920, CHEMBL69629, CHEMBL304271)
   Core SMARTS: pyrrole–CH=C attached to a 4-anilino-quinazolin-2(1H)-one

2. Triazene / hydrazone 4-anilinoquinazolines
   (e.g. CHEMBL136491, CHEMBL134312, CHEMBL133023)
   Core SMARTS: 4-anilinoquinazoline bearing -N=N-N(C)- or -N/N=N/-

3. 2,2'-Dithio-bis-indole-3-carboxamide dimers (irreversible)
   (e.g. CHEMBL292323, CHEMBL293913, CHEMBL304414)
   Core SMARTS: bis-indole linked via S–S bond

4. 4-Anilino-pyrrolo[2,3-d]pyrimidine / 7-deazapurine covalent inhibitors
   (the "Neratinib / EKB" series; e.g. CHEMBL285063, CHEMBL51659, CHEMBL334697)
   Core SMARTS: 4-anilino-7H-pyrrolo[2,3-d]pyrimidine with acrylamide warhead

5. Pyrido[3,4-d]pyrimidine / CK-636 scaffold series
   (e.g. CHEMBL52829, CHEMBL52591, CHEMBL298931)
   Core: 4-anilinopyrido[3,4-d]pyrimidine / 4-aminoquinazoline variants

6. Flavonoid / chalcone natural product-like series
   (e.g. CHEMBL44, CHEMBL428690)

Reproduction code
-----------------
Run the function `fetch_egfr_data()` to download fresh data from ChEMBL REST API.

pIC50 formula
-------------
pIC50 = -log10(IC50_nM * 1e-9)  =  9 - log10(IC50_nM)
  e.g. IC50 = 1 nM  → pIC50 = 9.00
       IC50 = 100 nM → pIC50 = 7.00
       IC50 = 10000 nM → pIC50 = 5.00
"""

import math

# ---------------------------------------------------------------------------
# Raw data  (molecule_chembl_id, smiles, ic50_nm, pchembl_value_from_chembl)
# pIC50 column is recomputed below for full transparency.
# ---------------------------------------------------------------------------

RAW = [
    # --- Batch 0-49 (offset=0) ---
    ("CHEMBL68920",  "Cc1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",                             41.0,   7.39),
    ("CHEMBL68920",  "Cc1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",                             300.0,  6.52),
    ("CHEMBL68920",  "Cc1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",                             7820.0, 5.11),
    ("CHEMBL69960",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",                  170.0,  6.77),
    ("CHEMBL69960",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",                  40.0,   7.40),
    ("CHEMBL69960",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",                  440.0,  6.36),
    ("CHEMBL137635", "CN(c1ccccc1)c1ncnc2ccc(N/N=N/Cc3ccccn3)cc12",                                           9300.0, 5.03),
    ("CHEMBL77085",  "N#CC(C#N)=Cc1cc(O)ccc1[N+](=O)[O-]",                                                    96000.0,4.02),
    ("CHEMBL443268", "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)N(C)c2ncnc(Nc3ccc(F)c(Cl)c3)c21",            5310.0, 5.28),
    ("CHEMBL76589",  "N#CC(C#N)=C(N)/C(C#N)=C/c1ccc(O)cc1",                                                   125.0,  6.90),
    ("CHEMBL76904",  "N#CC(C#N)=Cc1ccc(O)c(O)c1",                                                             35000.0,4.46),
    ("CHEMBL304271", "CCN(CC)CC(O)CNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",          0.45,   9.35),
    ("CHEMBL296407", "N#C/C(=C\\c1ccc(O)c(O)c1)C(N)=O",                                                       10000.0,5.00),
    ("CHEMBL309625", "N#CC(C#N)=C(N)/C(C#N)=C/c1ccc(O)c([N+](=O)[O-])c1",                                    60000.0,4.22),
    ("CHEMBL310798", "N#CC(C#N)=Cc1cc(O)c(O)c(O)c1",                                                          3000.0, 5.52),
    ("CHEMBL135592", "C/N=N/Nc1ccc2ncnc(N(C)c3ccccc3)c2c1",                                                   8000.0, 5.10),
    ("CHEMBL136491", "CN(C)/N=N/c1ccc2ncnc(Nc3cccc(Cl)c3)c2c1",                                               14.0,   7.85),
    ("CHEMBL336113", "Cc1cccc(Nc2ncnc3ccc(N)cc23)c1",                                                         1000.0, 6.00),
    ("CHEMBL133024", "C/N=N/Nc1ccc2ncnc(Nc3cccc(Cl)c3)c2c1",                                                  100.0,  7.00),
    ("CHEMBL344652", "CN(C)/N=N/c1ccc2ncnc(N(C)c3ccccc3)c2c1",                                                13350.0,4.88),
    ("CHEMBL47986",  "O=CN/C=C/c1cc(O)ccc1O",                                                                 14000.0,4.85),
    ("CHEMBL138125", "Brc1cccc(Nc2ncnc3ccc(N/N=N/Cc4ccccn4)cc23)c1",                                          238.0,  6.62),
    ("CHEMBL137364", "Brc1cccc(Nc2ncnc3ccc(N/N=N/CCN4CCOCC4)cc23)c1",                                         64.0,   7.19),
    ("CHEMBL302552", "O=C1Nc2ncnc(Nc3ccc(F)c(Cl)c3)c2/C1=C/c1ccc(C(=O)NCCN2CCOCC2)[nH]1",                   22.0,   7.66),
    ("CHEMBL309598", "N#CC(C#N)=C(O)c1cc(O)c(O)c(O)c1",                                                       13500.0,4.87),
    ("CHEMBL67057",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(ccn4Cc4ccccc4)c3)c21",         4.5,    8.35),
    ("CHEMBL67057",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(ccn4Cc4ccccc4)c3)c21",         100.0,  7.00),
    ("CHEMBL332096", "Nc1ccc2ncnc(NCc3ccccc3)c2c1",                                                           508.0,  6.29),
    ("CHEMBL432903", "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)c1C",             200.0,  6.70),
    ("CHEMBL69964",  "Cc1cc(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]c1C(=O)N1CCN(C)CC1",               210.0,  6.68),
    ("CHEMBL65848",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(c3)CCC4)c21",                  1410.0, 5.85),
    ("CHEMBL65848",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(c3)CCC4)c21",                  5000.0, 5.30),
    ("CHEMBL65848",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(c3)CCC4)c21",                  4320.0, 5.37),
    ("CHEMBL76985",  "N#C/C(=C\\c1ccc(C=O)cc1)C(=O)O",                                                        47000.0,4.33),
    ("CHEMBL69629",  "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               6.5,    8.19),
    ("CHEMBL69629",  "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               40.0,   7.40),
    ("CHEMBL69629",  "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               2050.0, 5.69),
    ("CHEMBL136492", "CN(c1ccccc1)c1ncnc2ccc(N)cc12",                                                         8400.0, 5.08),
    ("CHEMBL66570",  "CCN1CCN(C(=O)c2cc(C)c(/C=C3\\C(=O)Nc4ncnc(Nc5ccc(F)c(Cl)c5)c43)[nH]2)CC1",            13.0,   7.89),
    ("CHEMBL66570",  "CCN1CCN(C(=O)c2cc(C)c(/C=C3\\C(=O)Nc4ncnc(Nc5ccc(F)c(Cl)c5)c43)[nH]2)CC1",            200.0,  6.70),
    ("CHEMBL66570",  "CCN1CCN(C(=O)c2cc(C)c(/C=C3\\C(=O)Nc4ncnc(Nc5ccc(F)c(Cl)c5)c43)[nH]2)CC1",            7500.0, 5.12),
    ("CHEMBL311119", "N#CC(C#N)=Cc1ccc(C=O)cc1",                                                              60000.0,4.22),
    ("CHEMBL76905",  "N#CC(C#N)=Cc1cc(O)cc(O)c1",                                                             37000.0,4.43),
    ("CHEMBL305194", "Cc1cc(C(=O)O)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",                         29.0,   7.54),
    ("CHEMBL80155",  "COc1cc(C=C(C#N)C#N)cc(O)c1O",                                                           6000.0, 5.22),
    ("CHEMBL343352", "Nc1ccc2ncnc(Nc3cccc(Cl)c3)c2c1",                                                        200.0,  6.70),
    ("CHEMBL52765",  "Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                        44.0,   7.36),
    ("CHEMBL76642",  "N#C/C(=C/c1cc(O)ccc1O)C(=O)O",                                                          75000.0,4.12),
    ("CHEMBL308645", "Cc1cc(CCC(=O)O)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",                       4.7,    8.33),
    ("CHEMBL67003",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(cnn4Cc4ccccc4)c3)c21",         1.8,    8.74),
    # --- Batch 50-99 (offset=50) ---
    ("CHEMBL67003",  "Cc1cc(C(=O)N2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(cnn4Cc4ccccc4)c3)c21",         100.0,  7.00),
    ("CHEMBL134312", "CN(C)/N=N/c1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                               11.0,   7.96),
    ("CHEMBL77243",  "CCCCc1cc(/C=C(\\C#N)C(N)=O)cc(CCCC)c1O",                                                55000.0,4.26),
    ("CHEMBL80745",  "COc1cc(/C=C(\\C#N)C(N)=C(C#N)C#N)cc(O)c1O",                                            1200.0, 5.92),
    ("CHEMBL305246", "Cc1cc(C(=O)N2CCN(C)CC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               2.8,    8.55),
    ("CHEMBL305246", "Cc1cc(C(=O)N2CCN(C)CC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               100.0,  7.00),
    ("CHEMBL305246", "Cc1cc(C(=O)N2CCN(C)CC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",               750.0,  6.12),
    ("CHEMBL69966",  "C#Cc1cccc(Nc2ncnc3c2/C(=C/c2[nH]c(C(=O)NCCN4CCOCC4)cc2C)C(=O)N3)c1",                  18.0,   7.75),
    ("CHEMBL69966",  "C#Cc1cccc(Nc2ncnc3c2/C(=C/c2[nH]c(C(=O)NCCN4CCOCC4)cc2C)C(=O)N3)c1",                  40.0,   7.40),
    ("CHEMBL69966",  "C#Cc1cccc(Nc2ncnc3c2/C(=C/c2[nH]c(C(=O)NCCN4CCOCC4)cc2C)C(=O)N3)c1",                  620.0,  6.21),
    ("CHEMBL47940",  "Nc1ncnc2c1c(-c1ccc(Oc3ccccc3)cc1)cn2C1CCCC1",                                           3200.0, 5.50),
    ("CHEMBL69358",  "CCN(CC)CCNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",              1.2,    8.92),
    ("CHEMBL69358",  "CCN(CC)CCNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",              200.0,  6.70),
    ("CHEMBL69358",  "CCN(CC)CCNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)[nH]1",              1290.0, 5.89),
    ("CHEMBL137189", "C/N=N/Nc1ccc2ncnc(Nc3cccc(C)c3)c2c1",                                                   200.0,  6.70),
    ("CHEMBL137617", "C/N=N/Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                  39.0,   7.41),
    ("CHEMBL137027", "CN(C)/N=N/c1ccc2ncnc(NCc3ccccc3)c2c1",                                                  482.0,  6.32),
    ("CHEMBL78257",  "N#CC(C#N)=C(N)/C(C#N)=C/c1cc(O)c(O)c(O)c1",                                            800.0,  6.10),
    ("CHEMBL309334", "N#CC(C#N)=C(N)/C(C#N)=C/c1ccc(O)c(O)c1",                                               2500.0, 5.60),
    ("CHEMBL304971", "Cc1c(C(=O)N2CCN(C)CC2)c[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",              410.0,  6.39),
    ("CHEMBL69071",  "Cc1cc(C(=O)N2CC(C)NC(C)C2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",           5.9,    8.23),
    ("CHEMBL69071",  "Cc1cc(C(=O)N2CC(C)NC(C)C2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc(F)c(Cl)c3)c21",           2580.0, 5.59),
    ("CHEMBL137534", "CN(C)CC/N=N/Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                            38.0,   7.42),
    ("CHEMBL137788", "C/N=N/Nc1ccc2ncnc(NCc3ccccc3)c2c1",                                                     335.0,  6.47),
    ("CHEMBL424093", "Brc1cccc(Nc2ncnc3ccc(N/N=N/Cc4ccccc4)cc23)c1",                                          25.0,   7.60),
    ("CHEMBL77030",  "N#C/C(=C\\c1ccc(O)c(O)c1)C(N)=S",                                                       2400.0, 5.62),
    ("CHEMBL264382", "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(ccn4Cc4ccccc4)c3)c21",      1.1,    8.96),
    ("CHEMBL264382", "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(ccn4Cc4ccccc4)c3)c21",      200.0,  6.70),
    ("CHEMBL264382", "Cc1cc(C(=O)NCCN2CCOCC2)[nH]c1/C=C1\\C(=O)Nc2ncnc(Nc3ccc4c(ccn4Cc4ccccc4)c3)c21",      1950.0, 5.71),
    ("CHEMBL305782", "Cc1[nH]c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(F)c(Cl)c4)c32)c(C)c1CCC(=O)O",                    180.0,  6.75),
    ("CHEMBL308672", "CCN(CC)CCNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(Cl)cc4F)c32)[nH]1",                450.0,  6.35),
    ("CHEMBL308672", "CCN(CC)CCNC(=O)c1cc(C)c(/C=C2\\C(=O)Nc3ncnc(Nc4ccc(Cl)cc4F)c32)[nH]1",                500.0,  6.30),
    ("CHEMBL133023", "Cc1cccc(Nc2ncnc3ccc(/N=N/N(C)C)cc23)c1",                                                24.0,   7.62),
    ("CHEMBL136674", "COCC/N=N/Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                               71.0,   7.15),
    ("CHEMBL73820",  "N#C/C(=C\\c1ccc(O)c(O)c1)C(=O)O",                                                       25000.0,4.60),
    ("CHEMBL77803",  "N#CC(C#N)=C(N)/C(C#N)=C/c1cc(O)c(O)c(Br)c1",                                           500.0,  6.30),
    ("CHEMBL432941", "COc1cccc(-c2ccc3c(c2)NC(=O)/C3=C\\c2[nH]c3c(c2CCC(=O)O)CCCC3)c1",                      8880.0, 5.05),
    ("CHEMBL322395", "COc1cc(Nc2c(C#N)cnc3cc(-c4ccc(CN5CCOCC5)cc4)ccc23)c(Cl)cc1Cl",                         3000.0, 5.52),
    ("CHEMBL431996", "O=C(O)CCc1c(/C=C2\\C(=O)Nc3ccc(C(=O)O)cc32)[nH]c2c1CCCC2",                             56500.0,4.25),
    ("CHEMBL86531",  "NS(=O)(=O)c1ccc2c(c1)/C(=C/c1[nH]c3c(c1CCC(=O)O)CCCC3)C(=O)N2",                       29100.0,4.54),
    ("CHEMBL138691", "C=CC(=O)N(CCN(C)C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                       4.2,    8.38),
    ("CHEMBL138691", "C=CC(=O)N(CCN(C)C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                       2282.0, 5.64),
    ("CHEMBL52829",  "Cc1cccc(Nc2ncnc3cnc(N(C)C)nc23)c1",                                                     4.0,    8.40),
    ("CHEMBL52823",  "Cc1cccc(Nc2ncnc3cnc(Cl)nc23)c1",                                                        380.0,  6.42),
    ("CHEMBL139307", "C=CC(=O)N(CCCN1CCOCC1)c1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                   3.3,    8.48),
    ("CHEMBL139307", "C=CC(=O)N(CCCN1CCOCC1)c1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                   194.0,  6.71),
    ("CHEMBL322298", "CCn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccccc3)nc21",                                      200.0,  6.70),
    ("CHEMBL111434", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(Br)c3)nc21",                                   670.0,  6.17),
    ("CHEMBL111197", "COc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",                                   30.0,   7.52),
    ("CHEMBL136058", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCN1CCOCC1",                            0.81,   9.09),
    # --- Batch 100-149 (offset=100) ---
    ("CHEMBL136058", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCN1CCOCC1",                            8.8,    8.06),
    ("CHEMBL136058", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCN1CCOCC1",                            8.8,    8.06),
    ("CHEMBL51854",  "COc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       3.8,    8.42),
    ("CHEMBL939",    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",                                       515.0,  6.29),
    ("CHEMBL95320",  "C#Cc1cccc(Nc2c(C#N)cnc3cc(OCCOC)c(OCCOC)cc23)c1",                                       850.0,  6.07),
    ("CHEMBL298931", "Cc1cccc(Nc2ncnc3cnc(N)nc23)c1",                                                         17.0,   7.77),
    ("CHEMBL292323", "COc1cc2ncc(C(=O)Nc3ccccc3)c(SSc3c(C(=O)Nc4ccccc4)c4cccc(OC)c4n3C)n(C)c12",               6500.0, 5.19),
    ("CHEMBL334334", "COC(=O)c1cccc(NCc2cc(O)ccc2O)c1",                                                       9000.0, 5.05),
    ("CHEMBL545315", "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1.Cl.Cl",                         1.5,    8.82),
    ("CHEMBL545315", "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1.Cl.Cl",                         7.4,    8.13),
    ("CHEMBL545315", "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1.Cl.Cl",                         7.4,    8.13),
    ("CHEMBL111365", "COc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cn1",                                   1300.0, 5.89),
    ("CHEMBL52591",  "CNc1ncc2ncnc(Nc3cccc(C)c3)c2n1",                                                        4.3,    8.37),
    ("CHEMBL132907", "COC(=O)c1cc(NCc2ccccc2O)ccc1O",                                                         20000.0,4.70),
    ("CHEMBL422292", "CCOC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                      1.5,    8.82),
    ("CHEMBL422292", "CCOC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                      51.0,   7.29),
    ("CHEMBL113023", "COc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1C",                                  250.0,  6.60),
    ("CHEMBL417420", "Brc1cccc(Nc2ncnc3cnc(NCCCN4CCOCC4)nc23)c1",                                             2.9,    8.54),
    ("CHEMBL304414", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3ccccc3n2C)c(C(=O)Nc2ccccc2)c2ccccc21",                     10000.0,5.00),
    ("CHEMBL106232", "CN1CCN(CCCCCNc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)CC1",                               4500.0, 5.35),
    ("CHEMBL327127", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(CO)c3)nc21",                                   90.0,   7.05),
    ("CHEMBL294475", "Brc1cccc(Nc2ncnc3cnc(NCCN4CCOCC4)nc23)c1",                                              0.81,   9.09),
    ("CHEMBL294475", "Brc1cccc(Nc2ncnc3cnc(NCCN4CCOCC4)nc23)c1",                                              3.1,    8.51),
    ("CHEMBL293261", "O=C(CCc1c(SSc2[nH]c3ccccc3c2CCC(=O)NCc2ccccc2)[nH]c2ccccc12)NCc1ccccc1",               850.0,  6.07),
    ("CHEMBL53310",  "CNc1ncc2ncnc(Nc3ccccc3)c2n1",                                                           13.0,   7.89),
    ("CHEMBL267019", "O=C1NC(=O)c2cc(Nc3ccc(O)cc3)c(Nc3ccc(O)cc3)cc21",                                      2500.0, 5.60),
    ("CHEMBL285063", "C=CC(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                0.7,    9.15),
    ("CHEMBL285063", "C=CC(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                2.7,    8.57),
    ("CHEMBL136404", "C/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             0.5,    9.30),
    ("CHEMBL136404", "C/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             7.7,    8.11),
    ("CHEMBL291979", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cc(C#N)ccc3n2C)c(C(=O)Nc2ccccc2)c2cc(C#N)ccc21",           54100.0,4.27),
    ("CHEMBL62843",  "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cccc(O)c3n2C)c(C(=O)Nc2ccccc2)c2cccc(O)c21",               100000.0,4.00),
    ("CHEMBL334697", "C=CC(=O)N(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             0.17,   9.77),
    ("CHEMBL334697", "C=CC(=O)N(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             13.0,   7.89),
    ("CHEMBL434827", "CCN(CC)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                             0.73,   9.14),
    ("CHEMBL434827", "CCN(CC)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                             21.0,   7.68),
    ("CHEMBL434827", "CCN(CC)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                             21.0,   7.68),
    ("CHEMBL64430",  "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cc([N+](=O)[O-])ccc3n2C)c(C(=O)Nc2ccccc2)c2cc([N+](=O)[O-])ccc21", 5000.0, 5.30),
    ("CHEMBL137082", "C=CC(=O)N(CCCN1CCOCC1)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                   2.7,    8.57),
    ("CHEMBL137082", "C=CC(=O)N(CCCN1CCOCC1)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                   156.0,  6.81),
    ("CHEMBL324926", "COc1cc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc(OC)c1",                               330.0,  6.48),
    ("CHEMBL53311",  "Clc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       82.0,   7.09),
    ("CHEMBL62176",  "CN1C(=S)C(C(=O)Nc2ccccc2)c2ccccc21",                                                    1000.0, 6.00),
    ("CHEMBL423159", "COC(=O)c1cc(NCc2cccc(O)c2)ccc1O",                                                       21000.0,4.68),
    ("CHEMBL299622", "CNc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       0.76,   9.12),
    ("CHEMBL299622", "CNc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       30.0,   7.52),
    ("CHEMBL138363", "C=C(C)C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             1.6,    8.80),
    ("CHEMBL138363", "C=C(C)C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             44.0,   7.36),
    ("CHEMBL93464",  "C/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                             0.55,   9.26),
    ("CHEMBL93464",  "C/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                             8.7,    8.06),
    # --- Batch 150-199 (offset=150) ---
    ("CHEMBL53831",  "Nc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                        1.5,    8.82),
    ("CHEMBL299881", "Cc1cccc(Nc2ncnc3cnc(NCCN4CCOCC4)nc23)c1",                                               2.3,    8.64),
    ("CHEMBL50027",  "Clc1ncc2ncnc(Nc3ccccc3)c2n1",                                                           2550.0, 5.59),
    ("CHEMBL31965",  "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1",                                74.0,   7.13),
    ("CHEMBL344486", "CCOC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                      2.7,    8.57),
    ("CHEMBL344486", "CCOC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                      64.0,   7.19),
    ("CHEMBL138940", "O=C(/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cn1)NCCCN1CCOCC1",                         0.61,   9.21),
    ("CHEMBL138940", "O=C(/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cn1)NCCCN1CCOCC1",                         14.0,   7.85),
    ("CHEMBL138940", "O=C(/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cn1)NCCCN1CCOCC1",                         14.0,   7.85),
    ("CHEMBL321193", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(O)c3)nc21",                                    150.0,  6.82),
    ("CHEMBL419022", "CCOC(=O)CCCc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",                         700.0,  6.16),
    ("CHEMBL53426",  "Brc1cccc(Nc2ncnc3cnc(NCCc4c[nH]cn4)nc23)c1",                                            0.25,   9.60),
    ("CHEMBL53426",  "Brc1cccc(Nc2ncnc3cnc(NCCc4c[nH]cn4)nc23)c1",                                            5.1,    8.29),
    ("CHEMBL53711",  "CN(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                    2.6,    8.59),
    ("CHEMBL53711",  "CN(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                    12.0,   7.92),
    ("CHEMBL344177", "O=C(/C=C/Cl)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            0.69,   9.16),
    ("CHEMBL344177", "O=C(/C=C/Cl)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            20.0,   7.70),
    ("CHEMBL136178", "CC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                        1.2,    8.92),
    ("CHEMBL136178", "CC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                        1039.0, 5.98),
    ("CHEMBL136511", "C=CS(=O)(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                            1.4,    8.85),
    ("CHEMBL136511", "C=CS(=O)(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                            2.7,    8.57),
    ("CHEMBL416307", "Cc1cccc(Nc2ncnc3cnc(NCCc4c[nH]cn4)nc23)c1",                                             3.0,    8.52),
    ("CHEMBL44",     "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12",                                                  1000.0, 6.00),
    ("CHEMBL136102", "CN(C)CCCNC(=O)/C=C/C(=O)N(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                            1.45,   8.84),
    ("CHEMBL136102", "CN(C)CCCNC(=O)/C=C/C(=O)N(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                            193.0,  6.71),
    ("CHEMBL336264", "C=CS(=O)(=O)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             0.43,   9.37),
    ("CHEMBL104153", "COc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",                                   220.0,  6.66),
    ("CHEMBL53156",  "CNc1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       3.1,    8.51),
    ("CHEMBL53156",  "CNc1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       20.0,   7.70),
    ("CHEMBL95774",  "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)c(C#N)cnc2cc1OCCCN1CCOCC1",                           80.0,   7.10),
    ("CHEMBL131653", "COC(=O)c1cc(NCc2cc(O)ccc2O)ccc1O",                                                      9000.0, 5.05),
    ("CHEMBL268868", "O=C1NC(=O)c2cc(Nc3ccccc3)c(Nc3ccccc3)cc21",                                             300.0,  6.52),
    ("CHEMBL7775",   "CN(c1ccccc1)c1cc2c(cc1Nc1ccccc1)C(=O)NC2=O",                                            8700.0, 5.06),
    ("CHEMBL53711",  "CN(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                    21.0,   7.68),
    ("CHEMBL61119",  "CC(=O)Oc1ccc2c(c1)c(C(=O)Nc1ccccc1)c(SSc1c(C(=O)Nc3ccccc3)c3cc(OC(C)=O)ccc3n1C)n2C",  5300.0, 5.28),
    ("CHEMBL335440", "O=S(=O)(CCO)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                             93.5,   7.03),
    ("CHEMBL50647",  "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(CC(=O)O)cc3)nc21",                             220.0,  6.66),
    ("CHEMBL320705", "COc1ccccc1Nc1ncc2cc(-c3c(Cl)cccc3Cl)c(=O)n(C)c2n1",                                     66.0,   7.18),
    ("CHEMBL51228",  "Brc1cccc(Nc2ncnc3cnc(NCCCn4ccnc4)nc23)c1",                                              2.3,    8.64),
    ("CHEMBL129930", "COC(=O)c1cc(NCc2cc(O)ccc2O)ccc1OC",                                                     9500.0, 5.02),
    ("CHEMBL113863", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(NCCCCCC(=O)O)nc21",                                    1000.0, 6.00),
    ("CHEMBL553",    "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",                                            1450.0, 5.84),
    ("CHEMBL92086",  "C=C(C)C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                             1.2,    8.92),
    ("CHEMBL92086",  "C=C(C)C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                             16.0,   7.80),
    ("CHEMBL109631", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(C(=O)O)c3)nc21",                               190.0,  6.72),
    ("CHEMBL52015",  "CN(C)c1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                    84.0,   7.08),
    ("CHEMBL53796",  "CN(C)c1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                    9.6,    8.02),
    ("CHEMBL53796",  "CN(C)c1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                    32.0,   7.50),
    ("CHEMBL51853",  "Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                        0.13,   9.89),
    ("CHEMBL51853",  "Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                        16.0,   7.80),
    # --- Batch 200-249 (offset=200) ---
    ("CHEMBL293913", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3ccc(O)cc3n2C)c(C(=O)Nc2ccccc2)c2ccc(O)cc21",              44000.0,4.36),
    ("CHEMBL133283", "COC(=O)c1ccc(NCc2cc(O)ccc2O)cc1O",                                                      10000.0,5.00),
    ("CHEMBL7914",   "Cc1ccc(Nc2cc3c(cc2Nc2ccc(C)cc2)C(=O)NC3=O)cc1",                                        2500.0, 5.60),
    ("CHEMBL7914",   "Cc1ccc(Nc2cc3c(cc2Nc2ccc(C)cc2)C(=O)NC3=O)cc1",                                        4500.0, 5.35),
    ("CHEMBL334801", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                               1.1,    8.96),
    ("CHEMBL334801", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                               57.0,   7.24),
    ("CHEMBL334801", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                               57.0,   7.24),
    ("CHEMBL325589", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccnc3)nc21",                                       510.0,  6.29),
    ("CHEMBL99024",  "COc1cc2ncc(C#N)c(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",                                  1040.0, 5.98),
    ("CHEMBL59812",  "CC(=O)Oc1cccc2c1c(C(=O)Nc1ccccc1)c(SSc1c(C(=O)Nc3ccccc3)c3c(OC(C)=O)cccc3n1C)n2C",   20000.0,4.70),
    ("CHEMBL434828", "O=C(/C=C/c1ccccc1)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                      9.1,    8.04),
    ("CHEMBL434828", "O=C(/C=C/c1ccccc1)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                      77.0,   7.11),
    ("CHEMBL49596",  "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(N)nc21",                                               5600.0, 5.25),
    ("CHEMBL50519",  "COc1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                       4.3,    8.37),
    ("CHEMBL293250", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cc(O)ccc3n2C)c(C(=O)Nc2ccccc2)c2cc(O)ccc21",              40500.0,4.39),
    ("CHEMBL62701",  "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cccnc3n2C)c(C(=O)Nc2ccccc2)c2cccnc21",                    22300.0,4.65),
    ("CHEMBL51659",  "C=CC(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                0.91,   9.04),
    ("CHEMBL51659",  "C=CC(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                3.4,    8.47),
    ("CHEMBL135142", "O=C(/C=C/C(F)(F)F)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                      1.75,   8.76),
    ("CHEMBL135142", "O=C(/C=C/C(F)(F)F)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                      35.0,   7.46),
    ("CHEMBL343722", "C=C=CC(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                              1.6,    8.80),
    ("CHEMBL343722", "C=C=CC(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                              120.0,  6.92),
    ("CHEMBL448154", "CCn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccncc3)nc21",                                      490.0,  6.31),
    ("CHEMBL50470",  "CCN(CC)CCCNc1ncc2cc(-c3c(Cl)cccc3Cl)c(=O)n(C)c2n1",                                    1400.0, 5.85),
    ("CHEMBL104244", "CN1CCN(CCCCNc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)CC1",                               210.0,  6.68),
    ("CHEMBL53777",  "CN(C)c1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                    0.95,   9.02),
    ("CHEMBL53777",  "CN(C)c1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                    22.0,   7.66),
    ("CHEMBL293889", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cc(Cl)ccc3n2C)c(C(=O)Nc2ccccc2)c2cc(Cl)ccc21",            4300.0, 5.37),
    ("CHEMBL140561", "CN(C)CCCOC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               2.4,    8.62),
    ("CHEMBL140561", "CN(C)CCCOC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               108.0,  6.97),
    ("CHEMBL140561", "CN(C)CCCOC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               108.0,  6.97),
    ("CHEMBL345077", "C=C[S+]([O-])c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            4.6,    8.34),
    ("CHEMBL345077", "C=C[S+]([O-])c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            340.0,  6.47),
    ("CHEMBL107472", "CNc1ncc2cc(-c3c(Cl)cccc3Cl)c(=O)n(C)c2n1",                                             7500.0, 5.12),
    ("CHEMBL139044", "C=C/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                           1.1,    8.96),
    ("CHEMBL139044", "C=C/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                           27.0,   7.57),
    ("CHEMBL139095", "C=CS(=O)(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            0.76,   9.12),
    ("CHEMBL139095", "C=CS(=O)(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                            2.4,    8.62),
    ("CHEMBL113070", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(CCCC(=O)O)cc3)nc21",                           80.0,   7.10),
    ("CHEMBL104779", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccncc3)nc21",                                       910.0,  6.04),
    ("CHEMBL319620", "O=C(O)c1cc(NCc2cc(O)ccc2O)ccc1O",                                                       92000.0,4.04),
    ("CHEMBL134371", "O=C(O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                        0.37,   9.43),
    ("CHEMBL301612", "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccccc3)nc21",                                       260.0,  6.58),
    ("CHEMBL52765",  "Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                        0.78,   9.11),
    ("CHEMBL50344",  "CNc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                       11.0,   7.96),
    ("CHEMBL293251", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3cc(Br)ccc3n2C)c(C(=O)Nc2ccccc2)c2cc(Br)ccc21",            11400.0,4.94),
    ("CHEMBL60472",  "COc1ccc2c(C(=O)Nc3ccccc3)c(SSc3c(C(=O)Nc4ccccc4)c4ccc(OC)cc4n3C)n(C)c2c1",            3600.0, 5.44),
    ("CHEMBL342828", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCn1ccnc1",                             0.56,   9.25),
    ("CHEMBL342828", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCn1ccnc1",                             14.0,   7.85),
    ("CHEMBL342828", "O=C(/C=C/C(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cn1)NCCCn1ccnc1",                             14.0,   7.85),
    # --- Batch 250-299 (offset=250) ---
    ("CHEMBL299194", "CN1CCN(CCCNc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)CC1",                                8700.0, 5.06),
    ("CHEMBL54067",  "CN(C)CCNc1ncc2ncnc(Nc3cccc(Br)c3)c2n1",                                                 35.0,   7.46),
    ("CHEMBL291359", "CC(=O)Oc1cccc2c(C(=O)Nc3ccccc3)c(SSc3c(C(=O)Nc4ccccc4)c4cccc(OC(C)=O)c4n3C)n(C)c12", 9900.0, 5.00),
    ("CHEMBL109296", "CCNc1ncc2cc(-c3c(Cl)cccc3Cl)c(=O)n(C)c2n1",                                            4500.0, 5.35),
    ("CHEMBL111339", "Cc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",                                   480.0,  6.32),
    ("CHEMBL52197",  "COc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                                                       30.0,   7.52),
    ("CHEMBL7939",   "O=C1NC(=O)c2cc(Nc3ccc(F)cc3)c(Nc3ccc(F)cc3)cc21",                                      700.0,  6.16),
    ("CHEMBL7939",   "O=C1NC(=O)c2cc(Nc3ccc(F)cc3)c(Nc3ccc(F)cc3)cc21",                                      2200.0, 5.66),
    ("CHEMBL8095",   "COc1ccc(Nc2cc3c(cc2Nc2ccccc2)C(=O)NC3=O)cc1",                                          2000.0, 5.70),
    ("CHEMBL8095",   "COc1ccc(Nc2cc3c(cc2Nc2ccccc2)C(=O)NC3=O)cc1",                                          18000.0,4.75),
    ("CHEMBL345109", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               0.44,   9.36),
    ("CHEMBL345109", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               59.0,   7.23),
    ("CHEMBL345109", "CN(C)CCCNC(=O)/C=C/C(=O)Nc1ccc2ncnc(Nc3cccc(Br)c3)c2c1",                               59.0,   7.23),
    ("CHEMBL607707", "CCOc1cc2ncc(C#N)c(Nc3ccc(F)c(Cl)c3)c2cc1NC(=O)/C=C/CN(C)C",                            83.0,   7.08),
    ("CHEMBL293986", "Cn1c(SSc2c(C(=O)Nc3ccccc3)c3c(Cl)cccc3n2C)c(C(=O)Nc2ccccc2)c2c(Cl)cccc21",            100000.0,4.00),
    ("CHEMBL8223",   "CCN(CC)c1ccc(Nc2cc3c(cc2Nc2ccc(N(CC)CC)cc2)C(=O)NC3=O)cc1",                            50000.0,4.30),
    ("CHEMBL332269", "CCCCNc1ncc2cc(-c3c(Cl)cccc3Cl)c(=O)n(C)c2n1",                                           4400.0, 5.36),
    ("CHEMBL300217", "Nc1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                        7.6,    8.12),
    ("CHEMBL300217", "Nc1ccc2ncnc(Nc3cccc(Br)c3)c2n1",                                                        53.0,   7.28),
    ("CHEMBL53753",  "CNc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",                                                       15.0,   7.82),
    ("CHEMBL428690", "CN1CC[C@H](c2c(O)cc(O)c3c(=O)cc(-c4ccccc4Cl)oc23)[C@H](O)C1",                          7600.0, 5.12),
    ("CHEMBL14932",  "Oc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1",                                                98.0,   7.01),
    ("CHEMBL38199",  "Brc1cccc(Nc2ncnc3cc4c(cnn4CCN4CCOCC4)cc23)c1",                                          3.7,    8.43),
    ("CHEMBL38199",  "Brc1cccc(Nc2ncnc3cc4c(cnn4CCN4CCOCC4)cc23)c1",                                          192.0,  6.72),
    ("CHEMBL553",    "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",                                            1.0,    9.00),
    ("CHEMBL939",    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",                                       1.0,    9.00),
    ("CHEMBL29197",  "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC",                                                     0.025,  10.60),
    ("CHEMBL287007", "CN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                          2.6,    8.59),
    ("CHEMBL287007", "CN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                          9.7,    8.01),
    ("CHEMBL37346",  "Brc1cccc(Nc2ccnc3ccccc23)c1",                                                           5500.0, 5.26),
    ("CHEMBL15346",  "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1",                              74.0,   7.13),
    ("CHEMBL15202",  "COc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1",                                               68.0,   7.17),
    ("CHEMBL554",    "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1",                 10.0,   8.00),
    ("CHEMBL36819",  "O=C(O)Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             53.0,   7.28),
    ("CHEMBL36819",  "O=C(O)Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             1196.0, 5.92),
    ("CHEMBL40130",  "Brc1cccc(Nc2ncnc3cc4c(ccn4CCN4CCOCC4)cc23)c1",                                          3.7,    8.43),
    ("CHEMBL40130",  "Brc1cccc(Nc2ncnc3cc4c(ccn4CCN4CCOCC4)cc23)c1",                                          16.0,   7.80),
    ("CHEMBL14627",  "O=C(NCCCCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1)C(F)(F)F",                             74.0,   7.13),
    ("CHEMBL290096", "Brc1cccc(Nc2ncnc3ccccc23)c1",                                                           27.0,   7.57),
    ("CHEMBL460731", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc(F)c(Cl)c3)c2c1",                                   20.0,   7.70),
    ("CHEMBL460731", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc(F)c(Cl)c3)c2c1",                                   71.0,   7.15),
    ("CHEMBL37543",  "OCC(O)Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             12.0,   7.92),
    ("CHEMBL37543",  "OCC(O)Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             318.0,  6.50),
    ("CHEMBL287289", "CN(C)CCCn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                            21.0,   7.68),
    ("CHEMBL287289", "CN(C)CCCn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                            230.0,  6.64),
    ("CHEMBL514436", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc(OCc4cccc(F)c4)c(Cl)c3)c2c1",                      27.0,   7.57),
    ("CHEMBL39320",  "Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                                   0.37,   9.43),
    ("CHEMBL39320",  "Cn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                                   13.0,   7.89),
    ("CHEMBL289162", "O=C(O)Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             5.1,    8.29),
    ("CHEMBL289162", "O=C(O)Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             1780.0, 5.75),
    # --- Batch 300-349 (offset=300) ---
    ("CHEMBL14952",  "CS(=O)(=O)CCNCCCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1",                               104.0,  6.98),
    ("CHEMBL37373",  "CN(C)CCn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             41.0,   7.39),
    ("CHEMBL37373",  "CN(C)CCn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             145.0,  6.84),
    ("CHEMBL39715",  "COC(=O)CN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                   3.4,    8.47),
    ("CHEMBL39715",  "COC(=O)CN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                   9.6,    8.02),
    ("CHEMBL289213", "CN(CC(=O)O)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                    0.72,   9.14),
    ("CHEMBL289213", "CN(CC(=O)O)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                    584.0,  6.23),
    ("CHEMBL273789", "CS(=O)(=O)CCN(CCCCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1)C(=O)C(F)(F)F",              263.0,  6.58),
    ("CHEMBL517907", "C#Cc1cccc(Nc2ncnc3ccc(OCCCCNCCS(C)(=O)=O)cc23)c1",                                      11.0,   7.96),
    ("CHEMBL517907", "C#Cc1cccc(Nc2ncnc3ccc(OCCCCNCCS(C)(=O)=O)cc23)c1",                                      45.0,   7.35),
    ("CHEMBL460732", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc(S(=O)(=O)c4ccccc4)cc3)c2c1",                      93.0,   7.03),
    ("CHEMBL540701", "C=CCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1",                                            79.0,   7.10),
    ("CHEMBL287832", "CN(C)CCn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             40.0,   7.40),
    ("CHEMBL287832", "CN(C)CCn1ncc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             262.0,  6.58),
    ("CHEMBL39337",  "Brc1cccc(Nc2ncnc3cc4[nH]ncc4cc23)c1",                                                   0.44,   9.36),
    ("CHEMBL39337",  "Brc1cccc(Nc2ncnc3cc4[nH]ncc4cc23)c1",                                                   20.0,   7.70),
    ("CHEMBL39337",  "Brc1cccc(Nc2ncnc3cc4[nH]ncc4cc23)c1",                                                   0.44,   9.36),
    ("CHEMBL40734",  "Brc1cccc(Nc2ncnc3cc4[nH]ccc4cc23)c1",                                                   0.44,   9.36),
    ("CHEMBL40734",  "Brc1cccc(Nc2ncnc3cc4[nH]ccc4cc23)c1",                                                   22.0,   7.66),
    ("CHEMBL40734",  "Brc1cccc(Nc2ncnc3cc4[nH]ccc4cc23)c1",                                                   0.44,   9.36),
    ("CHEMBL39811",  "CN(C)CCN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                    7.5,    8.12),
    ("CHEMBL39811",  "CN(C)CCN(C)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                    81.0,   7.09),
    ("CHEMBL36727",  "Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                                   0.8,    9.10),
    ("CHEMBL36727",  "Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                                   22.0,   7.66),
    ("CHEMBL39355",  "Brc1cccc(Nc2ncnc3cc4c(cnn4CCCN4CCOCC4)cc23)c1",                                        8.8,    8.06),
    ("CHEMBL39355",  "Brc1cccc(Nc2ncnc3cc4c(cnn4CCCN4CCOCC4)cc23)c1",                                        77.0,   7.11),
    ("CHEMBL493428", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc4c(cnn4Cc4cccc(F)c4)c3)c2c1",                      24.0,   7.62),
    ("CHEMBL493428", "CS(=O)(=O)CCNCCCCOc1ccc2ncnc(Nc3ccc4c(cnn4Cc4cccc(F)c4)c3)c2c1",                      630.0,  6.20),
    ("CHEMBL289959", "c1ccc(Nc2ncnc3ccccc23)cc1",                                                             346.0,  6.46),
    ("CHEMBL47966",  "Oc1ccc(/N=N/c2ccc(O)cc2O)cc1",                                                          58000.0,4.24),
    ("CHEMBL14874",  "CN(CCCOc1ccc2ncnc(Nc3ccc(OCc4ccccc4)cc3)c2c1)CCS(C)(=O)=O",                            830.0,  6.08),
    ("CHEMBL416611", "OCC(O)Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             1.6,    8.80),
    ("CHEMBL416611", "OCC(O)Cn1ccc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",                                             49.0,   7.31),
    ("CHEMBL36967",  "Brc1cccc(Nc2ncnc3cc4[nH]cc(CN5CCOCC5)c4cc23)c1",                                       4.8,    8.32),
    ("CHEMBL36967",  "Brc1cccc(Nc2ncnc3cc4[nH]cc(CN5CCOCC5)c4cc23)c1",                                       9.6,    8.02),
    ("CHEMBL36164",  "OCCN(CCO)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                      3.5,    8.46),
    ("CHEMBL36164",  "OCCN(CCO)Cc1c[nH]c2cc3ncnc(Nc4cccc(Br)c4)c3cc12",                                      135.0,  6.87),
    ("CHEMBL428690", "CN1CC[C@H](c2c(O)cc(O)c3c(=O)cc(-c4ccccc4Cl)oc23)[C@H](O)C1",                          22000.0,4.66),
    ("CHEMBL338448", "CO[C@H]1[C@@H](N(C)C(=O)c2ccccc2)C[C@@H]2O[C@@]1(C)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 1900.0, 5.72),
    ("CHEMBL131715", "CO[C@H]1[C@@H](N(C)C(=O)c2ccc(Cl)cc2)C[C@@H]2O[C@@]1(C)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 7000.0, 5.16),
    ("CHEMBL435809", "CO[C@H]1[C@@H](N(C)C(=O)c2ccc(C(=O)O)cc2)C[C@@H]2O[C@@]1(C)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 110.0, 6.96),
    ("CHEMBL57323",  "Cc1cccc(C)c1-c1cc2cnc(N)nc2nc1NC(=O)NC(C)(C)C",                                        610.0,  6.21),
    ("CHEMBL57347",  "CN1CCN(CCCNc2ncc3cc(-c4c(Cl)cccc4Cl)c(NC(=O)NC(C)(C)C)nc3n2)CC1",                      150.0,  6.82),
    ("CHEMBL75188",  "Cc1cc(C)c(C)c(-c2cc3cnc(NCCCN4CCN(C)CC4)nc3nc2NC(=O)NC(C)(C)C)c1C",                   7000.0, 5.16),
    ("CHEMBL131020", "CO[C@H]1[C@@H](N(C)S(C)(=O)=O)C[C@@H]2O[C@@]1(C)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 1100.0, 5.96),
    ("CHEMBL131590", "CO[C@H]1[C@@H](N(C)C(=O)c2ccccc2Cl)C[C@@H]2O[C@@]1(C)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 4700.0, 5.33),
    ("CHEMBL160846", "CS(=O)(=O)CCNc1nc(-c2ccc3ncnc(Nc4ccc(OCc5ccccc5)c(Cl)c4)c3c2)cs1",                    8.0,    8.10),
    ("CHEMBL160846", "CS(=O)(=O)CCNc1nc(-c2ccc3ncnc(Nc4ccc(OCc5ccccc5)c(Cl)c4)c3c2)cs1",                    230.0,  6.64),
    ("CHEMBL444337", "COC(=O)CN(C)[C@H]1C[C@@H]2O[C@](C)([C@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", 1300.0, 5.89),
    ("CHEMBL422347", "CS(=O)(=O)CCNc1nc(-c2ccc3ncnc(Nc4ccc(OCc5ccccc5)cc4)c3c2)cs1",                        82.0,   7.09),
]


def build_dataframe():
    """
    Build a pandas DataFrame from the raw ChEMBL activity data.

    Returns
    -------
    pd.DataFrame with columns:
        molecule_chembl_id  : ChEMBL compound identifier
        smiles              : canonical SMILES string
        ic50_nm             : IC50 in nanomolar (standard_value)
        pIC50               : -log10(IC50_nM * 1e-9)  i.e. 9 - log10(IC50_nM)
        pchembl_value       : pChEMBL value as reported by ChEMBL (may differ
                              slightly due to rounding or assay-type adjustments)
    """
    import pandas as pd

    records = []
    for mol_id, smiles, ic50_nm, pchembl in RAW:
        # Recompute pIC50 from IC50 for transparency
        pic50 = round(9.0 - math.log10(ic50_nm), 2) if ic50_nm > 0 else None
        records.append({
            "molecule_chembl_id": mol_id,
            "smiles":             smiles,
            "ic50_nm":            ic50_nm,
            "pIC50":              pic50,
            "pchembl_value":      pchembl,
        })

    df = pd.DataFrame(records)
    return df


def build_unique_compounds():
    """
    Return a de-duplicated DataFrame keeping the most potent measurement
    per compound (lowest IC50 / highest pIC50).
    """
    df = build_dataframe()
    best = (
        df.sort_values("ic50_nm")
          .groupby("molecule_chembl_id", as_index=False)
          .first()
    )
    return best.reset_index(drop=True)


def fetch_egfr_data(limit_per_page=50, max_offset=300):
    """
    Reproduce this dataset live from the ChEMBL REST API.

    Parameters
    ----------
    limit_per_page : int
        Page size for each API request.
    max_offset : int
        Stop after fetching this many records (use 0 for all ~18 000).

    Returns
    -------
    pd.DataFrame with columns molecule_chembl_id, smiles, ic50_nm, pIC50.
    """
    import requests
    import pandas as pd

    base = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    params = dict(
        target_chembl_id="CHEMBL203",
        standard_type="IC50",
        standard_units="nM",
        pchembl_value__isnull="false",
        limit=limit_per_page,
        format="json",
    )

    all_records = []
    offset = 0
    while True:
        params["offset"] = offset
        resp = requests.get(base, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        activities = data.get("activities", [])
        if not activities:
            break
        for act in activities:
            ic50 = act.get("standard_value")
            smi  = act.get("canonical_smiles")
            mol  = act.get("molecule_chembl_id")
            pch  = act.get("pchembl_value")
            if ic50 and smi and mol:
                ic50 = float(ic50)
                all_records.append({
                    "molecule_chembl_id": mol,
                    "smiles":             smi,
                    "ic50_nm":            ic50,
                    "pIC50":              round(9.0 - math.log10(ic50), 2),
                    "pchembl_value":      float(pch) if pch else None,
                })
        offset += limit_per_page
        if max_offset and offset >= max_offset:
            break

    return pd.DataFrame(all_records)


# ---------------------------------------------------------------------------
# Quick summary when run as a script
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    df = build_dataframe()
    print(f"Total records (with duplicates across assays): {len(df)}")
    print(f"Unique compounds: {df['molecule_chembl_id'].nunique()}")
    print(f"\npIC50 range: {df['pIC50'].min():.2f} – {df['pIC50'].max():.2f}")
    print(f"IC50 range:  {df['ic50_nm'].min()} nM – {df['ic50_nm'].max()} nM\n")
    print("pIC50 distribution (binned):")
    bins = [4, 5, 6, 7, 8, 9, 10, 11]
    labels = ["4-5", "5-6", "6-7", "7-8", "8-9", "9-10", "10-11"]
    import pandas as pd
    df["pIC50_bin"] = pd.cut(df["pIC50"], bins=bins, labels=labels, right=False)
    print(df["pIC50_bin"].value_counts().sort_index().to_string())
    print()
    print("First 10 rows:")
    print(df[["molecule_chembl_id", "smiles", "ic50_nm", "pIC50"]].head(10).to_string())
