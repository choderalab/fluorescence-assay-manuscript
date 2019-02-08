# data
Raw data for fluorescence assay manuscript.

Table summary of the results we have:

| SINGLET  |  Src   | Abl   | Src GK  | Abl GK  | p38   | CA II  |
|:--------:|:------:|:-----:|:-------:|:-------:|:-----:|------- |
| Bos      | x      | x     | x       | x       | YES   |   x    |
| Bsi      | x      | x     | x       | x       | YES   |   x    |
| Erl      | x      | x     | x       | x       | YES   |   x    |
| Gef      | x      | x     | x       | x       | YES   |   x    |
| Pon      | x      | x     | x       | x       | YES   |   x    |
| Lap      | x      | x     | x       | x       | YES   |   x    |
| Sar      | x      |    x |x  |     x |YES    |    x |
| Van      | x      |    x |x  |     x |YES    |     x |
| Ima      | x      | x |x   |      x |x      |     x |
| Das      | x      |   x |x   |      x|x      |     x |
| DQA      | x      |    x |x   |     x |x      |    x |

|   SPECTRA     |  Src     | Abl   |  Src GK   | Abl GK  | p38   | CA II  |
|:-------------:|:--------:|:-----:|:---------:|:-------:|:-----:|------- |
| Bos      | YES | YES | YES  |     YES |YES    |     YES |
| Bsi      | YES      |  YES |YES   |      YES|YES    |    YES |
| Erl      | YES    | YES |YES   |     YES|YES    |    YES |
| Gef      | YES     | YES |YES   |     YES |YES    |    YES |
| Pon      | x | x |x   |      x |YES    |     YES |
| Lap      | x      |   x |x   |      x |x   |     x |
| Sar      | x      |   x|x   |      x |x   |    x|
| Van      | x      |   x |x   |      x |x   |    x |
| Ima      | x | x|x |x  |YES   |    YES |
| Das      |x      |   x |x |x  |YES    |    YES |
| DQA      | x      |    x |x   |      x |YES   |     YES |

Below is a description of each set of data.
* `singlet` - data from singlet assays (emission and excitation at a single wavelength)
  * `DMSO-backfill` - *p38_singlet1_20160420_153238.xml* and *p38_singlet2_20160420_154750.xml*
     * DMSO BACKFILL EXPERIMENT - first half of plate has DMSO backfill, second half does not
     * performed on April 20, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase 
     * ligands: Bosutinib [AB,IJ], Bosutinib Isomer [CD, KL], Erlotinib [EF, MN], and Gefitinib [GH, OP]
  * *p38_8lig1_20160426_120449.xml* and *p38_8lig2_20160426_122008.xml*
     * 8 LIGAND EXPERIMENT - first trial at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on April 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase 
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * *p38_0.5uM_8lig1_20161026_132449.xml* and *p38_0.5uM_8lig2_20161026_134006.xml*
     * 8 LIGAND EXPERIMENT - at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on October 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase at 0.5 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * *p38_0.25uM_8lig1_20161026_155648* and *p38_0.25uM_8lig2_20161026_161159.xml*
     * 8 LIGAND EXPERIMENT - at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on October 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase at 0.25 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * *p38_0.125uM_8lig1_20161026_171600* and *p38_0.125uM_8lig2_20161026_173114.xml*
     * 8 LIGAND EXPERIMENT - at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on October 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase at 0.125 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * `mixing-by-diffusion` - data from time-course fluorescence experiment to test if diffusion can help mixing in 384-well plates
     * EXTRA MIXING BY DIFFUSION EXPERIMENT - p38 binding assay with 8 fluorescent ligands (different than default set)
     * After standard protocol, plate was left for diffusion for indicated time and multiple fluorescence reads were taken in time.
     * Performed on June 15-16, 2017
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib HCl [EF], Gefitinib [GH], Lapatinib [IJ], Ponatinib [KL], Vandetanib (new stock) [MN], and Vandetanib (old stock) [OP] 
     * Without DMSO backfill
     * Time points: 1, 2, 18 and 26 hours after standard plate preparation. 

* `spectra` - data from spectra assays (excitation at a single wavelength, full emission spectra)
  * `Abl`
    * `2015-12-18`
        * *AblD382N_BosI_20151218_bw2020_gain12014-55-55_plate_1.xml*, *AblD382N_Bos_20151218_bw2020_gain120_ 14-41-14_plate_1.xml*, *AblD382N_Erl_20151218_bw2020_gain120 15-11-11_plate_1.xml*, and *AblD382N_Gef_20151218_bw2020_gain120 15-25-58_plate_1.xml*
      * performed on December 18, 2015
      * 96 well plate, assay volume 100 uL
      * protein: Abl kinase, note D382N mutant makes kinase catalytically inactive
      * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
    * `2016-03-11`
      * *Abl_D382N_BosI_20160311_135952.xml*, *Abl_D382N_Bos_20160311_132205.xml*, *Abl_D382N_Erl_20160311_143642.xml*, and *Abl_D382N_Gef_20160311_152340.xml*
      * performed on March 11, 2016
      * 96 well plate, assay volume 100 uL
      * protein: Abl kinase, note D382N mutant makes kinase catalytically inactive
      * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
  * `AblT334I`
    * *AblD382N-T334I_BosI_20151217_bw2020_gain120123031.xml*, *AblD382N-T334I_Bos_20151217_bw2020_gain120_120553.xml*, *AblD382N-T334I_Erl_20151217_bw2020_gain120_125515.xml*, and *AblD382N-T334I_Gef_20151217_bw2020_gain120_132641.xml*
    * performed on December 17, 2015
    * 96 well plate, assay volume 100 uL
    * protein: AblT334I kinase with gatekeeper mutant
    * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
  * `CAII`
     * `2016-03-15`
       * *CAII_BosI_20160315_152656.xml*, *CAII_Bos_20160315_144336.xml*, *CAII_Erl_20160315_160342.xml*, and *CAII_Gef_20160315_164029.xml*
       * performed on March 15, 2016
       * 96 well plate, assay volume 100 uL
       * protein: Carbonic Anhydrase II
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
     * `2016-03-11`
       * *CAII_DQA_20160331_153316.xml*, *CAII_Das_20160331_145548.xml*, *CAII_Ima_20160331_131533.xml*, and *CAII_Pon_20160331_135301.xml*
       * performed on March 11, 2016
       * 96 well plate, assay volume 100 uL
       * protein: Carbonic Anhydrase II
       * ligands: Imatinib [AB], Ponatinib [CD], Dasatinib [EF], and diamino-quinazoline [GH]
  * `Src`
     * `2015-12-15`
       * *Src_BosI_20151215_bw2020_gain120_163633.xml*, *Src_Bos_20151215_bw2020_gain120_161211.xml*, *Src_Erl_20151215_bw2020_gain120_170056.xml*, and *Src_Gef_20151215_bw2020_gain120_172518.xml*
       * performed on December 15, 2015
       * 96 well plate, assay volume 100 uL
       * protein: Src kinase
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
     * `2016-03-09`
       * *Src_BosI_20160309_151610.xml*, *Src_Bos_20160309_143920.xml*, *Src_Erl_20160309_155259.xml*, and *Src_Gef_20160309_163417.xml*
       * performed on March 9, 2016
       * 96 well plate, assay volume 100 uL
       * protein: Src kinase
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
     * `2019-02-06`
       * *WT_Src_4ti-0234_new_Bosutinib_20190206_153028.xml*, *WT_Src_4ti-0234_new_Bosutinib_Isomer_20190206_160.xml*, *WT_Src_4ti-0234_new_Erlotinib_20190206_164401.xml*, and *WT_Src_4ti-0234_new_Gefitinib_20190206_172124.xml*
       * performed on February 6, 2019
       * 96 well plate, assay volume 100 uL
       * 0.5 uM protein: Src kinase
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
  * `SrcT338I`
    * *SrcT338I_BosI_20151216_bw2020_gain120_154934.xml*, *SrcT338I_Bos_20151216_bw2020_gain120_152505.xml*, *SrcT338I_Erl_20151216_bw2020_gain120_161404.xml*, and *SrcT338I_Gef_bw2020_gain120_20151216_164154.xml*
    * performed on December 16, 2015
    * 96 well plate, assay volume 100 uL
    * protein: SrcT338I kinase with gatekeeper mutant
    * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
  * `p38`
     * `2016-01-26`
       * *p38_BosI_20160125_160254.xml*, *p38_Bos_20160125_153819.xml*, *p38_Erl_20160125_163232.xml*, and *p38_Gef_20160125_165832.xml*
       * performed on January 25, 2016
       * 96 well plate, assay volume 100 uL
       * protein: p38 kinase
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
     * `2016-03-07`
       * *p38_BosI_20160307_163847.xml*, *p38_Bos_20160307_160155.xml*, *p38_Erl_20160307_171539.xml*, and *p38_Gef_20160307_175343.xml*
       * performed on March 7, 2016
       * 96 well plate, assay volume 100 uL
       * protein: p38 kinase
       * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
     * `2016-03-30`
       * *p38_DQA_20160330_171906.xml*, *p38_Das_20160330_164208.xml*, *p38_Ima_20160330_152613.xml*, and *p38_Pon_20160330_160326.xml*
       * performed on March 30, 2016
       * 96 well plate, assay volume 100 uL
       * protein: p38 kinase
       * ligands: Imatinib [AB], Ponatinib [CD], Dasatinib [EF], and diamino-quinazoline [GH]
