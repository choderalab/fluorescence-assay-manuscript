# data
Raw data for fluorescence assay manuscript.

Table summary of the results we have:

| SINGLE WAVELENGTH  |  Src   | Abl   | Src GK  | Abl GK  | p38   | CA II  |
|:--------:|:------:|:-----:|:-------:|:-------:|:-----:|------- |
| Bos      | YES    | YES   | YES     | YES     | YES   |   x  |
| Bsi      | YES    | YES   | YES     | YES     | YES   |   x  |
| Erl      | YES    | YES   | YES     | YES     | YES   |   x  |
| Gef      | YES    | YES   | YES     | YES     | YES   |   x  |
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
| Axi      | x      |   x |x   |      x |YES   |     x |
| Lap      | x      |   x |x   |      x |YES   |     x |
| Pal      | x      |   x |x   |      x |YES   |     x |
| Paz      | x      |   x |x   |      x |YES   |     x |
| Sar      | x      |   x|x   |      x |x   |    x|
| Van      | x      |   x |x   |      x |x   |    x |
| Afa      | x      |    YES |x   |      x |x   |     x|
| Ner      | x      |    YES |x   |      x |x   |     x|
| Ima      | x | x|x |x  |YES   |    YES |
| Das      |x      |   x |x |x  |YES    |    YES |
| DQA      | x      |    x |x   |      x |YES   |     YES |
| Staur      | x      |   x |x   |      x |YES   |     x |

Below is a description of each set of data.
* `single-wavelength` - data from singlet assays (emission and excitation at a single wavelength)
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
  * `p38_0.5uM` - *p38_0.5uM_8lig1_20161026_132449.xml* and *p38_0.5uM_8lig2_20161026_134006.xml*
     * 8 LIGAND EXPERIMENT - at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on October 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase at 0.5 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * `p38_0.25uM` - *p38_0.25uM_8lig1_20161026_155648* and *p38_0.25uM_8lig2_20161026_161159.xml*
     * 8 LIGAND EXPERIMENT - at a standard 384-well plate fluorescent kinase inhibitor assay with 8 fluorescent ligands
     * performed on October 26, 2016
     * 384 well plate, assay volume 50 uL
     * protein: p38 kinase at 0.25 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH], Ponatinib [IJ], Lapatinib [KL], Saracatinib [MN], and Vandetanib [OP]
  * `p38_0.125uM` - *p38_0.125uM_8lig1_20161026_171600* and *p38_0.125uM_8lig2_20161026_173114.xml*
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
  * `p38_0.5uM` - *p38_Bos_NB_20190709_154049.xml*, *p38_Bos_Iso_NB_20190709_155414.xml*, *p38_Erl_NB_20190709_160747.xml*,  and *p38_Gef_NB_20190709_162114.xml*
     * 4 LIGAND EXPERIMENT - at a standard 96-well plate (Corning 3651) fluorescent kinase inhibitor assay with 4 fluorescent ligands
     * performed on July 9, 2019
     * 96 well plate, assay volume 100 uL
     * protein: p38 kinase at 0.5 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH]
   * `Abl_D382N_0.5uM` - *Abl_Bos_NB_20190710_111828.xml*, *Abl_Bos_Iso_NB_20190710_113155.xml*, *Abl_Erl_NB_20190710_114526.xml*,  and *Abl_Gef_NB_20190710_115855.xml*
     * 4 LIGAND EXPERIMENT - at a standard 96-well plate (Corning 3651) fluorescent kinase inhibitor assay with 4 fluorescent ligands
     * performed on July 10, 2019
     * 96 well plate, assay volume 100 uL
     * protein: Abl D382N kinase at 0.5 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH]
  * `Src_0.5uM` - *Src_Bos_NB_20190710_135049.xml*, *Src_Bos_Iso_NB_20190710_140416.xml*, *Src_Erl_NB_20190710_141748.xml*,  and *Src_Gef_NB_20190710_143120.xml*
     * 4 LIGAND EXPERIMENT - at a standard 96-well plate (Corning 3651) fluorescent kinase inhibitor assay with 4 fluorescent ligands
     * performed on July 10, 2019
     * 96 well plate, assay volume 100 uL
     * protein: Src kinase at 0.5 uM
     * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], Gefitinib [GH]
  * `Abl_D382N_T334I_0.5uM` - *Abl_D382N_T334I_Bos_NB_20190715_161412.xml*, *Abl_D382N_T334I_Bos_Iso_NB_20190715_162627.xml*, *Abl_D382N_T334I_Erl_NB_20190715_164047.xml*,  and *Abl_D382N_T334I_Gef_NB_20190715_165326.xml*
     * 4 LIGAND EXPERIMENT - at a standard 96-well plate (Corning 3651) fluorescent kinase inhibitor assay with 4 fluorescent ligands
     * performed on July 15, 2019
     * 96 well plate, assay volume 100 uL
     * protein: Abl D382N/T334I kinase at 0.5 uM
  * `Src_T338I_0.5uM` - *Src_T338I_Bos_NB_20190715_125826.xml*, *Src_T338I_Bos_Iso_NB_20190715_131151.xml*, *Src_T338I_Erl_NB_20190715_132525.xml*,  and *Src_T338I_Gef_NB_20190715_133852.xml*
     * 4 LIGAND EXPERIMENT - at a standard 96-well plate (Corning 3651) fluorescent kinase inhibitor assay with 4 fluorescent ligands
     * performed on July 15, 2019
     * 96 well plate, assay volume 100 uL
     * protein: Src T338I kinase at 0.5 uM

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
    * `2016-05-26`
      * *Abl_Afa_20160526_161533.xml*, *Abl_Ner_20160526_165224.xml*
      * performed on May 26, 2016
      * 96 well plate, assay volume 100 uL
      * protein: Abl kinase, note D382N mutant makes kinase catalytically inactive
      * ligands: Afatinib [AB], Neratinib [CD]
  * `AblT334I`
    * `2015-12-17` 
      * *AblD382N-T334I_BosI_20151217_bw2020_gain120123031.xml*, *AblD382N-T334I_Bos_20151217_bw2020_gain120_120553.xml*, *AblD382N-T334I_Erl_20151217_bw2020_gain120_125515.xml*, and *AblD382N-T334I_Gef_20151217_bw2020_gain120_132641.xml*
      * performed on December 17, 2015
      * 96 well plate, assay volume 100 uL
      * protein: AblT334I kinase with gatekeeper mutant
      * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib [EF], and Gefitinib [GH]
    * `2019-07-16`
      * *Abl_D382N_T334I_Bos_Iso_UV_20190716_143313.xml*, *Abl_D382N_T334I_Bos_UV_20190716_135721.xml*, *Abl_D382N_T334I_Erl_UV_20190716_151055*, and *Abl_D382N_T334I_Gef_UV_20190716_154656.xml*
      * performed on July 16, 2019
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
     * `2016-12-20`
       * *p38_Sta_ab_20161220_112406.xml*
       * performed on December 20, 2016
       * 96 well plate, assay volume 100 uL
       * protein: p38 kinase
       * ligands: Staurosporine[AB]
	   * Note: The fluorescence excitation and emission values are different than for other experiments, specific to staurosporine's properties.
     * `2017-01-19`
       * *p38_Axi_20170119_144258.xml*, *p38_Lap_20170119_152110.xml*, *p38_Pal_20170119_160546.xml*, and *p38_Paz_20170119_164846.xml*
       * performed on January 19, 2017
       * 96 well plate, assay volume 100 uL
       * protein: p38 kinase
       * ligands: Axitinib [AB], Lapatinib [CD], Palbociclib [EF], and Pazopanib [GH]
	   * Note: Palbociclib is at a different concentration range due to solibility limitations.
* `single-well-assay` - Data from single-well binding assays (could be single wavelength or spectra)  
   * `20171119_p38_single_well_assay`  
      * performed on November 20, 2017  
      * protein: p38 kinase, 1 uM, dialyzed [column 1-6] vs undialyzed [7-11]  
      * ligands: Bosutinib [row A], Bosutinib Isomer [row B], Gefitinib [row C], Erlotinib [row D], Ponatinib [row E], Lapatinib [row F],   Pazopanib [row G], Axitinib [row H]  
      * fluorescence measurement: single wavelenth measurement of the whole plate, spectra measurement of columns 5 and 6  

