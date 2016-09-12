Efficiency study package:

=====================
- Developed in D0 W mass analysis
- A general package to study lepton id / iso / trigger efficiency etc. for getting the data/mc efficiency correction.

Orignial Readme from D0 W mass analysis:

 ```
 README_Original_D0
 ```

Original Note:

 ```
 DataEff_v1.1.pdf
 ```

================
Example usage when apply in CMS:

0.) take muon id efficiency data/mc scale factors vs. sumEt as example.

1.) Make the package:

   ```
     gmake all
   ```

   If you updated the codes please:

   ```
   gmake clean
   gmake all
   ```

2.) prepare the root inputs:

 e.g. copy or link inputs files from: 
  ```
   /data/XZZ/76X_Ntuple/76X_effTree/
  ```
  to local working area.

3.) store the histograms for data efficiency extraction:


  - for data:
  ```
   ./storehist_mee.exe config/storehist_mee.config.npnm_data_mu_set 
  ```

  - for background shapes:
  ```
   ./storehist_mee.exe config/storehist_mee.config.npnm_bkgd_mu_set
  ```

  - for signal shapes:
  ```
   ./storehist_mee.exe config/storehist_mee.config.npnm_sgnl_mu_set
  ```


4.) fit Mll spectrum to extract the N pass and N fail:

   ```
    ./fiteff_mee_keys.exe  config/fiteff_mee.config.data_npnm_mu_set
   ```

   This will print you a set of Mll fitting plots in each sumEt bin for N pass and N fail, stored in plots folder.

5.) get the fullSim efficiency N pass and N fail histograms:

   ```
    ./storehist_fullmc.exe config/storehist_fullmc.config.npnm_full_mu_set 
   ```

6.) calculate the data/mc scale factors vs. sumET:

   ```
     root -l macros/make_effratio_npnm_mu_set.C 
   ```

    the example data eff , mc eff, and data/mc eff scale factor plot will be printed here:

    ```
     plots/effratio_npnm_mu_set.pdf
    ```
 
   ```

