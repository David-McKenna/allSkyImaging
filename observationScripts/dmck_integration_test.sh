
#!/bin/bash
swlevel 3

bash ./dmck_xst_generic_mode.sh 3 210 5 IntTest
bash ./dmck_xst_generic_mode.sh 3 210 10 IntTest
bash ./dmck_xst_generic_mode.sh 3 210 15 IntTest
bash ./dmck_xst_generic_mode.sh 3 210 20 IntTest

bash ./dmck_xst_generic_mode.sh 5 100 20 IntTest
bash ./dmck_xst_generic_mode.sh 5 100 30 IntTest
bash ./dmck_xst_generic_mode.sh 5 100 45 IntTest
bash ./dmck_xst_generic_mode.sh 5 100 60 IntTest
bash ./dmck_xst_generic_mode.sh 5 100 75 IntTest

bash ./dmck_xst_generic_mode.sh 7 100 20 IntTest
bash ./dmck_xst_generic_mode.sh 7 100 30 IntTest
bash ./dmck_xst_generic_mode.sh 7 100 45 IntTest
bash ./dmck_xst_generic_mode.sh 7 100 60 IntTest
bash ./dmck_xst_generic_mode.sh 7 100 75 IntTest

exit;


