#!/bin/bash
#SBATCH -c 5  ## Number of cores to be reserved
#SBATCH -t 0-06:00:00
#SBATCH --mem-per-cpu=4G

rm -rf Fuel/Linux_x86_64 Fuel/save Fuel/kinf Fuel/Saphyb
cd Fuel/data

# EvoNominale
sed -i 's/STRING Save := "TIH_UOX_.*/STRING Save := "TIH_UOX_255" ;/' EvoNominale.x2m
sed -i 's/INTEGER Teneur_I.*/INTEGER Teneur_I := 255    ;/' EvoNominale.x2m
sed -i 's/STRING  ConfigPyrex.*/STRING  ConfigPyrex := "Py8" ;/' EvoNominale.x2m
srun dragon EvoNominale.x2m v5bev1487 &
sleep 10
sed -i 's/STRING Save := "TIH_UOX_.*/STRING Save := "TIH_UOX_255" ;/' EvoNominale.x2m
sed -i 's/INTEGER Teneur_I.*/INTEGER Teneur_I := 255    ;/' EvoNominale.x2m
sed -i 's/STRING  ConfigPyrex.*/STRING  ConfigPyrex := "Py12" ;/' EvoNominale.x2m
srun dragon EvoNominale.x2m v5bev1487 &
sleep 10
sed -i 's/STRING Save := "TIH_UOX_.*/STRING Save := "TIH_UOX_195" ;/' EvoNominale.x2m
sed -i 's/INTEGER Teneur_I.*/INTEGER Teneur_I := 195    ;/' EvoNominale.x2m
sed -i 's/STRING  ConfigPyrex.*/STRING  ConfigPyrex := "None" ;/' EvoNominale.x2m
srun dragon EvoNominale.x2m v5bev1487 &
sleep 10
sed -i 's/STRING Save := "TIH_UOX_.*/STRING Save := "TIH_UOX_310" ;/' EvoNominale.x2m
sed -i 's/INTEGER Teneur_I.*/INTEGER Teneur_I := 310    ;/' EvoNominale.x2m
sed -i 's/STRING  ConfigPyrex.*/STRING  ConfigPyrex := "None" ;/' EvoNominale.x2m
srun dragon EvoNominale.x2m v5bev1487 &
sleep 10
sed -i 's/STRING Save := "TIH_UOX_.*/STRING Save := "TIH_UOX_310" ;/' EvoNominale.x2m
sed -i 's/INTEGER Teneur_I.*/INTEGER Teneur_I := 310    ;/' EvoNominale.x2m
sed -i 's/STRING  ConfigPyrex.*/STRING  ConfigPyrex := "Py12" ;/' EvoNominale.x2m
srun dragon EvoNominale.x2m v5bev1487
## MakeBib1BU
#sed -i 's/EvoName=.*/EvoName="TIH_UOX_195"/' MakeBib1BU.access
#srun dragon MakeBib1BU.x2m v5bev1487
#sleep 10
#sed -i 's/EvoName=.*/EvoName="TIH_UOX_255-Py8"/' MakeBib1BU.access
#srun dragon MakeBib1BU.x2m v5bev1487 &
#sleep 10
#sed -i 's/EvoName=.*/EvoName="TIH_UOX_255-Py12"/' MakeBib1BU.access
#srun dragon MakeBib1BU.x2m v5bev1487 &
#sleep 10
#sed -i 's/EvoName=.*/EvoName="TIH_UOX_310-Py12"/' MakeBib1BU.access
#srun dragon MakeBib1BU.x2m v5bev1487 &
#sleep 10
#sed -i 's/EvoName=.*/EvoName="TIH_UOX_310"/' MakeBib1BU.access
#srun dragon MakeBib1BU.x2m v5bev1487
