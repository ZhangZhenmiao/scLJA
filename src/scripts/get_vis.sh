stats=$1

# get color
python /datahome/datasets/ericteam/zmzhang/Microbe_Seq/scLJA/src/scripts/generate_color_for_bandage.py $stats $stats.color

# get plot
Bandage image /datahome/datasets/ericteam/zmzhang/Microbe_Seq/assemblies/LJA_24ZF07438/03_Polishing/mdbg.gfa $stats.jpg --colors $stats.color --width 2000

echo $stats.jpg