#!/bin/sh
# In MT3 (cross sections), retrieve
# (n,elastic) MF2 (' 3  2') and
# (n,gamma) MF102 (' 3102')

rm -rf strip
mkdir -p strip

for file in ../PyNjoy2016/output/SANDY/H1_H2O/1-H-1g.jeff33-*
do
    base=${file##*/}
    # Remove header (3 lines)
    # Remove right characters columns (beyond 67 characters)
    # Get them two-by-two on each line (energy and cross section)
    # Replace every "-" with the usual scientific notation "E-" /!\ It does not deal with E+ values
    # Remove empty lines
    grep -E '^.{70} 3  2' $file | tail -n +4 | cut -c 1-67 | awk '{print $1" "$2"\n"$3" "$4"\n"$5" "$6}' | sed 's/-/E-/g' | sed 's/+/E+/g' | sed '/^[[:space:]]*$/d' > strip/"$base"_nelastic
    grep -E '^.{70} 3102' $file | tail -n +4 | cut -c 1-67 | awk '{print $1" "$2"\n"$3" "$4"\n"$5" "$6}' | sed 's/-/E-/g' | sed 's/+/E+/g' | sed '/^[[:space:]]*$/d' > strip/"$base"_ngamma
done
