#!/bin/sh
rm -f DRGLIB draglib* # Remove any previous binary draglib that may persist
cd ../data
echo "Transforming ASCII-formatted draglib into binary-formatted draglib:"
ascii_draglibs=" ../../../Njoy/JEFF-3.1.1/draglib
                 ../../../Njoy/TENDL/*/draglib*   "
for ascii_draglib in $ascii_draglibs
do
    if [[ $ascii_draglib != *".bis.gz" ]]
    then
        ascii_draglib=$(readlink -f $ascii_draglib)
        echo $ascii_draglib
        # Modify access file : pointing toward current draglib file
        sed -i "s|ln.*|ln -s $ascii_draglib EXPORT|" convert_draglib.access
        # Run draglib_conversion.x2m
        ../runVersion5.sh convert_draglib.x2m
        # Rename output (binary-formatted) draglib
        if [[ $ascii_draglib == *"JEFF-3.1.1"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/draglibJeff3.1.1
        else # TENDL library
            mv ../draglib/DRGLIB ../draglib/$(basename $ascii_draglib)
        fi
    fi
done


