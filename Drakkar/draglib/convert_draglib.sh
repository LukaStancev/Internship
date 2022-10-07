#!/bin/bash
rm -f DRGLIB draglib* # Remove any previous binary draglib that may persist
cd ../data
echo "Transforming ASCII-formatted draglib into binary-formatted draglib:"
ascii_draglibs=" ../../PyNjoy2016/output/JEFF-3.3_shem295/draglib
                 ../../PyNjoy2016/output/JEFF-3.3/draglib
                 ../../PyNjoy2016/output/JEFF-3.1.1/draglib
                 ../../PyNjoy2016/output/TENDL-2019/*/draglib*
                 ../../PyNjoy2016/output/SANDY/*/draglib* "
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
        if [[ $ascii_draglib == *"JEFF-3.3_shem295"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/drglibJEFF-3.3_295
        elif [[ $ascii_draglib == *"JEFF-3.3"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/drglibJEFF-3.3
        elif [[ $ascii_draglib == *"JEFF-3.1.1"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/drglibJEFF-3.1.1
        elif [[ $ascii_draglib == *"TENDL"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/drglib$(basename $ascii_draglib | cut -c 8-)
        elif [[ $ascii_draglib == *"SANDY"* ]]
        then
            mv ../draglib/DRGLIB ../draglib/drglib$(basename $ascii_draglib | cut -c 8-)
        else
            echo $ascii_draglib " unknown."
            exit 1
        fi
    fi
done


