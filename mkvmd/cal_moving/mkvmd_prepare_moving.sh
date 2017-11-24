#!/bin/bash

if [ $# -eq 0 ]; then
    echo "make_XXX <rmsf>"
    echo "You better choose 1 only!"
    exit 1
fi

for name in $@; do
    rsync -a /home/kevin/Dropbox/QWD/scripts/Analysis/cal_moving/template-make-cal-moving-$name-sh make_cal_moving.sh
done

# Choose for moving mode
PS3='Please enter your choice: '
options=("Customize" "Backward" "Forward")
select name in "${options[@]}"
do
    case $name in
        "Customize")
            echo "you chose choice 1"
            sed -i 's/^customize=.*$/customize=true/g' make_cal_moving.sh

            read -p "How many ns per window? " WIN_ns
            if ! [[ $WIN_ns =~ ^-?[0-9]+$ ]]; then
                echo "Must be an integer!"
                exit 1
            fi
            sed -i "s/^WIN_ns=.*/WIN_ns=$WIN_ns/g" make_cal_moving.sh
            
            read -p "How many ns per move? " MOVE_ns
            if ! [[ $MOVE_ns =~ ^-?[0-9]+$ ]]; then
                echo "Must be an integer!"
                exit 1
            fi
            sed -i "s/^MOVE_ns=.*/MOVE_ns=$MOVE_ns/g" make_cal_moving.sh

            read -p "How many blocks? " BLOCK
            if ! [[ $BLOCK =~ ^-?[0-9]+$ ]]; then
                echo "Must be an integer!"
                exit 1
            fi
            sed -i "s/^BLOCK=.*/BLOCK=$BLOCK/g" make_cal_moving.sh

            break
            ;;
        "Backward")
            echo "you chose choice 2"
            sed -i 's/^backward=.*$/backward=true/g' make_cal_moving.sh

            read -p "How many ns per window? " WIN_ns
            if ! [[ $WIN_ns =~ ^-?[0-9]+$ ]]; then
                echo "Must be an integer!"
                exit 1
            fi
            sed -i "s/^WIN_ns=.*/WIN_ns=$WIN_ns/g" make_cal_moving.sh
            
            break
            ;;
        "Forward")
            echo "you chose choice 3"
            sed -i 's/^forward=.*$/forward=true/g' make_cal_moving.sh

            read -p "How many ns per window? " WIN_ns
            if ! [[ $WIN_ns =~ ^-?[0-9]+$ ]]; then
                echo "Must be an integer!"
                exit 1
            fi
            sed -i "s/^WIN_ns=.*/WIN_ns=$WIN_ns/g" make_cal_moving.sh

            break
            ;;
        *) echo invalid option;;
    esac
done

read -p "How many ps per frame? [2fs*1000=2] " psperframe
if ! [[ $psperframe =~ ^-?[0-9]+$ ]]; then
    echo "Must be an integer!"
    exit 1
fi
sed -i "s/^psperframe=.*/psperframe=$psperframe/g" make_cal_moving.sh
