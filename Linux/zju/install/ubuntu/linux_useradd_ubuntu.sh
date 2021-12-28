#!/bin/bash

adduser_and_data () {
    adduser $1 --disabled-password
    echo "$1:1st-time-$1" | chpasswd
    passwd -e $1
    mkdir /data/$1
    chown $1:$1 /data/$1
}

change_user_pw () {
    echo "$1:1st-time-$1" | chpasswd
    passwd -e $1
}

add_keys () {
    mkdir -p /home/$1/.ssh
    rm -f /home/$1/.ssh/authorized_keys
    touch /home/$1/.ssh/authorized_keys
    chmod 700 /home/$1/.ssh
    chown $1:$1 /home/$1/.ssh
    chmod 600 /home/$1/.ssh/authorized_keys
    chown $1:$1 /home/$1/.ssh/authorized_keys
    cat zju.pub > /home/$1/.ssh/authorized_keys
}

add_power () {
    cat utopia > /home/$1/.ssh/utopia
    chmod 600 /home/$1/.ssh/utopia
    chown $1:$1 /home/$1/.ssh/utopia
    sed -i 's/^alias.*//g' /home/$1/.bashrc
    echo alias helloworld=\"ssh -i \~/.ssh/utopia root@localhost\" >> /home/$1/.bashrc
}

if true; then
  mkdir -p .ssh
  cat utopia.pub > .ssh/authorized_keys
  chmod 700 .ssh
  chmod 600 .ssh/authorized_keys
fi

echo ">These are all the users:"
for user in $(cat ~/user-add.txt); do
  echo "$user"
  adduser_and_data $user
  add_keys $user
  change_user_pw $user
  ln -s /data/$user /home/$user/data
#  echo -e "\n#GCC 11.1.0\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/gcc/11.1.0/lib64/\n" >> /home/$user/.bashrc
  echo -e "#Gromacs Environment\nsource /opt/gromacs/2020.6/bin/GMXRC.bash\n" >> /home/$user/.bashrc
  echo -e "#NAMD Environment\nexport PATH=/opt/namd/3.0alpha9/:\$PATH\n" >> /home/$user/.bashrc
  echo -e "#VMD Environment\nexport PATH=/opt/vmd/1.9.4/:\$PATH\n" >> /home/$user/.bashrc
done

echo ">These are the power users:"
for user in $(cat ~/user-power.txt); do
  echo "$user"
  add_power $user
done
