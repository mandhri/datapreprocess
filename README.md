   35  $ lsblk
   36  cd /media/
   37  ls
   38  exit
   39  cd
   40  ls
   41  cd /mnt/
   42  ls
   43  mkdir vol1
   44  sudo mkdir vol1
   45  #sudo mount /dev/vdb vol1
   46  sudo fdisk -l
   47  sudo parted /dev/vdb
   48  lsblk
   49  cd /mnt/
   50  ls
   51  sudo mount /dev/vdb1 /mnt/
   52  sudo mount /dev/vdb /mnt/
   53  sudo mkfs.ext4 /dev/vdb1
   54  lsblk
   55  sudo mount /dev/vdb1 /mnt/
   56  ls
   57  sudo umount /dev/vdb1
   58  sudo mount /dev/vdb1 /mnt/vol1/
   59  cd vol1/
   60  ls
   61  mkdir projects
   62  sudo mkdir projects
   63  cd projects/
   64  mkdir project1
   65  cd ..
   66  ls -l
   67  sudo chown -R mandhri:mandhri projects/
   68  ls -l
   69  sudo su mandhri
   70  lsblk
   71  sudo nano /etc/fstab 
   72  blkid
   73  blkid /dev/vdb
   74  blkid /dev/vdb1 
   75  sudo blkid
   76  sudo nano /etc/fstab 
   77  reboot
   78  sudo reboot
   79  lsblk
