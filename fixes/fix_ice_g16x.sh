# fix for ice system connectivity issue

#### for g16x
# Removes the pattern: [space] 123 [space] 0.100
# from the EIn file

# USAGE: sh fix_ice.sh FILENAME.EIn

# Make a copy of .EIn file to a backup folder
#FNAME=`echo "$1" | awk -F'/' '{print $NF}'`
#cp $1 /home/b/b32708/s1s12_amoeba/ein/${FNAME}

# Remove partial connectivity from .EIn file (sed skips coordination section)
EIN = $1
NATOMS=`awk 'NR==1{print $1}' ${EIN}`
sed -E "1,${NATOMS}b;s/\s+[0-9]+\s+0.100//g" ${EIN}
