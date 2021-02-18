# fix for ice system connectivity issue

#### for g16x
# Removes the pattern: [space] 123 [space] 0.100
# from the EIn file

# USAGE: sh fix_ice.sh FILENAME.EIn

sed -i -E "s/\s+[0-9]+\s+0.100//g" $1
