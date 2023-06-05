#!/bin/sh
i=1
file_no_limit=21

while [ $i -lt $file_no_limit ]
do

reb_output_file="REBOUND_output/rebound${i}.out"
impact_parameters_file="impact_parameters/impact_parameters${i}.txt"
coll_report_raw_file="REBOUND_output/collision_report_raw${i}.txt"
ejections_file="ejections/ejections${i}.txt"
velocities_file="velocities/velocities${i}.txt"

#cat $reb_output_file | grep b/Rtarg >> $impact_parameters_file
#cat $reb_output_file | grep EJECTION: >> $ejections_file
#cat $reb_output_file | grep -E 'TIME|Target|Projectile|Mlr/Mt:|TYPE|FRAG' >> $coll_report_raw_file
cat $reb_output_file | grep Vimp/Vesc: >> $velocities_file

echo $i
echo $reb_output_file
i=`expr $i + 1`
done
