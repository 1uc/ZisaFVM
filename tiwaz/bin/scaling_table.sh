stems=(o1_ForwardEuler_constant_L4
       o1_ForwardEuler_isentropic_L4
       o2222_SSP2_constant_L4
       o2222_SSP2_isentropic_L4
       o3222_SSP3_constant_L4
       o3222_SSP3_isentropic_L4
       o4222_RK4_constant_L4
       o4222_RK4_isentropic_L4
       o5222_Fehlberg_constant_L4
       o5222_Fehlberg_isentropic_L4)

for s in ${stems[@]}
do
    echo $s
    for n in 2 4 12 36 72 120 240
    # for n in 2 4 12
    do
        f=scaling_experiment_${n}_${s}/run_time.txt
        if [[ -f ${f} ]]
        then
            t=$(cat $f)
            echo "${n} ${t}"
        fi
    done
done
