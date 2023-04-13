email=$1

n_countries=$(wc -l < "bundles.csv")

ncpu=8
mem_free=10

mkdir -p derived

for ((i=1; i<=n_countries; i++)); do
    iso3c=$(awk -F',' -v i=1 -v j=$i 'NR==j {print $i}' "bundles.csv")
    id=$(awk -F',' -v i=2 -v j=$i 'NR==j {print $i}' "bundles.csv")

    #check if output exists already
    if [ -f "derived/${id}" ]; then
        echo "${iso3c} already run"
    else
        #check if it's been unpacked
        id_nozip=$(echo ${id} | sed 's/.zip//')
        if [ -d "derived/${id_nozip}" ]; then
            rm -rf "derived/${id_nozip}"
        fi
        qsub -V -cwd -M ${email} -m ase -N excess_fit_${iso3c} -q parallel.q -pe smp ${ncpu} -l mem_free=${mem_free}G,h_vmem=${mem_free}.2G ./script.sh ${id}
    fi
done
