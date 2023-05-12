#!/usr/bin/env bash
# created by hzhou@ecust, 2020.9.28

BINDIR=$(realpath $0 | xargs dirname)
LIBDIR=${BINDIR}/../lib
chargemol=${BINDIR}/Chargemol_09_02_2017_linux_parallel
cm5pac=${BINDIR}/cm5pac.exe

# check files exist
for file in AECCAR0 AECCAR2 CHGCAR OUTCAR POSCAR POTCAR; do
    if [ ! -f "$file" ]; then
        echo -e "\n \033[1;31m \`$file\` not exist, please check! \033[0m \n"
        exit 1
    fi
done

NIONS=$(grep "NIONS" OUTCAR | awk -F= '{print $NF}')
Linenu=$(($NIONS / 10 + 1))
POSEND=$(($NIONS + 10 - 1))

cat >job_control.txt <<EOF
<net charge>
0.0 
</net charge>
<periodicity along A, B, and C vectors>
.true. 
.true. 
.false. 
</periodicity along A, B, and C vectors>
<atomic densities directory complete path>
/home/users/hzhou/soft/chargemol/atomic_densities/
</atomic densities directory complete path>
<charge type>
DDEC6 
</charge type>
EOF
atomic_densities_dir=${LIBDIR}/atomic_densities/ # Important: the end of "/" must exist!!! (see https://github.com/pzarabadip/chargemol-light/issues/1)
sed -i "/atomic_densities/c\\$atomic_densities_dir" job_control.txt

echo -e "\n\033[1;31m \tCHGCAR --> DDEC6 process, Please wait a moment...\033[0m\n"
$chargemol

sed -n '1,5p' POSCAR >input
sed -n '9p' POSCAR >>input
grep -A$NIONS 'Missing' VASP_DDEC_analysis.output | grep -v 'Missing' | awk '{print $2}' >temp_atom
sed -n "10,$POSEND p" POSCAR | awk '{print $1,$2,$3}' >temp_pos
paste temp_atom temp_pos >>input
echo "----" >>input
grep -A$Linenu 'computed CM5' VASP_DDEC_analysis.output | grep -v 'CM5' | xargs | sed -n 's/ /\n/gp' >>input

echo -e "\033[1;32m \t\t CM5 Charge Compution process...\033[0m\n"
$cm5pac 4 <input >charge_CM5
echo -e "\033[1;33m \t\t\tGood Luck! ^_^ \033[0m"

rm -rf temp_*
