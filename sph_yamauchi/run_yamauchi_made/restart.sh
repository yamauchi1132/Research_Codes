#ftype=data/r001k/init/s1.00_t000
#pori=1.5
#ftype=polytropic_initial_condition/data_n_$pori/poly.00_t000
#ftype=data/polytropic/r003k/poly.00_t000
ftype=snap_rp1.4e+11_05/t0020

frst=1

odir=snap
fexe=run

export OMP_NUM_THREADS=4
#export OMP_NUM_THREADS=1

if ! test -e $odir
then
    mkdir $odir
else
    echo "There is $odir directory!"
    exit
fi

cp $0 $odir/job.sh
cp Makefile $odir/
cp run $odir/

./"$fexe" "$frst" "$ftype"
