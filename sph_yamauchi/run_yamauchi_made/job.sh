#ftype=data/r001k/init/s1.00_t000
#pori=1.5
#ftype=polytropic_initial_condition/data_n_$pori/poly.00_t000
ftype=data/polytropic/r015k/poly.00_t000
#ftype=data/two_body/r003k/r_p_1.4e+11/two_body.00_t000

frst=0
ffle=0

odir=snap
fexe=run

export OMP_NUM_THREADS=4
#export OMP_NUM_THREADS=1

if ! test -e $odir
then
    mkdir $odir
fi

cp $0 $odir/job.sh
cp Makefile $odir/
cp run $odir/

./"$fexe" "$frst" "$ffle" "$ftype"
