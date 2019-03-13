力学的エネルギーのテェックプログラム

以下のようにして実行できます。
python energy_check.py data1 data2 time.log (timestep) #

data1は初期データでsph_t0000.datです。data2は最後のデータで初期位置と対称となる位置でのデータです。
time.logも必要です。
Timestepはオプションで基本的には指定しなくて大丈夫です。

energy_chech.pyのview==1とするとエネルギー誤差の時間変動のグラフが描かれます。view==0とすると描かれません。
