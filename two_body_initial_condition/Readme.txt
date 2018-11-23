Setting two body initail condition 

you need to edit user_difine.hpp to set parameters.
you can run this program as follow.

./two_body_initial.out data1 (if two of sphere are completely same : case1)   
./two_body_initial.out data1 data2  (if two of sphere are not same : case2)   


Be careful when you input two different argument files (case2).
Don't reverse data1 and data2.
Data1 must have the data of sphere1, and data2 have the data of sphere2.
So you must edit m1 r1 pori_num1 m2 r2 pori_num2 in user_define.hpp so that there is no contradiction.
That is,
./two_body_initial.out data2 data1
is not correct and will be a strange results.

Data file name means,
data _ ParticleNumbier _ MassOfM1 _ RadiOfR1_PoriNum1_MassOfM2 _ RadiOfR2 _ PoriNum2_Rp(C*(R1+R2))_VelRelinf


