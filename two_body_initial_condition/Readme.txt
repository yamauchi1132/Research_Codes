######################################
# setting two body initail condition #
######################################

you need to edit user_difine.hpp to set parameters.
you can run this program as follow.
#######################################################################
./two_body_initial.out data1 (if two of sphere are same)              #
./two_body_initial.out data1 data2  (if two of sphere are not same)   #
#######################################################################

Please be careful when you input 2 different files. Don't reverse data1 and data1 when you input.
data1 is inputed into m1, r1, pori_num1 and data2 is into m2, r2, pori_num2 respectively which are defined in user_define.hpp
That is
./two_body_initial.out data2 data1 is not collect.


Data file name means,
data _ ParticleNumber _ MassOfM1 _ RadiOfR1_PoriNum1_MassOfM2 _ RadiOfR2 _ PoriNum2_Rp(C*(R1+R2))_VelRelinf


