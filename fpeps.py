import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass
import fUD as UD
import time
import math
from scipy import optimize
L=1.0
N_x=4
N_y=4
La_S=(L*L)/((N_x+1.0)*(N_x+1.0))
La_S=1.0
#No-symmetry
D=[4]
chi_boundry=[20]
chi_single=[25]
chi_try=[10]
d_in=[2]

#Z2-symmetry
D=[2,2]
chi_boundry=[20]
chi_single=[60]
chi_try=[60]
d_in=[2,2]

#U1-symmetry
# D=[1,2,1]
# chi_boundry=[30]
# chi_single=[40]
# chi_try=[30]
# d_in=[2,2]

##########################################
interval=+1.0e-2
threshold=[1.0,1.0e+1]
accuracy=+1.0e-9
Model=[ "Fer_Z2", "on", "off"]               #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1, FFI_Z2, Fer_BOS, Fer_BOS_Z2#, "off, Sin, Tanh, SinShift, TanhShift"
N_iter_total=1
N_tebd=[100, "on", "on"]    #last one is printing (make code slow)

N_part=2.0
N_sz=N_part-2
RG_Int_Coupl_final=-1.0
E_fermi=math.pi*(N_part/(N_x*N_x))
E_unit=0.5*N_part*E_fermi

h_coupling=[ -1.0/La_S, 4.0/La_S, -3.0, 2.*RG_Int_Coupl_final, 3.50]
#Last ones: inteaction (g) and magnetic field (h)

Sys=['Fer','single', "QR", 15, "cg", "iden", 20, 140, "TEBD_SVD", "Z2" ]                 #"Swap:Fer, Bos", "Opt_energy:double, single, simple", "Opt_truncation=QR, Inv, Grad", QR_update_iteration,  "Inv=SVD, cg" , "Q_init=rand, iden, part", Q_update_iteration, Grad_update_iteration,  "Previous, TEBD_SVD"

#start_itebd=La_S*4.0
start_itebd=1.4
division_itebd=5.0
N_iter_SU=[50, "QR"]     #"full" or "QR"

#print Model, Sys, "h=", h_coupling, "D=", D, "chi_single", chi_single, "chi_boundry", chi_boundry

Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
E_iter_listRG=[]
E_iter_listRG1=[]
E_iter_list1=[]
N_list=[]
h_parameter=[]

time_list=[]
Fidel_val=1
E_coulmn=[]
E_row=[]
E_mag_coulmn=[]
E_mag_row=[]
E_0=1.0
E_1=1.0
E_00=1.0
E_11=1.0
E_list_try=[]
E_try1=[]
Sz_list=[]
Up_list=[]
Down_list=[]

PEPS_mps=[None]*(N_x)
PEPS_mps_left=[None]*(N_x)
PEPS_mps_right=[None]*(N_x)
PEPS_listten=[None]*(N_x)

PEPS_mps_leftU=[None]*(N_x)
PEPS_mps_rightU=[None]*(N_x)
PEPS_listtenU=[None]*(N_x)

mps_boundry_left=[None]*(N_x)
mps_boundry_right=[None]*(N_x)
mps_boundry_temp=[None]*(N_x)

mps_boundry_left1=[None]*(N_x)
mps_boundry_right1=[None]*(N_x)
mps_boundry_temp1=[None]*(N_x)

q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in=UD.full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_in)

print "q_d_in", q_d_in
print "q_D", q_D

for i in xrange(N_x):
 PEPS_listten[i] =UD.Init_PEPS( N_x, q_D, q_d_in, i)
 PEPS_listtenU[i]=UD.Init_PEPS( N_x, q_D, q_d_in, i)

#PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d_in, threshold, interval,Sys)
#print norm_val

for iter in xrange(N_iter_total):
 for i in xrange(N_x):
  PEPS_listten[i] =UD.Init_PEPS( N_x, q_D, q_d_in, i)
  PEPS_listtenU[i]=UD.Init_PEPS( N_x, q_D, q_d_in, i)

 h_coupling[2]=h_coupling[2]
 print h_coupling
 H_col=UD.make_H_col( N_x, q_d_in, h_coupling, Model, L)
 H_row=UD.make_H_row( N_x, q_d_in, h_coupling, Model, L)
 H_long=UD.make_H_long( N_x, q_d_in, h_coupling, Model)

 N_col=UD.make_N_col( N_x, q_d_in, h_coupling, Model)
 N_row=UD.make_N_row( N_x, q_d_in, h_coupling, Model)
 N_long=UD.make_N_long( N_x, q_d_in, h_coupling, Model)

 Sz_col=UD.make_Sz_col( N_x, q_d_in, h_coupling, Model)
 Sz_row=UD.make_Sz_row( N_x, q_d_in, h_coupling, Model)


 Up_col=UD.make_Up_col( N_x, q_d_in, h_coupling, Model)
 Up_row=UD.make_Up_row( N_x, q_d_in, h_coupling, Model)

 Down_col=UD.make_Down_col( N_x, q_d_in, h_coupling, Model)
 Down_row=UD.make_Down_row( N_x, q_d_in, h_coupling, Model)



 particle_col=UD.make_particle_col( N_x, q_d_in, h_coupling, Model)
 Magz_col=UD.make_Magz_col( N_x, q_d_in, h_coupling, Model)

 Magz_col_direct=UD.make_Magz_col_direct( N_x, q_d_in, h_coupling, Model)
 particle_col_direct=UD.make_particle_col_direct( N_x, q_d_in, h_coupling, Model)


 ############# Start: init guess ##################
 Landa_col=UD.Landa_f_col( q_D, N_x, "U") #Q, U, UQ, H, I
 Landa_row=UD.Landa_f_row( q_D, N_x, "U") #Q, U, UQ, H, I

 PEPS_listten=UD.Randomize_pepes(PEPS_listten, "U", N_x, N_y)  #Q, U, UQ, H, I

 #UD.Landa_f_row_rebonding(PEPS_listten, Landa_row, N_x)
 #UD.Landa_f_col_rebonding(PEPS_listten, Landa_col, N_x)
 #Landa_col=UD.Landa_f_col_iden( Landa_col,N_x/2)
 #Landa_row=UD.Landa_f_row_iden( Landa_row,N_x/2)
 #UD.Reload_fRG( PEPS_listten, N_x)

 # UD.Store_Gamma( PEPS_listten, N_x)
 # UD.Store_Landa_row( Landa_row, N_x)
 # UD.Store_Landa_col( Landa_col, N_x)

# UD.Reload_Gamma( PEPS_listten, N_x)
# UD.Reload_Landa_row( Landa_row, N_x)
# UD.Reload_Landa_col( Landa_col, N_x)

 PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, H_col, H_row, q_d_in, Model, N_x, q_D, q_chi_try, threshold, interval, Sys)
 UD.Store_Gamma( PEPS_listten, N_x)
 UD.Store_Landa_row( Landa_row, N_x)
 UD.Store_Landa_col( Landa_col, N_x)

 PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)

 ############# End: init guess ##################

 #UD.Reload_fRG( PEPS_listten, N_x)
 #UD.Store_fRG(PEPS_listten, N_x)

 #PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d_in, threshold, interval,Sys)
 #print 'norm_val', norm_val



 print iter
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 t0=time.time()



 if Sys[1] is "double":
  rho_row, rho_col=UD.make_density_matrix_double( PEPS_listten, N_x, q_chi_single, q_d_out, q_D,Sys)
 if Sys[1] is "single":
  rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x, q_chi_single, q_d_in, q_D,Sys)

 E_0=UD.Energy_from_Density(rho_row, rho_col, H_col, H_row, N_x)
 print "pure energy=", E_0

 N_0=UD.Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
 print "total particle=", N_0
 N_list.append(N_0)

 Sz_0=UD.Energy_from_Density(rho_row, rho_col, Sz_col, Sz_row, N_x)
 print "total mag=", Sz_0
 Up_0=UD.Energy_from_Density(rho_row, rho_col, Up_col, Up_row, N_x)
 print "up spins=", Up_0
 Down_0=UD.Energy_from_Density(rho_row, rho_col, Down_col, Down_row, N_x)
 print "Down spins=", Down_0
 print "Down_spins+up_spins", Down_0+Up_0
 Up_list.append(Up_0)
 Down_list.append(Down_0)
 
 print  "no chamical/mag energy",(E_0-h_coupling[2]*N_part-h_coupling[3]*N_sz)*0.5


 print "regularized energy", ((E_0-h_coupling[2]*N_part-h_coupling[3]*N_sz)*0.5)/E_unit
 E_iter_list1.append( E_0 )
 E_iter_list.append(  (E_0-h_coupling[2]*N_part-h_coupling[3]*N_sz)*0.5    )
 E_try1.append(   ((E_0-h_coupling[2]*N_part-h_coupling[3]*N_sz)*0.5) /E_unit    )


 Magz_col_direct_1=UD.Sz_from_Density_direct( rho_col, Magz_col_direct, N_x, Sys)
 particle_col_direct_1=UD.Sz_from_Density_direct( rho_col, particle_col_direct, N_x, Sys)

 sumZ=0
 file = open("Data/Zdirect.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str(index) + " "+str(index1)+" " + str(Magz_col_direct_1[index][index1])+" "+ "\n")
   sumZ=Magz_col_direct_1[index][index1]+sumZ
 file.close()
 #print "sumZ", sumZ
 Sz_list.append(sumZ)

 sum=0
 file = open("Data/particle.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str((index+1.0)*(L/(N_x+1))) + " "+str((index1+1.0)*(L/(N_x+1)))+" " + str(particle_col_direct_1[index][index1])+" "+ "\n")
   sum=particle_col_direct_1[index][index1]+sum
 file.close()
 #print "sumN", sum
 




 file = open("Data/Energy.txt", "w")
 for index in xrange(len(E_iter_list)):
  #print index, E_try1[index]
  file.write(str(index) + " " + str(E_iter_list1[index])+" " + str(E_iter_list[index])+" " + str(E_try1[index])+" " +str(N_list[index])+" "+str(Sz_list[index])+" " +str(Up_list[index])+" " +str(Down_list[index])+" "+ "\n")
 file.close()



 if Sys[1] is "single":
   UD.TEBD_Full_RG( H_col, H_row, N_x, PEPS_listten, q_D, accuracy, N_tebd, i, q_d_in, q_chi_single, q_chi_try, threshold, interval, Model, Sys, N_part, h_coupling, N_col, N_row, Sz_col, Sz_row, Up_col, Up_row, Down_col, Down_row, N_sz)





rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x, q_chi_single, q_d_in, q_D, Sys)
N_0=UD.Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
Energy_val=UD.Energy_cal( PEPS_listten, q_d_in, q_chi_single, N_x, q_D, H_col, H_row, Sys)
E_0=Energy_val
print "pure energy=", E_0
print "total particle=", N_0
print  "no chamical/mag energy",(E_0-h_coupling[2]*N_part-h_coupling[3]*N_sz)*0.5

