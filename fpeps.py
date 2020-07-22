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
N_x=6
N_y=6
La_S=(L*L)/((N_x+1.0)*(N_x+1.0))

#No-symmetry
D=[4]
chi_boundry=[20]
chi_single=[25]
chi_try=[10]
d_in=[2]

#Z2-symmetry
D=[2,2]
chi_boundry=[20]
chi_single=[80]
chi_try=[60]
d_in=[1,1]

#U1-symmetry
D=[1,2,2]
chi_boundry=[30]
chi_single=[70]
chi_try=[30]
d_in=[2,2]

##########################################
interval=+1.0e-2
threshold=[1.0,1.0e+1]
accuracy=+1.0e-9
Model=[ "Fer_U1", "on", "off"]               #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1, FFI_Z2, Fer_BOS, Fer_BOS_Z2#, "off, Sin, Tanh, SinShift, TanhShift"
N_iter_total=1
N_tebd=[100, "on"]

N_part=10.0
RG_Int_Coupl_final=-16.6720
La_S=1.0

h_coupling=[ -1.0/La_S, 4.0/La_S, +12.90, 2.*RG_Int_Coupl_final, 0.0]
#Last ones: inteaction (g) and magnetic field (h)

Sys=['Fer','single', "QR", 15, "cg", "iden", 20, 140, "TEBD_SVD", "U1" ]                 #"Swap:Fer, Bos", "Opt_energy:double, single, simple", "Opt_truncation=QR, Inv, Grad", QR_update_iteration,  "Inv=SVD, cg" , "Q_init=rand, iden, part", Q_update_iteration, Grad_update_iteration,  "Previous, TEBD_SVD"

#start_itebd=La_S*4.0
start_itebd=1.0
division_itebd=5.0
N_iter_SU=[60, "full"]     #"full" or "QR"

#print Model, Sys, "h=", h_coupling, "D=", D, "chi_single", chi_single, "chi_boundry", chi_boundry

Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
E_iter_listRG=[]
E_iter_listRG1=[]
E_iter_list1=[]
count_list=[]
h_parameter=[]
#E_iter_listQ=[]

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

 particle_col=UD.make_particle_col( N_x, q_d_in, h_coupling, Model)
 Magz_col=UD.make_Magz_col( N_x, q_d_in, h_coupling, Model)

 Magz_col_direct=UD.make_Magz_col_direct( N_x, q_d_in, h_coupling, Model)
 particle_col_direct=UD.make_particle_col_direct( N_x, q_d_in, h_coupling, Model)


 ############# Start: init guess ##################
 Landa_col=UD.Landa_f_col( q_D, N_x)
 Landa_row=UD.Landa_f_row( q_D, N_x)
 # UD.Store_Gamma( PEPS_listten, N_x)
 # UD.Store_Landa_row( Landa_row, N_x)
 # UD.Store_Landa_col( Landa_col, N_x)

 UD.Reload_Gamma( PEPS_listten, N_x)
 UD.Reload_Landa_row( Landa_row, N_x)
 UD.Reload_Landa_col( Landa_col, N_x)

# PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, H_col, H_row, q_d_in, Model, N_x, q_D, q_chi_try, threshold, interval, Sys)
 UD.Store_Gamma( PEPS_listten, N_x)
 UD.Store_Landa_row( Landa_row, N_x)
 UD.Store_Landa_col( Landa_col, N_x)

 #PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)

 ############# End: init guess ##################

 UD.Reload_fRG( PEPS_listten, N_x)
 #UD.Store_fRG(PEPS_listten, N_x)

 PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x, q_D, q_chi_try, q_d_in, threshold, interval,Sys)
 print 'norm_val', norm_val



 print iter
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 t0=time.time()



 if Sys[1] is "double" or Sys[1] is "simple":
  rho_row, rho_col=UD.make_density_matrix_double( PEPS_listten, N_x, q_chi_single, q_d_out, q_D,Sys)
 if Sys[1] is "single" or Sys[1] is "simple":
  rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x, q_chi_single, q_d_in, q_D,Sys)

 E_0=UD.Energy_from_Density(rho_row, rho_col, H_col, H_row, N_x)
 print "E_f",E_0

 N_0=UD.Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
 print "N_f",N_0
 count_list.append(N_0)

 Sz_0=UD.Energy_from_Density(rho_row, rho_col, Sz_col, Sz_row, N_x)
 print "Sz_f", Sz_0


 print "E_final",(E_0-h_coupling[2]*N_0)*0.5
 print "E_final",(E_0-h_coupling[2]*N_part)*0.5

 E_peps=(E_0-h_coupling[2]*N_part)*0.5
 E_fermi=math.pi*(N_part/(N_x*N_x))
 E_unit=0.5*N_part*E_fermi
 print "E_total",  E_unit, E_peps/E_unit
 E_iter_list1.append( E_0 )
 E_iter_list.append(  (E_0-h_coupling[2]*N_0)*0.5  )
 E_list_try.append(   (E_0-h_coupling[2]*N_part)*0.5    )
 E_try1.append(   E_peps/E_unit    )


 Magz_col_direct_1=UD.Sz_from_Density_direct( rho_col, Magz_col_direct, N_x, Sys)
 particle_col_direct_1=UD.Sz_from_Density_direct( rho_col, particle_col_direct, N_x, Sys)

 sumZ=0
 file = open("Zdirect.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str(index) + " "+str(index1)+" " + str(Magz_col_direct_1[index][index1])+" "+ "\n")
   sumZ=Magz_col_direct_1[index][index1]+sumZ
 file.close()
 print "sumZ", sumZ
 Sz_list.append(sumZ)

 sum=0
 file = open("particle.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str((index+1.0)*(L/(N_x+1))) + " "+str((index1+1.0)*(L/(N_x+1)))+" " + str(particle_col_direct_1[index][index1])+" "+ "\n")
   sum=particle_col_direct_1[index][index1]+sum
 file.close()
 print "sumN", sum
 


 h_parameter.append(h_coupling[2])
 #count_list.append(N_0)
 print h_parameter, count_list
 file = open("NumberQ.txt", "w")
 for index in range(len(count_list)):
   file.write( str(h_parameter[index])+ " "+ str(index) + " " + str(count_list[index])+ " "+ "\n" )
 file.close()


 time_list.append( time.time() - t0 )
 file = open("time.txt", "w")
 for index in range(len(time_list)):
   file.write(str(index) + " " + str(time_list[index])+" "+ "\n")
 file.close()

 file = open("Energy.txt", "w")
 for index in xrange(len(E_iter_list)):
  #print index, E_try1[index]
  file.write(str(h_parameter[index])+ " "+str(index) + " " + str(E_iter_list[index])+" " + str(E_try1[index])+" " + str(E_list_try[index])+" " + str(E_iter_list1[index])+" "+str(count_list[index])+" "+str(Sz_list[index])+" "+ "\n")
 file.close()



 if Sys[1] is "single":
   UD.TEBD_Full_RG( H_col, H_row, N_x, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d_in, q_chi_single, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model, count_list, E_iter_list1, Sys, N_part, h_coupling, N_col, N_row, Sz_col, Sz_row)





rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x, q_chi_single, q_d_in, q_D, Sys)
N_0=UD.Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
Energy_val=UD.Energy_cal( PEPS_listten, q_d_in, q_chi_single, N_x, q_D, H_col, H_row, Sys)
E_0=Energy_val
print "E_peps", E_0
print "N_f",N_0
print "E_final",(E_0-h_coupling[2]*N_0)*0.5
print "E_final",(E_0-h_coupling[2]*N_part)*0.5

