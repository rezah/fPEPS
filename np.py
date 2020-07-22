import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy import stats
from numpy.linalg import inv
N_list=[]
E_list=[]

L=6.0
E_list=[]
N_list=[]
Start=4
End=17
for N in xrange(Start,End,1):
 template = np.zeros([N, N])
 for idx in np.ndindex(template.shape):
    if idx[0]==idx[1]+1 or idx[0]==idx[1]-1:
      template[idx] = (-1.0)*(N+1)*(N+1)*(1/(L**2))
      #template[idx] = (-1.0)*(N)*(N)*(1/(L**2))
    if idx[0]==idx[1]:
      template[idx] = (2.0*(N+1)*(N+1)*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0
      #template[idx] = (2.0*(N)*(N)*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0

 #print N ,template

 w, v = la.eig(template)

 #print w

# N=10000
# val=0
# for i in xrange(N):
#  val=val+2.0*abs(np.cos( (2*np.pi*i)/N ))
#  print (2*np.pi*i)/N, np.cos( (2*np.pi*i)/N )
# print val, 1/np.pi


 #E_list.append(min(w))
 E_list.append(sorted(w)[1])
 N_list.append(N)
 #print min(w)
 #print  N, min(w)*3+sorted(w)[1]
# print  N,"2", min(w)*3+sorted(w)[1]
# print  N,"4", min(w)*4+sorted(w)[1]*4
# print  N,"6", min(w)*6+sorted(w)[1]*4+sorted(w)[2]*2
# print  N,"5", min(w)*5+sorted(w)[1]*4+sorted(w)[2]

# print  N,"8", min(w)*6+sorted(w)[1]*6+sorted(w)[2]*4
 print  N,"10", min(w)*8+sorted(w)[1]*6+sorted(w)[2]*4+sorted(w)[3]*2
# print  N,"14", min(w)*8+sorted(w)[1]*8+sorted(w)[2]*7+sorted(w)[3]*5
# print  N, min(w)*6+sorted(w)[1]*2

E_list1=[E_list[i]  for i in xrange(0,End-Start)]
N_list1=[(1.0/N_list[i])  for i in xrange(0,End-Start)]


#E_list11=[E_list[i]*N_list[i]*N_list[i]  for i in xrange(len(E_list))]
#N_list11=[(1.0/N_list[i])  for i in xrange(len(E_list))]

#print E_list1, "second"












L=6.0
Start=4
End=10
E_list=[]
N_list=[]
for N in xrange(Start,End):
 template = np.zeros([N, N])
 for idx in np.ndindex(template.shape):
    if idx[0]==idx[1]+1 or idx[0]==idx[1]-1:
      template[idx] = (-16.*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2))
    if idx[0]==idx[1]+2 or idx[0]==idx[1]-2:
      template[idx] = (+1.*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2))
    if idx[0]==idx[1] and idx[1]==0 :
      template[idx] = (((30.0-1)*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0
    elif idx[0]==idx[1] and idx[1]==N-1 :
      template[idx] = (((30.0-1)*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0
    elif idx[0]==idx[1]:
      template[idx] = ((30.0*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0

 w, v = la.eig(template)
 E_list.append(min(w))
 #E_list.append(sorted(w)[1])
 #E_list.append(sorted(w)[2])
# print  N,"2", min(w)*3+sorted(w)[1]
 #print  N, min(w)*3+sorted(w)[1]


# print  N,"8", min(w)*6+sorted(w)[1]*6+sorted(w)[2]*4
# print  N,"6", min(w)*6+sorted(w)[1]*4+sorted(w)[2]*2
# print  N,"14", min(w)*8+sorted(w)[1]*8+sorted(w)[2]*7+sorted(w)[3]*5

# print  N, min(w)*6+sorted(w)[1]*4+sorted(w)[2]*2


# print  N,"4", min(w)*4+sorted(w)[1]*4
 print  N,"10", min(w)*8+sorted(w)[1]*6+sorted(w)[2]*4+sorted(w)[3]*2

# print  N, min(w)*8+sorted(w)[1]*6+sorted(w)[2]*4+sorted(w)[3]*2

# print  N, min(w)*6+sorted(w)[1]*2

 N_list.append(N)
 #print template

E_list11=[E_list[i]  for i in xrange(0,End-Start)]
N_list11=[(1.0/N_list[i])  for i in xrange(0,End-Start)]

#print E_list11, "fourth"

E_list=[]
N_list=[]

N_list=[4,6,8,10,14,20,30,40,60,80]
E_list=[-51.741+80,-45.7603+80, -42.329+80, -40.11+80,-37.52+80, -35.48+80, -33.16+80, -32.77+80,-32.28436+80,-31.87436+80]

#N_list=[10,14,20,30,40,60,80]
#E_list=[-40.11+80,-37.52+80, -35.48+80, -33.16+80, -32.77+80,-32.28436+80,-31.87436+80]

N_list=[20,30,40,60,80]
E_list=[ -35.48+80, -33.16+80, -32.77+80,-32.28436+80,-31.87436+80]


#one partice in 2D box
#N_list=[4,8,12,16,20]
#E_list=[ +12.223, +15.482,+16.948,+17.6778,+18.315]


#N_list=[8,12,16,20]
#E_list=[  +15.482,+16.948,+17.6778,+18.315]


#E_list11=[E_list[i]  for i in xrange(len(N_list))]
#N_list11=[1.0/N_list[i]  for i in xrange(len(N_list))]

#E_list1=[E_list[i]  for i in xrange(len(N_list))]
#N_list1=[(1.0/N_list[i])  for i in xrange(len(N_list))]


#one partice in 2D box,  Ex=19.7392088022
#N_list=[4,6,8,10,12,14,16]
#E_list=[ +13.4723, +15.33482,+16.67948,+17.5778,+18.0315,18.2,18.6]

#one partice in 2D box, Fixed Ex=19.7392088022
#N_list=[4,6,8]
#E_list=[ 19.098276435, 19.4102838312,19.5386108703]


#one partice in 2D box, Fourth, Fixed Ex=19.7392088022, nu=25
#N_list=[4,6,8]
#E_list=[ 19.706209036464891, 19.73046928519, 19.738725]


#E_list11=[E_list[i]  for i in xrange(len(N_list))]
#N_list11=[(1.0/N_list[i])  for i in xrange(len(N_list))]

#two fermions in 2D box, fixed second order, Ex=69.08723
#N_list=[4,6,8]
#E_list=[63.19671496, 66.0308683954, 67.2855733]

#two fermions in 2D box,Fixed Fourth_order, Ex=69.08723
#N_list=[4,6,8,10,12]
#E_list=[68.08823, 68.81620, 69.020062195]



#two fermions in 2D box, Fourth_order, Ex=69.08723
#N_list=[4,6,8,10,12]
#E_list=[46.0274536036, 53.2528765615, 56.6956499265, 59.1175204027, 60.4940956636]


#N_list=[10,12]
#E_list=[ 59.1175204027,60.5883814552]

E_list11=[E_list[i]  for i in xrange(len(N_list))]
N_list11=[(1.0/N_list[i])  for i in xrange(len(N_list))]
E_list1=[E_list[i]  for i in xrange(len(N_list))]
N_list1=[(1.0/N_list[i])  for i in xrange(len(N_list))]


#print E_list11, N_list11




#print E_list, N_list
#plt.plot( N_list1,E_list1,'g+',label='2D box, second order' )
#plt.plot( N_list11, E_list11,'g.',label='2D box, fourth order' )
#plt.plot( Count_list3, variance_list3,'yh',label='CG-poly' )
#plt.plot( Count_list4, variance_list4,'c+',label='SteepestDescent-poly' )

# x = np.linspace( 0.0, 1.0/4, 100)
# 
# m1,b1 = np.polyfit(N_list1, E_list1, 1)
# m,b = np.polyfit(N_list11, E_list11, 1)
# 
# 
# y = m*x+b
# y1 = m1*x+b1
# 
# plt.plot(x, y, '-b', label='')
# #plt.plot(x, y1, '-g', label='')
# 
# print "fourth", m,b
# print "second", m1,b1
# 
# 
# plt.xlabel('$\epsilon$', fontsize=20)
# plt.ylabel('$E_{\epsilon}$', fontsize=20)
# plt.legend(loc='upper right')
# plt.savefig('Extra2D.pdf')
# #plt.show()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 




E_list=[]
N_list=[]
Start=100
End=101
for N in xrange(Start,End,1):
 template_A = np.zeros([N, N])
 for idx in np.ndindex(template_A.shape):
    if idx[0]==idx[1]+1 or idx[0]==idx[1]-1:
      template_A[idx] = (-1.0)
    if idx[0]==idx[1]:
      template_A[idx] = +2.0

 #print N ,template_A
 

 template_B = np.zeros([N, N])
 for idx in np.ndindex(template_B.shape):
    if idx[0]==idx[1]+1 or idx[0]==idx[1]-1:
      template_B[idx] = (1.0)
    if idx[0]==idx[1]:
      template_B[idx] = (10.0)


 template_V = np.zeros([N, N])
 for idx in np.ndindex(template_B.shape):
    if idx[0]==idx[1]:
      template_V[idx] = idx[0]



 #print N ,template_V

 template_B=template_B*(1.0/12.0)
 template_V=template_V*(0.0/N)
 #print template_B
 template_B_inv = inv(template_B)
 template=(np.dot(template_B_inv, template_A))
 template=(template*(N+1)*(N+1))+template_V
 w, v = la.eig(template)
 #print np.dot(template_B,template_B_inv)
 #print w
 
# N=10000
# val=0
# for i in xrange(N):
#  val=val+2.0*abs(np.cos( (2*np.pi*i)/N ))
#  print (2*np.pi*i)/N, np.cos( (2*np.pi*i)/N )
# print val, 1/np.pi


 E_list.append(min(w))
 #E_list.append(sorted(w)[2])
 N_list.append(N)
 #print min(w)
E_list1=[E_list[i]  for i in xrange(0,End-Start)]

#N_list1=[(1.0/N_list[i])  for i in xrange(0,End-Start)]





print E_list1

#m1,b1 = np.polyfit(N_list1, E_list1, 1)


#y = m*x+b

#plt.plot(x, y, '-b', label='')

#print "fourth", m,b













