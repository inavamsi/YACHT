import streamlit as st
import pandas as pd
import random
import sys
import pickle
import importlib.util
import os.path as osp
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np


pop=1000
dt=0.1
days=2000


st.sidebar.write("Epidemic parameters")
beta=st.sidebar.slider("Rate of infection", min_value=0.0 , max_value=1.0 , value=0.10 , step=0.01 , format=None , key=None )
gamma=st.sidebar.slider("Rate of recovery", min_value=0.0 , max_value=1.0 , value=0.05 , step=0.01 , format=None , key=None )
initial_infection=st.sidebar.slider("Select inital infection proportion", min_value=0.0 , max_value=1.0 , value=0.01 , step=1/pop , format=None , key=None )
st.sidebar.write("------------------------------------------------------------------------------------")

st.sidebar.write("Testing parameters")

fn=st.sidebar.slider("Select false negative rate", min_value=0.0 , max_value=1.0 , value=0.1 , step=0.1 , format=None , key=None )
fp=st.sidebar.slider("Select false positive rate", min_value=0.0 , max_value=1.0 , value=0.1 , step=0.1 , format=None , key=None )
n=st.sidebar.slider("Select number of pools", min_value=0 , max_value=100 , value=30 , step=1 , format=None , key=None )
frac_poolsize = st.sidebar.slider("Select Poolsize/Population", min_value=0.0 , max_value=1.0/n , value=0.01 , step=1/pop , format=None , key=None )
r = st.sidebar.slider("Select rate of testing per timestep", min_value=0.0 , max_value=1.0 , value=0.5 , step=0.01 , format=None , key=None)
poolsize=frac_poolsize*pop


st.sidebar.write("------------------------------------------------------------------------------------")

badges=['Green','Orange','Red']
states=['S','I','R']

st.sidebar.write("Badge parameters")
delta=st.sidebar.slider("Rate of faking a badge", min_value=0.0 , max_value=1.0 , value=0.01 , step=0.01 , format=None , key=None)

st.sidebar.write("Testing ratio for Green badge should be greater than Orange which in turn should be greater than Red.")
testing_ratio={}
testing_ratio['Green']=st.sidebar.slider("Select testing ratio for Green badge", min_value=1 , max_value=10 , value=3 , step=1 , format=None , key=None )
testing_ratio['Orange']=st.sidebar.slider("Select testing ratio for Orange badge", min_value=0 , max_value=10 , value=2 , step=1 , format=None , key=None )
testing_ratio['Red']=st.sidebar.slider("Select testing ratio for Red badge", min_value=0 , max_value=10 , value=1 , step=1 , format=None , key=None )

freedom={}
freedom['Green']=st.sidebar.slider("Select freedom for Green badge (recommended 1)", min_value=0.0 , max_value=1.0 , value=1.0 , step=0.1 , format=None , key=None )
freedom['Orange']=st.sidebar.slider("Select freedom for Orange badge", min_value=0.0 , max_value=1.0 , value=0.5 , step=0.1 , format=None , key=None )
freedom['Red']=st.sidebar.slider("Select freedom for Red badge (recommended 0)", min_value=0.0 , max_value=1.0 , value=0.0 , step=0.1 , format=None , key=None )




def simulate(n, frac_poolsize):

	def weighted_sum(d):
		wsum=0
		for state in states:
			for badge in badges:
				wsum+=testing_ratio[badge]*d[state][badge]
		return wsum

	def effective_pool_prevalence(d):
		numerator=0
		for badge in badges:
			numerator+=d['I'][badge]*testing_ratio[badge]
		denominator =weighted_sum(d)
		return numerator/denominator

	def effective_false_positive_rate(d):
		return fp*(1-effective_pool_prevalence(d))**(frac_poolsize*pop-1)+(1-fn)*(1 - (1-effective_pool_prevalence(d))**(frac_poolsize*pop-1))

	def effective_false_negative_rate(d):
		return fn

	def prop_tested(state,badge,d,n,frac_poolsize):
		return r*dt*n*frac_poolsize*testing_ratio[badge]*d[state][badge]/weighted_sum(d)

	def tau(state1,badge1,state2,badge2,d,n,frac_poolsize):
		if state1!=state2:
			print("Error!Tau requires both states to be the same.")
			return None

		state=state1
		badge_value={'Green':3,'Orange':2,'Red':1}
		value=prop_tested(state1,badge1,d,n,frac_poolsize)
		e_fp=effective_false_positive_rate(d)
		e_fn=effective_false_negative_rate(d)

		if state =='S' or state == 'R':
			if badge_value[badge1]>badge_value[badge2]:
				return e_fp*value
			else:
				return (1-e_fp)*value

		if state == 'I':
			if badge_value[badge1]>badge_value[badge2]:
				return (1-e_fn)*value
			else:
				return e_fn*value

	def update_cumulative_quarantine(d):
		cumulative_quarantine=0
		for state in states:
			for badge in badges:
				cumulative_quarantine+=(1-freedom[badge])*d[state][badge]
		return cumulative_quarantine

	def beta_summation(rep_badge,d):
		bsum=0
		for badge in badges:
			bsum+=freedom[rep_badge]*freedom[badge]*d['S'][rep_badge]*d['I'][badge]
		return bsum

	cumulative_quarantine=0
	compartment_value={}
	d_compartment_value={}
	for state in states:
		compartment_value[state]={}
		d_compartment_value[state]={}
		for badge in badges:
			compartment_value[state][badge]=0
			d_compartment_value[state][badge]=0
	compartment_value['I']['Green']=initial_infection
	compartment_value['S']['Green']=(1-initial_infection)

	cumulative_quarantine+=update_cumulative_quarantine(compartment_value)

	for i in range(days):
		d_compartment_value['S']['Green']=dt*(
			-tau('S','Green','S','Orange',compartment_value, n, frac_poolsize)
			+tau('S','Orange','S','Green',compartment_value, n, frac_poolsize)
			+tau('S','Red','S','Orange',compartment_value, n, frac_poolsize)
			-beta*beta_summation('Green',compartment_value)
			+delta*compartment_value['S']['Red']
			)

		d_compartment_value['S']['Orange']=dt*(
			-tau('S','Orange','S','Red',compartment_value, n, frac_poolsize)
			-tau('S','Orange','S','Green',compartment_value, n, frac_poolsize)
			+tau('S','Green','S','Orange',compartment_value, n, frac_poolsize)
			-beta*beta_summation('Orange',compartment_value)
			)

		d_compartment_value['S']['Red']=dt*(
			-tau('S','Red','S','Green',compartment_value, n, frac_poolsize)
			+tau('S','Orange','S','Red',compartment_value, n, frac_poolsize)
			-delta*compartment_value['S']['Red']
			)

		d_compartment_value['I']['Green']=dt*(
			-tau('I','Green','I','Orange',compartment_value, n, frac_poolsize)
			+tau('I','Orange','I','Green',compartment_value, n, frac_poolsize)
			+tau('I','Red','I','Green',compartment_value, n, frac_poolsize)
			+beta*beta_summation('Green',compartment_value)
			-gamma*compartment_value['I']['Green']
			+delta*compartment_value['I']['Red']
			)

		d_compartment_value['I']['Orange']=dt*(
			-tau('I','Orange','I','Red',compartment_value, n, frac_poolsize)
			-tau('I','Orange','I','Green',compartment_value, n, frac_poolsize)
			+tau('I','Green','I','Orange',compartment_value, n, frac_poolsize)
			+beta*beta_summation('Orange',compartment_value)
			-gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['I']['Red']=dt*(
			-tau('I','Red','I','Green',compartment_value, n, frac_poolsize)
			+tau('I','Orange','I','Red',compartment_value, n, frac_poolsize)
			-gamma*compartment_value['I']['Red']
			-delta*compartment_value['I']['Red']
			)

		d_compartment_value['R']['Green']=dt*(
			-tau('R','Green','R','Orange',compartment_value, n, frac_poolsize)
			+tau('R','Orange','R','Green',compartment_value, n, frac_poolsize)
			+tau('R','Red','R','Green',compartment_value, n, frac_poolsize)
			+gamma*compartment_value['I']['Green']
			+delta*compartment_value['R']['Red']
			)

		d_compartment_value['R']['Orange']=dt*(
			-tau('R','Orange','R','Green',compartment_value, n, frac_poolsize)
			-tau('R','Orange','R','Red',compartment_value, n, frac_poolsize)
			+tau('R','Green','R','Orange',compartment_value, n, frac_poolsize)
			+gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['R']['Red']=dt*(
			-tau('R','Red','R','Green',compartment_value, n, frac_poolsize)
			+tau('R','Orange','R','Red',compartment_value, n, frac_poolsize)
			+gamma*compartment_value['I']['Red']
			-delta*compartment_value['R']['Red']
			)

		for state in states:
			for badge in badges:
				compartment_value[state][badge]+=d_compartment_value[state][badge]

		cumulative_quarantine+=update_cumulative_quarantine(compartment_value)


	return compartment_value,cumulative_quarantine


x=np.arange(1, 11, 1)
y=np.arange(0.01, 0.11, 0.01)
X,Y = np.meshgrid(x,y)


data_list={'Infections':np.zeros((len(y),len(x))),'Cumulative Restriction':np.zeros((len(y),len(x)))}

for i,n in enumerate(x):
	for j,f in enumerate(y):
		ts,cq=simulate(n,f)
		data_list['Infections'][j][i]=1-ts['S']['Green']-ts['S']['Orange']-ts['S']['Red']
		data_list['Cumulative Restriction'][j][i]=cq


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, np.array(data_list['Infections']), cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
plt.xlabel("Number of Pools")
plt.ylabel("Fractional Poolsize")
plt.title("Pooling parameters vs total infections")
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.view_init(30, -45)
plt.show()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, np.array(data_list['Cumulative Restriction']), cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
plt.xlabel("Number of Pools")
plt.ylabel("Fractional Poolsize")
plt.title("Pooling parameters vs cumulative restriction")
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.view_init(30, -45)
plt.show()
