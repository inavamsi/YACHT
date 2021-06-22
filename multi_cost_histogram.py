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

badges=['Green','Orange','Red']
states=['S','I','R']

initial_infection=0.01
delta=0


def simulate(n, frac_poolsize, freedom, testing_ratio,fn,fp,beta,gamma,r,cost):

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

	inf=1-compartment_value['S']['Green']-compartment_value['S']['Orange']-compartment_value['S']['Red']
	return n*cost[0]+pop*inf*cost[1]+cumulative_quarantine*cost[2]

beta=0.1
gamma=0.05

fn=0.1
fp=0.1
n=30
frac_poolsize = 0.01
r = 0.5

testing_ratio={}
testing_ratio['Green']=5
testing_ratio['Orange']=2
testing_ratio['Red']=1

freedom={}
freedom['Green']=1
freedom['Orange']=0.5
freedom['Red']=0

costs=[(10,10,1),(30,5,1),(20,8,1)]

values={}
l={}
values['n']=[0,20,40,60,80,100]
for cost in costs:
	l[cost]={}

for cost in costs:
	l[cost]['n']=[]
	for frac_poolsize in [0.001,0.005,0.01]:
		temp=[]
		for n in values['n']:
			res=simulate(n, frac_poolsize, freedom, testing_ratio,fn,fp,beta,gamma,r,cost)
			temp.append(res)
		l[cost]['n'].append(temp)

def hist(ax,l,values,key,cost,title,xlab,ylab):
	x_axis=values[key]
	x = np.arange(len(x_axis))  # the label locations
	width = 0.2  # the width of the bars

	rects1 = ax.bar(x-width, l[cost][key][0], width, label='Poolsize 1')
	rects2 = ax.bar(x, l[cost][key][1], width, label='Poolsize 5')
	rects3 = ax.bar(x+width, l[cost][key][2], width, label='Poolsize 10')
	#rects4 = ax.bar(x+2*width, l[frac_poolsize][key][3], width, label=costs[3])

	#ax.set(ylim=[3000, 7000])
	ax.set_ylabel(ylab)
	ax.set_xlabel(xlab)
	ax.set_title(title)
	ax.set_xticks(x)
	ax.set_xticklabels(x_axis)
	ax.legend()


fig, ax0 = plt.subplots()
#fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)
#fig, (ax0, ax1,ax2) = plt.subplots(nrows=1, ncols=3)
hist(ax0,l,values,'n',costs[0],'Optimal number of pools with cost structure '+str(costs[0]),'Number of Pools (n)','Total Cost')
#hist(ax1,l,values,'n',costs[1],'Optimal number of pools with cost structure '+str(costs[1]),'Number of Pools (n)','Total Cost')
#hist(ax2,l,values,'n',costs[2],'Optimal number of pools with cost structure '+str(costs[2]),'Number of Pools (n)','Total Cost')

#hist(,l,0.1,values,'n',costs,',,,','..','.')
#hist(,l,0.5,values,'n',costs,',,,','..','.')
fig.tight_layout()
plt.show()

fig, ax1 = plt.subplots()
#fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)
#fig, (ax0, ax1,ax2) = plt.subplots(nrows=1, ncols=3)
#hist(ax0,l,values,'n',costs[0],'Optimal number of pools with cost structure '+str(costs[0]),'Number of Pools (n)','Total Cost')
hist(ax1,l,values,'n',costs[1],'Optimal number of pools with cost structure '+str(costs[1]),'Number of Pools (n)','Total Cost')
#hist(ax2,l,values,'n',costs[2],'Optimal number of pools with cost structure '+str(costs[2]),'Number of Pools (n)','Total Cost')

#hist(,l,0.1,values,'n',costs,',,,','..','.')
#hist(,l,0.5,values,'n',costs,',,,','..','.')
fig.tight_layout()
plt.show()
fig, ax2 = plt.subplots()
#fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)
#fig, (ax0, ax1,ax2) = plt.subplots(nrows=1, ncols=3)
#hist(ax0,l,values,'n',costs[0],'Optimal number of pools with cost structure '+str(costs[0]),'Number of Pools (n)','Total Cost')
#hist(ax1,l,values,'n',costs[1],'Optimal number of pools with cost structure '+str(costs[1]),'Number of Pools (n)','Total Cost')
hist(ax2,l,values,'n',costs[2],'Optimal number of pools with cost structure '+str(costs[2]),'Number of Pools (n)','Total Cost')

#hist(,l,0.1,values,'n',costs,',,,','..','.')
#hist(,l,0.5,values,'n',costs,',,,','..','.')
fig.tight_layout()
plt.show()
