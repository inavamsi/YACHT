import streamlit as st
import pandas as pd
import numpy as np
import random

st.write("""
# YACHT using a compartment model.
A badging cum pooling mechanism for effective epidemic control
""")
st.write("------------------------------------------------------------------------------------")

pop=100
dt=0.1
days=2000

cumulative_quarantine=0

st.write("Epidemic parameters")
beta=st.slider("Rate of infection", min_value=0.0 , max_value=1.0 , value=0.11 , step=0.01 , format=None , key=None )
gamma=st.slider("Rate of recovery", min_value=0.0 , max_value=1.0 , value=0.05 , step=0.01 , format=None , key=None )
initial_infection=st.slider("Select inital infection proportion", min_value=0.0 , max_value=1.0 , value=0.01 , step=1/pop , format=None , key=None )
st.write("------------------------------------------------------------------------------------")

st.write("YACHT parameters")
n=st.slider("Select number of pools", min_value=0 , max_value=100 , value=9 , step=1 , format=None , key=None )
frac_poolsize = st.slider("Select Poolsize/Population", min_value=0.0 , max_value=1.0/n , value=None , step=1/pop , format=None , key=None )
r = st.slider("Select rate of testing per timestep", min_value=0.0 , max_value=1.0 , value=0.5 , step=0.01 , format=None , key=None)
poolsize=frac_poolsize*pop
st.write("------------------------------------------------------------------------------------")

badges=['Green','Orange','Red']
states=['S','I','R']
compartment_ts={}
for state in states:
	compartment_ts[state]={}
	for badge in badges:
		compartment_ts[state][badge]=[]

st.write("Badge parameters")
st.write("Testing ratio for Green badge should be greater than Orange which in turn should be greater than Red.")
testing_ratio={}
testing_ratio['Green']=st.slider("Select testing ratio for Green badge", min_value=1 , max_value=10 , value=3 , step=1 , format=None , key=None )
testing_ratio['Orange']=st.slider("Select testing ratio for Orange badge", min_value=0 , max_value=10 , value=2 , step=1 , format=None , key=None )
testing_ratio['Red']=st.slider("Select testing ratio for Red badge", min_value=0 , max_value=10 , value=1 , step=1 , format=None , key=None )

freedom={}
freedom['Green']=1
freedom['Orange']=st.slider("Select freedom for Orange badge, given that it is 1 for Green and 0 for Red", min_value=0.0 , max_value=1.0 , value=0.5 , step=0.1 , format=None , key=None )
freedom['Red']=0

def weighted_sum(d):
	wsum=0
	for state in states:
		for badge in badges:
			wsum+=testing_ratio[badge]*d[state][badge]
	return wsum

def effective_pool_prevalence(d):
	numerator=0
	for badge in badges:
		numerator+=d['I'][badge]
	denominator =weighted_sum(d)
	return numerator/denominator

def effective_false_positive_rate(d):
	return 1 - (1-effective_pool_prevalence(d))**(frac_poolsize*pop-1)

def prop_tested(state,badge,d):
	return r*dt*n*frac_poolsize*testing_ratio[badge]*d[state][badge]/weighted_sum(d)

def tau(state1,badge1,state2,badge2,d):
	value=prop_tested(state1,badge1,d)
	if state1 =='I':
		return value
	elif badge2=='Green':
		return value*(1-effective_false_positive_rate(d))
	else :
		return value*effective_false_positive_rate(d)

def update_compartment_ts(d):
	for state in states:
		for badge in badges:
			compartment_ts[state][badge].append(d[state][badge])

def update_cumulative_quarantine(d):
	global cumulative_quarantine
	for state in states:
		for badge in badges:
			cumulative_quarantine+=(1-freedom[badge])*d[state][badge]
	

def beta_summation(rep_badge,d):
	bsum=0
	for badge in badges:
		bsum+=freedom[rep_badge]*freedom[badge]*d['S'][rep_badge]*d['I'][badge]
	return bsum

def simulate():
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

	update_compartment_ts(compartment_value)
	update_cumulative_quarantine(compartment_value)

	#st.write(tau('S','Green','S','Orange',compartment_value))
	#st.write(tau('S','Orange','S','Green',compartment_value))
	#st.write(tau('S','Red','S','Orange',compartment_value))
	#st.write(beta*beta_summation('Green',compartment_value))

	for i in range(days):
		d_compartment_value['S']['Green']=dt*(
			-tau('S','Green','S','Orange',compartment_value)
			+tau('S','Orange','S','Green',compartment_value)
			+tau('S','Red','S','Orange',compartment_value)
			-beta*beta_summation('Green',compartment_value)
			)

		d_compartment_value['S']['Orange']=dt*(
			-tau('S','Orange','S','Red',compartment_value)
			-tau('S','Orange','S','Green',compartment_value)
			+tau('S','Green','S','Orange',compartment_value)
			-beta*beta_summation('Orange',compartment_value)
			)

		d_compartment_value['S']['Red']=dt*(
			-tau('S','Red','S','Green',compartment_value)
			+tau('S','Orange','S','Red',compartment_value)
			)

		d_compartment_value['I']['Green']=dt*(
			-tau('I','Green','I','Orange',compartment_value)
			+beta*beta_summation('Green',compartment_value)
			-gamma*compartment_value['I']['Green']
			)

		d_compartment_value['I']['Orange']=dt*(
			-tau('I','Orange','I','Red',compartment_value)
			+tau('I','Green','I','Orange',compartment_value)
			+beta*beta_summation('Orange',compartment_value)
			-gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['I']['Red']=dt*(
			tau('I','Orange','I','Red',compartment_value)
			-gamma*compartment_value['I']['Red']
			)

		d_compartment_value['R']['Green']=dt*(
			-tau('R','Green','R','Orange',compartment_value)
			+tau('R','Orange','R','Green',compartment_value)
			+tau('R','Red','R','Green',compartment_value)
			+gamma*compartment_value['I']['Green']
			)

		d_compartment_value['R']['Orange']=dt*(
			-tau('R','Orange','R','Green',compartment_value)
			-tau('R','Orange','R','Red',compartment_value)
			+tau('R','Green','R','Orange',compartment_value)
			+gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['R']['Red']=dt*(
			-tau('R','Red','R','Green',compartment_value)
			+tau('R','Orange','R','Red',compartment_value)
			+gamma*compartment_value['I']['Red']
			)

		for state in states:
			for badge in badges:
				compartment_value[state][badge]+=d_compartment_value[state][badge]

		update_compartment_ts(compartment_value)
		update_cumulative_quarantine(compartment_value)


def plot(dict_ts):
	values=[]
	keys=[]

	for i in range(len(dict_ts['S']['Green'])):
		values.append([])
		for state in states:
			for badge in badges:
				values[-1].append(dict_ts[state][badge][i])
	
	for state in states:
		for badge in badges:
			keys.append(state+'_'+badge)


	chart_data = pd.DataFrame(
	values,
	columns=keys)
	st.line_chart(chart_data)
st.write("------------------------------------------------------------------------------------")
st.write("This is the simulation of for the above parameters.")

simulate()
plot(compartment_ts)

st.write("------------------------------------------------------------------------------------")

st.write("Goal is to minimise Cost function constrained by a maximum value for prevalence :")
st.write("Cost function =  a(Number of Pools) + b(Pandemic Size) + c(Economic Cost) ")#+ d(Logistic Cost)") 
st.write(" -- Number of pools refers to the total number of pooling station.")
st.write(" -- Pandmeic Size refers to percentage of total that got infected.")
st.write(" -- Economic cost refers to the economic loss due to quarantine and restrictions. This percentage represents the total time lost.")

total_infected_percentage=1-(compartment_ts['S']['Green'][-1]+compartment_ts['S']['Orange'][-1]+compartment_ts['S']['Red'][-1])
number_of_pools=n
cumulative_quarantine_percentage=cumulative_quarantine/days
st.write("------------------------------------------------------------------------------------")

st.write("The number of pools is "+str(n)+".")
a=st.slider("Select scaling parameter for 'Number of Pools'")
st.write("Total infected percentage is "+str(total_infected_percentage)+". Here 0 means no infection while 1 means no one is susceptible.")
b=st.slider("Select scaling parameter for 'Total Infected percentage'")
st.write("Cumulative Quarantine is "+str(cumulative_quarantine_percentage)+ ". Here 0 means nobody is quarantined at all, while 1 means everyone is fully quarantined for the entire duration.")
c=st.slider("Select scaling parameter for 'Cumulative Quarantine'")


cost = a*number_of_pools + b*total_infected_percentage + c*cumulative_quarantine_percentage
st.write("The cost is "+ str(cost))