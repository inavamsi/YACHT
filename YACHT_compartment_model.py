import streamlit as st
import pandas as pd
import numpy as np

st.write("""
# YACHT protocol
Yet Another Contagion Health Testing (YACHT) protocol. Also known as 'ABCDEFG Protocol' which stands for 'Adaptive Badging Control for a Decentralized Epidemic policy using Feasible Group testing Protocol'.
""")
st.write("This is a simple badging cum pooling mechanism for epidemic control. In this system we show improved testing rate while improving privacy at a cost of quarantining those that are infected. This is carried out through a badge which can be one of 3 colors : Green(G), Orange(O) and Red(R). A person testing negative in a pool is immediately givena green badge which entails complete freedom. An Orange badge obtained on testing ngetaive in a pool partially restricts the agent till the next test. A Red badge begotten on testing negative twice in a row mandates strict quarantine. We see that this simple system can control prevalence effectively without violating privacy of an individual.")
st.write("------------------------------------------------------------------------------------")

pop=1000
dt=0.1
days=2000

cumulative_quarantine=0

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
compartment_ts={}
for state in states:
	compartment_ts[state]={}
	for badge in badges:
		compartment_ts[state][badge]=[]

st.sidebar.write("Badge parameters")
delta=st.sidebar.slider("Rate of faking a badge", min_value=0.0 , max_value=1.0 , value=0.0 , step=0.01 , format=None , key=None)

st.sidebar.write("Testing ratio for Green badge should be greater than Orange which in turn should be greater than Red.")
testing_ratio={}
testing_ratio['Green']=st.sidebar.slider("Select testing ratio for Green badge", min_value=1 , max_value=10 , value=3 , step=1 , format=None , key=None )
testing_ratio['Orange']=st.sidebar.slider("Select testing ratio for Orange badge", min_value=0 , max_value=10 , value=2 , step=1 , format=None , key=None )
testing_ratio['Red']=st.sidebar.slider("Select testing ratio for Red badge", min_value=0 , max_value=10 , value=1 , step=1 , format=None , key=None )

freedom={}
freedom['Green']=st.sidebar.slider("Select freedom for Green badge (recommended 1)", min_value=0.0 , max_value=1.0 , value=1.0 , step=0.1 , format=None , key=None )
freedom['Orange']=st.sidebar.slider("Select freedom for Orange badge", min_value=0.0 , max_value=1.0 , value=0.5 , step=0.1 , format=None , key=None )
freedom['Red']=st.sidebar.slider("Select freedom for Red badge (recommended 0)", min_value=0.0 , max_value=1.0 , value=0.0 , step=0.1 , format=None , key=None )

dynamic_pooling=st.sidebar.checkbox("Dynamic Pooling",value=False)

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

def effective_false_positive_rate(d,frac_poolsize):
	return fp*(1-effective_pool_prevalence(d))**(frac_poolsize*pop-1)+(1-fn)*(1 - (1-effective_pool_prevalence(d))**(frac_poolsize*pop-1))

def effective_false_negative_rate(d):
	return fn

def prop_tested(state,badge,d,frac_poolsize):
	return r*dt*n*frac_poolsize*testing_ratio[badge]*d[state][badge]/weighted_sum(d)

def tau(state1,badge1,state2,badge2,d,frac_poolsize):
	if state1!=state2:
		print("Error!Tau requires both states to be the same.")
		return None

	state=state1
	badge_value={'Green':3,'Orange':2,'Red':1}
	value=prop_tested(state1,badge1,d,frac_poolsize)
	e_fp=effective_false_positive_rate(d,frac_poolsize)
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

def update_compartment_ts(d):
	for state in states:
		for badge in badges:
			compartment_ts[state][badge].append(d[state][badge])

def update_cumulative_quarantine(d):
	global cumulative_quarantine
	for state in states:
		for badge in badges:
			cumulative_quarantine+=(1-freedom[badge])*d[state][badge]


def beta_summation(rep_badge,compartment_value):
	bsum=0
	for badge in badges:
		bsum+=freedom[rep_badge]*freedom[badge]*compartment_value['S'][rep_badge]*compartment_value['I'][badge]
	return bsum

def simulate():
	global frac_poolsize
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


	for i in range(days):


		d_compartment_value['S']['Green']=dt*(
			-tau('S','Green','S','Orange',compartment_value,frac_poolsize)
			+tau('S','Orange','S','Green',compartment_value,frac_poolsize)
			+tau('S','Red','S','Orange',compartment_value,frac_poolsize)
			-beta*beta_summation('Green',compartment_value)
			+delta*compartment_value['S']['Red']
			)

		d_compartment_value['S']['Orange']=dt*(
			-tau('S','Orange','S','Red',compartment_value,frac_poolsize)
			-tau('S','Orange','S','Green',compartment_value,frac_poolsize)
			+tau('S','Green','S','Orange',compartment_value,frac_poolsize)
			-beta*beta_summation('Orange',compartment_value)
			)

		d_compartment_value['S']['Red']=dt*(
			-tau('S','Red','S','Green',compartment_value,frac_poolsize)
			+tau('S','Orange','S','Red',compartment_value,frac_poolsize)
			-delta*compartment_value['S']['Red']
			)

		d_compartment_value['I']['Green']=dt*(
			-tau('I','Green','I','Orange',compartment_value,frac_poolsize)
			+tau('I','Orange','I','Green',compartment_value,frac_poolsize)
			+tau('I','Red','I','Green',compartment_value,frac_poolsize)
			+beta*beta_summation('Green',compartment_value)
			-gamma*compartment_value['I']['Green']
			+delta*compartment_value['I']['Red']
			)

		d_compartment_value['I']['Orange']=dt*(
			-tau('I','Orange','I','Red',compartment_value,frac_poolsize)
			-tau('I','Orange','I','Green',compartment_value,frac_poolsize)
			+tau('I','Green','I','Orange',compartment_value,frac_poolsize)
			+beta*beta_summation('Orange',compartment_value)
			-gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['I']['Red']=dt*(
			-tau('I','Red','I','Green',compartment_value,frac_poolsize)
			+tau('I','Orange','I','Red',compartment_value,frac_poolsize)
			-gamma*compartment_value['I']['Red']
			-delta*compartment_value['I']['Red']
			)

		d_compartment_value['R']['Green']=dt*(
			-tau('R','Green','R','Orange',compartment_value,frac_poolsize)
			+tau('R','Orange','R','Green',compartment_value,frac_poolsize)
			+tau('R','Red','R','Green',compartment_value,frac_poolsize)
			+gamma*compartment_value['I']['Green']
			+delta*compartment_value['R']['Red']
			)

		d_compartment_value['R']['Orange']=dt*(
			-tau('R','Orange','R','Green',compartment_value,frac_poolsize)
			-tau('R','Orange','R','Red',compartment_value,frac_poolsize)
			+tau('R','Green','R','Orange',compartment_value,frac_poolsize)
			+gamma*compartment_value['I']['Orange']
			)

		d_compartment_value['R']['Red']=dt*(
			-tau('R','Red','R','Green',compartment_value,frac_poolsize)
			+tau('R','Orange','R','Red',compartment_value,frac_poolsize)
			+gamma*compartment_value['I']['Red']
			-delta*compartment_value['R']['Red']
			)

		if dynamic_pooling:
			if d_compartment_value['I']['Green']+d_compartment_value['I']['Orange']+d_compartment_value['I']['Red'] > 0:
				frac_poolsize*=1.001
			else:
				frac_poolsize/=1.001
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

def plot_SIR(dict_ts):
	values=[]
	keys=[]

	for i in range(len(dict_ts['S']['Green'])):
		values.append([])
		for state in states:
			values[-1].append(0)

	for i in range(len(dict_ts['S']['Green'])):
		for indx,state in enumerate(states):
			for badge in badges:
				values[i][indx]+=dict_ts[state][badge][i]

	for state in states:
			keys.append(state)


	chart_data = pd.DataFrame(
	values,
	columns=keys)
	st.line_chart(chart_data)





simulate()
st.write("This is the timeseries plot for the 9 compartments.")
plot(compartment_ts)
st.write("This is the timeseries plot for the underlying disease states {S,I,R}.")
plot_SIR(compartment_ts)



st.write("------------------------------------------------------------------------------------")

st.header("Cost Function")
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
a=st.slider("Select scaling parameter or cost for setting up a pool", min_value=0 , max_value=10000 , value=0 , step=1 , format=None , key=None )
st.write("Total infected percentage is "+str(total_infected_percentage)+". Here 0 means no infection while 1 means no one is susceptible.")
b=st.slider("Select scaling parameter for 'Total Infected percentage'", min_value=0 , max_value=10000 , value=0 , step=1 , format=None , key=None )
st.write("Cumulative Quarantine is "+str(cumulative_quarantine_percentage)+ ". Here 0 means nobody is quarantined at all, while 1 means everyone is fully quarantined for the entire duration.")
c=st.slider("Select scaling parameter for 'Cumulative Quarantine'", min_value=0 , max_value=10000 , value=0 , step=1 , format=None , key=None )

cost = a*number_of_pools + b*total_infected_percentage + c*cumulative_quarantine_percentage
st.write("The cost is "+ str(cost))
