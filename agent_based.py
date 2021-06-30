import random
import copy
from statistics import mean
import matplotlib.pyplot as plt
import matplotlib

infected_G2R=None
not_infected_R2G=None
false_positives=None
false_negatives=None

n=1000
poolsize=None

no_pools=30
beta=0.0003
gamma=0.08
infected_percentage=0.01
days=30

states=['Susceptible','Infected','Recovered']
badges=['Green','Orange','Red']
freedom={}
freedom['Green']=1
freedom['Orange']=0.5
freedom['Red']=0

testing_ratio={}
testing_ratio['Green']=1
testing_ratio['Orange']=0.4
testing_ratio['Red']=0.2

ts=None

class Agent():
    def __init__(self,infected_percentage):
        if random.random()<infected_percentage:
            self.state='Infected'
        else:
            self.state='Susceptible'

        self.badge='Green'
        self.results=[]

        self.inf_G2R=None
        self.sus_R2G=None


    def update_state(self,new_state):
        if self.state==new_state:
            return
        self.state=new_state


    def update_badge(self,new_badge):
        if self.badge==new_badge:
            return

        if self.state=='Infected' and self.inf_G2R is not None:
            if new_badge=='Red':
                infected_G2R.append(self.inf_G2R)

        if self.state !='Infected' and self.sus_R2G is not None:
            if new_badge=='Green':
                not_infected_R2G.append(self.sus_R2G)

        self.badge=new_badge

        if self.badge=='Green' and self.state=='Infected':
            self.inf_G2R=0
        if self.badge=='Red' and self.state !='Infected':
            self.sus_R2G=0


    def new_day(self):
        self.results=[]

        if self.inf_G2R is not None and self.badge=='Green' and self.state=='Infected':
            self.inf_G2R+=1
        if self.sus_R2G is not None and self.badge=='Red' and self.state !='Infected':
            self.sus_R2G+=1

    def next_badge(self):
        final_result='Negative'
        for result in self.results:
            if result=='Positive':
                final_result='Positive'
                break
        if final_result=='Negative':
            self.update_badge('Green')
        else:
            if self.badge=='Green':
                self.update_badge('Orange')
            else:
                self.update_badge('Red')

    def update_false_rates(self):
        if self.results==[]:
            return
        final_result='Negative'
        for result in self.results:
            if result=='Positive':
                final_result='Positive'
                break
        if final_result=='Negative' and self.state=='Infected':
            global false_negatives
            false_negatives+=1

        if final_result=='Positive' and self.state!='Infected':
            global false_positives
            false_positives+=1


    def add_result(self,result):
        self.results.append(result)


class Pool():
    def __init__(self):
        self.fn=0
        self.fp=0
        self.agents=[]

    def new_day(self):
        self.agents=[]

    def add_agent(self,poolsize,agent):
        if len(self.agents)>=poolsize:
            print('Error! Pool limit is crossed.')
            return
        self.agents.append(agent)


    def release_results(self):
        result='Negative'
        for agent in self.agents:
            if agent.state=='Infected':
                result='Positive'
                break

        if result=='Positive':
            if random.random()<self.fn:
                result='Negative'
        else:
            if random.random()<self.fp:
                result='Positive'

        for agent in self.agents:
            agent.add_result(result)

def beta_summation(agents):
    val=0
    for agent in agents:
        if agent.state =='Infected':
            val+=freedom[agent.badge]
    return val*beta

def state_change(agents):
    inf_prob=beta_summation(agents)
    rec_prob=gamma

    next_states=[]
    for agent in agents:
        if agent.state =='Susceptible':
            if random.random()<inf_prob*freedom[agent.badge]:
                next_states.append('Infected')
                continue
        elif agent.state=='Infected':
            if random.random()<rec_prob:
                next_states.append('Recovered')
                continue
        next_states.append(agent.state)

    for indx,agent in enumerate(agents):
        agent.update_state(next_states[indx])

def sample_agent(agents):
    while(True):
        r=random.random()
        chosen=random.choice(agents)
        if r <testing_ratio[chosen.badge]:
            return chosen


def test(pools,agents,poolsize):

    testing_agents=[]
    for i in range(len(pools)*poolsize):
        testing_agents.append(sample_agent(agents))

        j=0
    for pool in pools:
        for i in range(poolsize):
            pool.add_agent(poolsize,testing_agents[j])
            j+=1

    for pool in pools:
        pool.release_results()

def simulate_day( pools, agents):
    for badge in badges:
        for state in states:
            ts[badge][state].append(0)

    for agent in agents:
        agent.new_day()
    for pool in pools:
        pool.new_day()

    test(pools,agents,poolsize)
    for agent in agents:
        #agent.update_false_rates()
        agent.next_badge()
    state_change(agents)

    for agent in agents:
        ts[agent.badge][agent.state][-1]+=1

def simulate(days,pools,agents):
    for day in range(days):
        simulate_day(pools,agents)

def main():
    agents=[]
    pools=[]

    global infected_G2R
    infected_G2R=[]
    global not_infected_R2G
    not_infected_R2G=[]
    global false_positives
    false_positives=0
    global false_negatives
    false_negatives=0

    global ts
    ts={}
    for badge in badges:
        ts[badge]={}
        for state in states:
            ts[badge][state]=[]

    for i in range(n):
        agents.append(Agent(infected_percentage))
    for i in range(no_pools):
        pools.append(Pool())

    simulate(days,pools,agents)

    '''for state in states:
        temp=0
        for badge in badges:
            temp+=ts[badge][state][-1]
        print(state+' : '+str(temp))
        '''
    if infected_G2R==[]:
        avg_infected_G2R=0
    else:
        avg_infected_G2R=mean(infected_G2R)

    if not_infected_R2G==[]:
        avg_not_infected_R2G=0
    else:
        avg_not_infected_R2G=mean(not_infected_R2G)

    infections=n - ts['Green']['Susceptible'][-1]-ts['Orange']['Susceptible'][-1]-ts['Red']['Susceptible'][-1]

    return avg_infected_G2R,avg_not_infected_R2G,false_positives,false_negatives,infections



def worlds(no):
    infected_G2R=0
    not_infected_R2G=0
    false_positives=0
    false_negatives=0
    infections=0

    for i in range(no):
        x,y,a,b,c=main()
        infected_G2R+=x
        not_infected_R2G+=y
        false_positives+=a
        false_negatives+=b
        infections+=c

    return round(infected_G2R/no,2), round((not_infected_R2G)/no,2), round(infections/(n*no),2)

    print('Average time for an Infected with Green Badge to get a Red Badge : ' +str((infected_G2R)/no))
    print('Average time for a Non-Infected with Red Badge to get a Green Badge : '+str((not_infected_R2G)/no))
    #print('Total False Positives : '+str((false_positives)/no))
    #print('Total False Negatives : '+str((false_negatives)/no))
    print('Proportion of population Infected : '+str(infections/(n*no)))

def histogram():

    random.seed(42)

    global poolsize
    inf_G2R_list=[]
    not_inf_R2G_list=[]
    infection_proportion_list=[]

    poolsize_list=[0,1,2,3,4,5,10,15,20,25,30]

    for poolsize in poolsize_list:
        a,b,c=worlds(50)
        inf_G2R_list.append(a)
        not_inf_R2G_list.append(b)
        infection_proportion_list.append(c)

    plt.plot(poolsize_list,inf_G2R_list)
    plt.title('Poolsize vs Average time for an Infected \n agent with a Green Badge to get Red Badge')
    plt.xlabel('Poolsize')
    plt.ylabel('Average Time in Days')
    plt.show()

    plt.plot(poolsize_list,not_inf_R2G_list)
    plt.title('Poolsize vs Average time for a Non-Infected \n agent with a Red Badge to get Green Badge')
    plt.xlabel('Poolsize')
    plt.ylabel('Average Time in Days')
    plt.show()

    plt.plot(poolsize_list,infection_proportion_list)
    plt.title('Poolsize vs Total Proportion of population \n that gets infected in a duration of '+str(days)+' days')
    plt.xlabel('Poolsize')
    plt.ylabel('Proportion of population that got infected')
    plt.show()

    print(inf_G2R_list, not_inf_R2G_list, infection_proportion_list)

histogram()
