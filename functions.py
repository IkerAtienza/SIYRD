from scipy.integrate import ode, solve_ivp
import matplotlib.patches as mpatches ## Used by plot_population_pyramid function
import matplotlib.ticker as ticker ## To set labels as millions
import matplotlib.pyplot as plt
from itertools import islice ## This is used to split one list in parts of given length
import pandas as pd
import numpy as np
import os # operative system commands
import glob # takes all files from directory to be deleted
import csv
import sys

# NEW! IFR and contacts plot
def plot_ifr_contacts(x,country):
    # Ifr file is opened and information is saved
    ifr_file = 'ifr_'+str(x)+'.csv'
    reader2 = csv.reader(open('infiles/'+country+'/'+ifr_file, "r"), delimiter=",")

    y2 = list(reader2)
    ifr = np.array(y2).astype('float')

    # Contact matrix file is opened
    matrix_file = 'new_M_'+str(x)+'.csv'
    reader = csv.reader(open(('infiles/'+country+'/'+matrix_file), "r"), delimiter=",")
    y = list(reader)
    m_matrix = np.array(y).astype("float")
    # Contacts are calculated for each group as the row-elements sum. These values account for the contacts within the group and with the other group
    contacts =  [sum(m) for m in m_matrix]
    # Ifr is calculated as percentage
    ifr = [fr*100 for fr in ifr]
    ifr = ifr[0]

    # labels are created
    labels = []
    labels.append('G1 ('+r'$<$'+str(x)+')')
    labels.append('G2 ('+r'$\geq$'+str(x)+')')

    # IFR and contacts are plotted as a bar plot.
    fig, ax = plt.subplots(figsize=(8,8))
    ax2=ax.twinx()
    # x = np.arange(0,len(labels))  # label locations
    x = np.arange(0,len(labels))  # label locations
    width = 0.2  # width of the bars
    rects1 = ax.bar(x - width/1.8, contacts, width, label='Contacts',color='steelblue',edgecolor='#004c98',linewidth=2)
    rects2 = ax2.bar(x + width/1.8, ifr, width, label='ifr',color='maroon',edgecolor='#4c0026',linewidth=2)
    # rects1 = ax.bar(x , N_gr, width, color='steelblue', edgecolor='#004c98',linewidth=2)

    # ax.set_ylim(0,13)
    y1 = np.arange(0,16,2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize = 20)
    ax.set_yticks(y1)
    ax.set_yticklabels(y1, fontsize = 20, color='steelblue')
    ax.set_ylabel('Number of contacts',color='steelblue',fontsize=20)
    # ax2 information
    y2 = np.arange(0,27.5,2.5)
    # ax2.set_ylim(0,20)
    ax2.set_yticks(y2)
    ax2.set_yticklabels(y2, fontsize = 20, color='maroon')
    ax2.set_ylabel('% Infection Fatality Risk',color='maroon',fontsize=20)
    # ax2.spines['left'].set_color('steelblue')
    # ax2.spines['right'].set_color('maroon')
    plt.show()


# Groups density plot
def plot_groups_density(x,country):
    # print(country)
    # Age distribution input file is opened. Data are saved in two arrays, ages and amounts
    ages, amounts = [], []
    # with open('infiles/Spain/Spain_country_level_age_distribution_85.csv','r') as agefile:
    with open('infiles/'+country+'/'+country+'-2020_gr10.csv','r') as agefile:
        rows = csv.reader(agefile, delimiter=',')
        for row in rows:
            ages.append(row[0])
            amounts.append(float(row[1]))

    # Population density of each group is calculated given the age limit
    num_groups = 2 # two age grups
    age_limits = [] # age limits will always contain [0, limit1,len(ages)]
    age_limits.append(0)
    age_limits.append(int(x/10)) # we add the /10 because we have 9 groups an not 85 as before
    age_limits.append(len(ages))

    # Total population is splitted in two sublists. Then densitiy for each group is calculated.
    # groups = 2    age_limits = [0,80,85]      length_to_split = [80,5]     groups = [[80 elements],[5 elements]]
    length_to_split = []
    for i in range(1,len(age_limits)):
        length_to_split.append(age_limits[i]-age_limits[i-1])
    amounts_iter = iter(amounts) # iter object is needed
    N_gr = [list(islice(amounts_iter, elem)) for elem in length_to_split]
    # by adding all elements in each sublist we get the # of individuals of each group defined
    N_gr = [sum(element)/1E6 for element in N_gr]
    total_population = sum(N_gr)

    if int(total_population) in range(0,10):
        y_step = 2
    elif int(total_population) in range(10,50):
        y_step = 5
    elif int(total_population) in range(50,100):
        y_step = 10
    elif int(total_population) in range(100,300):
        y_step = 20
    elif int(total_population) in range(300,500):
        y_step = 50
    elif int(total_population) > 500:
        y_step = 100

    # labels are created
    labels = []
    labels.append('G1 ('+r'$<$'+str(x)+')')
    labels.append('G2 ('+r'$\geq$'+str(x)+')')


    # Population for each group is plotted as a barplot
    fig, ax = plt.subplots(figsize=(8,8))
    x = np.arange(0,len(labels))  # the label locations
    y = np.arange(0,total_population+y_step, y_step)
    width = 0.25  # the width of the bars
    rects1 = ax.bar(x , N_gr, width, color='steelblue', edgecolor='#004c98',linewidth=2)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize = 20)
    ax.set_yticks(y)
    ax.set_ylabel('Population in millions',fontsize=20)
    ax.tick_params(axis='y', labelsize=18)

    # plt.title('Population Density',fontsize=20)
    plt.show()

# Population pyramid plot
def plot_population_pyramid(country):
    # print(country)
    file = 'infiles/'+str(country)+'/'+str(country)+'-2020.csv'
    df = pd.read_csv(file)
    ages = df['Age']
    males = df['M']
    females = df['F']
    # From both peaks of males and females, the highest value is selected for being the limit of x-axis
    max_males = round(max(males)/1E6)
    max_females = round(max(females)/1E6)
    x_limit = max(max_males,max_females)
    # according to the x_limit the step for the x-axis is changed
    if x_limit in range(0,5):
        if x_limit < 1: # population peak below 1M. X-axis limit is 500k == 0.5M
            step = 0.25
            x_limit = 0.5
        else:
            step = 1
    elif x_limit in range(5,13):
        x_limit = 2 * round((x_limit+1)/2) # x_limit is rounded to the nearest higher 2 multiple
        step = 2
    elif x_limit in range(13,20):
        x_limit = 5 * round(x_limit/5) # x_limit is rounded to the nearest 5 multiple
        step = 5
    elif x_limit in range(20,100):
        x_limit = round(x_limit,-1) # x_limit is rounded to the nearest 10 multiple
        step = 10
    # x-axis values are created
    x = np.arange(-x_limit,x_limit+step,step)

    # Population pyramid is plotted
    fig, ax = plt.subplots(figsize=(9,12))
    ax.barh(ages, -males, color = 'steelblue', label = 'Male', edgecolor='#004c98',linewidth=0.5)
    ax.barh(ages, females, color = 'maroon', edgecolor='#4c0026',linewidth=0.5)
    # ax.set_title(country,fontsize=30)
    ax.set_xticks(x*1E6)
    ax.set_xticklabels([str(abs(v)) for v in x])
    ax.set_xlabel('Population in millions',fontsize=20)
    # two legends are displayed
    leg1 = ax.legend(loc='upper left', fontsize=18)
    red_patch = mpatches.Patch(color='maroon',label='Female')
    ax.legend(handles=[red_patch], loc = 'upper right', fontsize=18)
    ax.add_artist(leg1)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    plt.show()

# Save simulation results in output file
def saveResults(filename,t,sim_results):
    # print("Saving results for file %s" %(filename))
    fout = open(('outfiles/'+filename),"w+")
    fout.write("Time\tS1\tS2\tI1\tI2\tY1\tY2\tR1\tR2\tD1\tD2\tV1\tV2\tNewI\tNewY\n")

    for s in range(len(sim_results["t"])):
        fout.write("%f\t" % (sim_results["t"][s]))
        for c in range(len(sim_results["y"])):
            fout.write("%f\t" % (sim_results["y"][c][s]))
        fout.write('\n')
    fout.close()

# SIYRD no vaccination function
def noVacc(times,init,parms):
    N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M = parms
    S1, S2, I1, I2, Y1, Y2, R1 ,R2, D1, D2, V1, V2, new_I, new_Y =  init
    # ODEs
    dS1dt = -(bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - v[0]
    dS2dt = -(bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2 - v[1]
    dI1dt = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - r1*I1 - mui1*I1
    dI2dt = (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2 - r2*I2 - mui2*I2
    dY1dt = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 - r1*Y1 - muy1*Y1
    dY2dt = (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 - r2*Y2 - muy2*Y2
    dR1dt = r1*I1 + r1*Y1 - (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + v[0]
    dR2dt = r2*I2 + r2*Y2 - (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 + v[1]
    dD1dt = mui1*I1 + muy1*Y1
    dD2dt = mui2*I2 + muy2*Y2
    dV1dt = v[0]
    dV2dt = v[1]
    new_inf = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 + (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2
    new_reinf = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2

    return [dS1dt, dS2dt, dI1dt, dI2dt, dY1dt, dY2dt, dR1dt, dR2dt, dD1dt, dD2dt, dV1dt, dV2dt, new_inf, new_reinf]

# SIYRD simultaneous vaccination function
def simultVacc(times, init, parms):
    N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M = parms
    S1, S2, I1, I2, Y1, Y2, R1 ,R2, D1, D2, V1, V2, new_I, new_Y =  init
    if S1 < (0.3*N_gr[0]): # if 70% of G1 is recovered (S1<30%), vaccination is stopped for this group
        v[0] = 0
    if S2 < (0.3*N_gr[1]): # if 70% of G2 is recovered (S2<30%), vaccination is stopped for this group
        v[1] = 0
    # ODEs
    dS1dt = -(bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - v[0]
    dS2dt = -(bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2 - v[1]
    dI1dt = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - r1*I1 - mui1*I1
    dI2dt = (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2- r2*I2 - mui2*I2
    dY1dt = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 - r1*Y1 - muy1*Y1
    dY2dt = (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 - r2*Y2 - muy2*Y2
    dR1dt = r1*I1 + r1*Y1 - (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + v[0]
    dR2dt = r2*I2 + r2*Y2 - (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 + v[1]
    dD1dt = mui1*I1 + muy1*Y1
    dD2dt = mui2*I2 + muy2
    dV1dt = v[0]
    dV2dt = v[1]
    new_inf = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 + (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2
    new_reinf = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2

    return [dS1dt, dS2dt, dI1dt, dI2dt, dY1dt, dY2dt, dR1dt, dR2dt, dD1dt, dD2dt, dV1dt, dV2dt, new_inf, new_reinf]

# SIYRD G1 priority vaccination function
def priorG1Vacc(times,init,parms):
    N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M = parms
    S1, S2, I1, I2, Y1, Y2, R1 ,R2, D1, D2, V1, V2, new_I, new_Y =  init

    if S1 < (0.3*N_gr[0]): # if 70% of G1 is recovered (S1<30%), vaccination is stopped for this group
        temp = []
        temp.append(v[1])
        temp.append(v[0])
        v = temp
        if S2 < (0.3*N_gr[1]):# when vaccination is stopped for G1. We check if S2<30%. In that case, vaccination is also stopped
            v[1] = 0
    # ODEs
    dS1dt = -(bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - v[0]
    dS2dt = -(bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2 - v[1]
    dI1dt = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - r1*I1 - mui1*I1
    dI2dt = (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2- r2*I2 - mui2*I2
    dY1dt = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 - r1*Y1 - muy1*Y1
    dY2dt = (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 - r2*Y2 - muy2*Y2
    dR1dt = r1*I1 + r1*Y1 - (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + v[0]
    dR2dt = r2*I2 + r2*Y2 - (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 + v[1]
    dD1dt = mui1*I1 + muy1*Y1
    dD2dt = mui2*I2 + muy2*Y2
    dV1dt = v[0]
    dV2dt = v[1]
    new_inf = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 + (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2
    new_reinf = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2

    return [dS1dt, dS2dt, dI1dt, dI2dt, dY1dt, dY2dt, dR1dt, dR2dt, dD1dt, dD2dt, dV1dt, dV2dt, new_inf, new_reinf]

# SIYRD G2 priority vaccination function
def priorG2Vacc(times,init,parms):
    N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M = parms
    S1, S2, I1, I2, Y1, Y2, R1 ,R2, D1, D2, V1, V2, new_I, new_Y =  init

    if S2 < (0.3*N_gr[1]): # if 70% of G2 is recovered (S2<30%), vaccination is stopped for this group
        temp = []
        temp.append(v[1])
        temp.append(v[0])
        v = temp
        if S1 < (0.3*N_gr[0]): # when vaccination is stopped for G2. We check if S1<30%. In that case, vaccination is also stopped
            v[0] = 0
    # ODEs
    dS1dt = -(bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - v[0]
    dS2dt = -(bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2 - v[1]
    dI1dt = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 - r1*I1 - mui1*I1
    dI2dt = (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2- r2*I2 - mui2*I2
    dY1dt = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 - r1*Y1 - muy1*Y1
    dY2dt = (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 - r2*Y2 - muy2*Y2
    dR1dt = r1*I1 + r1*Y1 - (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + v[0]
    dR2dt = r2*I2 + r2*Y2 - (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2 + v[1]
    dD1dt = mui1*I1 + muy1*Y1
    dD2dt = mui2*I2 + muy2*Y2
    dV1dt = v[0]
    dV2dt = v[1]
    new_inf = (bsi*M[0][0]*I1/N_gr[0] + bsi*M[0][1]*I2/N_gr[1] + bsy*M[0][0]*Y1/N_gr[0] + bsy*M[0][1]*Y2/N_gr[1])*S1 + (bsi*M[1][0]*I1/N_gr[0] + bsi*M[1][1]*I2/N_gr[1] + bsy*M[1][0]*Y1/N_gr[0] + bsy*M[1][1]*Y2/N_gr[1])*S2
    new_reinf = (bri*M[0][0]*I1/N_gr[0] + bri*M[0][1]*I2/N_gr[1] + bry*M[0][0]*Y1/N_gr[0] + bry*M[0][1]*Y2/N_gr[1])*R1 + (bri*M[1][0]*I1/N_gr[0] + bri*M[1][1]*I2/N_gr[1] + bry*M[1][0]*Y1/N_gr[0] + bry*M[1][1]*Y2/N_gr[1])*R2

    return [dS1dt, dS2dt, dI1dt, dI2dt, dY1dt, dY2dt, dR1dt, dR2dt, dD1dt, dD2dt, dV1dt, dV2dt, new_inf, new_reinf]

# SIYRD model implementation
def agesiyrd(country,bsi,bri,bsy,bry,r1,r2,mui1,mui2,muy1,muy2,vacc_choice,vacc_rate,limit):
    # Age distribution input file is opened. Data are saved in two arrays, ages and amounts
    ages, amounts = [], []
    # with open('infiles/Spain/Spain_country_level_age_distribution_85.csv','r') as agefile:
    with open('infiles/'+country+'/'+country+'-2020_gr10.csv','r') as agefile:
        rows = csv.reader(agefile, delimiter=',')
        for row in rows:
            ages.append(row[0])
            amounts.append(float(row[1]))

    # Population density of each group is calculated given the age limit
    num_groups = 2 # two age grups
    age_limits = [] # age limits will always contain [0, limit1,len(ages)]
    age_limits.append(0)
    age_limits.append(int(limit/10)) # sl1 is the age limit /10 is added for the new distribution
    age_limits.append(len(ages))

    # Total population is splitted in two sublists. Then densitiy for each group is calculated.
    # groups = 2    age_limits = [0,80,85]      length_to_split = [80,5]     groups = [[80 elements],[5 elements]]
    length_to_split = []
    for i in range(1,len(age_limits)):
        length_to_split.append(age_limits[i]-age_limits[i-1])
    amounts_iter = iter(amounts) # iter object is needed
    N_gr = [list(islice(amounts_iter, elem)) for elem in length_to_split]
    # by adding all elements in each sublist we get the # of individuals of each group defined
    N_gr = [sum(element) for element in N_gr]
    # by adding sub groups # of individuals we get the total population
    N_total = sum(N_gr)
    # Labels for the plots are created according to age limit
    # groups = 2    age_limits = [0,80,85]----->[80]    numbers = '80'      labels = ['<80','≥ 80']
    labels = []
    # First and last elements are deleted. Corresponding to 0 and 85. We do not need that information from now.
    age_limits = age_limits[1:-1]
    # string containing the limits is created. Needed to open the correct input files regarding contact matrix and ifr
    # numbers = ''.join(map(str, age_limits))
    labels.append('<'+str(limit))
    labels.append('≥'+str(limit))

    # Contact matrix is opened and save in the M array with float number.
    M = []
    with open(('infiles/'+country+'/new_M_'+str(limit)+'.csv'),'r') as matrix:
        rows = csv.reader(matrix, delimiter=',')
        for row in rows:
            M.append(row)
    # convert string values to float values
    for i in range(len(M)):
        for j in range(len(M[i])):
            M[i][j] = float(M[i][j])

    # IFR file is opened, Info is saved in the array ifr as float numbers
    ifr = []
    with open(('infiles/'+country+'/ifr_'+str(limit)+'.csv'),'r') as ifr_file:
        rows = csv.reader(ifr_file, delimiter=',')
        for row in rows:
            ifr = row
    for i in range(len(ifr)):
        ifr[i] = float(ifr[i])

    # ---------------------------------------------------------------------------------------------- #

    # alpha, d, R0 = 0.5, 13, 1
    # bsi = R0 / d
    # bri = 0.01 * bsi
    # bry = alpha * bri
    # bsy = (bry/bri)*bsi
    v_rate = vacc_rate / 100
    # mui, muy, r , v = [], [], [],[]
    v = []
    S, I, Y , R , D = [], [], [] , [], []
    new_I, new_Y = [],[] # new infections and reinfections are saved at each time step
    V = [] # new vaccinations are saved for both groups at each time step

    # Vaccination order is created.
    vacc_order = []
    if v_rate == 0: # no vaccination, vacc_order = [] empty list: choice argument is not considered
        pass
    else:
        if vacc_choice == 1: # vacc_order = [1,2] G1 is first vaccinated
            vacc_order.append(vacc_choice)
            vacc_order.append(vacc_choice+1)
        elif vacc_choice == 2: # vacc_order = [2,1] G2 is first vaccinated
            vacc_order.append(vacc_choice)
            vacc_order.append(vacc_choice-1)

    # Initial conditions and parameters are defined
    for k in range(num_groups):
        I.append(1) # I_i(t=0) = 1
        S.append(N_gr[k]-I[k]) # S_i(t=0) = N_i - I_i
        Y.append(0)
        R.append(0)
        D.append(0)
        V.append(0)
        # mui.append(ifr[k]/d)
        # muy.append(0)
        # r.append((1-ifr[k])/d)
        if len(vacc_order) == 0: # simultaneous vaccination. Same rate for both groups
            v.append(v_rate*N_gr[k]) # daily # of vaccines is proportional to group size
        else:
            if k == (vacc_order[0]-1):
                v.append(v_rate*N_total)
            else:
                v.append(0)

    new_I.append(0)
    new_Y.append(0)

    # Initial conditions vector is created
    y0 = S+I+Y+R+D+V+new_I+new_Y
    # ---------------------------------------------------------------------------------------------- #

    # ------------------------ MODEL SIMULATION ------------------------ #
    # Time steps simulated
    times = np.linspace(0,365,365)
    # Initial conditions
    init = S+I+Y+R+D+V+new_I+new_Y


    # Different conditions for the vaccination lead to diferent model simulation running
    if v_rate == 0: # Execute noVacc function
        parms = N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M
        siyrd_sol = solve_ivp(fun=lambda t, y: noVacc(t,y,parms), t_span=[min(times),max(times)], y0=init, t_eval=times)
        siyrd_out = pd.DataFrame({"t":siyrd_sol["t"],"S1":siyrd_sol["y"][0],"S2":siyrd_sol["y"][1],"I1":siyrd_sol["y"][2],"I2":siyrd_sol["y"][3],"Y1":siyrd_sol["y"][4],"Y2":siyrd_sol["y"][5],"R1":siyrd_sol["y"][6],"R2":siyrd_sol["y"][7],"D1":siyrd_sol["y"][8],"D2":siyrd_sol["y"][9],"V1":siyrd_sol["y"][10],"V2":siyrd_sol["y"][11],"New_inf":siyrd_sol["y"][12],"New_reinf":siyrd_sol["y"][13]})

        graph_title = 'No vaccination'
        outfile_name = 'out_ageSIYRD_noVacc_limit%s.tsv' % (str(limit))

    else:
        if len(vacc_order) == 0: #  Execute simultVacc function
            parms = N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M
            siyrd_sol = solve_ivp(fun=lambda t, y: simultVacc(t,y,parms), t_span=[min(times),max(times)], y0=init, t_eval=times)
            siyrd_out = pd.DataFrame({"t":siyrd_sol["t"],"S1":siyrd_sol["y"][0],"S2":siyrd_sol["y"][1],"I1":siyrd_sol["y"][2],"I2":siyrd_sol["y"][3],"Y1":siyrd_sol["y"][4],"Y2":siyrd_sol["y"][5],"R1":siyrd_sol["y"][6],"R2":siyrd_sol["y"][7],"D1":siyrd_sol["y"][8],"D2":siyrd_sol["y"][9],"V1":siyrd_sol["y"][10],"V2":siyrd_sol["y"][11],"New_inf":siyrd_sol["y"][12],"New_reinf":siyrd_sol["y"][13]})

            graph_title = 'Simultaneous vaccination'
            outfile_name = 'out_ageSIYRD_simultVacc%s_limit%s.tsv' % (vacc_rate,str(limit))

        else: #  Execute priorVacc function
            if vacc_order[0] == 1: # title changes according to vaccination order
                parms = N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M
                siyrd_sol = solve_ivp(fun=lambda t, y: priorG1Vacc(t,y,parms), t_span=[min(times),max(times)], y0=init, t_eval=times)
                siyrd_out = pd.DataFrame({"t":siyrd_sol["t"],"S1":siyrd_sol["y"][0],"S2":siyrd_sol["y"][1],"I1":siyrd_sol["y"][2],"I2":siyrd_sol["y"][3],"Y1":siyrd_sol["y"][4],"Y2":siyrd_sol["y"][5],"R1":siyrd_sol["y"][6],"R2":siyrd_sol["y"][7],"D1":siyrd_sol["y"][8],"D2":siyrd_sol["y"][9],"V1":siyrd_sol["y"][10],"V2":siyrd_sol["y"][11],"New_inf":siyrd_sol["y"][12],"New_reinf":siyrd_sol["y"][13]})

                graph_title = 'G1-priority vaccination'
                outfile_name = 'out_ageSIYRD_priorG1Vacc%s_limit%s.tsv' % (vacc_rate,str(limit))

            else:
                parms = N_gr,bsi,bsy,bri,bry,mui1,mui2,muy1,muy2,r1,r2,v,M
                siyrd_sol = solve_ivp(fun=lambda t, y: priorG2Vacc(t,y,parms), t_span=[min(times),max(times)], y0=init, t_eval=times)
                siyrd_out = pd.DataFrame({"t":siyrd_sol["t"],"S1":siyrd_sol["y"][0],"S2":siyrd_sol["y"][1],"I1":siyrd_sol["y"][2],"I2":siyrd_sol["y"][3],"Y1":siyrd_sol["y"][4],"Y2":siyrd_sol["y"][5],"R1":siyrd_sol["y"][6],"R2":siyrd_sol["y"][7],"D1":siyrd_sol["y"][8],"D2":siyrd_sol["y"][9],"V1":siyrd_sol["y"][10],"V2":siyrd_sol["y"][11],"New_inf":siyrd_sol["y"][12],"New_reinf":siyrd_sol["y"][13]})

                graph_title = 'G2-priority vaccination'
                outfile_name = 'out_ageSIYRD_priorG2Vacc%s_limit%s.tsv' % (vacc_rate,str(limit))


    # Results of the simulation are saved in an output file
    saveResults(outfile_name, times, siyrd_sol)

    # ------------------------ SIMULATION EVOLUTION PLOTTING ------------------------ #
    # Two plots, one for each group of age
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(27.2,7.5))
    fig.suptitle('Scenario: '+graph_title,fontsize=24,fontweight="bold",style='italic',color='#004c98')
    # susceptible
    ax1.plot(siyrd_out["t"], siyrd_out["S1"], 'blue', alpha=0.5, lw=2, label='S1')
    ax2.plot(siyrd_out["t"], siyrd_out["S2"], 'blue', alpha=0.5, lw=2, label='S2')
    # first-infected
    ax1.plot(siyrd_out["t"], siyrd_out["I1"], 'red', alpha=0.5, lw=2, label='I1')
    ax2.plot(siyrd_out["t"], siyrd_out["I2"], 'red', alpha=0.5, lw=2, label='I2')
    # reinfected
    ax1.plot(siyrd_out["t"], siyrd_out["Y1"], 'orangered', alpha=0.5, lw=2, label='Y1')
    ax2.plot(siyrd_out["t"], siyrd_out["Y2"], 'orangered', alpha=0.5, lw=2, label='Y2')
    # recovered
    ax1.plot(siyrd_out["t"], siyrd_out["R1"], 'purple', alpha=0.5, lw=2, label='R1')
    ax2.plot(siyrd_out["t"], siyrd_out["R2"], 'purple', alpha=0.5, lw=2, label='R2')
    # deaths
    ax1.plot(siyrd_out["t"], siyrd_out["D1"], 'darkgreen', alpha=0.5, lw=2, label='D1')
    ax2.plot(siyrd_out["t"], siyrd_out["D2"], 'darkgreen', alpha=0.5, lw=2, label='D2')
    # to plot the population amounts in millions
    scale_y = 1e6
    ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
    # ax1 information
    ax1.set_xlabel('Days',fontsize=20)
    ax1.set_ylabel('Population in millions',fontsize=20)
    ax1.set_title('Group 1 ('+r'$<$'+'%s)' % (str(limit)),fontsize=20,fontweight="bold")
    ax1.legend(loc='center right',fontsize=20)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.yaxis.set_major_formatter(ticks_y)
    # ax2 information
    ax2.set_xlabel('Days',fontsize=20)
    ax2.set_ylabel('Population in millions',fontsize=20)
    ax2.set_title('Group 2 ('+r'$\geq$'+'%s)' % (str(limit)),fontsize=20,fontweight="bold")
    ax2.legend(loc='center right',fontsize=20)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.yaxis.set_major_formatter(ticks_y)

    plt.show()
    fig.tight_layout()

    return outfile_name

# Comparison of simulation results
def resultsAnalysis(age_limit,v_perc,datafiles):
    labels,colors = [], []
    for file in datafiles:
        if 'noVacc' in file:
            labels.append('No vaccination')
            colors.append('blue')
        elif 'G1' in file:
            labels.append('G1 priority')
            colors.append('red')
        elif 'G2' in file:
            labels.append('G2 priority')
            colors.append('orange')
        else:
            labels.append('Simultaneous')
            colors.append('purple')

    # cumulative ammounts and cumulative percentages are saved in two arrays for both: infections and deaths
    cumulative_inf, percentages_inf, cumulative_def, percentages_def = [],[],[],[]
    # total infections and deaths
    infections, deaths = [],[]
    # total vaccinated individuals are saved for each group
    vaccG1, vaccG2 = [],[]
    # total vaacinated individuals and coverages are saved for each group
    vaccs,vacc_coverages = [],[]
    prueba = []
    # we compute the analysis for each vaccination strategy file
    for file in datafiles:
        # each file is opened and information is saved in temporal arrays
        time,S1,S2,I1,I2,Y1,Y2,R1,R2,D1,D2,V1,V2,new_I, new_Y = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        with open(('outfiles/'+file), 'r') as outfile:
            rows = csv.reader(outfile, delimiter='\t')
            next(rows) # header is skipped
            for row in rows:
                time.append(float(row[0]))
                S1.append(float(row[1]))
                S2.append(float(row[2]))
                I1.append(float(row[3]))
                I2.append(float(row[4]))
                Y1.append(float(row[5]))
                Y2.append(float(row[6]))
                R1.append(float(row[7]))
                R2.append(float(row[8]))
                D1.append(float(row[9]))
                D2.append(float(row[10]))
                V1.append(float(row[11]))
                V2.append(float(row[12]))
                new_I.append(float(row[13]))
                new_Y.append(float(row[14]))

        # total number of individuals in the population is calculated
        N_total = S1[0]+S2[0]+I1[0]+I2[0]+Y1[0]+Y2[0]+R1[0]+R2[0]+D1[0]+D2[0]+V1[0]+V2[0]

        # total active infections are calculated at each time step as the sum of (I_i + Y_i)
        total_inf = [I1[t]+I2[t]+Y1[t]+Y2[t] for t in range(len(time))]
        # percentage of infections with respect total population is calculated
        perc_inf = [(inf*100)/N_total for inf in total_inf]

        # total deaths are calculated at each time step as the sum of (D1 + D2)
        total_def = [D1[t]+D2[t] for t in range(len(time))]
        # percentage of deaths with respect total population is calculated
        perc_def = [(df*100)/N_total for df in total_def]

        # total new infections ate calculated at each time step as the sum of new_I + new_Y
        new_inf = [new_I[t]+new_Y[t] for t in range(len(time))]
        # we add the information for each vaccination strategy to its corresponding array for plotting
        percentages_inf.append(perc_inf)
        cumulative_inf.append(total_inf)
        percentages_def.append(perc_def)
        cumulative_def.append(total_def)

        # we calculate the total number of infections, deaths and vaccinated (G1 and G2 separately)
        num_infections, num_deaths, num_vacc1, num_vacc2 = 0,0,0,0
        num_infections = max(new_inf)
        num_deaths = max(total_def)
        num_vacc1 = max(V1)
        num_vacc2 = max(V2)

        total_vacc = num_vacc1+num_vacc2
        vacc_coverage = (total_vacc/N_total)*100
        prueba.append(new_inf)

        # add the results to its corresponging array for future analysis
        infections.append(num_infections)
        deaths.append(num_deaths)
        vaccG1.append(num_vacc1)
        vaccG2.append(num_vacc2)
        vaccs.append(total_vacc)
        vacc_coverages.append(vacc_coverage)

    # Output file is created for the age limit and vaccination rate pre-defined.
    # We save the total infections, vaccinated, deaths and reductions for each vaccina strategy.
    outfile = 'summary.tsv'
    fout = open(('outfiles/'+outfile),"w+")
    fout.write('Vaccination strategy\tTotal infections\tTotal deaths\tInfection reduction(%)\tDeaths reduction(%)\tG1-vaccinated\tG2-vaccinated\tTotal vaccinated\tVaccination coverage (%)\n')
    for i in range(len(datafiles)):
        inf_reduction = 0
        def_reduction = 0
        if i != 0:
            inf_reduction = ((infections[0] - infections[i]) * 100) / infections[0]
            def_reduction = ((deaths[0] - deaths[i]) * 100) / deaths[0]
        fout.write('%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (labels[i],infections[i],deaths[i],inf_reduction,def_reduction, vaccG1[i], vaccG2[i], vaccs[i], vacc_coverages[i]))
    fout.close()
    return outfile

# Delete output files from previous simulations
def deleteOutfiles():
    outfiles = []
    outfiles = glob.glob('outfiles/*.tsv')
    for f in outfiles:
        os.remove(f)
