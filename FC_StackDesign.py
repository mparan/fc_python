import numpy as np
from math import sqrt
import pylab as plt

"""Input the single-cell results (stack Area and power) for the 15 A, 40 A,
100 A, and 200 A class cells."""

"First entry: Cell Area.  Second entry: Cell Power"
data_15A = np.array([7.81,9.49])
data_40A = np.array([25.31,28.23])
data_100A = np.array([61.74,65.74])
data_200A = np.array([137.81,149.23])

V_cell = 0.75

"Formatting for the figure:"
fs = 12     #font size for plots
lw = 2.0    #line width for plots
font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-2)

def plotStack(V_cell,Area,Power,Rating):
    t_MEA = 0.125
    t_Int = 0.04
    t_endplate = 0.4

    rho_MEA = 5. #g/cm3
    rho_SM = 7.7 #g/cm3


    Stack_power = np.empty([150,1])
    Stack_height = np.empty([150,1])
    Stack_mass = np.empty([150,1])
    Stack_volume = np.empty([150,1])
    Power_mass = np.empty([150,1])
    Power_volume = np.empty([150,1])
    Aspect_ratio = np.empty([150,1])

    cell_array = np.arange(150)

    for i in cell_array:
        n_cells = i+1


        Stack_width = sqrt(Area)
        Stack_power[i] = i*Power #cells*power=stackpower 
        Stack_height[i] = n_cells*t_MEA + (n_cells+1)*t_Int + 2*t_endplate
        Stack_volume[i] = (Area*(Stack_height[i])**2)*.001 #L
        Stack_mass[i] = (n_cells*Area*t_MEA*rho_MEA + (n_cells+1)*Area*t_Int*rho_SM + 2*t_endplate*Area*rho_SM)/1000 #kg
        Power_mass[i] =  Stack_power[i]/Stack_mass[i]  #W/kg
        Power_volume[i] =  Stack_power[i]/Stack_volume[i] #W/L
        Aspect_ratio[i] =  Stack_height[i]/Stack_width

        print('%d cells, Aspect ratio = %1.2f, Power = %1.2f W, Mass = %1.3f kg, Volume = %1.2f L, Power-to-mass = %1.3f W/kg, Power-to-volume = %1.3f W/L'
            %(n_cells, Aspect_ratio[i], Stack_power[i], Stack_mass[i], Stack_volume[i], Power_mass[i], Power_volume[i]))

    fig = plt.figure()
    ax = fig.add_axes([0.2,0.2,0.6,0.75])
    fig.set_size_inches((4.5,3.0))


    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')

    p1, = plt.plot(cell_array,Power_mass,linewidth=lw,color='k')
    p2, = plt.plot(cell_array,Power_volume,linewidth=lw,color='r')

    plt.ylabel('Specific Power (W/kg, W/L)',fontname='Times new Roman',fontsize=fs+2,labelpad=5.0)

    ax2 = ax.twinx()
    p3, = ax2.plot(cell_array,Aspect_ratio,linewidth=lw,color='b')

    "Add a line to denote an aspect ratio of 1.0:"
    plt.plot([0,150],[1.0,1.0],'k--',linewidth=lw-1)

    for label in ax2.get_yticklabels():
        label.set_fontsize(fs)
        label.set_fontname('Times New Roman')

    plt.xlabel('Number of cells',fontname='Times new Roman',fontsize=fs+2,labelpad=5.0)
    plt.ylabel('Aspect Ratio',fontname='Times new Roman',fontsize=fs+2,labelpad=5.0)

    plt.legend([p1,p2,p3],['Power-to-Mass (W/kg)', 'Power-to-Volume (W/L)', 'Asepct Ratio'],handletextpad=0.1,prop=font,frameon=False,loc=7,borderaxespad=0.5)
    plt.savefig('Figure_%s.pdf'%Rating,dpi=350)

"""Create the plots for each current.  It is probably easier to do them one at
a time (comment and uncomment as necessary), to locate the aspect ratio data"""
#plotStack(V_cell,data_15A[0],data_15A[1],'15 A')
#plotStack(V_cell,data_40A[0],data_40A[1],'40 A')
#plotStack(V_cell,data_100A[0],data_100A[1],'100 A')
plotStack(V_cell,data_200A[0],data_200A[1],'200 A')
