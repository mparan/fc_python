import numpy as np

"""The file 'PlanarSOFC.py' has to be in the same folder as this file.  Our
SOFC Model is saved as a function called 'Model'.  Once we load this function,
we can now call it."""

from HW12_PlanarSOFC import Model

"""The function requires two inputs: V_cell and Q_anode.  You can either load
them in that order:"""

V_cell = 0.75
w_ch = 0.5

max_length = 21.0
n_steps = int(max_length/w_ch) + 1

"A range of Q values, to find that which gives E_fuel = 66.7%"
Q_range = np.arange(0.7,0.05,-0.005)

"Make some flags to specify whether the size for a given current has been set:"
flag_15 = 0
flag_40 = 0
flag_100 = 0
flag_200 = 0

#Empty arrays for the power and current.
Power = np.empty([n_steps,1])
Current = np.empty([n_steps,1])

"Vary the number/length of channels:"
for i in np.arange(1,n_steps):

    "For each length, find the flow rate that gives E_fuel = 66.7%"
    for Q_value in Q_range:
        "Calculate the performance of a single channel SOFC:"
        [E_fuel, E_agg, channel_Power] = Model(V_cell, Q_value,w_ch,i)

        if E_fuel >= 2./3.:
            break

    """Calculate the Power (W) and current (A) from the model output.  Don't
    forget that you have 'i' channels, while the PlanarSOFC results are for
    just one channel."""
    Power[i] = channel_Power*i
    Current[i] = Power[i]/V_cell

    """Check to see if the current size exceeds the lowest rating (15, 40, 100,
    or 200 A) that has not yet been specified:"""
    if not flag_15:
        print('%d channels (length = %1.2f cm), I = %1.2f A'
            %(i,i*w_ch,Current[i]))
        if Current[i] > 15.: #Power exceeds target.  Use previous value.
            # We have found the largest cell where I <= 15 A.  Set flag = 1:
            flag_15 = 1
            print('For 15 A-rated cell, Active Area = %1.2f cm2, Cell Area = %1.2f cm2, Power = %1.2f W, Current = %1.2f A'
                %(((i-1)*w_ch)**2, 1.25*((i-1)*w_ch)**2, Power[i-1],Current[i-1]))
    elif not flag_40:
        print('%d channels (length = %1.2f cm), I = %1.2f A'
            %(i,i*w_ch,Current[i]))
        if Current[i] > 40.:  #Power exceeds target.  Use previous value.
            # We have found the largest cell where I <= 40 A.  Set flag = 1:
            flag_40 = 1
            print('For 40 A-rated cell, Active Area = %1.2f cm2, Cell Area = %1.2f cm2, Power = %1.2f W, Current = %1.2f A'
                %(((i-1)*w_ch)**2, 1.25*((i-1)*w_ch)**2, Power[i-1],Current[i-1]))
    elif not flag_100:
        print('%d channels (length = %1.2f cm), I = %1.2f A'
            %(i,i*w_ch,Current[i]))
        if Current[i] > 100.:  #Power exceeds target.  Use previous value.
            # We have found the largest cell where I <= 100 A.  Set flag = 1:
            flag_100 = 1
            print('For 100 A-rated cell, Active Area = %1.2f cm2, Cell Area = %1.2f cm2, Power = %1.2f W, Current = %1.2f A'
                %(((i-1)*w_ch)**2, 1.25*((i-1)*w_ch)**2, Power[i-1],Current[i-1]))
    elif not flag_200:
        print('%d channels (length = %1.2f cm), I = %1.2f A'%(i,i*w_ch,Current[i]))
        if Current[i] > 200.:  #Power exceeds target.  Use previous value.
            # We have found the largest cell where I <= 200 A.  Set flag = 1:
            flag_200 = 1
            print('For 200 A-rated cell, Active Area = %1.2f cm2, Cell Area = %1.2f cm2, Power = %1.2f W, Current = %1.2f A'
            %(((i-1)*w_ch)**2, 1.25*((i-1)*w_ch)**2, Power[i-1],Current[i-1]))
            break
    else:
        # All cell sizes found.  Exit the loop.
        break
