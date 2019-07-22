def Model(V_cell, Q_anode, w_ch, n_ch):
    import numpy as np
    import math

    #Note: our fuel cell reaciton is H2 + 1/2 O2 -> H2O
    #  So all of our reaction terms (delta h_rxn, delta g_rxn) will be per mole of h2

    "Inputs from the assignment:"
    T_cell = 900    # cell temperature (K)
    #V_cell = 0.75   # cell voltage (V)
    #w_ch = 1.5      # channel width (cm)
    h_ch = 0.5      # channel height (cm)
    l_cell = n_ch*w_ch   # channel length  - equal to channel width times
                        # number of channels.

    dz = 0.05       # axial step size (cm)

    x_h2_inlet = 0.97     # mole fraction of hydrogen fed to the anode
    x_h2o_inlet = 0.03    # mole fraction of steam fed to the anode
    x_o2 = 0.21     # mole fraction of oxygen fed to the cathode

    #Q_anode = 0.75   # anode flow rate, SLPM
    rho_STP = 0.1021 # density of anode flow at STP (kg per standard m3) - needed to convert from SLPM

    # We only need the j-V data for the two points surrounding our oparting voltage.  We will use this data
    #   to calcaulte the ASR (Ohm-cm2)
    j_data = np.array([1.333, 1.667, 2.000])
    V_data = np.array([0.80, 0.75, 0.70])

    # ASR is the negative slope of the j-V curve:
    dV = V_data[2] - V_data[0]
    dj = j_data[2] - j_data[0]
    ASR = -dV/dj

    # Reaction info:
    n_h2 = -1  #number of h2 produced per reaction (it is consumed, so we put a negative number)
    n_h2o = 1  #number of h2o produced per reaction
    n_o2 = -0.5 #number of o2 produced per reaction
    n_elec = 2  #number of electrons transferred per reaction

    "Constants"
    pi = math.pi
    R = 8.3145*1e-3      # universal gas constant (kJ/mol)
    F = 96485      # Faraday's constant (mol of electron per Coulomb of charge)
    MW_h2 = 0.002    # molecular weight of hydrogen (kg/mol)
    MW_h2o = 0.018  # molecular weight of water/steam (kg/mol)
    MW_avg_loc = x_h2_inlet*MW_h2 + x_h2o_inlet*MW_h2o

    "Table Lookups"
    #Standard Gibbs energies (kJ/mol) from your textbook
    g_o_h2o = -425.51
    g_o_h2 = -129.07
    g_o_o2 = -196.73

    #Standard enthalpies (kJ/mol) from your textbook
    h_o_h2o =-219.89
    h_o_h2 = 17.68
    h_o_o2 = 19.29


    "Efficiency calculation terms that will not vary with current density:"
    # delta_g_o_rxn - the standard-state gibbs free energy of reaction:
    Delta_g_o_rxn = n_h2o*g_o_h2o + n_h2*g_o_h2 + n_o2*g_o_o2

    #Equations for delta_h_rxn and delta_g_rxn:
    Delta_h_rxn = n_h2o*h_o_h2o + n_h2*h_o_h2 + n_o2*h_o_o2

    # Total molar flow rate to anode, in mol/s - PAY ATTENTION TO UNITS! 1 m3 = 1000 L
    Q_unitConvert  = Q_anode/60.0/1000.0 # convert to standard m3 per s
    ndot_anode_inlet = Q_unitConvert*rho_STP/MW_avg_loc

    # Molar flow rate of hydrogen (mol/s) sent to the anode:
    ndot_h2_provided = ndot_anode_inlet*x_h2_inlet

    # Total power provided via fuel flow to anode:
    Firing_rate = ndot_h2_provided*(-1)*Delta_h_rxn*1000

    # Numger of steps to get down the channel (we add 1 because the exit is an extra step - but we only
    #   need/want the flows at the exit):
    nsteps = int(l_cell/dz) + 1

    # Create arrays for all of the terms we want to calculate. They should be the same size as nsteps
    #   We will leave them empty, for now.
    z = np.empty([nsteps,1])                  # location of the local 'slice'.
    x_h2 = np.empty([nsteps,1])               # hydrogen mole fraction
    x_h2o = np.empty([nsteps,1])              # water/steam mole fraction
    j_local = np.empty([nsteps,1])            # local current density (A/cm2) for each 'slice'.
    Delta_g_local = np.empty([nsteps,1])        # local Gibbs free energy of reaction
    E_rev_local = np.empty([nsteps,1])        # local reversible cell potential (V)
    PowerDensity_local = np.empty([nsteps,1]) # local power density (W/cm2) produced by each slice.
    Power_agg = np.empty([nsteps,1])        # Aggregate power (W) produced by the SOFC from th inlet, up through
                                            #       the local slice.
    E_fuel = np.empty([nsteps,1])             # Aggregate fuel utilization efficiency - % of fuel consumed up through
                                            #       the local slice.
    E_agg = np.empty([nsteps,1])              # Aggregate efficiency up through the local slice.  We could use eq. 1
                                            #       from HW 8, but it is perhaps clearer to use eq. 2.
    mdot_h2 = np.empty([nsteps,1])            # Mass flux of H2 *into* cell i (kg/cm2/s)
    mdot_h2o = np.empty([nsteps,1])           # Mass flux of H2O *into* cell i (kg/cm2/s)

    "Load in the inlet conditions (the zeroeth step)"
    x_h2[0] = x_h2_inlet
    x_h2o[0] = x_h2o_inlet
    MW_avg_loc = x_h2[0]*MW_h2 + x_h2o[0]*MW_h2o
    mdot_tot_loc = MW_avg_loc*ndot_anode_inlet/h_ch/w_ch
    Y_h2_loc = x_h2[0]*MW_h2/MW_avg_loc
    Y_h2o_loc = x_h2o[0]*MW_h2o/MW_avg_loc
    mdot_h2[0] = mdot_tot_loc*Y_h2_loc
    mdot_h2o[0] = mdot_tot_loc*Y_h2o_loc
    z[0] = dz
    Delta_g_local[0] = Delta_g_o_rxn + R*T_cell*math.log(x_h2o[0]**n_h2o*x_h2[0]**n_h2*x_o2**n_o2)

    # Reversible cell voltage - PAY ATTENTION TO UNITS.
    E_rev_local[0] = -Delta_g_local[0]*1000./n_elec/F

    # Overpotential - net driving force (V)
    Eta = E_rev_local[0] - V_cell

    # Local current density
    j_local[0] = Eta/ASR

    J_h2 = j_local[0]/n_elec/F
    J_h2o = -J_h2

    # local power denisty:
    PowerDensity_local[0] = V_cell*j_local[0]

    Power_agg[0] = w_ch*dz*PowerDensity_local[0]

    mdot_h2[1] = mdot_h2[0] - J_h2*MW_h2*dz/h_ch
    mdot_h2o[1] = mdot_h2o[0] - J_h2o*MW_h2o*dz/h_ch

    mdot_tot = mdot_h2[1] + mdot_h2o[1]
    Y_h2_loc = mdot_h2[1]/mdot_tot
    Y_h2o_loc= mdot_h2o[1]/mdot_tot

    MW_avg_loc = 1/(Y_h2_loc/MW_h2 + Y_h2o_loc/MW_h2o)

    x_h2[1] = Y_h2_loc*MW_avg_loc/MW_h2
    x_h2o[1] = Y_h2o_loc*MW_avg_loc/MW_h2o

    E_fuel[0] = (mdot_h2[0] - mdot_h2[1])/mdot_h2[0]
    E_agg[0] = Power_agg[0]/Firing_rate

    "Loop through each slice and calculate the required quantities:"
    for i in np.arange(1,nsteps-1):

        # distance from inlet:
        z[i] = i*dz

        Delta_g_local[i] = Delta_g_o_rxn + R*T_cell*math.log(x_h2o[i]**n_h2o*x_h2[i]**n_h2*x_o2**n_o2)

        # Reversible cell voltage - PAY ATTENTION TO UNITS.
        E_rev_local[i] = -Delta_g_local[i]*1000./n_elec/F

        # Overpotential - net driving force (V)
        Eta = E_rev_local[i] - V_cell

        # Local current density
        j_local[i] = Eta/ASR

        # local power denisty:
        PowerDensity_local[i] = V_cell*j_local[i]

        Power_agg[i] = Power_agg[i-1] + w_ch*dz*PowerDensity_local[i]

        # Molar rate of h2 and h2o consumption (mol/cm2/s):
        J_h2 = j_local[i]/n_elec/F
        J_h2o = -J_h2

        mdot_h2[i+1] = mdot_h2[i] - J_h2*MW_h2*dz/h_ch
        mdot_h2o[i+1] = mdot_h2o[i] - J_h2o*MW_h2o*dz/h_ch

        mdot_tot = mdot_h2[i+1] + mdot_h2o[i+1]
        Y_h2_loc = mdot_h2[i+1]/mdot_tot
        Y_h2o_loc= mdot_h2o[i+1]/mdot_tot

        MW_avg_loc = 1/(Y_h2_loc/MW_h2 + Y_h2o_loc/MW_h2o)

        x_h2[i+1] = Y_h2_loc*MW_avg_loc/MW_h2
        x_h2o[i+1] = Y_h2o_loc*MW_avg_loc/MW_h2o

        E_fuel[i] = (mdot_h2[0]-mdot_h2[i])/mdot_h2[0]
        E_agg[i] = Power_agg[i]/Firing_rate


    "Fill in the final values for E_fuel and z:"
    E_fuel[-1] = (mdot_h2[0]-mdot_h2[-1])/mdot_h2[0]
    z[-1] = nsteps*dz
    A_cell = w_ch*l_cell

    return np.hstack([E_fuel[-1], E_agg[-2], Power_agg[-2]]);
