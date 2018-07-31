from sesame import builder, solvers, utils, analyzer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes
# creating the system
x = np.concatenate((np.linspace(0,.75e-4,5,False),np.linspace(.75e-4,1e-4,20,False),
                    np.linspace(1e-4,2e-4,120,False),np.linspace(2e-4,2.25e-4,20,False),
                    np.linspace(2.25e-4,3e-4,5)))
y = np.concatenate((np.linspace(0,1.25e-4,60,False),np.linspace(1.25e-4,1.75e-4,120,False),np.linspace(1.75e-4,3e-4,60)))
# 0, 3e-4

material = {'Nc':8e17, 'Nv':1.8e19, 'Eg':1.5, 'affinity':3.9, 'epsilon':9.4,
        'mu_e':100, 'mu_h':100, 'tau_e':10e-9, 'tau_h':10e-9, 'Et':0} # arbitrary values

transitions = {'donor' : (1, 0), 'acceptor' : (0, -1), 'symmetric': (1, -1)}

junction = 1.5e-4 # extent of the junction from the left contact [cm]
def n_region(pos):
    x = pos[0]
    return x < junction

def p_region(pos):
    x = pos[0]
    return x >= junction

# Add the donors
nD = 1e16 # [cm^-3]


# Add the acceptors
nA = 1e16 # [cm^-3]


# Define contacts: CdS and CdTe contacts are Ohmic

# Define the surface recombination velocities for electrons and holes [cm/s]
Sn_left, Sp_left, Sn_right, Sp_right = 1e7, 0, 0, 1e7
# This function specifies the simulation contact recombination velocity

# GB defect state properties
rho_GB = 1e14               # defect density [1/cm^2]
S_GB = 1e-14                # capture cross section [cm^2]
# Specify the two points that make the line containing additional charges
p1 = (1e-4, 1.5e-4)      # [cm]
p2 = (2e-4, 1.5e-4)     # [cm]

def makesystem(transition):
    sys = builder.Builder(x, y)
    sys.add_material(material)
    sys.add_donor(nD, n_region)
    sys.add_acceptor(nA, p_region)
    sys.contact_type('Ohmic', 'Ohmic')
    sys.contact_S(Sn_left, Sp_left, Sn_right, Sp_right)
    sys.add_line_defects([p1, p2], rho_GB, S_GB, E=E_GB, transition=transition)

    return(sys)
def noGB():
    sys = builder.Builder(x, y)
    sys.add_material(material)
    sys.add_donor(nD, n_region)
    sys.add_acceptor(nA, p_region)
    sys.contact_type('Ohmic', 'Ohmic')
    sys.contact_S(Sn_left, Sp_left, Sn_right, Sp_right)

    sys.generation(lambda x, y: 1e21)
    equilibrium = solvers.solve_equilibrium(sys)
    solution = solvers.solve(sys, equilibrium)
    analyze = analyzer.Analyzer(sys, solution)
    utils.save_sim(sys, solution, "noGB.gzip")
    return float(analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4)

E_GB_to_test = [-.75, -.7, -.6, -.5, -.4, -.35, -.3, -.2, -.1, 0, .1, .2, .3, .35, .4, .5, .6, .7, .75]
G_for_guess = [1e15,1e17,1e19,1e20,1e21]

failed = [[],[],[]]
success = [[],[],[]]
currents = [[],[],[]]

count = 0
for t in transitions:

    file = f"GB_test_{t}_"
    for E_GB in E_GB_to_test:  # energy of gap state from intrinsic level [eV]
        try:
            sys, solution = utils.load_sim(file+f"{E_GB}.gzip")
            print(f"loading {E_GB} from file")
            analyze = analyzer.Analyzer(sys, solution)
            currents[count].append(analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4)
            success[count].append(E_GB)
        except FileNotFoundError:
            print(f"calculating {E_GB}")
            #building system
            sys = makesystem(transitions[t])
            sys.generation(lambda x, y: 1e15)
             # find base solution
            equilibrium = solvers.solve_equilibrium(sys)
            solution = solvers.solve(sys, equilibrium)
            # iterate up to 1e21
            for g in G_for_guess:
                sys.generation(lambda x, y: g)
                solution = solvers.solve(sys, solution)
                if solution is None:
                    print("Failed at " + str(g))
                    break
            try:
                analyze = analyzer.Analyzer(sys, solution)
                currents[count].append(analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4)  # units A/cm^2
                print(str(E_GB) + ": " + str(analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4) + " Amp/cm^2")
                # .03 Amp/cm^2 is usual

                name = file + f"{E_GB}"
                # add some system settings to the saved results
                filename = "%s.gzip" % name
                utils.save_sim(sys, solution, filename)
                success[count].append(E_GB)
            except TypeError:
                print("Failed to converge for: " + str(E_GB))
                failed[count].append(E_GB)
                pass
    if not failed[count] == []:
        print("Failed to converge for" + str(failed))
    count += 1

try:
    sys, solution = utils.load_sim("noGB.gzip")
    analyze = analyzer.Analyzer(sys,solution)
    noGBsol = analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4
except FileNotFoundError:
    noGBsol = noGB()

donarplot,     = plt.plot(success[0],currents[0],'ro-')
acceptorplot,  = plt.plot(success[1],currents[1],'go-')
symmetricplot, = plt.plot(success[2],currents[2],'bo-')

plt.title("EGB vs Jsc")
plt.xlabel("E_GB [eV]")
plt.ylabel("Jsc [A/cm^2]")

# line = mlines.Line2D([E_GB_to_test[0],[noGBsol]],[E_GB_to_test[-1],noGBsol],color='r')
plt.hlines(noGBsol,E_GB_to_test[0],E_GB_to_test[-1],color = 'y')
# fig, ax = plt.subplots()
# line = axes.Axes.axhline(noGBsol/ax.get_ybound())
# ax.add_line(line)
plt.show()