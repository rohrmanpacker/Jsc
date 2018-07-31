from sesame import utils, analyzer
import matplotlib.pyplot as plt
currents = []
E_GB_to_test = [-.7, -.6, -.5, -.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7]
for E_GB in E_GB_to_test:
    sys, solution = utils.load_sim("GB_test_{0}.gzip".format(E_GB))
    analyze = analyzer.Analyzer(sys, solution)
    currents.append(analyze.full_current() * sys.scaling.current * sys.scaling.length / 3e-4)
plt.plot(E_GB_to_test,currents,'bo-')
plt.title("E_GB vs Jsc")
plt.xlabel("E_GB [eV]")
plt.ylabel("Jsc [A/cm^2]")
plt.show()