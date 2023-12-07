import numpy as np
import sys
sys.path.append('utils')
from plotsUtil import *
import fileio as io
import os

def computeFieldErr(ind, A, Aref):
    interpField = np.interp(A[:,0], Aref[:,0], Aref[:,ind])
    if A[-1,0] < 0.099:
        error = np.nan
    else:
        error = np.linalg.norm(A[:,ind] - interpField) / A.shape[0]
    return error

def computeError(file, fileRef_np=None, fileRef= 'ref/PPreaction.txt'):
    if fileRef_np is not None:
        Aref = np.load(fileRef_np)
    else:
        Aref = io.readMultiColFile(fileRef)
    A = io.readMultiColFile(file)
    return computeFieldErr(1, A, Aref)


def read_sol(folder):
    try:
        timing=io.readFile(folder, "timing.txt")
    except FileNotFoundError:
        return {"err": np.nan, "timing": np.nan}
    err=io.readFile(folder, "err.txt")
    minTime = np.amin(timing)
    meanTime = np.mean(timing)
    stdTime = np.std(timing)
    #if stdTime/np.sqrt(len(timing)) > 0.01* meanTime:
    #    print(timing)
    #    sys.exit(f"ERROR : high uncertainty on timing for {folder}")
    #err = computeError(2, os.path.join(folder, 'PPreaction.txt'))
    return {"err": np.mean(err), "timing": minTime}
 
def get_sol(ndt_list, tol_list, sol_list):
    allSolDict = {}
    for sol in sol_list:
        err = np.zeros((len(tol_list), len(ndt_list)))
        timing = np.zeros((len(tol_list), len(ndt_list)))
        ndt_tot = []
        tol_tot = []
        for indt, ndt in enumerate(ndt_list):
            for itol, tol in enumerate(tol_list):
                folder = f"{sol}_{tol}_{ndt}"
                folder = folder.replace("1e-0", "1e-")
                tmpDict = read_sol(folder)
                err[itol, indt] = tmpDict["err"]
                timing[itol, indt] = tmpDict["timing"]
                ndt_tot.append(ndt)
                tol_tot.append(tol)
        allSolDict[sol] = {}
        allSolDict[sol]["err"] = err
        allSolDict[sol]["timing"] = timing
    allSolDict["tol"] = np.array(tol_tot)
    allSolDict["ndt"] = np.array(ndt_tot)
    return allSolDict


#ndt_list = [1000, 10000, 400000]
#tol_list = [1e-6, 1e-8, 1e-10]
ndt_list=[100, 500, 1000, 5000, 10000, 50000, 100000, 400000, 800000, 3200000, 10000000]
tol_list=[1e-6, 1e-8, 1e-10, 1e-12, 1e-14]
ndt_list.sort()
tol_list.sort()
sol_list = ["aj", "denseDir", "gmres", "gmresTyp"]
sol_list_jac = ["aj", "denseDir"]
sol_list_gmres = ["gmres", "gmresTyp"]
#sol_list = ["aj", "denseDir", "gmresTyp"]
sol_title = {"aj" : "Analytical Jac.", "denseDir": "Numerical Jac.", "gmres": "GMRES", "gmresTyp": "GMRES + Typical Val"}

allSolDict = get_sol(ndt_list, tol_list, sol_list)

fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(allSolDict["gmres"]["timing"].flatten(), allSolDict["gmres"]["err"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(allSolDict["gmresTyp"]["timing"].flatten(), allSolDict["gmresTyp"]["err"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
plt.plot(allSolDict["aj"]["timing"].flatten(), allSolDict["aj"]["err"].flatten(), 'o', color='r', label=sol_title["aj"])
plt.plot(allSolDict["denseDir"]["timing"].flatten(), allSolDict["denseDir"]["err"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax = plt.gca()
ax.set_yscale("log")
plotLegend()
prettyLabels("Exec Time [s]", "Error [K]", 14)
plt.savefig("errorDrop.png")
plt.savefig("errorDrop.eps")


fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(0.1/allSolDict["ndt"], allSolDict["gmres"]["err"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["gmresTyp"]["err"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["aj"]["err"].flatten(), 'o', color='r', label=sol_title["aj"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["denseDir"]["err"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
plotLegend()
prettyLabels("Timestep [s]", "Error [K]", 14)


fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(allSolDict["tol"], allSolDict["gmres"]["err"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(allSolDict["tol"], allSolDict["gmresTyp"]["err"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
plt.plot(allSolDict["tol"], allSolDict["aj"]["err"].flatten(), 'o', color='r', label=sol_title["aj"])
plt.plot(allSolDict["tol"], allSolDict["denseDir"]["err"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
plotLegend()
prettyLabels("Tolerance", "Error [K]", 14)

fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(0.1/allSolDict["ndt"], allSolDict["gmres"]["timing"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["gmresTyp"]["timing"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["aj"]["timing"].flatten(), 'o', color='r', label=sol_title["aj"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["denseDir"]["timing"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
plotLegend()
prettyLabels("Timestep [s]", "Exec time [s]", 14)

fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(0.1/allSolDict["ndt"], allSolDict["gmres"]["timing"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(0.1/allSolDict["ndt"], allSolDict["gmresTyp"]["timing"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
plotLegend()
prettyLabels("Timestep [s]", "Exec time [s]", 14)

fig = plt.figure()
if "gmres" in sol_list:
    plt.plot(allSolDict["tol"], allSolDict["gmres"]["timing"].flatten(), 'o', color='k', label=sol_title["gmres"])
plt.plot(allSolDict["tol"], allSolDict["gmresTyp"]["timing"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
plt.plot(allSolDict["tol"], allSolDict["aj"]["timing"].flatten(), 'o', color='r', label=sol_title["aj"])
plt.plot(allSolDict["tol"], allSolDict["denseDir"]["timing"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax = plt.gca()
ax.set_yscale("log")
ax.set_xscale("log")
plotLegend()
prettyLabels("Tolerance", "Exec time [s]", 14)


fig = plt.figure()
ax1=plt.gca()
if "gmres" in sol_list:
    ax1.plot(0.1/allSolDict["ndt"], allSolDict["gmres"]["timing"].flatten(), 'o', color='k', label=sol_title["gmres"])
ax1.plot(0.1/allSolDict["ndt"], allSolDict["gmresTyp"]["timing"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
ax1.plot(0.1/allSolDict["ndt"], allSolDict["aj"]["timing"].flatten(), 'o', color='r', label=sol_title["aj"])
ax1.plot(0.1/allSolDict["ndt"], allSolDict["denseDir"]["timing"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax2 = ax1.twinx()
if "gmres" in sol_list:
    ax2.plot(0.1/allSolDict["ndt"], allSolDict["gmres"]["err"].flatten(), 'x', color='k')
ax2.plot(0.1/allSolDict["ndt"], allSolDict["gmresTyp"]["err"].flatten(), 'x', color='b')
ax2.plot(0.1/allSolDict["ndt"], allSolDict["aj"]["err"].flatten(), 'x', color='r')
ax2.plot(0.1/allSolDict["ndt"], allSolDict["denseDir"]["err"].flatten(), 'x', color='g')
ax1.set_yscale("log")
ax1.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xscale("log")
#axplotLegend(ax1)
#axplotLegend(ax2)
axprettyLabels(ax1,"Timestep [s]", "Exec time [s]", 14)
axprettyLabels(ax2,"", "Error [K]", 14)


fig = plt.figure()
ax1=plt.gca()
if "gmres" in sol_list:
    ax1.plot(allSolDict["tol"], allSolDict["gmres"]["timing"].flatten(), 'o', color='k', label=sol_title["gmres"])
ax1.plot(allSolDict["tol"], allSolDict["gmresTyp"]["timing"].flatten(), 'o', color='b', label=sol_title["gmresTyp"])
ax1.plot(allSolDict["tol"], allSolDict["aj"]["timing"].flatten(), 'o', color='r', label=sol_title["aj"])
ax1.plot(allSolDict["tol"], allSolDict["denseDir"]["timing"].flatten(), 'o', color='g', label=sol_title["denseDir"])
ax2 = ax1.twinx()
if "gmres" in sol_list:
    ax2.plot(allSolDict["tol"], allSolDict["gmres"]["err"].flatten(), 'x', color='k')
ax2.plot(allSolDict["tol"], allSolDict["gmresTyp"]["err"].flatten(), 'x', color='b')
ax2.plot(allSolDict["tol"], allSolDict["aj"]["err"].flatten(), 'x', color='r')
ax2.plot(allSolDict["tol"], allSolDict["denseDir"]["err"].flatten(), 'x', color='g')
ax1.set_yscale("log")
ax1.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xscale("log")
#axplotLegend(ax1)
#axplotLegend(ax2)
axprettyLabels(ax1,"Tolerance", "Exec time [s]", 14)
axprettyLabels(ax2,"", "Error [K]", 14)


#plt.show()

qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list],
    [substep_arr for _ in sol_list],
    [allSolDict[sol][qoi] for sol in sol_list],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list)-1)] + ["Exec time [s]"],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list],
    listYAxisName=["Tolerance" for _ in sol_list],
    listXAxisName=["Timestep [s]" for _ in sol_list],
    interp="nearest",
    vminList=[1e-1, 1e-1, 1e-1, 1e-1],
    vmaxList=[1e3, 1e3, 1e3, 1e3], 
)
#plt.show()
#stop
plt.savefig("exec_sol.png")
plt.savefig("exec_sol.eps")

qoi="err"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list],
    [substep_arr for _ in sol_list],
    [allSolDict[sol][qoi] for sol in sol_list],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list)-1)] + ["Error [K]"],
    #listCBLabel=["Error [K]" for _ in sol_list],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list],
    listYAxisName=["Tolerance" for _ in sol_list],
    listXAxisName=["Timestep [s]" for _ in sol_list],
    interp="nearest",
)
plt.savefig("err_sol.png")
plt.savefig("err_sol.eps")

qoi="err"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list_jac],
    [substep_arr for _ in sol_list_jac],
    [allSolDict[sol][qoi] for sol in sol_list_jac],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list_jac)-1)] + ["Error [K]"],
    #listCBLabel=["Error [K]" for _ in sol_list],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list_jac],
    listYAxisName=["Tolerance" for _ in sol_list_jac],
    listXAxisName=["Timestep [s]" for _ in sol_list_jac],
    interp="nearest",
)
plt.savefig("err_sol_jac.png")
plt.savefig("err_sol_jac.eps")


qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list_jac],
    [substep_arr for _ in sol_list_jac],
    [allSolDict[sol][qoi] for sol in sol_list_jac],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list_jac)-1)] + ["Execution time [s]"],
    #listCBLabel=["Error [K]" for _ in sol_list],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list_jac],
    listYAxisName=["Tolerance" for _ in sol_list_jac],
    listXAxisName=["Timestep [s]" for _ in sol_list_jac],
    interp="nearest",
    vminList=[1e-1, 1e-1],
    vmaxList=[1e3, 1e3], 
)
plt.savefig("exec_sol_jac.png")
plt.savefig("exec_sol_jac.eps")



qoi="err"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list_gmres],
    [substep_arr for _ in sol_list_gmres],
    [allSolDict[sol][qoi] for sol in sol_list_gmres],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list_gmres)-1)] + ["Error [K]"],
    #listCBLabel=["Error [K]" for _ in sol_list],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list_gmres],
    listYAxisName=["Tolerance" for _ in sol_list_gmres],
    listXAxisName=["Timestep [s]" for _ in sol_list_gmres],
    interp="nearest",
)
plt.savefig("err_sol_gmres.png")
plt.savefig("err_sol_gmres.eps")

qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr for _ in sol_list_gmres],
    [substep_arr for _ in sol_list_gmres],
    [allSolDict[sol][qoi] for sol in sol_list_gmres],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    listCBLabel=["" for _ in range(len(sol_list_gmres)-1)] + ["Execution time [s]"],
    #listCBLabel=["Error [K]" for _ in sol_list],
    #listTitle=[sol_title[sol] for sol in sol_list],
    listTitle=["" for sol in sol_list_gmres],
    listYAxisName=["Tolerance" for _ in sol_list_gmres],
    listXAxisName=["Timestep [s]" for _ in sol_list_gmres],
    interp="nearest",
    vminList=[1e-1, 1e-1],
    vmaxList=[1e3, 1e3],
)
plt.savefig("exec_sol_gmres.png")
plt.savefig("exec_sol_gmres.eps")


qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
plot2dContour_log(
    [tol_arr],
    [substep_arr],
    [allSolDict["denseDir"][qoi]/allSolDict["aj"][qoi]],
    [np.amin(tol_arr), np.amax(tol_arr)],
    [np.amin(substep_arr), np.amax(substep_arr)],
    #listCBLabel=["Exec time numerical Jac. / Exec time analytical Jac."],
    listCBLabel=[""],
    listTitle=[" "],
    listYAxisName=["Tolerance"],
    listXAxisName=["Timestep [s]"],
    interp="nearest",
)
plt.savefig("speedupAJ.png")
plt.savefig("speedupAJ.eps")


qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
fig = plt.figure()
plt.plot(allSolDict["tol"], allSolDict["denseDir"]["timing"].flatten()/allSolDict["aj"]["timing"].flatten(), 'o', color='k')
ax=plt.gca()
ax.set_xscale("log")
#prettyLabels("Tolerance", "Exec time numerical Jac. / Exec time analytical Jac.", 14)
prettyLabels("Tolerance", "Speedup analyical Jac.", 14)


qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
fig = plt.figure()
plt.plot(0.1/allSolDict["ndt"], allSolDict["denseDir"]["timing"].flatten()/allSolDict["aj"]["timing"].flatten(), 'o', color='k')
ax=plt.gca()
ax.set_xscale("log")
#prettyLabels("Substep [s]", "Exec time numerical Jac. / Exec time analytical Jac.", 14)    
prettyLabels("Timestep [s]", "Speedup analytical Jac.", 14)    

tmp_gmres = allSolDict["gmres"]["timing"].flatten()
tmp_gmres[np.argwhere(0.1/allSolDict["ndt"].flatten() > 6e-6)] = np.nan

qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
fig = plt.figure()
plt.plot(allSolDict["tol"], tmp_gmres/allSolDict["aj"]["timing"].flatten(), 'o', color='k')
ax=plt.gca()
ax.set_xscale("log")
prettyLabels("Tolerance", "Exec time GMRES / Exec time analytical Jac.", 14)

qoi="timing"
tol_arr = np.array(tol_list)
substep_arr = 0.1/np.array(ndt_list)
fig = plt.figure()
plt.plot(0.1/allSolDict["ndt"], tmp_gmres/allSolDict["aj"]["timing"].flatten(), 'o', color='k')
ax=plt.gca()
ax.set_xscale("log")
prettyLabels("Timestep [s]", "Exec GMRES / Exec time analytical Jac.", 14)    







