import numpy as np
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib import cm
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

def prettyLabels(xlabel,ylabel,fontsize,title=None, grid=True):
    plt.xlabel(xlabel, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    plt.ylabel(ylabel, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    if not title==None:
        plt.title(title, fontsize=fontsize, fontweight='bold', fontname="Times New Roman")
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color("black")
    if grid:
        plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.tight_layout()

def axprettyLabels(ax, xlabel, ylabel, fontsize, title=None, grid=True):
    ax.set_xlabel(
        xlabel,
        fontsize=fontsize,
        fontweight="bold",
        fontname="Times New Roman",
    )   
    ax.set_ylabel(
        ylabel,
        fontsize=fontsize,
        fontweight="bold",
        fontname="Times New Roman",
    )
    if not title == None:
        ax.set_title(
            title,
            fontsize=fontsize,
            fontweight="bold",
            fontname="Times New Roman",
        )
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight("bold")
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontname("Times New Roman")
        tick.label1.set_fontweight("bold")
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color("black")
    if grid:
        ax.grid(color="k", linestyle="-", linewidth=0.5)
    try:
        plt.tight_layout()
    except:
        print("Could not call tight_layout")
        pass

def plotLegend():
    fontsize = 16
    plt.legend()
    leg=plt.legend(prop={'family':'Times New Roman','size': fontsize-3,'weight':'bold' })
    leg.get_frame().set_linewidth(2.0)
    leg.get_frame().set_edgecolor('k')    

def axplotLegend(ax):
    fontsize = 16
    ax.legend()
    leg = ax.legend(
        prop={
            "family": "Times New Roman",
            "size": fontsize - 3,
            "weight": "bold",
        }
    )
    leg.get_frame().set_linewidth(2.0)
    leg.get_frame().set_edgecolor("k")

def plot2dContour(
    listDatax,
    listDatay,
    listData,
    xbound,
    ybound,
    listCBLabel,
    listTitle,
    listXAxisName=None,
    listYAxisName=None,
    vminList=None,
    vmaxList=None,
    interp="bicubic",
    xticks=None,
    yticks=None,
):
    lim = -1
    lim_vmax_t = -1
    lim_vmax_x = -1
    lim_plot = -1
    fig, axs = plt.subplots(1, len(listData), figsize=(len(listData) * 3, 4))
    if len(listData) == 1:
        i_dat = 0
        data = listData[i_dat]
        data_x = np.squeeze(listDatax[i_dat])
        if vminList == None:
            vmin = np.nanmin(data[:lim, :])
        else:
            vmin = vminList[i_dat]
        if vmaxList == None:
            vmax = np.nanmax(data[:lim, :])
        else:
            vmax = vmaxList[i_dat]
        im = axs.imshow(
            data[:lim, :],
            cmap=cm.viridis,
            interpolation=interp,
            vmin=vmin,
            vmax=vmax,
            extent=[xbound[0], xbound[1], ybound[1], ybound[0]],
            aspect="auto",
        )  
        
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="10%", pad=0.2)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(listCBLabel[i_dat])
        ax = cbar.ax
        text = ax.yaxis.label
        font = matplotlib.font_manager.FontProperties(
            family="times new roman", weight="bold", size=14
        )
        text.set_font_properties(font)
        if i_dat == 0:
            axprettyLabels(
                axs,
                listXAxisName[i_dat],
                listYAxisName[i_dat],
                12,
                listTitle[i_dat],
            )
        else:
            axprettyLabels(axs, listXAxisName[i_dat], "", 12, listTitle[i_dat])

        #axs.set_xticks([])  # values
        #axs.set_xticklabels([])  # labels
        #axs.set_xscale("log")
        #axs.set_yscale("log")
        if xticks is not None:
            axs.set_xticks(xticks)
        if yticks is not None:
            axs.set_yticks(yticks)
        if not i_dat == 0:
            axs.set_yticks([])  # values
            axs.set_yticklabels([])  # labels
        for l in cbar.ax.yaxis.get_ticklabels():
            l.set_weight("bold")
            l.set_family("serif")
            l.set_fontsize(12)
    else:
        for i_dat in range(len(listData)):
            data = listData[i_dat]
            data_x = np.squeeze(listDatax[i_dat])
            if vminList == None:
                vmin = np.nanmin(data[:lim, :])
            else:
                vmin = vminList[i_dat]
            if vmaxList == None:
                vmax = np.nanmax(data[:lim, :])
            else:
                vmax = vmaxList[i_dat]
            im = axs[i_dat].imshow(
                data[:lim, :],
                cmap=cm.viridis,
                interpolation=interp,
                vmin=vmin,
                vmax=vmax,
                extent=[xbound[0], xbound[1], ybound[1], ybound[0]],
                aspect="auto",
            )
            divider = make_axes_locatable(axs[i_dat])
            cax = divider.append_axes("right", size="10%", pad=0.2)
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_label(listCBLabel[i_dat])
            ax = cbar.ax
            text = ax.yaxis.label
            font = matplotlib.font_manager.FontProperties(
                family="times new roman", weight="bold", size=14
            )
            text.set_font_properties(font)
            if i_dat == 0:
                axprettyLabels(
                    axs[i_dat],
                    listXAxisName[i_dat],
                    listYAxisName[i_dat],
                    12,
                    listTitle[i_dat],
                )
            else:
                axprettyLabels(
                    axs[i_dat],
                    listXAxisName[i_dat],
                    "",
                    12,
                    listTitle[i_dat],
                )
            #axs[i_dat].set_xticks([])  # values
            #axs[i_dat].set_xticklabels([])  # labels
            #if not i_dat == 0:
            #    axs[i_dat].set_yticks([])  # values
            #    axs[i_dat].set_yticklabels([])  # labels
            for l in cbar.ax.yaxis.get_ticklabels():
                l.set_weight("bold")
                l.set_family("serif")
                l.set_fontsize(12)
            #axs[i_dat].set_xscale("log")
            #axs[i_dat].set_yscale("log")
            if xticks is not None:
                axs[i_dat].set_xticks(xticks)
            if yticks is not None:
                axs[i_dat].set_yticks(yticks)


def plot2dContour_log(
    listDatax,
    listDatay,
    listData,
    xbound,
    ybound,
    listCBLabel,
    listTitle,
    listXAxisName=None,
    listYAxisName=None,
    vminList=None,
    vmaxList=None,
    interp="bicubic",
    xticks=None,
    yticks=None,
):
    lim = -1
    lim_vmax_t = -1
    lim_vmax_x = -1
    lim_plot = -1
    if len(listData)==1:
        fig, axs = plt.subplots(1, len(listData), figsize=(len(listData) * 4, 3)) 
    else:
        fig, axs = plt.subplots(1, len(listData), figsize=(len(listData) * 3, 3))
    if len(listData) == 1:
        i_dat = 0
        data = listData[i_dat]
        data_x = np.squeeze(listDatax[i_dat])
        if vminList == None:
            vmin = np.nanmin(data[:lim, :])
        else:
            vmin = vminList[i_dat]
        if vmaxList == None:
            vmax = np.nanmax(data[:lim, :])
        else:
            vmax = vmaxList[i_dat]
        X, Y = np.meshgrid(listDatay[i_dat], listDatax[i_dat])
        Z = listData[i_dat]
        Zmid = np.zeros((Z.shape[0]-1, Z.shape[1]-1))
        for i in range(Z.shape[0]-1):
            for j in range(Z.shape[1]-1):
                Zmid[i,j] = (Z[i,j] + Z[i,j+1] + Z[i+1,j] + Z[i+1,j+1])/4
        im = axs.pcolor(
            X, Y, Zmid,
            cmap=cm.viridis,
            #interpolation=interp,
            #vmin=vmin,
            #vmax=vmax,
            #norm=LogNorm(),
            shading='auto',
            #extent=[xbound[0], xbound[1], ybound[1], ybound[0]],
            #aspect="auto",
        )  
        divider = make_axes_locatable(axs)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(listCBLabel[i_dat])
        ax = cbar.ax
        text = ax.yaxis.label
        font = matplotlib.font_manager.FontProperties(
            family="times new roman", weight="bold", size=20
        )
        text.set_font_properties(font)
        if i_dat == 0:
            axprettyLabels(
                axs,
                listXAxisName[i_dat],
                listYAxisName[i_dat],
                20,
                listTitle[i_dat],
                grid=False
            )
        else:
            axprettyLabels(axs, listXAxisName[i_dat], "", 20, listTitle[i_dat], grid=False)

        #axs.set_xticks([])  # values
        #axs.set_xticklabels([])  # labels
        axs.set_xscale("log")
        axs.set_yscale("log")
        if xticks is not None:
            axs.set_xticks(xticks)
        if yticks is not None:
            axs.set_yticks(yticks)
        if not i_dat == 0:
            axs.set_yticks([])  # values
            axs.set_yticklabels([])  # labels
        for l in cbar.ax.yaxis.get_ticklabels():
            l.set_weight("bold")
            l.set_family("serif")
            l.set_fontsize(15)
    else:
        for i_dat in range(len(listData)):
            data = listData[i_dat]
            data_x = np.squeeze(listDatax[i_dat])
            if vminList == None:
                vmin = np.nanmin(data[:lim, :])
            else:
                vmin = vminList[i_dat]
            if vmaxList == None:
                vmax = np.nanmax(data[:lim, :])
            else:
                vmax = vmaxList[i_dat]
            X, Y = np.meshgrid(listDatay[i_dat], listDatax[i_dat])
            Z = listData[i_dat]
            Zmid = np.zeros((Z.shape[0]-1, Z.shape[1]-1))
            for i in range(Z.shape[0]-1):
                for j in range(Z.shape[1]-1):
                    Zmid[i,j] = (Z[i,j] + Z[i,j+1] + Z[i+1,j] + Z[i+1,j+1])/4
            
            im = axs[i_dat].pcolor(
                X, Y, Zmid,
                cmap=cm.viridis,
                #interpolation=interp,
                #vmin=vmin,
                #vmax=vmax,
                norm=LogNorm(vmin, vmax),
                shading='auto',
                #extent=[xbound[0], xbound[1], ybound[1], ybound[0]],
                #aspect="auto",
            )  
            if i_dat == len(listData)-1:
                divider = make_axes_locatable(axs[i_dat])
                cax = divider.append_axes("right", size="5%", pad=0.1)
                cbar = fig.colorbar(im, cax=cax)
                cbar.set_label(listCBLabel[i_dat])
                ax = cbar.ax
                text = ax.yaxis.label
                font = matplotlib.font_manager.FontProperties(
                    family="times new roman", weight="bold", size=20
                )
                text.set_font_properties(font)
            if i_dat == 0:
                axprettyLabels(
                    axs[i_dat],
                    listXAxisName[i_dat],
                    listYAxisName[i_dat],
                    20,
                    listTitle[i_dat],
                    grid=False
                )
            else:
                axprettyLabels(
                    axs[i_dat],
                    listXAxisName[i_dat],
                    "",
                    20,
                    listTitle[i_dat],
                    grid=False
                )
            #axs[i_dat].set_xticks([])  # values
            #axs[i_dat].set_xticklabels([])  # labels
            #if not i_dat == 0:
            #    axs[i_dat].set_yticks([])  # values
            #    axs[i_dat].set_yticklabels([])  # labels
            if i_dat == len(listData)-1:
                for l in cbar.ax.yaxis.get_ticklabels():
                    l.set_weight("bold")
                    l.set_family("serif")
                    l.set_fontsize(15)
            axs[i_dat].set_xscale("log")
            axs[i_dat].set_yscale("log")
            axs[i_dat].plot([5e-8, 5e-8], [np.amin(Y), np.amax(Y)], '--',  lw=3, color='r')
            axs[i_dat].plot([5e-5, 5e-5], [np.amin(Y), np.amax(Y)], '--',  lw=3, color='b')
            if xticks is not None:
                axs[i_dat].set_xticks(xticks)
            if yticks is not None:
                axs[i_dat].set_yticks(yticks)
            if not i_dat == 0:
                axs[i_dat].set_yticks([])  # values
                axs[i_dat].set_yticklabels([])  # labels
    fig.tight_layout()
    #plt.subplots_adjust(wspace=0.5, hspace=0)

