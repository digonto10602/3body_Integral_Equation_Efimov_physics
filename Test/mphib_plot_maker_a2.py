#this is the working code cell. The function takes the input file name as a string
#makes the contour plot and save it in a pdf file

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate



def print_mphib_files():
    plt.rcParams.update({'font.size': 12})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    a = 1.16

    filename1 = "Mphib_a_" + str(a) + "_N_900.dat"
    filename2 = "Mphib_a_" + str(a) + "_N_920.dat"
    filename3 = "Mphib_a_" + str(a) + "_N_940.dat"
    filename4 = "Mphib_a_" + str(a) + "_N_960.dat"
    filename5 = "Mphib_a_" + str(a) + "_N_980.dat"
    filename6 = "Mphib_a_" + str(a) + "_N_1000.dat"
    filename7 = "bs1st_a_2_N_1000_smoothcutoff.dat"


    res1,ims1,reM1,imM1,red1,imd1,reM21,imM21,reMdenom1,imMdenom1,a1,N1,gL1, gR1, pth1, eps1 = np.genfromtxt(filename1,unpack=True)
    res2,ims2,reM2,imM2,red2,imd2,reM22,imM22,reMdenom2,imMdenom2,a2,N2,gL2, gR2, pth2, eps2 = np.genfromtxt(filename2,unpack=True)
    res3,ims3,reM3,imM3,red3,imd3,reM23,imM23,reMdenom3,imMdenom3,a3,N3,gL3, gR3, pth3, eps3 = np.genfromtxt(filename3,unpack=True)
    res4,ims4,reM4,imM4,red4,imd4,reM24,imM24,reMdenom4,imMdenom4,a4,N4,gL4, gR4, pth4, eps4 = np.genfromtxt(filename4,unpack=True)
    res5,ims5,reM5,imM5,red5,imd5,reM25,imM25,reMdenom5,imMdenom5,a5,N5,gL5, gR5, pth5, eps5 = np.genfromtxt(filename5,unpack=True)
    res6,ims6,reM6,imM6,red6,imd6,reM26,imM26,reMdenom6,imMdenom6,a6,N6,gL6, gR6, pth6, eps6 = np.genfromtxt(filename6,unpack=True)

    am1,sb1,en1,someN1 = np.genfromtxt(filename7,unpack=True)

    fernando_pole_a2 = (2.6931)**2
    fernando_pole_a16 = (2.9636)**2
    fernando_pole2_a16 = (2.99598)**2

    print("file loaded: ",filename1)
    print("file loaded: ",filename2)
    print("file loaded: ",filename3)
    print("file loaded: ",filename4)
    print("file loaded: ",filename5)
    print("file loaded: ",filename6)

    totreM = (reM1 + reM2 + reM3 + reM4 + reM5 + reM6)/6.0
    totimM = (imM1 + imM2 + imM3 + imM4 + imM5 + imM6)/6.0

    indexpos = 0
    indexneg = 0
    for i in range(0,len(totreM),1):
        if(totreM[i]>0.0 and totreM[i+1]<0.0):
            indexpos = i
            indexneg = i+1
            break
    
    new_reM1 = []
    new_reM2 = []
    new_res1 = []
    new_res2 = []

    for i in range(0,indexpos,1):
        new_reM1.append(totreM[i])
        new_res1.append(res1[i])
    for i in range(indexneg,len(totreM),1):
        new_reM2.append(totreM[i])
        new_res2.append(res1[i])


    fig, ax = plt.subplots(figsize = (12,5))



    #title_str = "am = " + str(a[0])



    #fig.suptitle(title_str,fontsize=20)
    #ax[1].set_title(fig2_str,fontsize=20)

    ax.set_xlim([res1.min(),res1.max()+0.031])
    #ax[1].set_xlim([res.min(),res.max()])
    #ax[2].set_xlim([res.min(),res.max()])

    min_ylim=-8000

    ax.set_ylim([-8000,8000])
    plt.yticks(np.arange(-8000, 8100, 2000.0))
    #ax[1].set_ylim([-20000,20000])
    #ax[2].set_ylim([-500,500])


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    ax.set_xlabel("Re ($s/m^2$)",fontsize=20,horizontalalignment='right')
    ax.xaxis.set_label_coords(.99, -.12)
    
    ax.set_ylabel("$\mathcal{M}_{\phi b}^{HS}$",fontsize=20)
    #ax[1].set_ylabel("$\mathcal{M}_{\phi b}^{II}$",fontsize=20)
    #ax[2].set_ylabel("$1 + 2i \\rho_{\phi b} \mathcal{M}_{\phi b}^{I}$",fontsize=20)


    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].tick_params(axis='both', which='major', labelsize=20)
    #ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)
    

    ax.axhline(y = 0.0, color = 'black', linestyle = '-',linewidth=1)
    ax.plot(new_res1,new_reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    ax.plot(new_res2,new_reM2,linestyle="solid",linewidth=3,color="darkred")
    
    ax.plot(res1,totimM,linestyle="solid",linewidth=3,color="teal",label="Imag")
    ax.scatter(gR1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False, zorder=3)
    ax.scatter(gL1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    ax.scatter(pth1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.scatter(sb1,0,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    ax.axvline(x = fernando_pole_a2, color = 'darkorange', linestyle = 'solid',linewidth=3,label="Romero-L{\\'{o}}pez et al.")
    
    #ax[1].plot(res,reM2,linestyle="solid",linewidth=2,color="brown")
    #ax[2].plot(res,reMdenom,linestyle="solid",linewidth=2,color="brown")
    #ax[2].axhline(y = 0.0, color = 'black', linestyle = '-')

    outputfile_str = "Mphib_a_" + str(a) + ".pdf"

    print("output file: ",outputfile_str)

    plt.legend(loc='upper center')
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()

def print_mphib_singlefile_firstsheet():
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    a = 2


    filename1 = "Mphib_a_2.000000_N_500_imspos_newzoomed.dat"
    #filename1 = "Mphib_a_" + str(a) + "_N_500_imspos.dat"
    
    filename7 = "bs1st_a_2_N_500_smoothcutoff.dat"

    FVfile = "datakcot2_all_s3.txt"
    effAmp = "effRange_amplitude_a2.dat"


    res1,ims1,reM1,imM1,red1,imd1,reM21,imM21,reMdenom1,imMdenom1,a1,N1,gL1, gR1, pth1, eps1 = np.genfromtxt(filename1,unpack=True)
    
    am1,sb1,en1,someN1 = np.genfromtxt(filename7,unpack=True)

    fv_s, fv_ims, fv_qsq, fv_qcot, fv_reM, fv_imM = np.genfromtxt(FVfile, unpack=True)

    effR_s, effR_ims, effR_reM, effR_imM = np.genfromtxt(effAmp, unpack=True)

    fernando_pole_a2 = (2.6931)**2
    fernando_pole_a16 = (2.9636)**2
    fernando_pole2_a16 = (2.99598)**2

    print("file loaded: ",filename1)
    
    totreM = reM1
    totimM = imM1

    indexpos = 0
    indexneg = 0
    for i in range(0,len(totreM),1):
        if(totreM[i]>0.0 and totreM[i+1]<0.0):
            indexpos = i
            indexneg = i+1
            break
    
    new_reM1 = []
    new_reM2 = []
    new_res1 = []
    new_res2 = []

    new_effR_s1 = []
    new_effR_reM1 = []
    new_effR_s2 = []
    new_effR_reM2 = []
    effRindexpos = 0
    effRindexneg = 0

    for i in range(0,len(effR_s),1):
        if(effR_s[i]>6.5):
            if(effR_reM[i]>0.0 and effR_reM[i+1]<0.0):
                effRindexpos = i
                effRindexneg = i+1
                break

    for i in range(0,indexpos,1):
        new_reM1.append(totreM[i])
        new_res1.append(res1[i])
    for i in range(indexneg,len(totreM),1):
        new_reM2.append(totreM[i])
        new_res2.append(res1[i])

    for i in range(0,effRindexpos,1):
        new_effR_reM1.append(effR_reM[i])
        new_effR_s1.append(effR_s[i])
    for i in range(effRindexneg,len(effR_s),1):
        new_effR_reM2.append(effR_reM[i])
        new_effR_s2.append(effR_s[i])

    


    fig, ax = plt.subplots(figsize = (8,5))



    #title_str = "am = " + str(a[0])



    #fig.suptitle(title_str,fontsize=20)
    #ax[1].set_title(fig2_str,fontsize=20)

    ax.set_xlim([6.5,res1.max()+0.01])
    #ax[1].set_xlim([res.min(),res.max()])
    #ax[2].set_xlim([res.min(),res.max()])

    min_ylim=-8000

    
    ax.set_ylim([-8000,8000])
    plt.yticks(np.arange(-8000, 8100, 4000.0))

    labelsy = [item.get_text() for item in ax.get_yticklabels()]
    some_n=0
    for i in range(0,len(labelsy),1):
        item = -0.80 + some_n*0.4
        item = round(item,1)
        labelsy[i] = item
        some_n = some_n + 1

    ax.set_yticklabels(labelsy)
    #ax[1].set_ylim([-20000,20000])
    #ax[2].set_ylim([-500,500])


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    x_ticks = ax.xaxis.get_major_ticks()
    x_ticks[0].label1.set_visible(False)
    x_ticks[len(x_ticks)-1].label1.set_visible(False)

    #ax.set_xlabel("Re ($s/m^2$)",fontsize=20,horizontalalignment='right')
    ax.xaxis.set_label_coords(.99, -.12)
    
    #ax.set_ylabel("$\mathcal{M}_{\phi b}^{HS} \\times 10^{-3}$",fontsize=20)
    #ax[1].set_ylabel("$\mathcal{M}_{\phi b}^{II}$",fontsize=20)
    #ax[2].set_ylabel("$1 + 2i \\rho_{\phi b} \mathcal{M}_{\phi b}^{I}$",fontsize=20)


    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].tick_params(axis='both', which='major', labelsize=20)
    #ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)
    

    ax.axhline(y = 0.0, color = 'black', linestyle = '-',linewidth=2.5,zorder=1)
    #ax.plot(new_res1,new_reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    #ax.plot(new_res2,new_reM2,linestyle="solid",linewidth=3,color="darkred")
    ##ax.plot(res1,reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    
    

    ax.plot(new_res1,new_reM1,linestyle="solid",linewidth=5,color="darkred",label="Real",zorder=2)
    ax.plot(new_res2,new_reM2,linestyle="solid",linewidth=5,color="darkred",zorder=2)
    
    #ax.plot(res1,imM1,linestyle="solid",linewidth=3,color="teal",label="Imag")
    ax.plot(res1,imM1,linestyle="solid",linewidth=5,color="teal",label="Imag",zorder=2)

    #ax.plot(new_effR_s1,new_effR_reM1,linestyle="solid",linewidth=3,color="darkgrey",label="ERE Linear Fit",zorder=10)
    #ax.plot(new_effR_s2,new_effR_reM2,linestyle="solid",linewidth=3,color="darkgrey",zorder=10)
    
    #ax.scatter(fv_s,fv_reM,s=100, marker='o', facecolor = 'white', edgecolor='darkorange', linestyle = 'solid',linewidth=1.5,clip_on=True, zorder=10,label="Romero-L{\\'{o}}pez et al.")
    
    #ax.scatter(gR1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False, zorder=3)
    #ax.scatter(gL1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.scatter(pth1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    
    ax.axvline(x = gR1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    ax.axvline(x = gL1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    ax.axvline(x = pth1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    
    
    #ax.scatter(sb1,0,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    ax.axvline(x = fernando_pole_a2, color = 'darkorange', linestyle = 'solid',linewidth=5,label="FV w/o Cut",zorder=1)
    
    #ax[1].plot(res,reM2,linestyle="solid",linewidth=2,color="brown")
    #ax[2].plot(res,reMdenom,linestyle="solid",linewidth=2,color="brown")
    #ax[2].axhline(y = 0.0, color = 'black', linestyle = '-')

    outputfile_str = "Mphib_firstsheet_a_" + str(a) + "_5.pdf"

    print("output file: ",outputfile_str)

    #plt.legend(loc='lower left')
    #fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()

def print_dS_singlefile_firstsheet(filename):
    plt.rcParams.update({'font.size': 20})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    a = 2


    filename1 = "Mphib_a_2.000000_N_500_imspos_newzoomed.dat"
    #filename1 = "Mphib_a_" + str(a) + "_N_500_imspos.dat"
    
    filename7 = "bs1st_a_2_N_500_smoothcutoff.dat"

    FVfile = "datakcot2_all_s3.txt"
    effAmp = "effRange_amplitude_a2.dat"


    res1,ims1,reM1,imM1,red1,imd1,reM21,imM21,reMdenom1,imMdenom1,a1,N1,gL1, gR1, pth1, eps1,s2kcutoff,scutoffline = np.genfromtxt(filename,unpack=True)
    
    #am1,sb1,en1,someN1 = np.genfromtxt(filename7,unpack=True)

    #fv_s, fv_ims, fv_qsq, fv_qcot, fv_reM, fv_imM = np.genfromtxt(FVfile, unpack=True)

    #effR_s, effR_ims, effR_reM, effR_imM = np.genfromtxt(effAmp, unpack=True)

    fernando_pole_a2 = (2.6931)**2
    fernando_pole_a16 = (2.9636)**2
    fernando_pole2_a16 = (2.99598)**2

    print("file loaded: ",filename)
    
    totreM = reM1
    totimM = imM1

    indexpos = 0
    indexneg = 0
    


    


    fig, ax = plt.subplots(figsize = (12,5))



    title_str = "$am = $" + str(a1[0]) + ", $s_{2,k}^{min} = $" + str(s2kcutoff[0])



    fig.suptitle(title_str,fontsize=20)
    #ax[1].set_title(fig2_str,fontsize=20)

    ax.set_xlim([res1.min(),res1.max()])
    #ax[1].set_xlim([res.min(),res.max()])
    #ax[2].set_xlim([res.min(),res.max()])

    min_ylim=-8000

    
    ax.set_ylim([-25,25])
    #plt.yticks(np.arange(-8000, 8100, 4000.0))

    #labelsy = [item.get_text() for item in ax.get_yticklabels()]
    #some_n=0
    #for i in range(0,len(labelsy),1):
    #    item = -0.80 + some_n*0.4
    #    item = round(item,1)
    #    labelsy[i] = item
    #    some_n = some_n + 1

    #ax.set_yticklabels(labelsy)
    #ax[1].set_ylim([-20000,20000])
    #ax[2].set_ylim([-500,500])


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    #x_ticks = ax.xaxis.get_major_ticks()
    #x_ticks[0].label1.set_visible(False)
    #x_ticks[len(x_ticks)-1].label1.set_visible(False)

    ax.set_xlabel("Re ($s/m^2$)",fontsize=20,horizontalalignment='right')
    #ax.xaxis.set_label_coords(.99, -.12)
    
    ax.set_ylabel("$d_S$",fontsize=20)
    #ax[1].set_ylabel("$\mathcal{M}_{\phi b}^{II}$",fontsize=20)
    #ax[2].set_ylabel("$1 + 2i \\rho_{\phi b} \mathcal{M}_{\phi b}^{I}$",fontsize=20)


    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].tick_params(axis='both', which='major', labelsize=20)
    #ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)
    

    ax.axhline(y = 0.0, color = 'black', linestyle = '-',linewidth=2.5,zorder=1)
    #ax.plot(new_res1,new_reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    #ax.plot(new_res2,new_reM2,linestyle="solid",linewidth=3,color="darkred")
    ##ax.plot(res1,reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    
    

    ax.plot(res1,red1,linestyle="solid",linewidth=5,color="darkred",label="Real",zorder=2)
    ax.plot(res1,imd1,linestyle="solid",linewidth=5,color="teal",label="Imag",zorder=2)
    
    #ax.plot(res1,imM1,linestyle="solid",linewidth=3,color="teal",label="Imag")
    #ax.plot(res1,imM1,linestyle="solid",linewidth=5,color="teal",label="Imag",zorder=2)

    #ax.plot(new_effR_s1,new_effR_reM1,linestyle="solid",linewidth=3,color="darkgrey",label="ERE Linear Fit",zorder=10)
    #ax.plot(new_effR_s2,new_effR_reM2,linestyle="solid",linewidth=3,color="darkgrey",zorder=10)
    
    #ax.scatter(fv_s,fv_reM,s=100, marker='o', facecolor = 'white', edgecolor='darkorange', linestyle = 'solid',linewidth=1.5,clip_on=True, zorder=10,label="Romero-L{\\'{o}}pez et al.")
    
    #ax.scatter(gR1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False, zorder=3)
    #ax.scatter(gL1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.scatter(pth1[0],min_ylim,s=200, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    
    ax.axvline(x = gR1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    ax.axvline(x = gL1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    ax.axvline(x = pth1[0], color = 'black', linestyle = (0,(5,10)),linewidth=2.5,zorder=0)
    ax.axvline(x = scutoffline[0], color = 'orange', linestyle = 'solid',linewidth=2.5,zorder=0)
    
    
    #ax.scatter(sb1,0,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.axvline(x = fernando_pole_a2, color = 'darkorange', linestyle = 'solid',linewidth=5,label="FV w/o Cut",zorder=1)
    
    #ax[1].plot(res,reM2,linestyle="solid",linewidth=2,color="brown")
    #ax[2].plot(res,reMdenom,linestyle="solid",linewidth=2,color="brown")
    #ax[2].axhline(y = 0.0, color = 'black', linestyle = '-')

    outputfile_str = filename + ".pdf"

    print("output file: ",outputfile_str)

    #plt.legend(loc='lower left')
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()


def print_mphib_singlefile_secondsheet():
    plt.rcParams.update({'font.size': 12})
    plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)

    a = 16


    filename1 = "Mphib_a_" + str(a) + "_N_500_imspos.dat"
    
    filename7 = "bs1st_a_2_N_500_smoothcutoff.dat"


    res1,ims1,reM1,imM1,red1,imd1,reM21,imM21,reMdenom1,imMdenom1,a1,N1,gL1, gR1, pth1, eps1 = np.genfromtxt(filename1,unpack=True)
    
    am1,sb1,en1,someN1 = np.genfromtxt(filename7,unpack=True)

    fernando_pole_a2 = (2.6931)**2
    fernando_pole_a16 = (2.9636)**2
    fernando_pole2_a16 = (2.99598)**2

    print("file loaded: ",filename1)
    
    totreM = reM21
    totimM = imM21

    indexpos = 0
    indexneg = 0
    #for i in range(0,len(totreM),1):
    #   if(totreM[i]<0.0 and totreM[i+1]>0.0):
    #        indexpos = i
    #        indexneg = i+1
    #        break
    
    new_reM1 = []
    new_reM2 = []
    new_res1 = []
    new_res2 = []

    for i in range(0,indexpos,1):
        new_reM1.append(totreM[i])
        new_res1.append(res1[i])
    for i in range(indexneg,len(totreM),1):
        new_reM2.append(totreM[i])
        new_res2.append(res1[i])


    fig, ax = plt.subplots(figsize = (14,5))



    #title_str = "am = " + str(a[0])



    #fig.suptitle(title_str,fontsize=20)
    #ax[1].set_title(fig2_str,fontsize=20)

    ax.set_xlim([res1.min(),res1.max()+0.0031])
    #ax[1].set_xlim([res.min(),res.max()])
    #ax[2].set_xlim([res.min(),res.max()])

    min_ylim=-500

    
    ax.set_ylim([-500,500])
    plt.yticks(np.arange(-500, 510, 100.0))
    #ax[1].set_ylim([-20000,20000])
    #ax[2].set_ylim([-500,500])


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    ax.set_xlabel("Re ($s/m^2$)",fontsize=20,horizontalalignment='right')
    ax.xaxis.set_label_coords(.99, -.12)
    
    ax.set_ylabel("$\mathcal{M}_{\phi b}^{HS,II}$",fontsize=20)
    #ax[1].set_ylabel("$\mathcal{M}_{\phi b}^{II}$",fontsize=20)
    #ax[2].set_ylabel("$1 + 2i \\rho_{\phi b} \mathcal{M}_{\phi b}^{I}$",fontsize=20)


    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='both', which='minor', labelsize=20)
    #ax[1].tick_params(axis='both', which='major', labelsize=20)
    #ax[1].tick_params(axis='both', which='minor', labelsize=20)
    #ax[2].tick_params(axis='both', which='major', labelsize=20)
    #ax[2].tick_params(axis='both', which='minor', labelsize=20)
    

    ax.axhline(y = 0.0, color = 'black', linestyle = '-',linewidth=1)
    #ax.plot(new_res1,new_reM1,linestyle="solid",linewidth=3,color="darkred",label="Real")
    #ax.plot(new_res2,new_reM2,linestyle="solid",linewidth=3,color="darkred")
    #ax.plot(res1,reM21,linestyle="solid",linewidth=3,color="darkred",label="Real")
    ax.plot(res1,reM21,linestyle="solid",linewidth=1,color="darkred",label="Real")
    
    #ax.plot(res1,imM21,linestyle="solid",linewidth=3,color="teal",label="Imag")
    ax.plot(res1,imM21,linestyle="solid",linewidth=1,color="teal",label="Imag")
    
    ax.scatter(gR1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False, zorder=3)
    ax.scatter(gL1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    ax.scatter(pth1[0],min_ylim,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.scatter(sb1,0,s=100, marker='o', facecolor = 'white', edgecolor='black', linestyle = 'solid',linewidth=1.5,clip_on=False,zorder=3)
    #ax.axvline(x = fernando_pole_a2, color = 'darkorange', linestyle = 'solid',linewidth=3,label="Romero-L{\\'{o}}pez et al.")
    
    #ax[1].plot(res,reM2,linestyle="solid",linewidth=2,color="brown")
    #ax[2].plot(res,reMdenom,linestyle="solid",linewidth=2,color="brown")
    #ax[2].axhline(y = 0.0, color = 'black', linestyle = '-')

    outputfile_str = "Mphib_secondsheet_a_" + str(a) + ".pdf"

    print("output file: ",outputfile_str)

    plt.legend(loc='upper center')
    fig.tight_layout()
    plt.savefig(outputfile_str)
    plt.close()

scountinitial = 0
scountfinal = 1
pcountinitial = 0
pcountfinal = 204




#print_mphib_files()

#print_mphib_singlefile_firstsheet()
#print_mphib_singlefile_secondsheet()

for i in range(20):
    filename = "Mphib_a_16.000000_N_500_s2kcutoff_"+str(i)+".dat"
    if(os.path.exists(filename)):
        print_dS_singlefile_firstsheet(filename)

#for i in range(0,3337,1):
#    filename = "Mphib_secondsheet_acount_"+str(i)+".dat"
#    if(os.path.exists(filename)):
#        #print_contour_files(opefile,qvecfile,pcutfile,thresholdfile,opecuttracerfile,0,N)
#        print_mphib_files(filename)
