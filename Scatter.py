import matplotlib.pylab as plt
import numpy as np
from scipy import norm
import stat as st

N = 100


def Scatter(dim,A,B,C,D):
    fig1 = plt.figure(1)
    fig1.suptitle('Position Scatter Plots')
    fig1.set_size_inches(6,6)
    W = dim[0]
    H = dim[1]
    d1 = dim[2]
    d2 = dim[3]
    T = dim[4]
    ax1 = fig1.add_subplot(2,2,1)
    ax2 = fig1.add_subplot(2,2,2)
    ax3 = fig1.add_subplot(2,2,3)
    ax4 = fig1.add_subplot(2,2,4)
    ax1.scatter(A[0],A[1],s = 3,c='lime')
    ax2.scatter(A[0],A[1],s = 3,c='lime')
    ax3.scatter(A[1],A[2],s = 3,c='lime')
    ax4.scatter(A[0],A[2],s = 3,c='lime')
    ax1.scatter(B[0],B[1],s = 3,c='red')
    ax3.scatter(B[1],B[2],s = 3,c='red')
    ax4.scatter(B[0],B[2],s = 3,c='red')
    ax2.scatter(C[0],C[1],s = 3,c='blue')
    ax3.scatter(C[1],C[2],s = 3,c='blue')
    ax4.scatter(C[0],C[2],s = 3,c='blue')
    ax1.scatter(D[0],D[1],s = 3,c='purple')
    ax3.scatter(D[1],D[2],s = 3,c='purple')
    ax4.scatter(D[0],D[2],s = 3,c='purple')
    ax1.set_title('Plate 1')
    ax1.set_ylabel('H (mm)')
    ax1.set_xlabel('W (mm)')
    ax1.vlines(x = -W/2, ymin = -H/2, ymax = H/2, linewidth = 1, color ='k')
    ax1.vlines(x = W/2, ymin = -H/2, ymax = H/2, linewidth = 1, color ='k')
    ax1.hlines(y = H/2, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k')
    ax1.hlines(y = -H/2, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k')  
    ax2.set_title('Plate 2')
    ax2.set_ylabel('H (mm)')
    ax2.set_xlabel('W (mm)')
    ax2.vlines(x = -W/2, ymin = -H/2, ymax = H/2, linewidth = 1, color ='k')
    ax2.vlines(x = W/2, ymin = -H/2, ymax = H/2, linewidth = 1, color ='k')
    ax2.hlines(y = H/2, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k')
    ax2.hlines(y = -H/2, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k') 
    ax3.set_title('d x H')
    ax3.set_ylabel('d (mm)')
    ax3.set_xlabel('H (mm)')
    ax3.hlines(y = d2+T, xmin = -H/2, xmax = H/2, linewidth = 1, color ='k')
    ax3.hlines(y = d2, xmin = -H/2, xmax = H/2, linewidth = 1, color ='k')
    ax3.hlines(y = d1, xmin = -H/2, xmax = H/2, linewidth = 1, color ='k')  
    ax3.hlines(y = 0, xmin = -H/2, xmax = H/2, linewidth = 1, color ='lime')  
    ax3.hlines(y = -(T-d1), xmin = -H/2, xmax = H/2, linewidth = 1, color ='k') 
    ax4.set_title('d x W')
    ax4.set_ylabel('d (mm)')
    ax4.set_xlabel('W (mm)')
    ax4.hlines(y = d2+T, xmin = -H/2, xmax = H/2, linewidth = 1, color ='k')
    ax4.hlines(y = d2, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k')
    ax4.hlines(y = d1, xmin = -W/2, xmax = W/2, linewidth = 1, color ='k')
    ax4.hlines(y = 0, xmin = -H/2, xmax = H/2, linewidth = 1, color ='lime') 
    ax4.hlines(y = -(T-d1), xmin = -H/2, xmax = H/2, linewidth = 1, color ='k') 
    a = W/10
    b = H/10
    D = d2-d1
    c = D/10
    ax1.set_xlim([-W/2-a,W/2+a])
    ax1.set_ylim([-H/2-b,H/2+b])
    ax2.set_xlim([-W/2-a,W/2+a])
    ax2.set_ylim([-H/2-b,H/2+b])
    ax3.set_xlim([-W/2-a,W/2+a])
    ax3.set_ylim([-(T-d1)-c,d2+T+c])
    ax4.set_xlim([-H/2-b,H/2+b])
    ax4.set_ylim([-(T-d1)-c,d2+T+c])
    plt.tight_layout()

def HistE(A,B,C,fig,string):
    D = np.concatenate((B,C))
    ax1 = fig.add_subplot(241)
    ax2 = fig.add_subplot(242)
    ax3 = fig.add_subplot(243)
    ax4 = fig.add_subplot(244)
    n_bins = int(10)
    ax1.hist(A, bins = n_bins, facecolor = 'lime')
    ax2.hist(B, bins = n_bins, facecolor = 'blue')
    ax3.hist(C, bins = n_bins, facecolor = 'red')
    ax4.hist(D, bins = n_bins, facecolor = 'purple')
    ax1.set_ylabel(f'N = {N}')
    ax1.set_title('Source')
    ax2.set_title('Plate 1')
    ax3.set_title('Plate 2')
    ax4.set_title('Total (P1 + P2)')


def PDF_E(A,B,C,fig,string):
    A = np.sort(A)
    B = np.sort(B)
    C = np.sort(C)
    D = np.concatenate((B,C))
    D = np.sort(D)
    m1 = st.mean(A)
    sd1 = st.stdev(A)
    m2 = st.mean(B)
    sd2 = st.stdev(B)
    m3 = st.mean(C)
    sd3 = st.stdev(C)
    m4 = st.mean(D)
    sd4 = st.stdev(D)
    ax5 = fig.add_subplot(245)
    ax6 = fig.add_subplot(246)
    ax7 = fig.add_subplot(247)
    ax8 = fig.add_subplot(248)
    ax5.plot(A, norm.pdf(A,m1,sd1),color = 'lime')
    ax6.plot(B, norm.pdf(B,m2,sd2),color = 'blue')
    ax7.plot(C, norm.pdf(C,m3,sd3),color = 'red')
    ax8.plot(D, norm.pdf(D,m4,sd4),color = 'purple')
    ax5.set_ylabel('PDF')
    ax5.set_xlabel('Initial Energy (keV)')
    ax6.set_xlabel(string+' '+'(keV)')
    ax7.set_xlabel(string+' '+'(keV)')
    ax8.set_xlabel(string+' '+'(keV)')

def PositionGraphs(dim,A,B,C,D):
    Scatter(dim,A,B,C,D)
    #HistPos(A,B,C)
    #PDF_Pos(A,B,C)


def EnergyGraphs(A,B,C,D,E):
    # Measured
    st1 = 'Measured Energy'
    fig4 = plt.figure(4)
    fig4.suptitle(st1)
    fig4.set_size_inches(10,4)
    HistE(A,B,C,fig4,st1)
    PDF_E(A,B,C,fig4,st1)
    plt.tight_layout()
    # True
    st2 = 'True Energy'
    fig5 = plt.figure(5)
    fig5.suptitle(st2)
    fig5.set_size_inches(10,4)
    HistE(A,D,E,fig5,st2)
    PDF_E(A,D,E,fig5,st2)
    plt.tight_layout()

def Graphs(A1,B1,C1,D1,E1,A2,B2,C2,D2,E2):
    plt.clf
    PositionGraphs(A1,B1,C1,D1,E1)
    #EnergyGraphs(A2,B2,C2,D2,E2)
    plt.show()