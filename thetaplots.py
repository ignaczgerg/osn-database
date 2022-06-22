import numpy as np
import pandas as pd
import process as pro
import base64
from io import BytesIO
from matplotlib.figure import Figure

def calculate_values(user_r_rejection, user_s_rejection, user_r_racemate):
    P = 4e-6; 
    A = 300; # cm2
    dp = 20; # bar
    alfa = 0.01;
    F0 = 0.05; # L/min
    V_loop = 0.05

    R1 = float(user_r_rejection)/100.0
    R2 = float(user_s_rejection)/100.0
    ratio = float(user_r_racemate)/100.0

    # try:
    #     R1 = float(input('Rejection of R compound: '))/100
    #     R2 = float(input('Rejection of S compound: '))/100
    #     ratio = float(input("Ratio of R solute in the feed (50 for racemic): "))/100
    # except ValueError:
    #     print("Invalid input. Please enter rejection values equal to or lower than 100%")
    #     exit(1)

    # if R1 > 1 or R2 > 1:
    #     print("Invalid input. Please enter rejection values equal to or lower than 100%")
    #     exit(2)

    # if ratio > 1 or ratio < 0:
    #     print("Invalid input. Please enter a ratio value between 0 and 100% (excluding 0 and 100%)")
    #     exit(3)

    # R1 = 0.95
    # R2 = 0.9

    c0 = [2*ratio,2*(1-ratio)]

    theta = np.linspace(0,1,100)
    ee_Rr = np.array([float(0)]*100)
    ee_Rp = np.array([float(0)]*100)
    ee_Sr = np.array([float(0)]*100)
    ee_Sp = np.array([float(0)]*100)
    eta_Rr = np.array([float(0)]*100)
    eta_Rp = np.array([float(0)]*100)
    eta_Sr = np.array([float(0)]*100)
    eta_Sp = np.array([float(0)]*100)

    j = 0
    for i in theta:
        if i == 1 and R1 == 1:
            ee_Rr[j] = 100
        else:
            cc_Rr = (c0[0]/c0[1])*(1-i*R2)/(1-i*R1) # problem: i = 1, R1 = 1
            ee_Rr[j] = 100*(cc_Rr-1)/(cc_Rr+1)
        
        if R2 == 1:
            ee_Rp[j] = -100
        else:
            cc_Rp = cc_Rr*(1-R1)/(1-R2) # problem: R2 = 1
            ee_Rp[j] = 100*(cc_Rp-1)/(cc_Rp+1)
        
        ee_Sr[j] = -ee_Rr[j]
        ee_Sp[j] = -ee_Rp[j]

        if i == 1 and R1 == 1:
            eta_Rr[j] = 100
            eta_Rp[j] = 0
        else:
            eta_Rr[j] = 100*(1-i)/(1-R1*i) # problem: i = 1, R1 = 1
            eta_Rp[j] = 100*i*(1-R1)/(1-i*R1) # problem: i = 1, R1 = 1

        if i == 1 and R2 == 1:
            eta_Sr[j] = 100
            eta_Sp[j] = 0
        else:
            eta_Sr[j] = 100*(1-i)/(1-R2*i) # problem: i = 1, R2 = 1
            eta_Sp[j] = 100*i*(1-R2)/(1-i*R2) # problem: i = 1, R2 = 1
        j += 1

    R = [R1,R2]
    tspan = np.linspace(0,10000,1000)
    z = 0.99

    cel = pro.cell(P)
    system = pro.sys(cel,True)
    system.per_cells.clear()
    system.ractype = 'srrz'
    system.z = z

    theta2 = np.linspace(0.01,1,99)
    ee2R = np.array([float(0)]*99)
    eta2R = np.array([float(0)]*99)
    ee2S = np.array([float(0)]*99)
    eta2S = np.array([float(0)]*99)

    j = 0
    try:
        for i in theta2:
            # system.set_flow(F0,alfa)
            system.set_flow_cut(F0,i)
            system.initial_cell.R = R
            system.initial_cell.c0 = c0
            system.initial_cell.V_loop = V_loop
            system.run_sys(tspan)
            ee2S[j] = 100*((system.res[-1,1]-system.res[-1,0])/(system.res[-1,0]+system.res[-1,1]))
            ee2R[j] = -ee2S[j]
            eta2S[j] = 100*(system.res[-1,1]*system.initial_cell.Fp/(F0*c0[1]))
            eta2R[j] = 100*(system.res[-1,0]*system.initial_cell.Fp/(F0*c0[0]))
            j += 1
    except IndexError:
        print('INDEX_ERROR SRR')

    initial_cel = pro.cell(P)
    system_cas = pro.sys(initial_cel,False)
    per_cel = pro.cell(P)
    system_cas.per_cells.clear()
    system_cas.per_cells.append(per_cel)

    system_cas.initial_cell.R = R
    system_cas.initial_cell.c0 = c0
    system_cas.initial_cell.V_loop = V_loop

    theta3 = np.linspace(0.01,1,99)
    ee3Rr = np.array([float(0)]*99)
    ee3Sr = np.array([float(0)]*99)
    ee3Rp = np.array([float(0)]*99)
    ee3Sp = np.array([float(0)]*99)
    eta3Rr = np.array([float(0)]*99)
    eta3Sr = np.array([float(0)]*99)
    eta3Rp = np.array([float(0)]*99)
    eta3Sp = np.array([float(0)]*99)

    k = 0
    try:
        for i in theta3:
            system_cas.initial_cell.F0 = F0
            system_cas.initial_cell.Fp = F0
            system_cas.initial_cell.Fr = F0*(1-i)

            system_cas.per_cells[0].Fp = F0*i
            system_cas.per_cells[0].Fr = F0*(1-i)
            system_cas.per_cells[0].R = R
            system_cas.run_sys(tspan)

            ee3Rr[k] = 100*((system_cas.res[-1,0]-system_cas.res[-1,4])/(system_cas.res[-1,0]+system_cas.res[-1,4]))
            ee3Sr[k] = -ee3Rr[k]
            ee3Sp[k] = 100*((system_cas.res[-1,7]-system_cas.res[-1,3])/(system_cas.res[-1,7]+system_cas.res[-1,3]))
            ee3Rp[k] = -ee3Sp[k]

            eta3Sp[k] = 100*(system_cas.res[-1,7]*system_cas.per_cells[0].Fp/(F0*c0[1]))
            eta3Sr[k] = 100-eta3Sp[k]
            eta3Rr[k] = 100*(system_cas.res[-1,0]*system_cas.initial_cell.Fr/(F0*c0[0]))
            eta3Rp[k] = 100-eta3Rr[k]
            k+=1
    except IndexError:
        print('INDEX_ERROR CAS')

    return theta, theta2, theta3, ee_Rr, ee_Rp, ee_Sr, ee_Sp, eta_Rr, eta_Rp, eta_Sr, eta_Sp, ee2R, ee2S, eta2R, eta2S, ee3Rr, ee3Rp, ee3Sr, ee3Sp, eta3Rr, eta3Rp, eta3Sr, eta3Sp

def plotting(theta, theta2, theta3, ee_Rr, ee_Rp, ee_Sr, ee_Sp, eta_Rr, eta_Rp, eta_Sr, eta_Sp, ee2R, ee2S, eta2R, eta2S, ee3Rr, ee3Rp, ee3Sr, ee3Sp, eta3Rr, eta3Rp, eta3Sr, eta3Sp):
    # fig = Figure(figsize=(10,5.5), tight_layout=True)
    fig_ss = Figure(figsize=(3.5,5.5), tight_layout=True)
    axs_ss = fig_ss.subplots(2, 1)
    fig_srr = Figure(figsize=(3.5,5.5), tight_layout=True)
    axs_srr = fig_srr.subplots(2, 1)
    fig_c01 = Figure(figsize=(3.5,5.5), tight_layout=True)
    axs_c01 = fig_c01.subplots(2, 1)

    axs_ss[0].plot(theta, ee_Rr, label='Retentate - R', color='navy')
    axs_ss[0].plot(theta, ee_Rp, label='Permeate - R', color='navy', linestyle='dashed')
    axs_ss[0].plot(theta, ee_Sr, label='Retentate - S', color='coral')
    axs_ss[0].plot(theta, ee_Sp, label='Permeate - S', color='coral', linestyle='dashed')
    axs_ss[0].legend(loc='upper left')
    axs_ss[0].set_title('ee - Single stage')
    axs_ss[1].plot(theta, eta_Rr, label='Retentate - R', color='navy')
    axs_ss[1].plot(theta, eta_Rp, label='Permeate - R', color='navy', linestyle='dashed')
    axs_ss[1].plot(theta, eta_Sr, label='Retentate - S', color='coral')
    axs_ss[1].plot(theta, eta_Sp, label='Permeate - S', color='coral', linestyle='dashed')
    axs_ss[1].legend(loc='upper left')
    axs_ss[1].set_title('Recovery (η) - Single stage')

    axs_srr[0].plot(theta2, ee2R, label='R', color='navy')
    axs_srr[0].plot(theta2, ee2S, label='S', color='coral')
    axs_srr[0].legend(loc='upper left')
    axs_srr[0].set_title('ee - SRR')
    axs_srr[1].plot(theta2, eta2R, label='R', color='navy')
    axs_srr[1].plot(theta2, eta2S, label='S', color='coral')
    axs_srr[1].legend(loc='upper left')
    axs_srr[1].set_title('Recovery (η) - SRR')

    axs_c01[0].plot(theta3, ee3Rr, label='Retentate - R', color='navy')
    axs_c01[0].plot(theta3, ee3Rp, label='Permeate - R', color='navy', linestyle='dashed')
    axs_c01[0].plot(theta3, ee3Sr, label='Retentate - S', color='coral')
    axs_c01[0].plot(theta3, ee3Sp, label='Permeate - S', color='coral', linestyle='dashed')
    axs_c01[0].legend(loc='upper left')
    axs_c01[0].set_title('ee - Permeate cascade')
    axs_c01[1].plot(theta3, eta3Rr, label='Retentate - R', color='navy')
    axs_c01[1].plot(theta3, eta3Rp, label='Permeate - R', color='navy', linestyle='dashed')
    axs_c01[1].plot(theta3, eta3Sr, label='Retentate - S', color='coral')
    axs_c01[1].plot(theta3, eta3Sp, label='Permeate - S', color='coral', linestyle='dashed')
    axs_c01[1].legend(loc='upper left')
    axs_c01[1].set_title('Recovery (η) - Permeate cascade')

    axs_ss[0].set(xlabel='', ylabel='ee (%)')
    axs_ss[1].set(xlabel='θ (-)', ylabel='η (%)')
    axs_srr[0].set(xlabel='', ylabel='ee (%)')
    axs_srr[1].set(xlabel='θ (-)', ylabel='η (%)')
    axs_c01[0].set(xlabel='', ylabel='ee (%)')
    axs_c01[1].set(xlabel='θ (-)', ylabel='η (%)')

    axs_ss[0].set_xlim([0, 1])
    axs_ss[1].set_xlim([0, 1])
    axs_ss[0].set_ylim([0, 100])
    axs_ss[1].set_ylim([0, 100])

    axs_srr[0].set_xlim([0, 1])
    axs_srr[1].set_xlim([0, 1])
    axs_srr[0].set_ylim([0, 100])
    axs_srr[1].set_ylim([0, max(max(eta2R),max(eta2S),200)])


    axs_c01[0].set_xlim([0, 1])
    axs_c01[1].set_xlim([0, 1])
    axs_c01[0].set_ylim([0, 100])
    axs_c01[1].set_ylim([0, 100])


    #for ax in axs.flat:
    #    ax.label_outer()

    # export data
    # notes and information
    # 3 process diagram
    # changes in article (based on Gergo's article) Highlights, abstract, methods, res., concl.

    # export data
    # RES = np.transpose(np.concatenate([[theta2], [ee_Rr[1:]],[ee_Rp[1:]],[ee_Sr[1:]],[ee_Sp[1:]],
    #     [eta_Rr[1:]],[eta_Rp[1:]],[eta_Sr[1:]],[eta_Sp[1:]],[ee2R],[ee2S],[eta2R],[eta2S],[ee3Rr],[ee3Rp],[ee3Sr],[ee3Sp],
    #     [eta3Rr],[eta3Rp],[eta3Sr],[eta3Sp]]))

    # headers = ['stage_cut','ss_ee_Rr','ss_ee_Rp','ss_ee_Sr','ss_ee_Sp','ss_rec_Rr','ss_rec_Rp','ss_rec_Sr','ss_rec_Sp',
    #     'srr_ee_R','srr_ee_S','srr_rec_R','srr_rec_S','c01_ee_Rr','c01_ee_Rp','c01_ee_Sr','c01_ee_Sp',
    #     'c01_rec_Rr','c01_rec_Rp','c01_rec_Sr','c01_rec_Sp']

    # df = pd.DataFrame(RES)
    # df.columns = headers
    # df.to_csv('./static/results.csv', index=False)

    #plt.savefig('thetaplot_2.png', dpi=300)

    
    filepath_ss = 'static/theta_ss.png'
    fig_ss.savefig(filepath_ss, dpi=500, bbox_inches='tight')

    filepath_srr = 'static/theta_srr.png'
    fig_srr.savefig(filepath_srr, dpi=500, bbox_inches='tight')

    filepath_c01 = 'static/theta_c01.png'
    fig_c01.savefig(filepath_c01, dpi=500, bbox_inches='tight')

    # buf = BytesIO()
    # fig.savefig(buf, format='png')
    # data = base64.b64encode(buf.getbuffer()).decode("ascii")
    ''' do not change data else the code breaks and I do not know why.'''


# if __name__ == "__main__":
#     calculate_diagrams(95, 90, 50)