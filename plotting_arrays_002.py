import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import os
import numpy as np

path = os.getcwd()
#path = 'g:\\My Drive\\Data\\PyCharmProjects\\New FRAP analysis software\\20190116 - testing new vs old supercharged data\\eco fast\\0\\20160310\\'
print('Processing directory:\n%s ...' % path)

f_d = os.listdir(path)


cell_names = []
cell_ids = []

try:
    Ds_table = pd.read_csv(os.path.join(path, 'diffusion-coefficients.csv')).astype({'cell_id': str})
except FileNotFoundError:
    Ds_table = ''
    print \
        ('#############################################\n####diffusion-coefficients.csv not found!####\n#############################################')
#print(Ds_table['cell_id'])
#print(Ds_table['diffusion_coefficient'])


for f in f_d:
    if f.endswith('.csv') and f[0].isdigit():
        # print(f)
        cell_name = f.split('-')[0]
        cell_names.append(cell_name)

        for t in cell_name.split('_'):
            # print(int(t))
            try:
                cell_ids.append(int(t))
            except ValueError:
                pass


cell_ids = list(set(cell_ids))
cell_ids.sort()




print('There are %s cells to analyze.' % len(cell_ids))

# print(cell_ids)
pdf = PdfPages('allcells.pdf')
for cid_int in cell_ids:

    for name in cell_names:
        if name.startswith(str(cid_int)):
            cid = name
    #print(cid)
    #print(type(cid))
    #print(Ds_table['cell_id'])
    #print(type(Ds_table['cell_id']))
    print('Processing %s...' % cid)
    #print(Ds_table['cell_id'] == cid)
    #print(Ds_table.loc[Ds_table['cell_id'] == cid, 'diffusion_coefficient'].iloc[0])


    file_not_found = 0
    not_in_table = 0
    try:

        D_coe = np.format_float_scientific \
            ((10 ** 12) * (Ds_table.loc[Ds_table['cell_id'] == cid, 'diffusion_coefficient'].iloc[0]), precision=2)
        D_err = np.format_float_scientific \
            ((10 ** 12) * (Ds_table.loc[Ds_table['cell_id'] == cid, 'diffusion_coefficient_error'].iloc[0]), precision=2)
    except IndexError:
        not_in_table = 1
        D_coe = 'not found'
        D_err = 'not found'
        print('##\n## %s not in table! \n##' % cid)
    except AttributeError:
        D_coe = 'no file'
        D_err = 'no file'


    try:
        g_val = pd.read_csv(os.path.join(path, '%s-values.csv' % cid))
        g_sim = pd.read_csv(os.path.join(path, '%s-simulated.csv' % cid))
        g_res = pd.read_csv(os.path.join(path, '%s-residuals.csv' % cid))
        g_sim_hr = pd.read_csv(os.path.join(path, '%s-simulated-high-res.csv' % cid))
        g_res_hr = pd.read_csv(os.path.join(path, '%s-residuals-high-res.csv' % cid))
    except FileNotFoundError:
        file_not_found = 1




    if file_not_found == 0:


        z_g_val = np.asarray(g_val.iloc[0:, 1:], dtype=np.float)
        x_g_val = np.asarray(g_val.columns.values[1:], dtype=np.float)
        y_g_val = np.asarray(g_val.iloc[:, 0].values, dtype=np.float)

        x_g_val, y_g_val = np.meshgrid(x_g_val, y_g_val)

        z_g_sim = np.asarray(g_sim.iloc[0:, 1:], dtype=np.float)
        x_g_sim = np.asarray(g_sim.columns.values[1:], dtype=np.float)
        y_g_sim = np.asarray(g_sim.iloc[:, 0].values, dtype=np.float)

        x_g_sim, y_g_sim = np.meshgrid(x_g_sim, y_g_sim)

        z_g_res = np.asarray(g_res.iloc[0:, 1:], dtype=np.float)
        x_g_res = np.asarray(g_res.columns.values[1:], dtype=np.float)
        y_g_res = np.asarray(g_res.iloc[:, 0].values, dtype=np.float)

        x_g_res, y_g_res = np.meshgrid(x_g_res, y_g_res)

        z_g_sim_hr = np.asarray(g_sim_hr.iloc[0:, 1:], dtype=np.float)
        x_g_sim_hr = np.asarray(g_sim_hr.columns.values[1:], dtype=np.float)
        y_g_sim_hr = np.asarray(g_sim_hr.iloc[:, 0].values, dtype=np.float)

        x_g_sim_hr, y_g_sim_hr = np.meshgrid(x_g_sim_hr, y_g_sim_hr)

        z_g_res_hr = np.asarray(g_res_hr.iloc[0:, 1:], dtype=np.float)
        x_g_res_hr = np.asarray(g_res_hr.columns.values[1:], dtype=np.float)
        y_g_res_hr = np.asarray(g_res_hr.iloc[:, 0].values, dtype=np.float)

        x_g_res_hr, y_g_res_hr = np.meshgrid(x_g_res_hr, y_g_res_hr)

        x_ax_label = 'Length [µm]'
        y_ax_label = 'Time [s]'

        fig, ax = plt.subplots(figsize=(7, 8))
        # ax = fig.add_subplot(1,1,1)
        # ax.set_ylabel('Time [s]')
        # ax.set_xlabel('Length [um]')
        plt.set_cmap(cm.magma)
        fig.gca()

        plt.subplot(3, 2, 1)
        plt.pcolor(x_g_val, y_g_val, z_g_val)
        plt.colorbar()
        plt.ylabel(y_ax_label)
        plt.xlabel(x_ax_label)
        plt.title('Experimental values')

        plt.subplot(3, 2, 3)
        plt.pcolor(x_g_sim, y_g_sim, z_g_sim)
        plt.colorbar()
        plt.ylabel(y_ax_label)
        plt.xlabel(x_ax_label)
        plt.title('1D Heat Equation Simulation')

        plt.subplot(3, 2, 5)
        plt.pcolor(x_g_res, y_g_res, z_g_res, cmap=cm.PuOr_r)
        plt.colorbar()
        plt.ylabel(y_ax_label)
        plt.xlabel(x_ax_label)
        plt.title('Residuals')

        plt.subplot(3, 2, 4)
        plt.pcolor(x_g_sim_hr, y_g_sim_hr, z_g_sim_hr)
        plt.colorbar()
        plt.ylabel(y_ax_label)
        plt.xlabel(x_ax_label)
        plt.title('HR 1D Heat Equation Simulation')

        plt.subplot(3, 2, 6)
        plt.pcolor(x_g_res_hr, y_g_res_hr, z_g_res_hr, cmap=cm.PuOr_r)
        plt.colorbar()
        plt.ylabel(y_ax_label)
        plt.xlabel(x_ax_label)
        plt.title('HR Residuals')

        plt.subplot(3, 2, 2)
        plt.text(0, 0.6, 'Cell name - %s\n' r'$D_t$ - %s $(µm^2)/s$' '\n' r'Error - %s' % (cid, D_coe, D_err), fontsize=14)
        plt.text(0, 0, '\n\nValues make sense only\nif the model fits the data!\nCheck the residuals!\n\nIf the model fits,\nerror is standard deviation.', fontsize=10)
        plt.axis('off')


        plt.tight_layout()
        pdf.savefig()
        plt.savefig(os.path.join(path, '%s-graphs.png' % cid), dpi=300)

        print('%s done!' % cid)

    else:

        print('Files for %s missing!' % cid)

pdf.close()
print('all done!')