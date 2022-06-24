import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import os
import math
import re

initial_diffusion_coef = 1e-13
#path = 'g:\\My Drive\\Data\\PyCharmProjects\\New FRAP analysis software\\20181212_002\\test metabolic\\20181206\\'
path = os.getcwd()
# https://pycav.readthedocs.io/en/latest/api/pde/crank_nicolson.html
# https://www.quantstart.com/articles/Crank-Nicholson-Implicit-Scheme
def crank_nicolson(initial_heat: np.array, delta_x: float, delta_t: float, diffusion_coef: float):

    A = np.zeros([len(initial_heat), len(initial_heat)])

    r = diffusion_coef * delta_t / (delta_x * delta_x)
    a = -r / 2.0
    b = 1.0 + r
    c = -r / 2.0

    np.fill_diagonal(A[:,1:], a)
    np.fill_diagonal(A, b)
    np.fill_diagonal(A[1:], c)

    A[0, 1] *= 2
    A[-1, -2] *= 2

    D = np.zeros(len(initial_heat))
    for i in range(1, len(initial_heat) - 1):
        D[i] = r / 2.0 * initial_heat[i - 1] + (1.0 - r) * initial_heat[i] + r / 2.0 * initial_heat[i + 1]
    D[0] = (1.0 - r) * initial_heat[0] + r * initial_heat[1]
    D[-1] = (1.0 - r) * initial_heat[-1] + r * initial_heat[-2]

    return np.linalg.solve(A, D)


def simulation(initial_heat, delta_x, timestamps, diffusion_coef):

    after = initial_heat.copy()
    simulated_values = after

    # calculate delta t's
    delta_ts = [t1 - t0 for t0, t1 in zip(timestamps[:-1], timestamps[1:])]

    for delta_t in delta_ts:
        # delta_x from micro meter to meter
        after = crank_nicolson(after, delta_x / 1e6, delta_t, diffusion_coef)
        simulated_values.extend(after.tolist())

    return simulated_values



def open_csv(path):
    with open(path, 'r') as f:
        lines = f.readlines()
        delta_x = float(lines[0].split(',')[2])
        timestamps = [float(line.split(',')[0]) for line in lines[1:]]
        values = [[float(col) for col in line.strip().split(',')[1:]] for line in lines[1:]]

        # timestamps start from 0
        timestamps = [t - timestamps[0] for t in timestamps]
        return delta_x, timestamps, values


def fit(values, delta_x, timestamps):

    # flatten all values (necessary for curve_fit function)
    all_values = [item for sublist in values for item in sublist]

    f = lambda xdata, *params: simulation(values[0], delta_x, xdata, params[0])
    popt, pcov = curve_fit(f, timestamps, all_values, (initial_diffusion_coef,))
    perr = np.sqrt(np.diag(pcov))

    return popt[0], perr


# grows a list to double its size
# [1,2,3,6] => [1,1.5,2,2.5,3,4.5,6]
def grow(l):
    for a, b in zip(l[:-1], l[1:]):
        yield a
        yield (a + b) / 2
    yield l[-1]


def determine_residuals(values1, values2):
    residuals = list()

    for row1, row2 in (zip(values1, values2)):
        residuals.append([value1 - value2 for value1, value2 in zip(row1, row2)])

    return residuals


def write_image(path, values, delta_x, timestamps):
    table = [['time'] + [delta_x * i for i in range(len(values[0]))]]

    for timestamp, row in zip(timestamps, values):
        table.append([timestamp] + row)
    write_values(path, table)


def write_values(path, values):
    with open(path, 'w') as f:
        table = '\n'.join([','.join([str(value) for value in row]) for row in values])
        f.write(table)


def flatten_list(l):
    return [item for sublist in l for item in sublist]


def to_matrix(l, w):
    return [l[offset:offset + w] for offset in range(0, len(l), w)]


def process(path):
    print('process %s' % path)

    delta_x, timestamps, values = open_csv(path)
    diffusion_coef, diffusion_coef_error = fit(values, delta_x, timestamps)

    simulated = simulation(values[0], delta_x, timestamps, diffusion_coef)
    residuals = determine_residuals(values, to_matrix(simulated, len(values[0])))

    # convert flattend list back to matrix
    simulated = to_matrix(simulated, len(values[0]))

    # create high resolution simulation
    high_res_values = list(grow(values[0]))
    range_t = timestamps[-1] - timestamps[0]

    # take all delta t values
    high_res_delta_ts = [t1 - t0 for t0, t1 in zip(timestamps[:-1], timestamps[1:])]

    # discard the first 5 values
    high_res_delta_ts = high_res_delta_ts[5:]
    # sort all delta_t's in ascending order
    high_res_delta_ts = sorted(high_res_delta_ts)
    # take the middle value (median)
    high_res_delta_t = high_res_delta_ts[len(high_res_delta_ts) // 2]

    steps = int(math.ceil(range_t / high_res_delta_t))
    high_res_timestamps = [i * high_res_delta_t for i in range(len(values))]

    high_res_delta_x = delta_x / 2                              # we doubled the resolution (see grow function)
    high_res = simulation(high_res_values, high_res_delta_x, high_res_timestamps, diffusion_coef)

    high_res = to_matrix(high_res, len(high_res_values))

    # upscaled residuals
    high_res_values = [list(grow(row)) for row in values]
    high_res_residuals = determine_residuals(high_res_values, high_res)


 

    # save to csv files
    base_path = path[:-11]
    cell_id = re.split(r'\\', os.path.basename(base_path))[-1]
    simulation_path = os.path.join(path, base_path + '-simulated.csv')
    residuals_path = os.path.join(path, base_path + '-residuals.csv')
    high_res_path = os.path.join(path, base_path + '-simulated-high-res.csv')
    high_res_residuals_path = os.path.join(path, base_path + '-residuals-high-res.csv')

    write_image(simulation_path, simulated, delta_x, timestamps)
    write_image(residuals_path, residuals, delta_x, timestamps)
    write_image(high_res_path, high_res, high_res_delta_x, high_res_timestamps)
    write_image(high_res_residuals_path, high_res_residuals, high_res_delta_x, high_res_timestamps)

    return diffusion_coef, diffusion_coef_error, simulation_path, residuals_path, high_res_path, high_res_delta_x, high_res_delta_t, cell_id


def process_all(path):
    results = []

    for filename in os.listdir(path):
        full_path = os.path.join(path, filename)

        if os.path.isdir(full_path):
            process_all(full_path)
        elif filename.endswith('-values.csv'):
            results.append(process(full_path))

    table = [['cell_id', 'diffusion_coefficient', 'diffusion_coefficient_error', 'simulation', 'residuals', 'high_res', 'high_res_delta_x', 'high_res_delta_t']]
    for diffusion_coef, diffusion_coef_error, simulation_path, residuals_path, high_res_path, delta_x, delta_t, cell_id in results:
            table.append([cell_id, '%g' % diffusion_coef, '%g' % diffusion_coef_error, simulation_path, residuals_path, high_res_path, str(delta_x), str(delta_t)])

    output_path = os.path.join(path, 'diffusion-coefficients.csv')
    write_values(output_path, table)


process_all(path)
