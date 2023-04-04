from pyprojroot.here import here
from scipy.io import savemat
from tests import generate_sine_data, plot_input_data

data = generate_sine_data()
plot_input_data(data)

mat_file = here('data/sine_data.mat')
mat_dict = {'data': data}
savemat(mat_file, mat_dict)
