import numpy as np
from scipy.interpolate import splrep, splev

def interpolate(interp_voltage):
    #Get the raw data and access columns of data by transposing
    lakeshore_data = np.loadtxt('lakeshore.txt')
    lakeshore_data = np.transpose(lakeshore_data)

    #Find spline representation for temperature vs voltage  
    spl = splrep(lakeshore_data[1][::-1], lakeshore_data[0][::-1])

    #Return temperature interpolated
    return splev(interp_voltage, spl)

#Enter you input voltage here
input_voltage = 1

#Call the interpolate function to get the temperature and its error.
output_temperature = interpolate(input_voltage)
print('The temperature interpolated for a voltage of ', input_voltage, ' is ', output_temperature, '.')