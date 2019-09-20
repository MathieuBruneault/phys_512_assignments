import numpy as np
from scipy.interpolate import splrep, splev

def interpolate(interp_voltage):
    #Get the raw data and access columns of data by transposing
    lakeshore_data = np.loadtxt('lakeshore.txt')
    lakeshore_data = np.transpose(lakeshore_data)

    #Find spline representation for temperature vs voltage  
    spl = splrep(lakeshore_data[1][::-1], lakeshore_data[0][::-1])

    #Interpolate temperature
    output_temperature = splev(interp_voltage, spl)

    #Find the nearest voltage to our desired voltage for which we know the temperature, its index and corresponding temperature
    closest_voltage_idx = (np.abs(lakeshore_data[1] - interp_voltage)).argmin()
    closest_voltage = lakeshore_data[1][closest_voltage_idx]
    actual_temp = lakeshore_data[0][closest_voltage_idx]

    #Perform a new interpolation for this voltage value, after removing that point from the data
    new_volt_data = np.delete(lakeshore_data[1],closest_voltage_idx)
    new_temp_data = np.delete(lakeshore_data[0],closest_voltage_idx)
    spl = splrep(new_volt_data[::-1], new_temp_data[::-1])
    err_temperature = splev(closest_voltage, spl)

    #Estimate the error as being the actual error from that second interpolation
    err = np.abs(actual_temp - err_temperature)

    return output_temperature, err

#Enter you input voltage here
input_voltage = 1

#Call the interpolate function to get the temperature and its error.
output_temperature, err = interpolate(input_voltage)
print('The temperature interpolated for a voltage of ', input_voltage, ' is ', output_temperature, ' with an error of ', err)