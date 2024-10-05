import netCDF4 as nc

# Open NetCDF file
ncfile = nc.Dataset("gebco_2023_n27.5_s12.5_w125.5_e140.5.nc", 'r')

# Get list of variable names
variable_names = ncfile.variables.keys()

# Print variable names
for var_name in variable_names:
    print(var_name)

lat_variable = ncfile.variables['lat'][:]
lon_variable = ncfile.variables['lon'][:]
elevation_variable = ncfile.variables['elevation'][:]

for x in range(10):
    print(lat_variable[x], lon_variable[x], elevation_variable[x].shape)
    x += 1

# Close the file
ncfile.close()
