import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Proj, transform

# Sample data (replace this with your actual data)
data = {
    'lat': [31.0, 32.0, 30.0, 29.0],
    'long': [-100.0, -99.0, -101.0, -102.0]
}
df = pd.DataFrame(data)

# Define the Lambert Conformal Conic projection parameters for Texas
lambert_conformal_conic = Proj(
    proj='lcc',
    lat_1=33.0,  # First standard parallel (near the northern boundary of Texas)
    lat_2=27.0,  # Second standard parallel (near the southern boundary of Texas)
    lat_0=31.0,  # Latitude of origin (central latitude of Texas)
    lon_0=-100.0,  # Central meridian (central longitude of Texas)
    x_0=0,
    y_0=0,
    datum='WGS84'
)

# Function to project latitude and longitude to x, y
def project_coords(row):
    x, y = lambert_conformal_conic(row['long'], row['lat'])
    return pd.Series({'x': x, 'y': y})

# Apply the projection function to each row in the DataFrame
df[['x', 'y']] = df.apply(project_coords, axis=1)

# Adjust the coordinates so that the minimum x and y values start at zero
#df['x'] = df['x'] - df['x'].min()
#df['y'] = df['y'] - df['y'].min()

# Inverse transformation function
def inverse_project_coords(row):
    lon, lat = lambert_conformal_conic(row['x'], row['y'], inverse=True)
    return pd.Series({'lat': lat, 'long': lon})

# Apply the inverse transformation to each row in the DataFrame
df[['inv_lat', 'inv_long']] = df.apply(inverse_project_coords, axis=1)

# Plot the original and inverse-transformed coordinates
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.scatter(df['long'], df['lat'], color='blue', label='Original')
plt.title('Original Coordinates')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.legend()

plt.subplot(1, 2, 2)
plt.scatter(df['inv_long'], df['inv_lat'], color='red', label='Inverse Transformed')
plt.title('Inverse Transformed Coordinates')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.legend()

plt.tight_layout()
plt.show()
