# The lines 350 and 357 both have been modified (.T) this for changes in the teach solution, but in the compute RPY the teach solution not *-1 the euler angles

import numpy as np
from math import sin, cos, tan, radians, sqrt, atan2, degrees, asin
from numpy.linalg import svd
import pyproj

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class GeoSensors():
    """
    Class for georeferencing sensors.
    """

    def __init__(self, eframe_geod_coords:dict, bframe_coords:dict, ins_orientation:dict, angles_bframe:dict, ellipsoid:str= "GRS80", origin:str="REC2", magnetic_declination:float=0.0, geoid_undulation:float=0.0):
        # Initialization of attributes
        self.total_rec = 0
        self.eframe_geod_coords = eframe_geod_coords
        self.bframe_coords = bframe_coords
        self.ins_orientation = ins_orientation
        self.angles_bframe = angles_bframe
        self.ellipsoid = ellipsoid

        # Origin validation
        if origin not in self.eframe_geod_coords.keys() and origin != "barycenter":
            raise ValueError("Invalid origin, please enter a valid receiver id or 'barycenter'")
        
        self.origin = origin
        self.magnetic_declination = magnetic_declination
        self.geoid_undulation = geoid_undulation

        # Initialization of dictionaries to store different coordinate representations
        self.rec_bframe_coords = {}
        self.eframe_geod_coords_dec = {} # e-frame coordinates in decimal
        self.eframe_geod_coords_rad = {} # e-frame coordinates in radians
        self.eframe_ecef_coords = {} # e-frame coordinates in ecef
        self.eframe_local_coords = {} # e-frame coordinates in local, includes origin and subtractions
        self.eframe_nframe_coords = {} # e-frame coordinates in n-frame
        self.reference_point_coords = None
        self.barycenter = {} # barycenter coordinates in radians
        self.eframe_nframe_rotation_matrix = None # rotation matrix e-frame to n-frame
        self.bn_rotation_matrix = None # rotation matrix n-frame to b-frame
        self.ins_rotation_matrix = None # rotation matrix ins to n-frame
        self.en_euler_angles = {} # Rotation angles in e-frame
        self.bn_euler_angles = {} # Rotation angles in b-frame
        self.platform_roll_pitch_yaw = {} # Rotation angles
        self.platform_center = None  # Geodetic coordinates of the platform center

        self.sensors_eframe_coords = {}

        if not set(self.eframe_geod_coords.keys()).issubset(self.bframe_coords.keys()):
            raise ValueError("The e-frame coordinates of the receivers must be provided in the b-frame coordinates.")
        else:
            for rec in self.eframe_geod_coords.keys():
                self.rec_bframe_coords[rec] = self.bframe_coords[rec]

        # Ellipsoid constants
        a = 6378137.0  # Equatorial radius of WGS84 ellipsoid
        b = 6356752.314140356 # Polar radius of WGS84 ellipsoid

        # Calculation of ellipsoid parameters
        f = (a-b)/a  # flattening
        inverse_flattening = 1/f # inverse flattening
        first_eccentricity = ((a**2)-(b**2))**(1/2)/a # first eccentricity
        second_eccentricity = ((a**2)-(b**2))**(1/2)/b # second eccentricity

        # Storage of ellipsoid data
        self.ellipsoid_data = {
            "ellipsoid": ellipsoid,
            "semi-major_axis": a,
            "semi-minor_axis": b,
            "inverse_flattening (1/f)": inverse_flattening,
            "flattening (f)": f,
            "first_eccentricity": first_eccentricity,
            "second_eccentricity": second_eccentricity,
            "eccentricity_squared": 2 * f - f ** 2
            }
        
        self._handle_process() # executes the complete process

    def __repr__(self):
        # String representation of the class
        return f'''
        - GeoSensors(eframe_geod_coords=\n{self.eframe_geod_coords}, \n
        - bframe_coords=\n{self.bframe_coords}, \n
        - ins_orientation=\n{self.ins_orientation}, \n
        - ellipsoid={self.ellipsoid}, \n
        - origin={self.origin}, \n
        - magnetic_declination={self.magnetic_declination}, \n
        - geoid_undulation={self.geoid_undulation})
        '''
        

    def _dms_to_decimal(self, deg, min, sec, direction):
        """
        Converts degrees, minutes, seconds to decimal degrees.
        :param deg: Degrees
        :param min: Minutes
        :param sec: Seconds
        :param direction: Direction ('N', 'S', 'E', 'W')
        :return: Value in decimal degrees, with correct sign.
        """
        # Conversion to decimal
        decimal = deg + min / 60 + sec / 3600
        
        # Apply sign according to direction
        if direction == 'S' or direction == 'W':
            decimal = -decimal
        
        return decimal

    def _geodetic_to_ecef(self):
        """
        Converts geodetic coordinates from DMS to decimal degrees.
        :param rec_geod_coords: Dictionary with coordinates in DMS format
        :return: Dictionary with coordinates in [latitude, longitude, height] format
        """
        # Convert geodetic coordinates to decimal coordinates
        rec_geod_coords = self.eframe_geod_coords
        count = 0
        coord_radians= []
        
        # Iteration over the coordinates of each receiver
        for rec, coords in rec_geod_coords.items():
            lat_deg, lat_min, lat_sec, lat_dir = coords[:4]  # Latitude part
            lon_deg, lon_min, lon_sec, lon_dir = coords[4:8]  # Longitude part
            height = coords[8]  # Ellipsoidal height
            
            # Convert latitude and longitude to decimal degrees
            latitude = self._dms_to_decimal(lat_deg, lat_min, lat_sec, lat_dir)
            longitude = self._dms_to_decimal(lon_deg, lon_min, lon_sec, lon_dir)
            
            # Save the result in decimal and radians
            self.eframe_geod_coords_dec[rec] = [latitude, longitude, height]
            count += 1

            coord_radians = [radians(latitude), radians(longitude)]
            self.eframe_geod_coords_rad[rec]={ 'coord': coord_radians, 'height': height }

        # Convert geodetic coordinates to ECEF coordinates
        a = self.ellipsoid_data["semi-major_axis"]
        b = self.ellipsoid_data["semi-minor_axis"]
        eccentricity = sqrt(self.ellipsoid_data["eccentricity_squared"])
        for sensor, coords in self.eframe_geod_coords_rad.items():
            lat_rad = coords['coord'][0]
            long_rad = coords['coord'][1]
            h = coords['height']
            
            # Radius of curvature of the first vertical
            v = a/sqrt(1-(eccentricity**2*sin(lat_rad)**2))
            
            # Conversion from geodetic to cartesian
            X = (v+h)*cos(lat_rad)*cos(long_rad)
            Y = (v+h)*cos(lat_rad)*sin(long_rad)
            Z = (v*(1-eccentricity**2)+h)*sin(lat_rad)

            # Storage of cartesian coordinates
            self.eframe_ecef_coords[sensor] = [X,Y,Z]

        self.total_rec = count

    def _decimal_to_dms(self, decimal):
        """
        Converts decimal coordinates to degrees, minutes and seconds.
        """
        # Extraction of degrees, minutes and seconds
        degrees = int(decimal)
        minutes_dec = (decimal - degrees) * 60
        minutes = int(minutes_dec)
        seconds = (minutes_dec - minutes) * 60
        return [degrees, minutes, seconds]

    def _dec_to_geodetic(self, platform_coords_dec, ellipsoidal_height):
        """
        Converts ECEF coordinates to geodetic coordinates for all receivers. (used in the calculation of the platform center)
        """
        lat_dec, lon_dec = platform_coords_dec
        alt = ellipsoidal_height

        # Conversion of decimal coordinates to degrees, minutes and seconds
        lat_dms = self._decimal_to_dms(abs(lat_dec))
        lon_dms = self._decimal_to_dms(abs(lon_dec))

        # Determination of direction (N/S, E/W)
        lat_dir = 'N' if lat_dec >= 0 else 'S'
        lon_dir = 'E' if lon_dec >= 0 else 'W'

        return [lat_dms, lon_dms, alt]

    def _set_origin(self):
        """
        Defines the origin (0,0) of the e-frame at a point in the b-frame.
        It can be the barycenter or a specific receiver.
        The coordinates of the receivers take value with respect to the origin.
        """
        if self.origin != "barycenter":
            # If the origin is not the barycenter, a specific receiver is used
            if self.origin not in self.eframe_ecef_coords.keys():
                raise ValueError("Invalid origin")
            else:
                rec_ref = self.eframe_ecef_coords[self.origin]
                self.reference_point_coords = {'name':self.origin, 'ecef':rec_ref, 'geodetic':self.eframe_geod_coords[self.origin], 'decimal':self.eframe_geod_coords_dec[self.origin]}
                adjusted_coords = {}
                # Adjustment of coordinates with respect to the origin
                for rec, ecef in self.eframe_ecef_coords.items():
                    adjusted_coords[rec] = []
                    for i in range(3):
                        adjusted_coords[rec].append(ecef[i] - rec_ref[i])
                
                self.eframe_local_coords['coords'] = adjusted_coords
                self.eframe_local_coords['origin'] = self.origin

        elif self.origin == "barycenter":
            # If the origin is the barycenter, it is calculated
            coord_bar_ecef, coord_bar_rad, coord_bar_dec = self._calculate_barycenter()

            self.reference_point_coords = {'name':'barycenter', 'ecef':coord_bar_ecef, 'geodetic':[], 'decimal':coord_bar_dec}
        else:
            raise ValueError("Invalid origin")
        
    def _calculate_barycenter(self):
        """
        Calculates the barycenter of the receivers.
        """
        # Calculation of the barycenter in ECEF coordinates
        x_mean, y_mean, z_mean = 0, 0, 0
        for rec, ecef in self.eframe_ecef_coords.items():
            x_mean += ecef[0]
            y_mean += ecef[1]
            z_mean += ecef[2]
        x_mean /= self.total_rec
        y_mean /= self.total_rec
        z_mean /= self.total_rec

        coord_bar = [x_mean, y_mean, z_mean]

        # Adjustment of coordinates with respect to the barycenter
        adjusted_coords = {}
        for rec, ecef in self.eframe_ecef_coords.items():
            adjusted_coords[rec] = [ecef[0] - coord_bar[0], ecef[1] - coord_bar[1], ecef[2] - coord_bar[2]]

        # Calculation of the barycenter in geodetic coordinates
        x_rad_mean, y_rad_mean, z_rad_mean = 0, 0, 0
        for rec, coord in self.eframe_geod_coords_rad.items():
            x_rad, y_rad = coord['coord']
            z_rad = coord['height']
            x_rad_mean += x_rad
            y_rad_mean += y_rad
            z_rad_mean += z_rad
        x_rad_mean /= len(self.eframe_geod_coords_rad)
        y_rad_mean /= len(self.eframe_geod_coords_rad)
        z_rad_mean /= len(self.eframe_geod_coords_rad)

        coord_bar_rad = [x_rad_mean, y_rad_mean]
        coord_bar_dec = [degrees(x_rad_mean), degrees(y_rad_mean)]
        ellipsoidal_height = z_rad_mean
        

        # Conversion to degrees, minutes and seconds
        dms = self._dec_to_geodetic(coord_bar_dec, ellipsoidal_height)
        utm = self._geodetic_to_utm(lat=coord_bar_dec[0], lon=coord_bar_dec[1])

        coord = {'radians': coord_bar_rad, 'decimals': coord_bar_dec, 'dms': dms, 'utm': utm}

        self.barycenter = {'coord': coord, 'height': ellipsoidal_height}

        self.eframe_local_coords['coords'] = adjusted_coords
        self.eframe_local_coords['origin'] = self.origin

        return [coord_bar, coord_bar_rad, coord_bar_dec]
        
    def _eframe_to_nframe(self):
        """
        Applies the Jekeli transformation to calculate the e-frame coordinates expressed in n-frame.
        """
        # Obtaining latitude and longitude of the origin
        if self.origin != "barycenter":
            lat_rad, lon_rad = self.eframe_geod_coords_rad[self.origin]['coord']

        else:
            lat_rad = self.barycenter['coord']['radians'][0]
            lon_rad = self.barycenter['coord']['radians'][1]

        # Calculation of the rotation matrix (transposed)
        R = np.array([[-sin(lat_rad) * cos(lon_rad), -sin(lat_rad) * sin(lon_rad), cos(lat_rad)],
                        [-sin(lon_rad), cos(lon_rad), 0],
                        [-cos(lat_rad) * cos(lon_rad), -cos(lat_rad) * sin(lon_rad), -sin(lat_rad)]])
        
        self.eframe_nframe_rotation_matrix = R

        # Application of the matrix to the coordinates of the receivers
        coord_eframe = []
        for rec, coord in self.eframe_local_coords['coords'].items():
            coord_eframe.append(coord)
        coord_eframe = np.array(coord_eframe)
        coord_nframe = np.dot(R, np.transpose(coord_eframe))
        self.eframe_nframe_coords = np.transpose(coord_nframe)

        # Calculation of Euler angles
        self.en_euler_angles = self._calculate_euler_angles(R)

    def _calculate_euler_angles(self, matrix):
        matrix = matrix
                
        # Calculation of Euler angles
        a1 = atan2(-matrix[2][1],matrix[2][2])
        a2 = asin(matrix[2][0])
        a3 = atan2(-matrix[1][0],matrix[0][0])

        angles = {'coord': {'radians': [a1, a2, a3], 'decimals': [degrees(a1), degrees(a2), degrees(a3)]}}

        return angles     

    def _calculate_rotation_matrix_bn(self):
        """
        Calculates the rotation matrix for the transformation of coordinates from b-frame to n-frame starting from the knowledge of the b (local) and n (ecef) coordinates -> In local, solving the matrix using SVD = Singular Value Decomposition.
        """
        # Formation of coordinate matrices of b and n
        matrix_b = [] 
        for rec, coord in self.rec_bframe_coords.items():
            matrix_b.append(coord)
        matrix_b = np.array(matrix_b)

        matrix_n = []
        for coord in self.eframe_nframe_coords:
            matrix_n.append(coord)
        matrix_n = np.array(matrix_n)
        
        # Calculation of rotation matrix using SVD
        # 1. Calculation of centroids of b and n
        b_centroid = np.mean(matrix_b, axis=0)
        n_centroid = np.mean(matrix_n, axis=0)

        # Subtract the centroids from the coordinates of b and n
        coor_n_centered = n_centroid - matrix_n
        coor_b_centered = b_centroid - matrix_b

        # 2. Correlation matrix
        H = np.dot(np.transpose(coor_n_centered), coor_b_centered)

        # 3. Singular Value Decomposition (SVD)
        U, S, Vt = svd(H)

        # 4. Rotation matrix
        R = np.dot(U, Vt)
        
        # Storage of the rotation matrix
        self.bn_rotation_matrix = R.T # modified (.T)

    def _calculate_rpy(self):
        """
        Calculates the rotation angles roll, pitch and yaw.
        """
        matrix = self.bn_rotation_matrix
        self.bn_euler_angles = self._calculate_euler_angles(matrix.T) # modified (.T)

        # Calculation of roll, pitch and yaw
        roll = self.bn_euler_angles['coord']['decimals'][0]
        pitch = self.bn_euler_angles['coord']['decimals'][1]
        yaw = self.bn_euler_angles['coord']['decimals'][2]

        self.platform_roll_pitch_yaw = {'roll': roll, 'pitch': pitch, 'yaw': yaw}

    def _calculate_sensor_coordinates(self, angle_source="calculate"):
        """
        Calculate coordinates with angles source defined in the parameter.
        
        :param angle_source: "calculate" or "input"
        """

        angle_source = angle_source

        if angle_source == "calculate":
            # Use calculated angles
            alpha, beta, gamma = self.angles_bframe['PLA']

            alpha_rad = radians(alpha)
            beta_rad = radians(beta)
            gamma_rad = radians(gamma)

            R_sb = self._calculate_rotation_matrix(alpha_rad, beta_rad, gamma_rad)

            R_bn = self.bn_rotation_matrix
            R_bn = R_bn.T

        else:  # angle_source == "input"
            # Use alpha, beta, gamma angles provided by the INS
            alpha,  beta, gamma = self.angles_bframe['INS']

            alpha_rad = radians(alpha)
            beta_rad = radians(beta)
            gamma_rad = radians(gamma)

            R_sb = self._calculate_rotation_matrix(alpha_rad, beta_rad, gamma_rad)

            # Use angles from input
            roll = self.ins_orientation['roll']
            pitch = self.ins_orientation['pitch']
            yaw = self.ins_orientation['yaw']

            # Adjust yaw with magnetic declination
            yaw += self.magnetic_declination

            # Convert angles to radians
            roll_rad, pitch_rad, yaw_rad = map(radians, [roll, pitch, yaw])

            # Calculate rotation matrix from b-frame to n-frame
            R_bn = self._calculate_rotation_matrix(roll_rad, pitch_rad, yaw_rad)
            self.ins_rotation_matrix = R_bn.T        

        # Transpose of the e-frame to n-frame rotation matrix
        R_ne = self.eframe_nframe_rotation_matrix.T

        # Combine rotation matrices
        R_be = np.dot(R_ne, np.dot(R_bn.T, R_sb))

        # Calculate sensor coordinates
        for sensor, coord_bframe in self.bframe_coords.items():
            coord_bframe_array = np.array(coord_bframe)
            coord_eframe_local = np.dot(R_be, coord_bframe_array)
            coord_eframe_ecef = coord_eframe_local + np.array(self.reference_point_coords['ecef'])  
            self.sensors_eframe_coords[sensor] = coord_eframe_ecef.tolist()

        # Convert ECEF coordinates to geodetic UTM
        for sensor, coord_ecef in self.sensors_eframe_coords.items():
            lat_dec, lon_dec, alt = self._ecef_to_geodetic(coord_ecef)
            lat_dms = self._decimal_to_dms(lat_dec)
            lon_dms = self._decimal_to_dms(lon_dec)
            utm = self._geodetic_to_utm(lat_dec, lon_dec)
            h_orto = alt - self.geoid_undulation
            
            self.sensors_eframe_coords[sensor] = {
                'ecef': coord_ecef,
                'geodetic': [lat_dms, lon_dms, alt],
                'utm': utm,
                'ortometric_height': h_orto
            }

        print(f"Calculated sensor coordinates using {angle_source} angles.")

        if angle_source == "calculate":
            self._plot_euler_angles(title="Euler Angles with Calculated Angles")
        else:
            self._plot_euler_angles(title="Euler Angles with Input Angles")

    def _calculate_rotation_matrix(self, roll, pitch, yaw):
        """
        Rotation matrix from b-frame to n-frame using Euler angles.
        """
        # Rotation matrices for each axis
        R_x = np.array([[1, 0, 0],
                        [0, cos(roll), -sin(roll)],
                        [0, sin(roll), cos(roll)]])
        
        R_y = np.array([[cos(pitch), 0, sin(pitch)],
                        [0, 1, 0],
                        [-sin(pitch), 0, cos(pitch)]])
        
        R_z = np.array([[cos(yaw), -sin(yaw), 0],
                        [sin(yaw), cos(yaw), 0],
                        [0, 0, 1]])
        
        # Complete rotation matrix
        R_bn = np.dot(R_z, np.dot(R_y, R_x))
        
        return R_bn

    def _ecef_to_geodetic(self, ecef):
        """
        Convert ECEF coordinates to geodetic coordinates (latitude, longitude, height).
        """
        # Implement the conversion from ECEF to geodetic here
        # This is a simplified function; it is recommended to use a library like pyproj for greater accuracy
        x, y, z = ecef
        a = self.ellipsoid_data['semi-major_axis']
        e = self.ellipsoid_data['first_eccentricity']
        
        lon = atan2(y, x)
        p = sqrt(x**2 + y**2)
        lat = atan2(z, p * (1 - e**2))
        
        for _ in range(5):  # Iterations to improve accuracy
            N = a / sqrt(1 - e**2 * sin(lat)**2)
            h = p / cos(lat) - N
            lat = atan2(z, p * (1 - e**2 * N / (N + h)))
        
        return degrees(lat), degrees(lon), h

    def _geodetic_to_utm(self, lat, lon):
        """
        Convert geodetic coordinates to UTM coordinates.
        """
        # Implement the conversion from geodetic to UTM here
        # It is recommended to use a library like pyproj for this conversion
        utm_proj = pyproj.Proj(proj='utm', zone=self._get_utm_zone(lon), ellps='GRS80')
        easting, northing = utm_proj(lon, lat)
        return easting, northing

    def _get_utm_zone(self, lon):
        """
        Determine the UTM zone based on longitude.
        """
        return int((lon + 180) / 6) + 1

    def _handle_process(self):
        """
        Handles automatic transformations.
        """
        # Sequential execution of transformation processes
        print("Converting geodetic coordinates to ECEF...")
        self._geodetic_to_ecef()
        print("Setting origin...")
        self._set_origin()
        print("Transforming from e-frame to n-frame...")
        self._eframe_to_nframe()
        print("Calculating barycenter...")
        self._calculate_barycenter()
        print("Calculating rotation matrix from b-frame to n-frame...")
        self._calculate_rotation_matrix_bn()
        print("Calculating roll, pitch and yaw...")
        self._calculate_rpy()
        print("Calculating sensor coordinates...")
        self._calculate_sensor_coordinates()
        self._calculate_sensor_coordinates(angle_source='input')
        

    def generate_report(self, filename = "report.txt"):
        """
        Generates a report with the coordinates of the receivers in geodetic and cartesian coordinates.
        """
        with open(filename, "w") as file:
            file.write("Student: Cristhian Andres Quiza - MUIGG\n\n Practice 1\nSteps\n")

            file.write("Sensor Georeferencing Report\n")
            file.write("=======================================\n\n")
            
            file.write("0. Input Data\n")
            file.write("-------------------\n")
            file.write("E-frame geodetic coordinates:\n")
            for key, value in self.eframe_geod_coords.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")
            
            file.write("B-frame coordinates:\n")
            for key, value in self.rec_bframe_coords.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")

            file.write("B-frame angles:\n")
            for key, value in self.angles_bframe.items():
                file.write(f"{key}:{value}\n")
            file.write("\n")
            
            file.write("INS orientation:\n")
            for key, value in self.ins_orientation.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")
            file.write(f"Ellipsoid: {self.ellipsoid}\n")
            file.write(f"Origin: {self.origin}\n")
            file.write(f"Magnetic declination: {self.magnetic_declination}\n")
            file.write(f"Geoid undulation: {self.geoid_undulation}\n\n")

            file.write("1. Ellipsoid Data\n")
            file.write("-----------------------\n")
            for key, value in self.ellipsoid_data.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")

            file.write("2. Decimal Geodetic Coordinates\n")
            file.write("-----------------------------------\n")
            for receptor, coords in self.eframe_geod_coords_dec.items():
                file.write(f"{receptor}: {coords}\n")
            file.write("\n")

            file.write("3. ECEF Coordinates\n")
            file.write("-------------------\n")
            for receptor, coords in self.eframe_ecef_coords.items():
                file.write(f"{receptor}: {coords}\n")
            file.write("\n")

            file.write("4. Local Coordinates\n")
            file.write("-----------------------\n")
            file.write(f"Origin: {self.eframe_local_coords['origin']}\n")
            for receptor, coords in self.eframe_local_coords['coords'].items():
                file.write(f"{receptor}: {coords}\n")
            file.write("\n")

            file.write("5. E-frame to N-frame Rotation Matrix\n")
            file.write("---------------------------------------\n")
            file.write(str(self.eframe_nframe_rotation_matrix))
            file.write("\n\n")

            file.write("6. N-frame Coordinates\n")
            file.write("-------------------------\n")
            for coords in self.eframe_nframe_coords:
                file.write(f"{coords}\n")
            file.write("\n")

            file.write("7. E-N Euler Angles\n")
            file.write("-----------------------\n")
            file.write(f"Radians: {self.en_euler_angles['coord']['radians']}\n")
            file.write(f"Decimals: {self.en_euler_angles['coord']['decimals']}\n")
            file.write("\n")

            file.write("8. Barycenter -----------> Practice 1.1\n")
            file.write("--------------\n")
            file.write(f"Radians: {self.barycenter['coord']['radians']}\n")
            file.write(f"Decimals: {self.barycenter['coord']['decimals']}\n")
            file.write(f"DMS: {self.barycenter['coord']['dms']}\n")
            file.write(f"utm: {self.barycenter['coord']['utm']}\n")
            file.write(f"Height: {self.barycenter['height']}\n")            
            file.write("\n")

            file.write("9. B-N Rotation Matrix\n")
            file.write("-------------------------\n")
            file.write(str(self.bn_rotation_matrix.T))
            file.write("\n\n")

            file.write("10. B-N Euler Angles\n")
            file.write("------------------------\n")
            file.write(f"Radians: {self.bn_euler_angles['coord']['radians']}\n")
            file.write(f"Decimals: {self.bn_euler_angles['coord']['decimals']}\n")
            file.write("\n")

            file.write("11. Platform Roll, Pitch, Yaw ------------> Practice 1.2\n")
            file.write("-------------------------------------\n")
            file.write(f"Roll: {self.platform_roll_pitch_yaw['roll']}\n")
            file.write(f"Pitch: {self.platform_roll_pitch_yaw['pitch']}\n")
            file.write(f"Yaw: {self.platform_roll_pitch_yaw['yaw']}\n\n")

            file.write("12. Coords of sensors in e-frame ---------------> Practice 2.1 y 2.2\n")
            file.write("------------------------------------------\n")
            file.write("With calculated angles:\nR_bn:\n")
            file.write(str(self.bn_rotation_matrix.T))
            file.write("\n\nR_ne:\n")
            file.write(str(self.eframe_nframe_rotation_matrix.T))
            file.write("\n\n")
            for sensor, coords in self.sensors_eframe_coords.items():
                file.write(f"{sensor}:\n")
                for coord_type, values in coords.items():
                    file.write(f"  {coord_type}: {values}\n")
            file.write("\n\n")
            
            file.write("With input angles:\nR_bn:\n")
            file.write(str(self.ins_rotation_matrix.T))
            file.write("\n\nR_ne:\n")
            file.write(str(self.eframe_nframe_rotation_matrix))
            file.write("\n\n")
            for sensor, coords in self.sensors_eframe_coords.items():
                file.write(f"{sensor}:\n")
                for coord_type, values in coords.items():
                    file.write(f"  {coord_type}: {values}\n")
            file.write("\n\n")

            file.write("\n****************************************************\n")

        print("\n\nReport generated successfully.")

    def _plot_euler_angles(self, title):
        """
        Generates a 3D plot of the Euler angles (roll, pitch, yaw) and saves it as a PDF.
        """
        title = title
        # Get the Euler angles
        if title == "Euler Angles with Calculated Angles":
            roll, pitch, yaw = self.platform_roll_pitch_yaw.values()
            pdf_filename = 'euler_angles_plot_calculated.pdf'
        else:
            roll = self.ins_orientation['roll']
            pitch = self.ins_orientation['pitch']
            yaw = self.ins_orientation['yaw']
            yaw += self.magnetic_declination
            pdf_filename = 'euler_angles_plot_input.pdf'

        # Convert to radians
        roll, pitch, yaw = map(np.radians, [roll, pitch, yaw])

        # Create the figure and 3D axis
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Create vectors for the axes
        x = np.array([1, 0, 0])
        y = np.array([0, 1, 0])
        z = np.array([0, 0, 1])

        # Rotation matrices
        Rx = np.array([[1, 0, 0],
                    [0, np.cos(roll), -np.sin(roll)],
                    [0, np.sin(roll), np.cos(roll)]])

        Ry = np.array([[np.cos(pitch), 0, np.sin(pitch)],
                    [0, 1, 0],
                    [-np.sin(pitch), 0, np.cos(pitch)]])

        Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                    [np.sin(yaw), np.cos(yaw), 0],
                    [0, 0, 1]])

        # Apply rotations
        x_rotated = Rz.dot(Ry.dot(Rx.dot(x)))
        y_rotated = Rz.dot(Ry.dot(Rx.dot(y)))
        z_rotated = Rz.dot(Ry.dot(Rx.dot(z)))

        # Plot the original axes
        ax.quiver(0, 0, 0, x[0], x[1], x[2], color='r', label='Original X')
        ax.quiver(0, 0, 0, y[0], y[1], y[2], color='g', label='Original Y')
        ax.quiver(0, 0, 0, z[0], z[1], z[2], color='b', label='Original Z')

        # Plot the rotated axes
        ax.quiver(0, 0, 0, x_rotated[0], x_rotated[1], x_rotated[2], color='m', label='Rotated X')
        ax.quiver(0, 0, 0, y_rotated[0], y_rotated[1], y_rotated[2], color='y', label='Rotated Y')
        ax.quiver(0, 0, 0, z_rotated[0], z_rotated[1], z_rotated[2], color='c', label='Rotated Z')

        # Configure the plot
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'{title}\nEuler Angles: Roll={np.degrees(roll):.2f}°, Pitch={np.degrees(pitch):.2f}°, Yaw={np.degrees(yaw):.2f}°')
        ax.legend()

        plt.tight_layout()
        
        # Save the plot as a PDF
        plt.savefig(pdf_filename, format='pdf')
        plt.close(fig)
        
        print(f"Saved graphic of Euler angles as: {pdf_filename}")


def practice1():

    bframe_coords = {
        "PLA":[0.0000,0.0000,0.0000],
        "INS":[0.0000,0.0000,-0.0250],
        "REC2":[0.4000,-0.4000,-0.0150],
        "REC3":[0.4000,0.4000,-0.0150],
        "REC4":[-0.4000,-0.4000,-0.0150],
        "REC5":[-0.4000, 0.4000,-0.0150],
        "CAM1":[ 0.7000,-0.1500,-0.0800],
        "CAM2":[ 0.7000, 0.1500,-0.0800]
    } 

    angles_bframe = {
        "PLA":[0.0000,0.0000,0.0000],
        "INS":[0.263568,-0.437974,2.869214],
        "REC2":[0.0000,0.0000,0.0000],
        "REC3":[0.0000,0.0000,0.0000],
        "REC4":[0.0000,0.0000,0.0000],
        "REC5":[0.0000,0.0000,0.0000],
        "CAM1":[0.196543,0.541674,-1.355847],
        "CAM2":[-0.287495,0.426588,1.128792]
    }

    eframe_geod_coords = {
        "REC2":[39,28,53.133637,"N",0,20,7.650964,"W", 55.8670],
        "REC3":[39,28,53.109885,"N",0,20,7.663129,"W", 55.8290],
        "REC4":[39,28,53.143535,"N",0,20,7.680928,"W", 55.8830],
        "REC5":[39,28,53.119474,"N",0,20,7.693816,"W", 55.8660]
    }

    ins_orientation = {'roll': 1.195278, 'pitch': -1.189671, 'yaw': 111.269874}

    magnetic_declination = -0.700
    geoid_undulation = 50.019

    geoReference = GeoSensors(eframe_geod_coords, bframe_coords, ins_orientation, angles_bframe = angles_bframe, origin="barycenter", magnetic_declination=magnetic_declination, geoid_undulation=geoid_undulation)

    geoReference.generate_report()

if __name__ == __name__:
    practice1()