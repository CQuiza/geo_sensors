# geo_sensors

## Clase utilizada para la transformación de coordenadas a traves de sistemas de referencia con el uso de matrices de rotación, implementación para sensores de imagen, angulares y gps.

## Generalidades

La clase `GeoSensors` se utiliza para la georreferenciación de sensores mediante la transformación de coordenadas a través de sistemas de referencia. Esta implementación es adecuada para sensores de imagen, angulares y GPS, y utiliza matrices de rotación para realizar las conversiones necesarias.

## Modo de Uso

**Advertencia:** Asegúrate de contar con las dependencias de python necesarias para el funcionamiento de la clase `GeoSensors`. Puedes instalarlas utilizando pip con el siguiente comando:
-  `pip install numpy`
- `pip install pyproj`
- `pip install matplotlib`

Para utilizar la clase `GeoSensors`, primero debes crear una instancia de la misma proporcionando los siguientes parámetros:

- `eframe_geod_coords`: Un diccionario que contiene las coordenadas geodésicas en el marco de referencia E.
- `bframe_coords`: Un diccionario que contiene las coordenadas en el marco de referencia B.
- `ins_orientation`: Un diccionario que describe la orientación del sistema de navegación inercial (INS).
- `angles_bframe`: Un diccionario que contiene los ángulos en el marco de referencia B.
- `ellipsoid`: (opcional) El tipo de elipsoide a utilizar, por defecto es "GRS80".
- `origin`: (opcional) El origen de las coordenadas, por defecto es "REC2".
- `magnetic_declination`: (opcional) La declinación magnética en grados, por defecto es 0.0.
- `geoid_undulation`: (opcional) La undulación del geoide en metros, por defecto es 0.0.

### Productos Esperados

Al utilizar la clase `GeoSensors`, se espera obtener:

- Coordenadas en diferentes representaciones, incluyendo:
  - Coordenadas en el marco de referencia E en decimal y radianes.
  - Coordenadas ECEF (Earth-Centered, Earth-Fixed).
  - Coordenadas locales que incluyen el origen y las restas necesarias.
  - Coordenadas en el marco de referencia N.
- Matrices de rotación que permiten transformar entre los diferentes marcos de referencia.
- Ángulos de rotación en los marcos de referencia E y B.


