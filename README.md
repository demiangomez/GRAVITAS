# gravitas
A Matlab package to compute the gravity differences between stations of multiple gravity lines, and combine them into a network adjustment

Authors: Kevin Ahlgren, Demián D. Gómez, Michael G. Bevis, Franco S. Sobrero, Jacob Heck, Arturo Echalar, Dana J. Caccamise II, Eric Kendrick, Paola Montenegro, Ariele Batistti, Lizeth Contreras Choque, Juan Carlos Catari, Roger Tinta Sallico, and Hernan Guerra Trigo

GRAVITAS is a Matlab graphic user interface and adjustment software developed by the Geodesy and Geodynamics group at the Division of Geodetic Science, The Ohio State University. It can be used to adjust relative gravity lines (surveyed with multiple relative gravimeters) and absolute gravity measurements, and perform a least-squares adjustment of large surface gravity networks.

GRAVITAS has three main modules:
    • Instruments: this module handles the calibration curves of the relative gravimeters used to survey the gravity lines. It allows to load the calibration curve of gravimeters, and remove or modify existing ones. It also provides statistics regarding the overall performance of the instruments, based on their adjustment residuals.
    • Lines: this module loads the relative gravity measurements collected by multiple gravimeters along a gravity line, and conducts individual line adjustments for each instrument. The module computes the average relative gravity reading at each station and applies the calibration and earth tides corrections to produce a reduced gravity measurement. The gravity line is then adjusted using the reduced measurements to estimate the gravity differences between consecutive stations and the gravimeter's drift rate. 
The user can compare the results from different instruments to evaluate their agreement, contrast the measurements timelines, and plot the raw observations to detect blunders. 
Every adjusted gravity line is stored as an individual Matlab structure that contains the list of instruments used in the line survey, the raw observations, the measurement timestamps, the adjusted gravity differences, and the drift rates estimates. This module allows to create, edit, or remove benchmarks, modify their properties, and compute their coordinates with PPP by directly importing a RINEX file.
    • Gravity: this module ingests the results from all stored gravity lines, builds the design matrix, introduces absolute gravity values as stochastic condition equations, performs the network adjustment, and computes free-air and Bouguer anomalies. It allows the user to: 
        1. Add/remove absolute gravity stations 
        2. Assign an absolute gravity value to an existing benchmark
        3. Visualize the network geometry on a map
        4. Include/exclude specific lines from the network adjustment
        5. Export the adjustment results 

The line adjustment is performed individually for each gravimeter within the “Lines” module. Each adjustment produces the adjusted gravity differences between consecutive stations, and the gravimeter drift rate.

The network least-squares adjustment is performed within the “Gravity” module. The adjustment combines all the relative gravity pseudoobservations from the line adjustments, and the absolute gravity values.

The geoid undulation of each gravity station is computed using the Earth Gravitational Model EGM08, and exported into the gravity solution file. The gravimetric anomalies are computed using the normal gravity from Somigliana Formula and the values from the WGS84 ellipsoid. 

The gravity solution contains the list of adjusted stations with their attributes:
    • Station ID
    • Geodetic coordinates latitude, longitude, and height (WGS84 ellipsoid)
    • Gravity value (mGal)
    • Uncertainty (mGal)
    • Free-Air anomaly (mGal)
    • Bouguer anomaly (mGal)
    • Geoid undulation (m)

Some of the abilities of GRAVITAS are:
    • Adjust gravity networks with thousands of stations
    • Handle gravity lines missing forward or reverse observations
    • Plot gravity line residuals, and flag those larger than 3 sigmas
    • Visualize the geometry of the gravity network on a map
    • Visualize and select individual gravity lines on a map
    • Export the network database (benchmarks, absolute gravity stations, gravity lines)
    • Export the gravity solution in .txt and .kml formats
    • Exclude stations without precise coordinates from the solution file
    • Flag gravity lines with large uncertainties
    • Produce histograms of residuals for the entire network 
    • Produce histograms of residuals per gravimeter

GRAVITAS requires the following dependencies:
    • m_map v1.4j or above (https://www.eoas.ubc.ca/~rich/map.html)
    • Earth Gravitational Model EGM2008 (https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84) 
    • NRCan-CGS GPSPACE and its dependencies to obtain GNSS Precise Point Positioning (PPP) solutions. (https://github.com/pomath/GPSPACE)
