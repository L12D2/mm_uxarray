Supported Datasets
==================

Supported Models and Observations are below. Please see
:ref:`Adding New Datasets <develop/datasets:Adding New Datasets>`
for advice on how to add new model and observational datasets to MELODIES MONET.

Supported Models
----------------

.. list-table:: Currently Connected Capabilities for Model Readers
   :widths: 20 20 30 30
   :header-rows: 1

   * - Model
     - Surface
     - Aircraft
     - Satellite
   * - `MERRA2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_
     - Yes
     - Needs testing
     - MODIS
   * - `WRF-Chem <https://www2.acom.ucar.edu/wrf-chem>`_
     - Yes
     - Yes
     - TROPOMI, TEMPO
   * - `CESM/CAM-chem FV <https://www2.acom.ucar.edu/gcm/cam-chem>`_
     - Yes
     - Needs testing
     - Needs testing
   * - `CESM/CAM-chem SE <https://www2.acom.ucar.edu/gcm/cam-chem>`_
     - Yes
     - | Needs testing & to 
       | add unstructured 
       | grid capabilities
     - | Needs testing & to 
       | add unstructured 
       | grid capabilities
   * - `CMAQ <https://www.epa.gov/cmaq/>`_
     - Yes
     - Needs testing
     - Needs testing
   * - `UFS-AQM (RRFS) <https://github.com/ufs-community/ufs-srweather-app/wiki/Air-Quality-Modeling>`_
     - Yes
     - Yes
     - Needs testing
   * - `CAMx <https://www.camx.com/>`_
     - Yes
     - Needs testing
     - TROPOMI, TEMPO
   * - `RAQMS <http://raqms-ops.ssec.wisc.edu/>`_
     - Yes
     - Needs testing
     - MOPITT, OMPS

In general, processing requires input to be in netCDF format. For the above 
models, scripts to configure the model data into a standard format for 
MELODIES MONET are available. If input datasets are in netCDF format and  
define latitude, longitude, altitude, and a datetime object, MELODIES MONET may be able 
to directly read the input files.

Please refer to the
`MELODIES MONET project board <https://github.com/orgs/NOAA-CSL/projects/6>`__ 
to learn more about our current and future development plans.

Supported Observations
----------------------

Surface
^^^^^^^
To use these surface datasets in MELODIES MONET specify "obs_type" equal to "pt_sfc" in your YAML file.

Available now:
   * `AirNow <https://www.airnow.gov/>`_
   * `AERONET <https://aeronet.gsfc.nasa.gov/>`_
   * `AQS <https://www.epa.gov/aqs/>`_
   * `ISH and ISH-Lite <https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database>`_
   * `OpenAQ <https://openaq.org/>`_

Under Development:
   * `IMPROVE <http://vista.cira.colostate.edu/Improve/>`_
   * `CRN <https://www.ncdc.noaa.gov/crn/>`_
   * `TOLNet <https://www-air.larc.nasa.gov/missions/TOLNet/>`_
   * `CEMS <https://www.epa.gov/emc/emc-continuous-emission-monitoring-systems/>`_

Please refer to the
`MELODIES MONET project board <https://github.com/orgs/NOAA-CSL/projects/6>`__ 
under milestone "Surface and Aircraft Evaluation Version 2" to learn more about our current and future development plans.

.. note::

   The :doc:`/cli` can be used to download and create MELODIES MONET-ready datasets for:
   AirNow, AERONET, AQS, ISH, ISH-Lite, and OpenAQ.

Aircraft, Sonde, Mobile, and Ground Campaign Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pairing capabilities include time, horizontal, and vertical interpolation. Users can evaluate aircraft data, 
ozonesonde data, mobile or walking data, and single ground site data. To use these options in MELODIES MONET 
specify "obs_type" equal to "aircraft", "ozone_sonder", "mobile", or "ground" in your YAML file. The table 
below describes these options in more detail. Available datafile formats include NetCDF, ICARTT, and CSV.

.. list-table:: Description of YAML File Options for "obs_type" For Campaign Data
   :widths: 70 30
   :header-rows: 1

   * - "obs_type"
     - Description
   * - "aircraft"
     - Aircraft - time, horizontal, and vertical interpolation
   * - "ozone_sonder"
     - Ozonesonde - time and vertical interpolation at a fixed horizontal location
   * - "mobile"
     - Mobile - time and horizontal interpolation at the surface
   * - "ground"
     - Ground - time interpolation at a fixed horizontal location at the surface

Tested datasets include the following: 
   * `FIREX-AQ <https://csl.noaa.gov/projects/firex-aq/>`_
   * `ATom <https://espo.nasa.gov/atom/content/ATom>`_
   * `SUNVEx <https://csl.noaa.gov/projects/sunvex/>`_
   * `AEROMMA <https://csl.noaa.gov/projects/aeromma/>`_
   * `ASIA-AQ <https://espo.nasa.gov/asia-aq>`_
   * `GML ozonesondes <https://gml.noaa.gov/ozwv/ozsondes/>`_

MELODIES MONET is written generally enough that other field campaign datasets should work well in the tool with 
minimal adjustments.

Please refer to the
`MELODIES MONET project board <https://github.com/orgs/NOAA-CSL/projects/6>`__ 
under milestone "Surface and Aircraft Evaluation Version 2" to learn more about our current and future development plans.

Satellite (under development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the
`MELODIES MONET project board <https://github.com/orgs/NOAA-CSL/projects/6>`__ 
under milestone "Satellite Evaluation Version 2" to learn more about our current and future development plans.
