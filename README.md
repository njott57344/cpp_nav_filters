# cpp_nav_filters
A collection of C-Plus-Plus INS filters to perform fusion of GNSS and INS nav solutions
To-Do:
- [] Loosely Coupled INS EKF
- [] Tightly Coupled INS EKF
- [] GPS Receiver PV Solution
- [] Utility
    - [] SV State from Ephemeris calculator
        - [] Calculate Satellite Pos and Vel
        - [] Calculate Satellite clock corretions
    - [] Intake receiver PVC and Satellite PV and return matrix of unit vectors
- [] Common 
    - [] Struct of ephemeris for each satellite
        - [] Struct contains vectors of variable dimension for each ephemeride
        - [] A vector of same variable dimension that contains ID for each SV     
