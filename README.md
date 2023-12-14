# cpp_nav_filters
A collection of C-Plus-Plus INS filters to perform fusion of GNSS and INS nav solutions

## To-Do
- [ ] Loosely Coupled INS EKF
- [ ] Tightly Coupled INS EKF
- [ ] GPS Receiver PV Solution
- [ ] Utility
    - [ ] SV State from Ephemeris calculator
        - [ ] Calculate Satellite Pos and Vel
        - [ ] Calculate Satellite clock corretions
    - [ ] Intake receiver PVC and Satellite PV and return matrix of unit vectors
- [ ] Common 
    - [ ] Struct of ephemeris for each satellite
        - [ ] Struct contains vectors of variable dimension for each ephemeride
        - [ ] A vector of same variable dimension that contains ID for each SV
- [ ] Class Data to .csv parser
## Test Data
Test data is pulled from Auburn University Fundamentals of GPS Navigation lab data as it is the data I am most familiar with:
[Fundamentals of GPS Navigation Class Website](https://eng.auburn.edu/~dmbevly/fund_gps/)
